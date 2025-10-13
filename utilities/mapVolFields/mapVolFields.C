/* ***************************************
Coarse to Fine Mapping - Volume Weighted
(c) 2025 Ganeshkumar V
**************************************** */

#include "fvCFD.H"
#include "meshToMesh.H"
#include "GeometricField.H"
#include "IOobjectList.H"

template<class Type> 
void MapVolFields
(
    const fvMesh& meshSource, 
    const fvMesh& meshTarget,
    const meshToMesh& mapper
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    IOobjectList objects(meshSource, meshSource.time().timeName());
    IOobjectList fields = objects.lookupClass(fieldType::typeName);

    forAllConstIters(fields, fieldIter)
    {
        Info << "    Mapping " << (*fieldIter)->name() << endl;

        fieldType fieldSource(*fieldIter(), meshSource, false);

        IOobject fieldTargetIOobject
        (
            (*fieldIter)->name(),
            meshTarget.time().timeName(),
            meshTarget,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        );

        if (fieldTargetIOobject.typeHeaderOk<fieldType>(true))
        {
            
            fieldType fieldTarget
            (
                fieldTargetIOobject, 
                mapper.mapSrcToTgt<Type>(fieldSource)
            );

            // Write field
            fieldTarget.write();
        }
        else
        {
            fieldTargetIOobject.readOpt(IOobject::NO_READ);

            // Interpolate field
            fieldType fieldTarget
            (
                fieldTargetIOobject,
                mapper.mapSrcToTgt<Type>(fieldSource)
            );

            // Write field
            fieldTarget.write();
        }

    }

    Info << endl;
}

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Map volume fields from one mesh to another"
    );
    argList::noParallel();
    argList::addArgument("sourceCase");

    argList::addOption
    (
        "sourceTime",
        "scalar|'latestTime'",
        "Specify the source time"
    );

    argList::addBoolOption
    (
        "d",
        "Enable debug output"
    );

    argList args(argc, argv);
    if (!args.check())
    {
        FatalError.exit();
    }
    #include "foamDlOpenLibs.H"

    fileName rootDirTarget(args.rootPath());
    fileName caseDirTarget(args.globalCaseName());

    const auto casePath = args.get<fileName>(1);
    const fileName rootDirSource = casePath.path().toAbsolute();
    const fileName caseDirSource = casePath.name();
    const bool debug = args.found("d");

    Info<< "Source: " << rootDirSource << " " << caseDirSource << endl;
    word sourceRegion = polyMesh::defaultRegion;
    if (args.found("sourceRegion"))
    {
        sourceRegion = args["sourceRegion"];
        Info<< "Source region: " << sourceRegion << endl;
    }

    Info<< "Target: " << rootDirTarget << " " << caseDirTarget << endl;
    word targetRegion = polyMesh::defaultRegion;
    if (args.found("targetRegion"))
    {
        targetRegion = args["targetRegion"];
        Info<< "Target region: " << targetRegion << endl;
    }

    #include "createTimes.H"

    HashTable<word> patchMap;
    wordList cuttingPatches;

    #include "setTimeIndex.H"

    Info<< "Create meshes\n" << endl;

    fvMesh meshSource
    (
        IOobject
        (
            sourceRegion,
            runTimeSource.timeName(),
            runTimeSource
        )
    );

    fvMesh meshTarget
    (
        IOobject
        (
            targetRegion,
            runTimeTarget.timeName(),
            runTimeTarget
        )
    );

    Info<< "Source mesh size: " << meshSource.nCells() << tab
        << "Target mesh size: " << meshTarget.nCells() << nl << endl;

    
    // Create Mesh-To-Mesh mapper
    meshToMesh mapper(meshSource, meshTarget, meshToMesh::interpolationMethod::imCellVolumeWeight);

    // --- Retrieve addressing and weights -- In Debug Mode 
    if (debug)
    {
        const labelListList& tgtToSrcCellAddr = mapper.tgtToSrcCellAddr();
        const scalarListList& tgtToSrcCellWght = mapper.tgtToSrcCellWght();

        Info<< "\n--- Mapping Information ---" << endl;
        forAll(tgtToSrcCellAddr, i)
        {
            Info<< "Target cell " << i << " receives contribution from: ";

            const labelList& addrList = tgtToSrcCellAddr[i];
            const scalarList& weightList = tgtToSrcCellWght[i];

            forAll(addrList, j)
            {
                Info<< " (src=" << addrList[j]
                    << ", w=" << weightList[j] << ")";
            }
            Info << endl;
        }
        Info<< "\nTotal target cells: " << tgtToSrcCellAddr.size() << endl;
    }

    // Map the fields
    Info<< nl
        << "Mapping scalar and vector volFields for time "
        << meshSource.time().timeName() << nl << endl;
    MapVolFields<scalar>(meshSource, meshTarget, mapper);
    MapVolFields<vector>(meshSource, meshTarget, mapper);

    Info<< "\nEnd\n" << endl;

    return 0;
}
