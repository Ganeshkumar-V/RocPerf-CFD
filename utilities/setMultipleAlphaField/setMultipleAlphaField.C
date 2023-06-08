/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 DHI
    Copyright (C) 2017-2020 OpenCFD Ltd.
    Copyright (C) 2017-2020 German Aerospace Center (DLR)
    Copyright (C) 2020 Johan Roenby
    Copyright (C) 2023 Ganeshkumar
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    setMultipleAlphaField

Description
    Uses cutCellIso to create a volume fraction field from either a cylinder,
    a sphere or a plane.

    Original code supplied by Johan Roenby, DHI (2016)
    Modification Henning Scheufler, DLR (2019)
    Modification Ganeshkumar, IITGN (2023)

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "triSurface.H"
#include "triSurfaceTools.H"

#include "implicitFunction.H"

#include "cutCellIso.H"
#include "cutFaceIso.H"
#include "OBJstream.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void isoFacesToFile
(
    const DynamicList<List<point>>& faces,
    const word& fileName
)
{
    // Writing isofaces to OBJ file for inspection in paraview
    OBJstream os(fileName + ".obj");

    if (Pstream::parRun())
    {
        // Collect points from all the processors
        List<DynamicList<List<point>>> allProcFaces(Pstream::nProcs());
        allProcFaces[Pstream::myProcNo()] = faces;
        Pstream::gatherList(allProcFaces);

        if (Pstream::master())
        {
            Info<< "Writing file: " << fileName << endl;

            for (const DynamicList<List<point>>& procFaces : allProcFaces)
            {
                for (const List<point>& facePts : procFaces)
                {
                    os.write(face(identity(facePts.size())), facePts);
                }
            }
        }
    }
    else
    {
        Info<< "Writing file: " << fileName << endl;

        for (const List<point>& facePts : faces)
        {
            os.write(face(identity(facePts.size())), facePts);
        }
    }
}

void setAlpha
(
    volScalarField& alpha1,
    DynamicList<List<point>>& facePts,
    scalarField& f,
    const bool writeOBJ
)
{
    const fvMesh& mesh = alpha1.mesh();
    cutCellIso cutCell(mesh, f);
    cutFaceIso cutFace(mesh, f);

    forAll(alpha1, cellI)
    {
        cutCell.calcSubCell(cellI, 0.0);

        alpha1[cellI] = max(min(cutCell.VolumeOfFluid(), 1), 0);

        if (writeOBJ && (mag(cutCell.faceArea()) >= 1e-14))
        {
            facePts.append(cutCell.facePoints());
        }
    }

    // Setting boundary alpha1 values
    forAll(mesh.boundary(), patchi)
    {
        if (mesh.boundary()[patchi].size() > 0)
        {
            const label start = mesh.boundary()[patchi].patch().start();
            scalarField& alphap = alpha1.boundaryFieldRef()[patchi];
            const scalarField& magSfp = mesh.magSf().boundaryField()[patchi];

            forAll(alphap, patchFacei)
            {
                const label facei = patchFacei + start;
                cutFace.calcSubFace(facei, 0.0);
                alphap[patchFacei] =
                    mag(cutFace.subFaceArea())/magSfp[patchFacei];
            }
        }
    }
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Uses cutCellIso to create a volume fraction field from an "
        "implicit function."
    );

    #include "addDictOption.H"
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    const word dictName("setMultipleAlphaFieldDict");
    #include "setSystemMeshDictionaryIO.H"

    IOdictionary setMultipleAlphaFieldDict(dictIO);

    Info<< "Reading " << setMultipleAlphaFieldDict.name() << endl;

    const word fieldName = setMultipleAlphaFieldDict.get<word>("field");
    const bool invert = setMultipleAlphaFieldDict.getOrDefault("invert", false);
    const bool writeOBJ = setMultipleAlphaFieldDict.getOrDefault("writeOBJ", false);
    const word invertFieldName = setMultipleAlphaFieldDict.get<word>("otherField");

    Info<< "Reading field " << fieldName << nl << endl;
    volScalarField alpha1
    (
        IOobject
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    autoPtr<implicitFunction> func = implicitFunction::New
    (
        setMultipleAlphaFieldDict.get<word>("type"),
        setMultipleAlphaFieldDict
    );

    scalarField f(mesh.nPoints(), Zero);

    forAll(f, pi)
    {
        f[pi] = func->value(mesh.points()[pi]);
    };

    DynamicList<List<point>> facePts;

    setAlpha(alpha1, facePts, f, writeOBJ);

    if (writeOBJ)
    {
        isoFacesToFile(facePts, fieldName + "0");
    }

    ISstream::defaultPrecision(18);

    if (invert)
    {
        alpha1 = scalar(1) - alpha1;
    }

    Info<< "Writing new alpha field " << alpha1.name() << endl;
    alpha1.write();

    Info<< "Creating other field " << invertFieldName << nl << endl;
    volScalarField alpha2
    (
        IOobject
        (
            invertFieldName,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        scalar(1) - alpha1
    );
    alpha2.write();

    const scalarField& alpha = alpha1.internalField();

    Info<< "sum(alpha*V):" << gSum(mesh.V()*alpha)
        << ", 1-max(alpha1): " << 1 - gMax(alpha)
        << " min(alpha1): " << gMin(alpha) << endl;

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
