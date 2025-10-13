/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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
    postProcess

Description
    Execute the set of functionObjects specified in the selected dictionary
    (which defaults to system/controlDict) or on the command-line for the
    selected set of times on the selected set of fields.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "profiling.H"
#include "timeSelector.H"
#include "ReadFields.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "uniformDimensionedFields.H"
#include "fileFieldSelection.H"
#include "mapPolyMesh.H"
#include "fvCFD.H"
#include "multiPhaseSystem.H"
#include "pimpleControl.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "phasePair.H"
#include "GeometricField.H"
#include "scalarIOField.H"
#include "processorFvPatch.H"
#include <stdio.h>
#include <omp.h>
#include <sstream>

using namespace Foam;

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Execute the set of functionObjects specified in the selected"
        " dictionary or on the command-line for the"
        " selected set of times on the selected set of fields"
    );

    timeSelector::addOptions();
    #include "addProfilingOption.H"
    #include "addRegionOption.H"
    #include "addFunctionObjectOptions.H"

    // Set functionObject post-processing mode
    functionObject::postProcess = true;

    #include "setRootCase.H"

    if (args.found("list"))
    {
        functionObjectList::list();
        return 0;
    }

    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createNamedMesh.H"
    #include "readGravitationalAcceleration.H"
    #include "readhRef.H"

    Info<< "Creating phaseSystem\n" << endl;

    autoPtr<multiPhaseSystem> fluidPtr
    (
        multiPhaseSystem::New(mesh)
    );
    multiPhaseSystem& fluid = fluidPtr();
    multiPhaseSystem::phaseModelList& phases = fluid.phases();

    #include "gh.H"

    // File to write time evolution data
    std::string fileName = "performance.csv";
    std::remove(fileName.c_str()); // delete if already present
    std::ofstream file;
    file.open(fileName, std::ios_base::app);
    file << "t, " << "mdotG, " << "mdotP, " << "Fg, " << "Fp, " << "Fpressure, " << "Pc\n";
    file.close();
    std::stringstream output;

    Info << "timeDirs: " << timeDirs.size() << endl;
    // #pragma omp parallel for ordered
    for (label timei = 0; timei < timeDirs.size(); ++timei)
    // forAll(timeDirs, timei)
    {

        // If outside start and end time ignore
        if (runTime.startTime().value() > timeDirs[timei].value()) continue;
        if (runTime.endTime().value() < timeDirs[timei].value()) continue;

        Info  << "Time: " << timeDirs[timei].value() << endl;
        runTime.setTime(timeDirs[timei], timei);

        #include "readFields.H"

        // OutletPatch Fields
        label outIndex = -1;
        forAll(mesh.boundary(), bFi)
        {
            if (mesh.boundary()[bFi].name() == "outlet")
            {
                outIndex = mesh.boundary()[bFi].index();
                Info << "Outlet Index: " << outIndex << endl;
            }
        }
        const scalarField& pF(p.boundaryField()[outIndex]);

        const scalarField& alphaParticlesF(alphaParticles.boundaryField()[outIndex]);
        // const scalarField& TparticlesF(Tparticles.boundaryField()[outIndex]);
        const scalarField& phiParticlesF(phiParticles.boundaryField()[outIndex]);
        const vectorField& UparticlesF(Uparticles.boundaryField()[outIndex]);

        const scalarField alphaGasF(1.0 - alphaParticlesF);
        const scalarField& TgasF(Tgas.boundaryField()[outIndex]);
        const scalarField& phiGasF(phiGas.boundaryField()[outIndex]);
        const vectorField& UgasF(Ugas.boundaryField()[outIndex]);

        const volScalarField W(phases[0].thermo().W());
        const scalarField rhoGas(pF*W.boundaryField()[outIndex]/(8314*TgasF));
        const scalarField rhoParticles(phases[1].thermo().rho()().boundaryField()[outIndex]);

        const vectorField& Sf(mesh.boundary()[outIndex].Sf());

        // Mass flow rates and momentum flow rates
        const scalarField mdotGas(alphaGasF*rhoGas*phiGasF);
        const scalarField mdotparticles(alphaParticlesF*rhoParticles*phiParticlesF);
        const vectorField Pgas(mdotGas*UgasF);
        const vectorField Pparticles(mdotparticles*UparticlesF);

        // Total flow rates and thrust
        scalar tmdotGas(sum(mdotGas));
        scalar tmdotparticles(sum(mdotparticles));
        scalar tFgas(sum(Pgas.component(vector::Z)));
        scalar tFparticles(sum(Pparticles.component(vector::Z)));
        scalar tFpressure(sum((pF - 101325)*mag(Sf)));

        // write scalar values to the table
        if (timei%100 == 0)
        {
            file.open(fileName, std::ios_base::app);
            file << output.str();
            output.str(std::string());
            // #pragma omp ordered
            {
              file << runTime.timeName() << ", "
                   << tmdotGas << ", "
                   << tmdotparticles << ", "
                   << tFgas << ", "
                   << tFparticles << ", "
                   << tFpressure << ", "
                   << p[0] << ", " << "\n";
            }
            file.close();
        }
        else
        {
            output << runTime.timeName() << ", "
                   << tmdotGas << ", "
                   << tmdotparticles << ", "
                   << tFgas << ", "
                   << tFparticles << ", "
                   << tFpressure << ", "
                   << p[0] << ", " << "\n";
        }

    }
    // store remaining data
    file.open(fileName, std::ios_base::app);
    file << output.str();
    file.close();

    // Write cell dimensions
    #include "cellDim.H"
    Info<< "End\n" << endl;

    return 0;
}
