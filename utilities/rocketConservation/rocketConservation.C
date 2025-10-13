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

    std::ofstream file;

    // File to write time evolution outlet data
    std::string outfileName = "outlet.csv";
    std::remove(outfileName.c_str()); // delete if already present
    file.open(outfileName, std::ios_base::app);
    file << "t, " << "mdotG, " << "mdotP, " << "Mgx, " << "Mgy, " << "Mgz, " << "magMg, " << "Mpx, " << "Mpy, " << "Mpz, " << "magMp, " << "Hg, " << "Hp, " << "Kg, " << "Kp\n";
    file.close();
    std::stringstream outString;

    // File to write time evolution inlet data
    std::string infileName = "inlet.csv";
    std::remove(infileName.c_str()); // delete if already present
    file.open(infileName, std::ios_base::app);
    file << "t, " << "mdotG, " << "mdotP, " << "Mgx, " << "Mgy, " << "Mgz, " << "magMg, " << "Mpx, " << "Mpy, " << "Mpz, " << "magMp, " << "Hg, " << "Hp, " << "Kg, " << "Kp\n";
    file.close();
    std::stringstream inString;

    // File to write time evolution total forces (Pressure + shear stress)
    std::string forcefileName = "forces.csv";
    std::remove(forcefileName.c_str()); // delete if already present
    file.open(forcefileName, std::ios_base::app);
    file << "t, " << "Fpx, " << "Fpy, " << "Fpz, "<< "magFp, "<< "Fsx, " << "Fsy, "<< "Fsz, "<< "magFs, "<< "Q, " << "Ws\n";
    file.close();
    std::stringstream forceString;

    Info << "timeDirs: " << timeDirs.size() << endl;
    for (label timei = 0; timei < timeDirs.size(); ++timei)
    {

        // If outside start and end time ignore
        if (runTime.startTime().value() > timeDirs[timei].value()) continue;
        if (runTime.endTime().value() < timeDirs[timei].value()) continue;

        Info  << "Time: " << timeDirs[timei].value() << endl;
        runTime.setTime(timeDirs[timei], timei);

        #include "readFields.H"
        volScalarField hgas(phases[0].thermo().he(p, Tgas));
        volScalarField hparticles(phases[1].thermo().he(p, Tparticles));

        volScalarField kgas(phases[0].thermo().kappa());
        volScalarField mugas(phases[0].thermo().mu());
        volTensorField tau((1.0 - alphaParticles)*mugas*(fvc::grad(Ugas) + dev2(T(fvc::grad(Ugas)))));
        volVectorField qgas(-(1.0 - alphaParticles)*kgas*fvc::grad(Tgas));

        #include "findOutletData.H"
        #include "findInletData.H"
        #include "findTotalForce.H"

        // print mass conservation
        forAll(mesh.boundary(), bFi)
        {
            Info << mesh.boundary()[bFi].name();
            const scalarField& pF(p.boundaryField()[bFi]);
            const scalarField& alphaParticlesF(alphaParticles.boundaryField()[bFi]);
            const scalarField& phiParticlesF(phiParticles.boundaryField()[bFi]);
            vectorField UparticlesF(Uparticles.boundaryField()[bFi]);
            const scalarField& hparticlesF(hparticles.boundaryField()[bFi]);

            const scalarField alphaGasF(1.0 - alphaParticlesF);
            const scalarField& TgasF(Tgas.boundaryField()[bFi]);
            const scalarField& phiGasF(phiGas.boundaryField()[bFi]);
            const vectorField& UgasF(Ugas.boundaryField()[bFi]);
            const scalarField& hgasF(hgas.boundaryField()[bFi]);

            const volScalarField W(phases[0].thermo().W());
            const scalarField rhoGas(pF*W.boundaryField()[bFi]/(8314*TgasF));
            const scalarField rhoParticles(phases[1].thermo().rho()().boundaryField()[bFi]);
            const scalarField mdotGas(alphaGasF*rhoGas*phiGasF);
            const scalarField mdotparticles(alphaParticlesF*rhoParticles*phiParticlesF);

            scalar tmdotGas(sum(mdotGas));
            scalar tmdotparticles(sum(mdotparticles));
            Info << "Gas: " << tmdotGas << " Particles: " << tmdotparticles << " Total: " << tmdotGas+tmdotparticles << endl;
        }
    }

    // store remaining outlet data
    file.open(outfileName, std::ios_base::app);
    file << outString.str();
    file.close();

    // store remaining inlet data
    file.open(infileName, std::ios_base::app);
    file << inString.str();
    file.close();

    // store remaining force data
    file.open(forcefileName, std::ios_base::app);
    file << forceString.str();
    file.close();

    Info<< "End\n" << endl;

    return 0;
}
