/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    reactingMultiphaseEulerFoam

Description
    Solver for a system of any number of compressible fluid phases with a
    common pressure, but otherwise separate properties. The type of phase model
    is run time selectable and can optionally represent multiple species and
    in-phase reactions. The phase system is also run time selectable and can
    optionally represent different types of momentum, heat and mass transfer.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "multiPhaseSystem.H"
#include "phaseCompressibleTurbulenceModel.H"
#include "pimpleControl.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "AdiabaticWall.H"

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "createFields.H"
    #include "createFieldRefs.H"

    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    Switch partialElimination
    (
        pimple.dict().getOrDefault<Switch>("partialElimination", false)
    );

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    // Correcting Phase Volume Fractions
    forAll(fluid.phases(), phasei)
    {
      fluid.phases()[phasei].clip(SMALL, 1 - SMALL);
    }

    // Setting variables to initialize
    labelList purePropellantCells(0);
    scalarField setTemp(0, 0);
    scalarField setPressure(0, 0);
    vectorField setVelocity(0, vector(0, 0, 0));
    label propellantIndex = fluid.get<label>("propellantIndex");

    while (runTime.run())
    {
        #include "readTimeControls.H"

        int nEnergyCorrectors
        (
            pimple.dict().getOrDefault<int>("nEnergyCorrectors", 1)
        );

        #include "CourantNo.H"
        #include "setDeltaTFactor.H"

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Store Old Time
        fluid.store();

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            fluid.solve();
            fluid.correct();

            //***********  Start Find Propellant size ***********//
            if (propellantIndex != -1)
            {
              label purePropellantSize = 0;
              scalar cutoff = 0.999; //(1.0 - SMALL);
              {
                const volScalarField& propellant = phases[propellantIndex];

                forAll(propellant, i)
                {
                  if (propellant[i] >= cutoff)
                  {
                    purePropellantSize++;
                  }
                }
              }
              purePropellantCells = labelList(purePropellantSize);
              setTemp = scalarField
              (
                purePropellantSize,
                fluid.getOrDefault<scalar>("Tset", 2000)
              );
              setPressure = scalarField(purePropellantSize, 101325);
              setVelocity = vectorField(purePropellantSize, vector(0, 0, 0));
              {
                const volScalarField& propellant = phases[propellantIndex];

                label j = 0;
                forAll(propellant, i)
                {
                  if (propellant[i] >= cutoff)
                  {
                    purePropellantCells[j] = i;
                    j++;
                  }
                }
              }
            }
            //***********  End Find Propellant size ***********//

            #include "YEqns.H"

            #include "pU/UEqns.H"
            #include "EEqns.H"
            #include "pU/pEqn.H"

            fluid.correctKinematics();

            if (pimple.turbCorr())
            {
                fluid.correctTurbulence();
            }
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
