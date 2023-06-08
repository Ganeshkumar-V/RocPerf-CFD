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
#include "phasePair.H"
#include "GeometricField.H"

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


        Mach = mag(phases[0].U())/sqrt(phases[0].thermo().gamma()*p/phases[0].thermo().rho());

        runTime.write();

        runTime.printExecutionTime(Info);

        // Check Mass and Energy Conservation
        // {
        //   // 3 - Inlet , 4 - Outlet
        //   const volScalarField& alphag(phases[0]);
        //   const volScalarField& rhog(phases[0].thermo().rho());
        //   const volScalarField& Tg(phases[0].thermo().T());
        //   const scalarField& Cpg(phases[0].thermo().Cpv()().boundaryField()[3]);
        //   const vectorField& Ugin(phases[0].U()().boundaryField()[3]);
        //   const vectorField& Ugout(phases[0].U()().boundaryField()[4]);
        //
        //   const volScalarField& alphap(phases[1]);
        //   const volScalarField& rhop(phases[1].thermo().rho());
        //   const volScalarField& Tp(phases[1].thermo().T());
        //   const scalarField& Cpp(phases[1].thermo().Cpv()().boundaryField()[3]);
        //   const vectorField& Upin(phases[1].U()().boundaryField()[3]);
        //   const vectorField& Upout(phases[1].U()().boundaryField()[4]);
        //
        //   const surfaceScalarField& magSf(mesh.magSf());
        //
        //   const Field<scalar> mgin
        //   (
        //     alphag.boundaryField()[3]*rhog.boundaryField()[3]*magSf.boundaryField()[3]*mag(Ugin)
        //   );
        //
        //   const Field<scalar> mpin
        //   (
        //     alphap.boundaryField()[3]*rhop.boundaryField()[3]*magSf.boundaryField()[3]*mag(Upin)
        //   );
        //
        //   const Field<scalar> mgout
        //   (
        //     alphag.boundaryField()[4]*rhog.boundaryField()[4]*magSf.boundaryField()[4]*mag(Ugout)
        //   );
        //
        //   const Field<scalar> mpout
        //   (
        //     alphap.boundaryField()[4]*rhop.boundaryField()[4]*magSf.boundaryField()[4]*mag(Upout)
        //   );
        //   const Field<scalar> Egin
        //
        //   (
        //     mgin*Cpg*Tg.boundaryField()[3] + 0.5*mgin*magSqr(Ugin)
        //   );
        //
        //   const Field<scalar> Epin
        //   (
        //     mpin*Cpp*Tp.boundaryField()[3] + 0.5*mpin*magSqr(Upin)
        //   );
        //
        //   const Field<scalar> Egout
        //   (
        //     mgout*Cpg*Tg.boundaryField()[4] + 0.5*mgout*magSqr(Ugout)
        //   );
        //
        //   const Field<scalar> Epout
        //   (
        //     mpout*Cpp*Tp.boundaryField()[4] + 0.5*mpout*magSqr(Upout)
        //   );
        //
        //   Info << "Mass Conservation Statistics: " << endl;
        //   Info << "  mgin = " << sum(mgin) << "   mpin = " << sum(mpin) << "   massIn = " << sum(mgin + mpin) << endl;
        //   Info << "  mgout = " << sum(mgout) << "   mpout = " << sum(mpout) << "   massOut = " << sum(mgout + mpout) << endl;
        //   Info << "  Difference = " << sum(mgin + mpin) - sum(mgout + mpout) << endl;
        //   Info << endl;
        //   Info << "Energy Conservation Statistics: " << endl;
        //   Info << "  Egin = " << sum(Egin) << "   Epin = " << sum(Epin) << "   EnergyIn = " << sum(Egin + Epin) << endl;
        //   Info << "  Egout = " << sum(Egout) << "   Epout = " << sum(Epout) << "   EnergyOut = " << sum(Egout + Epout) << endl;
        //   Info << "  Difference = " << sum(Egin + Epin) - sum(Egout + Epout) << endl;
        //   Info << endl;
        // }
        // {
        //   // 3 - Inlet , 4 - Outlet
        //   const volScalarField& alphag(phases[0]);
        //
        //   const tmp<volScalarField> trhog(phases[0].thermo().rho());
        //   const volScalarField& rhog(trhog());
        //
        //   const tmp<volScalarField> tTg(phases[0].thermo().T());
        //   const volScalarField& Tg(tTg());
        //
        //   const tmp<volScalarField> tCpg(phases[0].thermo().Cpv());
        //   const volScalarField& CpgF(tCpg());
        //   const scalarField& Cpg(CpgF.boundaryField()[3]);
        //
        //   const tmp<volVectorField> tUg(phases[0].U());
        //   const volVectorField& UgF(tUg());
        //   const vectorField& Ugin(UgF.boundaryField()[3]);
        //   const vectorField& Ugout(UgF.boundaryField()[4]);
        //
        //   const surfaceVectorField& Sf(mesh.Sf());
        //
        //   const Field<scalar> mgin
        //   (
        //     alphag.boundaryField()[3]*rhog.boundaryField()[3]*(Sf.boundaryField()[3]&Ugin)
        //   );
        //
        //   const Field<scalar> mgout
        //   (
        //     alphag.boundaryField()[4]*rhog.boundaryField()[4]*(Sf.boundaryField()[4]&Ugout)
        //   );
        //
        //   const Field<scalar> mgwall
        //   (
        //     alphag.boundaryField()[2]*rhog.boundaryField()[2]*(Sf.boundaryField()[2]&UgF.boundaryField()[2])
        //   );
        //
        //   const Field<scalar> mgfront
        //   (
        //     alphag.boundaryField()[0]*rhog.boundaryField()[0]*(Sf.boundaryField()[0]&UgF.boundaryField()[0])
        //   );
        //
        //   const Field<scalar> mgback
        //   (
        //     alphag.boundaryField()[1]*rhog.boundaryField()[1]*(Sf.boundaryField()[1]&UgF.boundaryField()[1])
        //   );
        //
        //   const Field<scalar> Egin
        //
        //   (
        //     mgin*Cpg*Tg.boundaryField()[3] + 0.5*mgin*magSqr(Ugin)
        //   );
        //
        //   const Field<scalar> Egout
        //   (
        //     mgout*Cpg*Tg.boundaryField()[4] + 0.5*mgout*magSqr(Ugout)
        //   );
        //
        //   Info << "Mass Conservation Statistics: " << endl;
        //   Info << "  mgin = " << sum(mgin) <<  endl;
        //   Info << "  mgout = " << sum(mgout) <<  endl;
        //   Info << "  mgwall = " << sum(mgwall) <<  endl;
        //   Info << "  mgfront = " << sum(mgfront) <<  endl;
        //   Info << "  mgback = " << sum(mgback) <<  endl;
        //   Info << "  Difference = " << sum(mgin) + sum(mgout) + + sum(mgwall) + sum(mgfront) + sum(mgback) << endl;
        //   Info << endl;
        //   Info << "Energy Conservation Statistics: " << endl;
        //   Info << "  Egin = " << sum(Egin) << endl;
        //   Info << "  Egout = " << sum(Egout) << endl;
        //   Info << "  Difference = " << sum(Egin) + sum(Egout) << endl;
        //   Info << endl;
        // }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
