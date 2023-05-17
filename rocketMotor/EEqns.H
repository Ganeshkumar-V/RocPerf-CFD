for (int Ecorr=0; Ecorr<nEnergyCorrectors; Ecorr++)
{
    fluid.correctEnergyTransport();

    autoPtr<phaseSystem::heatTransferTable>
        heatTransferPtr(fluid.heatTransfer());

    phaseSystem::heatTransferTable& heatTransfer = heatTransferPtr();

    forAll(fluid.anisothermalPhases(), anisothermalPhasei)
    {
        phaseModel& phase = fluid.anisothermalPhases()[anisothermalPhasei];

        const volScalarField& alpha = phase;
        const volScalarField& rho = phase.rho();
        const volVectorField& U = phase.U();

        fvScalarMatrix EEqn
        (
            phase.heEqn()
         ==
           *heatTransfer[phase.name()]
          + alpha*rho*(U&g)
          + fvOptions(alpha, rho, phase.thermoRef().he())
        );

        EEqn.relax();
        fvOptions.constrain(EEqn);

        // set Propellant phase Values
        if (propellantIndex != -1)
        {
          EEqn.setValues
          (
            purePropellantCells,
            phase.thermo().he(setTemp, setTemp, purePropellantCells)
          );
        }

        EEqn.solve();
        fvOptions.correct(phase.thermoRef().he());

        // Correct for Adiabatic Walls
        if (propellantIndex != -1)
        {
        ImposeWall(phase.thermoRef().he(), phases[propellantIndex]);
        }
    }

    forAll(phases, phasei)
    {
        phaseModel& phase = phases[phasei];

        // if (phasei == 1)
        // {
        //   scalar minT(min(phase.thermo().T()).value());
        //   volScalarField T(phase.thermo().T());
        //   forAll(T, celli)
        //   {
        //     if (T[celli] == minT)
        //     {
        //       Info << "Cell: " << celli << " has minT: " << minT << endl;
        //     }
        //   }
        //   // forAll(T.boundaryField(), bfi)
        //   // {
        //   //   forAll(T.boundaryField()[bfi], pi)
        //   //   {
        //   //     if (T.boundaryField()[bfi][pi] == minT)
        //   //     {
        //   //       Info << "boundaryField: " << T.mesh().boundary()[bfi].name() << " patchFace: " << pi << " has minT: " << minT << endl;
        //   //     }
        //   //   }
        //   // }
        //   // Info << T.boundaryField()[0] << endl;
        // }
        Info<< phase.name() << " min/max T "
            << min(phase.thermo().T()).value()
            << " - "
            << max(phase.thermo().T()).value()
            << endl;
    }

    fluid.correctThermo();
    fluid.correct();
}


forAll(phases, phasei)
{
    phaseModel& phase = phases[phasei];

    // if (phasei == 0)
    // {
    //   scalar minT(min(phase.thermo().T()).value());
    //   volScalarField T(phase.thermo().T());
    //   forAll(T, celli)
    //   {
    //     if (T[celli] == minT)
    //     {
    //       Info << "Cell: " << celli << " has minT: " << minT << endl;
    //     }
    //   }
    //   forAll(T.boundaryField(), bfi)
    //   {
    //     forAll(T.boundaryField()[bfi], pi)
    //     {
    //       if (T.boundaryField()[bfi][pi] == minT)
    //       {
    //         Info << "boundaryField: " << T.mesh().boundary()[bfi].name() << " patchFace: " << pi << " has minT: " << minT << endl;
    //       }
    //     }
    //   }
    //   Info << T.boundaryField()[0] << endl;
    // }
    Info<< phase.name() << " min/max T "
        << min(phase.thermo().T()).value()
        << " - "
        << max(phase.thermo().T()).value()
        << endl;
}
