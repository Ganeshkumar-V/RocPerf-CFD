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
        if (phase.name() == "particles")
        {
          EEqn.setValues
          (
            particleFreeCells,
            phase.thermo().he(setParticleTemp, setParticleTemp, particleFreeCells)
          );
        }

        EEqn.solve();
        fvOptions.correct(phase.thermoRef().he());

        // Correct for Adiabatic Walls
        if (propellantIndex != -1)
        {
          ImposeWall(phase.thermoRef().he(), phases[propellantIndex]);
        }

        // Limit Temperature
        if (limitTemperature)
	{
          labelList cells(1);
          scalarField Tmax(cells.size(), 3000);
          scalarField Tmin(cells.size(), 300);
          scalarField heMax(phase.thermo().he(Tmax, Tmax, cells));
          scalarField heMin(phase.thermo().he(Tmin, Tmin, cells));
          limitT(phase.thermoRef().he(), heMax, heMin);
	}
    }

    fluid.correctThermo();
    fluid.correct();
}


forAll(phases, phasei)
{
    phaseModel& phase = phases[phasei];

    // // Limit Temperature
    // labelList cells(1);
    // scalarField Tmax(cells.size(), 2740);
    // scalarField Tmin(cells.size(), 500);
    // scalarField heMax(phase.thermo().he(Tmax, Tmax, cells));
    // scalarField heMin(phase.thermo().he(Tmin, Tmin, cells));
    // limitT(phase.thermoRef().he(), heMax, heMin);

    // fluid.correctThermo();
    // fluid.correct();

    Info<< phase.name() << " min/max T "
        << min(phase.thermo().T()).value()
        << " - "
        << max(phase.thermo().T()).value()
        << endl;
    /***
    scalar minT = max(phase.thermo().T()).value();
    if (true)
    {
    const volScalarField& T(phase.thermo().T());
    forAll(T, i)
    {
	    if (T[i] == minT)
	    {
		Pout << "minT at: " << i << " , ";
	    }
    }
    Pout << endl;
    }
    ***/
}
