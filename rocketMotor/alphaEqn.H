{
  PtrList<volScalarField> dmdts(fluid.dmdts());

  word particlePhase = fluid.get<word>("particlePhase");
  word gasPhase = fluid.get<word>("gasPhase");

  volScalarField dalpha(0*fluid.phases()[0]);

  forAll(fluid.phases(), phasei)
  {
    phaseModel& phase = fluid.phases()[phasei];
    volScalarField& alpha = phase;

    if (phase.name() == gasPhase) continue;
    // Alpha Equation solved only for dispersed Phase
    else if (phase.name() == particlePhase)
    {
      const volScalarField& rho = phase.rho();
      const surfaceScalarField rhoPhi(fvc::interpolate(rho)*phase.phi());

      fvScalarMatrix alphaEqn
      (
        fvm::ddt(rho, alpha)
        + fvm::div(rhoPhi, alpha)
        == dmdts[phasei]
      );
      alphaEqn.solve();
    }
    else {}

    // Summation of dispersed phase to calculate gas alpha
    dalpha += alpha;
  }

  // Volume fraction of continuous phase
  forAll(fluid.phases(), phasei)
  {
    phaseModel& phase = fluid.phases()[phasei];
    if(phase.stationary()) continue;
    
    volScalarField& alpha = phase;

    if (phase.name() == gasPhase)
    {
      alpha = 1.0 - dalpha;
    }

    phase.alphaRhoPhiRef() = fvc::interpolate(phase.rho())*phase.alphaPhi();
    phase.clip(SMALL, 1 - SMALL);
  }

}
