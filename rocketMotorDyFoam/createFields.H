#include "createRDeltaT.H"
#include "readGravitationalAcceleration.H"
#include "readhRef.H"

Info<< "Creating phaseSystem\n" << endl;

autoPtr<multiPhaseSystem> fluidPtr
(
    multiPhaseSystem::New(mesh)
);
multiPhaseSystem& fluid = fluidPtr();
multiPhaseSystem::phaseModelList& phases = fluid.phases();

dimensionedScalar pMin
(
    "pMin",
    dimPressure,
    fluid
);

#include "gh.H"

volScalarField& p = phases[0].thermoRef().p();

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);
mesh.setFluxRequired(p_rgh.name());

word particlePhase = fluid.get<word>("particlePhase");
word gasPhase = fluid.get<word>("gasPhase");

label particleIndex = -1;
label gasIndex = -1;
label propellantIndex = -1;
forAll(fluid.phases(), phasei)
{
  phaseModel& phase = fluid.phases()[phasei];
  if(phase.name() == gasPhase)
      gasIndex = phasei;
  else if (phase.name() == particlePhase)
      particleIndex = phasei;
  else
      propellantIndex = phasei;
}
