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

if (runTime.value() == 0)
{
  // Write Surface Area of all faces
  surfaceVectorField Sf
  (
      IOobject
      (
          "Sf",
          runTime.timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
      ),
      mesh.Sf()
  );
  Sf.write();

  // Write Magnitude of surface Areas
  surfaceScalarField magSf
  (
      IOobject
      (
          "magSf",
          runTime.timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
      ),
      mesh.magSf()
  );
  magSf.write();

  // Write Control Volumes
  volScalarField::Internal V(mesh.V());
  V.write();

  // Write face centers
  surfaceVectorField Cf
  (
      IOobject
      (
          "Cf",
          runTime.timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
      ),
      mesh.Cf()
  );
  Cf.write();
}

// calculate Mach Number
volScalarField Mach
(
    IOobject
    (
      "Mach",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
    ),
    mag(phases[0].U())/sqrt(phases[0].thermo().gamma()*p/phases[0].thermo().rho())
);
