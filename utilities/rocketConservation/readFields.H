volScalarField p
(
    IOobject
    (
      "p_rgh",
      runTime.timeName(),
      mesh,
      IOobject::MUST_READ,
      IOobject::NO_WRITE
    ),
    mesh
);

volScalarField Tgas
(
    IOobject
    (
      "T.gas",
      runTime.timeName(),
      mesh,
      IOobject::MUST_READ,
      IOobject::NO_WRITE
    ),
    mesh
);

volScalarField Tparticles
(
    IOobject
    (
      "T.particles",
      runTime.timeName(),
      mesh,
      IOobject::MUST_READ,
      IOobject::NO_WRITE
    ),
    mesh
);

volScalarField alphaParticles
(
    IOobject
    (
      "alpha.particles",
      runTime.timeName(),
      mesh,
      IOobject::MUST_READ,
      IOobject::NO_WRITE
    ),
    mesh
);

volVectorField Ugas
(
    IOobject
    (
      "U.gas",
      runTime.timeName(),
      mesh,
      IOobject::MUST_READ,
      IOobject::NO_WRITE
    ),
    mesh
);
// Ugas.correctBoundaryConditions();

volVectorField Uparticles
(
    IOobject
    (
      "U.particles",
      runTime.timeName(),
      mesh,
      IOobject::MUST_READ,
      IOobject::NO_WRITE
    ),
    mesh
);

surfaceScalarField phiGas
(
    IOobject
    (
      "phi.gas",
      runTime.timeName(),
      mesh,
      IOobject::READ_IF_PRESENT,
      IOobject::NO_WRITE
    ),
    fvc::flux(Ugas)
);

surfaceScalarField phiParticles
(
    IOobject
    (
      "phi.particles",
      runTime.timeName(),
      mesh,
      IOobject::READ_IF_PRESENT,
      IOobject::NO_WRITE
    ),
    fvc::flux(Uparticles)
);
