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

// // Write Surface Area of all faces
// surfaceVectorField Sf
// (
//     IOobject
//     (
//         "Sf",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     mesh.Sf()
// );
// if (runTime.writeNow())
// {
//   Info << "Surface Area is written" << exit(FatalError);
// }

// write gas Density
volScalarField rhoG
(
  IOobject
  (
    "rhoG",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::AUTO_WRITE
  ),
  phases[0].thermo().rho()
);
if (runTime.writeNow())
{
  Info << "Density is written " << exit(FatalError);
}
