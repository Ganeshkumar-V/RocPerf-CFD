{
// OutletPatch Fields
label outIndex = -1;
forAll(mesh.boundary(), bFi)
{
    if (mesh.boundary()[bFi].name() == "outlet")
    {
        outIndex = mesh.boundary()[bFi].index();
    }
}
const scalarField& pF(p.boundaryField()[outIndex]);

const scalarField& alphaParticlesF(alphaParticles.boundaryField()[outIndex]);
// const scalarField& TparticlesF(Tparticles.boundaryField()[outIndex]);
const scalarField& phiParticlesF(phiParticles.boundaryField()[outIndex]);
const vectorField& UparticlesF(Uparticles.boundaryField()[outIndex]);
const scalarField& hparticlesF(hparticles.boundaryField()[outIndex]);

const scalarField alphaGasF(1.0 - alphaParticlesF);
const scalarField& TgasF(Tgas.boundaryField()[outIndex]);
const scalarField& phiGasF(phiGas.boundaryField()[outIndex]);
const vectorField& UgasF(Ugas.boundaryField()[outIndex]);
const scalarField& hgasF(hgas.boundaryField()[outIndex]);

const volScalarField W(phases[0].thermo().W());
const scalarField rhoGas(pF*W.boundaryField()[outIndex]/(8314*TgasF));
const scalarField rhoParticles(phases[1].thermo().rho()().boundaryField()[outIndex]);

// const vectorField& Sf(mesh.boundary()[outIndex].Sf());

// Mass flow rates and momentum flow rates
const scalarField mdotGas(alphaGasF*rhoGas*phiGasF);
const scalarField mdotparticles(alphaParticlesF*rhoParticles*phiParticlesF);
const vectorField Pgas(mdotGas*UgasF);
const vectorField Pparticles(mdotparticles*UparticlesF);
const scalarField Hgas(mdotGas*hgasF);
const scalarField Hparticles(mdotparticles*hparticlesF);
const scalarField Kgas(0.5*mdotGas*(UgasF & UgasF));
const scalarField Kparticles(0.5*mdotparticles*(UparticlesF & UparticlesF));

// Total flow rates, momentum and energy
scalar tmdotGas(sum(mdotGas));
scalar tmdotparticles(sum(mdotparticles));
vector tMgas(sum(Pgas));
vector tMparticles(sum(Pparticles));
scalar tHgas(sum(Hgas));
scalar tHparticles(sum(Hparticles));
scalar tKgas(sum(Kgas));
scalar tKparticles(sum(Kparticles));

// write scalar values to the table
if (timei%100 == 0)
{
    file.open(outfileName, std::ios_base::app);
    file << outString.str();
    outString.str(std::string());
    {
      file << runTime.timeName() << ", "
           << tmdotGas << ", "
           << tmdotparticles << ", "
           << tMgas.x() << ", "
           << tMgas.y() << ", "
           << tMgas.z() << ", "
           << mag(tMgas) << ", "
           << tMparticles.x() << ", "
           << tMparticles.y() << ", "
           << tMparticles.z() << ", "
           << mag(tMparticles) << ", "
           << tHgas << ", "
           << tHparticles << ", "
           << tKgas << ", "
           << tKparticles << "\n";
    }
    file.close();
}
else
{
    outString << runTime.timeName() << ", "
           << tmdotGas << ", "
           << tmdotparticles << ", "
           << tMgas.x() << ", "
           << tMgas.y() << ", "
           << tMgas.z() << ", "
           << mag(tMgas) << ", "
           << tMparticles.x() << ", "
           << tMparticles.y() << ", "
           << tMparticles.z() << ", "
           << mag(tMparticles) << ", "
           << tHgas << ", "
           << tHparticles << ", "
           << tKgas << ", "
           << tKparticles << "\n";
}
}
