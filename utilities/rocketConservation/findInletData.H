{
// OutletPatch Fields
label inIndex = -1;
forAll(mesh.boundary(), bFi)
{
    if (mesh.boundary()[bFi].name() == "inlet")
    {
        inIndex = mesh.boundary()[bFi].index();
    }
}
const scalarField& pF(p.boundaryField()[inIndex]);

const scalarField& alphaParticlesF(alphaParticles.boundaryField()[inIndex]);
// const scalarField& TparticlesF(Tparticles.boundaryField()[inIndex]);
const scalarField& phiParticlesF(phiParticles.boundaryField()[inIndex]);
vectorField UparticlesF(Uparticles.boundaryField()[inIndex]);
const scalarField& hparticlesF(hparticles.boundaryField()[inIndex]);

const scalarField alphaGasF(1.0 - alphaParticlesF);
const scalarField& TgasF(Tgas.boundaryField()[inIndex]);
const scalarField& phiGasF(phiGas.boundaryField()[inIndex]);
const vectorField& UgasF(Ugas.boundaryField()[inIndex]);
const scalarField& hgasF(hgas.boundaryField()[inIndex]);

const volScalarField W(phases[0].thermo().W());
const scalarField rhoGas(pF*W.boundaryField()[inIndex]/(8314*TgasF));
const scalarField rhoParticles(phases[1].thermo().rho()().boundaryField()[inIndex]);

// const vectorField& Sf(mesh.boundary()[inIndex].Sf());

// Mass flow rates and momentum flow rates
const scalarField mdotGas(alphaGasF*rhoGas*phiGasF);
const scalarField mdotparticles(alphaParticlesF*rhoParticles*phiParticlesF);
const vectorField Pgas(mdotGas*UgasF);
UparticlesF = UgasF;
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
    file.open(infileName, std::ios_base::app);
    file << inString.str();
    inString.str(std::string());
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
    inString << runTime.timeName() << ", "
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
