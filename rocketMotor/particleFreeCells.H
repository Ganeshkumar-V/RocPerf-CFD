// Find particle free cells (alpha < 1e-10) <- temperature is set to default 300 to those cells
label particleFreeCellSize = 0;
{
    const volScalarField alphaP(phases[1]);
    forAll(alphaP, i)
    {
        if(alphaP[i] < 1e-10) particleFreeCellSize++;
    } 
}
labelList particleFreeCells(particleFreeCellSize);
const volScalarField& gasTemp(gasPhase.thermo().T());
scalarField setParticleTemp(particleFreeCellSize, 300);
{
    const volScalarField alphaP(phases[1]);
    label j = 0;
    forAll(alphaP, i)
    {
        if(alphaP[i] < 1e-10) 
        {
            particleFreeCells[j] = i; 
            setParticleTemp[j] = gasTemp[i]; 
            j++; 
        }
    }
}