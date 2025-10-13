// Find propellant cells -> velocity and temperature set to a default value at the propellant cells
if (propellantIndex != -1)
{
    label purePropellantSize = 0;
    scalar cutoff = 0.999; // (1.0 - SMALL) is another alternative;
    {
        const volScalarField& propellant = phases[propellantIndex];
        forAll(propellant, i)
        {
            if (propellant[i] >= cutoff)
            {
                purePropellantSize++;
            }
        }
    }
    purePropellantCells = labelList(purePropellantSize);
    setTemp = scalarField
    (
        purePropellantSize,
        fluid.getOrDefault<scalar>("Tset", 2000)
    );
    setVelocity = vectorField(purePropellantSize, vector(0, 0, 0));
    {
        const volScalarField& propellant = phases[propellantIndex];
        label j = 0;
        forAll(propellant, i)
        {
            if (propellant[i] >= cutoff)
            {
                purePropellantCells[j] = i;
                j++;
            }
        }
    }
}