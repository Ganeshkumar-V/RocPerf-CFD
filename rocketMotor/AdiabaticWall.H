inline void ImposeWall(volScalarField& psi, const volScalarField& alpha)
{
  const labelList& Own(psi.mesh().owner());
  const labelList& Nei(psi.mesh().neighbour());

  scalar One(0.999); //1.0 - SMALL);

  forAll(Own, celli)
  {
    if ((alpha[Own[celli]] < One) && (alpha[Nei[celli]] >= One))
    {
      psi[Nei[celli]] = psi[Own[celli]];
    }
    else if ((alpha[Nei[celli]] < One) && (alpha[Own[celli]] >= One))
    {
      psi[Own[celli]] = psi[Nei[celli]];
    }
  }

}

inline void correctFlux(surfaceScalarField& flux, const volScalarField& propellant)
{
  /***/
    const fvMesh& mesh(propellant.mesh());
    const labelList& Own(mesh.owner());
    const labelList& Nei(mesh.neighbour());
    scalar One(1 - SMALL);

    forAll(Own, i)
    {
      if(propellant[Own[i]] == One && propellant[Nei[i]] == One)
      {
         flux[i] = 0;
      }
      else continue;
    }
    /***/
}

inline void correctU(volVectorField& U, const volScalarField& propellant)
{
  U = neg(propellant - (1 - SMALL))*U;
}
