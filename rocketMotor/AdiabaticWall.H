void findWalls
(
  const volScalarField& alpha,
  labelList& Wall,
  labelList& WallPatch
)
{
  const labelList& Own(alpha.mesh().owner());
  const labelList& Nei(alpha.mesh().neighbour());

  scalar One(1.0 - SMALL);
  label cells = 0;
  Info << "---------------------------Finding Wall------------------------" << endl;
  forAll(Own, celli)
  {
    if
    (
      ((alpha[Own[celli]] < One) && (alpha[Nei[celli]] == One)) ||
      ((alpha[Nei[celli]] < One) && (alpha[Own[celli]] == One))
    )
    {
      cells++;
    }
  }

  Wall = labelList(cells, -1);
  WallPatch = labelList(cells, -1);
  label j = 0;
  forAll(Own, celli)
  {
    if ((alpha[Own[celli]] < One) && (alpha[Nei[celli]] == One))
    {
      Wall[j] = Nei[celli];
      WallPatch[j] = Own[celli];
      j++;
    }
    else if ((alpha[Nei[celli]] < One) && (alpha[Own[celli]] == One))
    {
      WallPatch[j] = Nei[celli];
      Wall[j] = Own[celli];
      j++;
    }
    else continue;
  }
}

// inline void ImposeWall
// (
//   volScalarField& psi,
//   const volScalarField& alpha,
//   labelList& Wall,
//   labelList& WallPatch
// )
// {
//   scalar One(1.0 - SMALL);
//
//   bool doAgain = false;
//   forAll(Wall, celli)
//   {
//     if (alpha[Wall[celli]] < One)
//     {
//       findWalls(alpha, Wall, WallPatch);
//       doAgain = true;
//       break;
//     }
//     else
//     {
//       psi[Wall[celli]] = psi[WallPatch[celli]];
//     }
//   }
//
//   // If changed Impose zeroGradient BC again
//   if (doAgain)
//   {
//     forAll(Wall, celli)
//     {
//       psi[Wall[celli]] = psi[WallPatch[celli]];
//     }
//   }
//
// }

inline void ImposeWall(volScalarField& psi, const volScalarField& alpha)
{
  const labelList& Own(psi.mesh().owner());
  const labelList& Nei(psi.mesh().neighbour());

  scalar One(1.0 - SMALL);

  forAll(Own, celli)
  {
    if ((alpha[Own[celli]] < One) && (alpha[Nei[celli]] == One))
    {
      psi[Nei[celli]] = psi[Own[celli]];
    }
    else if ((alpha[Nei[celli]] < One) && (alpha[Own[celli]] == One))
    {
      psi[Own[celli]] = psi[Nei[celli]];
    }
  }

}
