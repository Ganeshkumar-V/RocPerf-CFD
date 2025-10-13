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


  forAll(psi.mesh().boundary(),patchi)
  {
    if (isType<processorFvPatch>(psi.mesh().boundary()[patchi]))
    {
      const processorPolyPatch& pp
          = refCast<const processorPolyPatch>(psi.mesh().boundaryMesh()[patchi]);
      if (pp.owner())
      {
        const scalarField& psibF(psi.boundaryField()[patchi].patchNeighbourField());
        const scalarField& alphaF(alpha.boundaryField()[patchi].patchInternalField());
        const scalarField& alphabF(alpha.boundaryField()[patchi].patchNeighbourField());
        const labelList& fC(psi.mesh().boundary()[patchi].faceCells());
        forAll(fC, i)
        {
          if ((alphabF[i] < One) && (alphaF[i] >= One))
          {
            psi[fC[i]] = psibF[i];
          }
        }
      }
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
      if(propellant[Own[i]] == One || propellant[Nei[i]] == One)
      {
         flux[i] = 0;
      }
      else continue;
    }

    forAll(mesh.boundary(),patchi)
    {
      if (isType<processorFvPatch>(mesh.boundary()[patchi]))
      {
        const processorPolyPatch& pp
            = refCast<const processorPolyPatch>(mesh.boundaryMesh()[patchi]);
        if (pp.neighbour())
        {
          const scalarField& alphaF(propellant.boundaryField()[patchi].patchInternalField());
          const scalarField& alphabF(propellant.boundaryField()[patchi].patchNeighbourField());
          scalarField& psi(flux.boundaryFieldRef()[patchi]);
          forAll(psi, i)
          {
            if ((alphabF[i] == One) || (alphaF[i] == One))
            {
              psi[i] = 0;
            }
          }
        }
      }
    }
    /***/
}

inline void correctU(volVectorField& U, const volScalarField& propellant)
{
  scalar One(1 - SMALL);

  forAll(propellant, i)
  {
    if(propellant[i] == One)
    {
       U[i] = vector(0, 0, 0);
    }
    else continue;
  }


  const fvMesh& mesh(propellant.mesh());
  forAll(mesh.boundary(),patchi)
  {
    if (isType<processorFvPatch>(mesh.boundary()[patchi]))
    {
      const processorPolyPatch& pp
          = refCast<const processorPolyPatch>(mesh.boundaryMesh()[patchi]);
      if (pp.neighbour())
      {
        const scalarField& alphaF(propellant.boundaryField()[patchi].patchInternalField());
        const scalarField& alphabF(propellant.boundaryField()[patchi].patchNeighbourField());
        vectorField& Uf(U.boundaryFieldRef()[patchi]);
        forAll(Uf, i)
        {
          if ((alphabF[i] == One) || (alphaF[i] == One))
          {
              Uf[i] = vector(0, 0, 0);
          }
        }
      }
    }
  }
  // U = neg(propellant - (1 - SMALL))*U;
}

void findYplus(const phaseModel& phase)
{
    const volVectorField& U(phase.U());
    const volVectorField::Boundary& UBf(U.boundaryField());

    const tmp<volScalarField> tnut(phase.nut());
    const volScalarField& nut(tnut());
    const volScalarField::Boundary& nutBf(nut.boundaryField());

    const tmp<volScalarField> tnuEff(phase.nuEff());
    const volScalarField& nuEff(tnuEff());
    const volScalarField::Boundary& nuEffBf(nuEff.boundaryField());

    const fvMesh& mesh(U.mesh());

    const fvPatchList& patches = mesh.boundary();
    const nearWallDist nwd(mesh);
    const volScalarField::Boundary& d = nwd.y();

    forAll(patches, i)
    {
        const fvPatch& patch = patches[i];
        if
        (
            isA<nutWallFunctionFvPatchScalarField>(nutBf[i])
        )
        {
            const nutWallFunctionFvPatchScalarField& nutPf =
                dynamic_cast<const nutWallFunctionFvPatchScalarField&>
                (
                    nutBf[i]
                );

            scalarField yPlus(nutPf.yPlus());
            Info << "yPlus in patch: " << patch.name() << " -> min/max: " << min(yPlus) << " - " << max(yPlus) << endl;
        }
        // else if (isA<wallFvPatch>(patch))
        if (isA<wallFvPatch>(patch))
        {
            Info  << "Not using nut wall functions! " << endl;
            scalarField yPlus(d[i]*sqrt((nuEffBf[i]*mag(UBf[i].snGrad()))())/(nuEffBf[i] - nutBf[i]));
            Info << "yPlus in patch: " << patch.name() << " -> min/max: " << min(yPlus) << " - " << max(yPlus) << endl;
            Info << "nu:  min/max: " << min((nuEffBf[i] - nutBf[i])()) << " - " << max((nuEffBf[i] - nutBf[i])()) << endl;
            Info << "nutBf:  min/max: " << min(nutBf[i]) << " - " << max(nutBf[i]) << endl;
            Info << "d:  min/max: " << min(d[i]) << " - " << max(d[i]) << endl;
            Info << "du/dn:  min/max: " << min(mag(UBf[i].snGrad())) << " - " << max(mag(UBf[i].snGrad())) << endl;
        }
    }
    Info << exit(FatalError);
}

inline void correctUY(volVectorField& U)
{
  forAll(U, i)
  {
    U[i].y() = 0;
  }
}

inline void limitT(volScalarField& he, const scalarField heMax, const scalarField heMin)
{
  forAll(he, i)
  {
    if (he[i] > heMax[0])
    {
      he[i] = heMax[0];
    }
    if (he[i] < heMin[0])
    {
      he[i] = heMin[0];
    }
  }
  // Info << "Limit Temperature Internal cells" << endl;
}
