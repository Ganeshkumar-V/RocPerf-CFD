/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "Loth.H"
#include "mathematicalConstants.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace particleDragModels
{
    defineTypeNameAndDebug(Loth, 0);
    addToRunTimeSelectionTable(particleDragModel, Loth, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::particleDragModels::Loth::Loth
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject
)
:
    particleDragModel(dict, pair, registerObject),
    residualRe_("residualRe", dimless, dict),
    R_
    (
      "R", dimVelocity*dimVelocity/dimTemperature,
      dict.get<scalar>("R")
    ),
    gamma_(dict.get<scalar>("gamma"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::particleDragModels::Loth::~Loth()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::tmp<Foam::volScalarField> Foam::particleDragModels::Loth::erf
(const volScalarField& gf) const
{
  volScalarField erfF
  (
    IOobject
    (
      "erf(" + gf.name() + ")",
      gf.instance(),
      gf.db(),
      IOobject::NO_READ,
      IOobject::NO_WRITE
    ),
    gf.mesh(),
    gf.dimensions()
  );

  Internal::Field<scalar>& pF(erfF.primitiveFieldRef());
  // volScalarField& pF(terf.ref());
  forAll(pF, i)
  {
    pF[i] = std::erf(gf[i]);
  }

  Boundary& bFs(erfF.boundaryFieldRef());
  const Boundary& gfbFs = gf.boundaryField();
  forAll(bFs, i)
  {
    Field<scalar> bFpF = bFs[i];
    const Field<scalar> gfbFpF = gfbFs[i];

    forAll(bFpF, j)
    {
      bFpF[j] = std::erf(gfbFpF[j]);
    }
  }

  return Foam::tmp<Foam::volScalarField>(new volScalarField ("terf", erfF));
}

Foam::tmp<Foam::volScalarField> Foam::particleDragModels::Loth::tanh
(const volScalarField& gf) const
{
  volScalarField tanhF
  (
    IOobject
    (
        "tanh(" + gf.name() + ")",
        gf.instance(),
        gf.db(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    gf.mesh(),
    dimensionedScalar("", dimless, 0)
  );

  Internal::Field<scalar>& pF = tanhF.primitiveFieldRef();
  forAll(pF, i)
  {
    pF[i] = std::tanh(gf[i]);
  }
  // tanhF.correctBoundaryConditions();
  Boundary& bFs = tanhF.boundaryFieldRef();
  const Boundary& gfbFs = gf.boundaryField();
  forAll(bFs, i)
  {
    Field<scalar> bFpF = bFs[i];
    const Field<scalar> gfbFpF = gfbFs[i];
    forAll(bFpF, j)
    {
      bFpF[j] = std::tanh(gfbFpF[j]);
    }
  }

  return Foam::tmp<Foam::volScalarField>(new volScalarField ("ttanhF", tanhF));
}
// *************************** Rarefied Flow ****************************//
Foam::tmp<Foam::volScalarField> Foam::particleDragModels::Loth::fKn
(const Foam::volScalarField& Kn) const
{
  return 1/(1 + Kn*(2.514 + 0.8*exp(-0.55/max(Kn, SMALL))));
}

Foam::tmp<Foam::volScalarField> Foam::particleDragModels::Loth::JM
(const volScalarField& Ma) const
{
  volScalarField MaMax(max(Ma, SMALL));
  return
      neg0(Ma - 1)*(2.26 - 0.1/MaMax + 0.14/pow(MaMax, 3))
      + pos(Ma - 1)*(1.6 + 0.25/MaMax + 0.11/sqr(MaMax) + 0.44/pow(MaMax, 3));
}

Foam::tmp<Foam::volScalarField> Foam::particleDragModels::Loth::Cdfm
(const volScalarField& S) const
{
  volScalarField Smax(max(S, SMALL));
  volScalarField sqrS(sqr(Smax));
  volScalarField S4(pow(Smax, 4));

  return
    (1 + 2*sqrS)*exp(-sqrS)/max((sqrt(constant::mathematical::pi)*pow(Smax, 3)), SMALL)
    + 2*sqrt(constant::mathematical::pi)/max((3*Smax), SMALL)
    + (4*S4 + 4*sqrS - 1)*erf(Smax)/max((2*S4), SMALL);
}

Foam::tmp<Foam::volScalarField> Foam::particleDragModels::Loth::CdfmRe
(const volScalarField& Re, const volScalarField& Ma, const volScalarField& S) const
{
  volScalarField CdfmF(Cdfm(S));
  volScalarField JMF(JM(Ma));

  return CdfmF/(1 + sqrt(Re/45)*(CdfmF/JMF - 1));
}

Foam::tmp<Foam::volScalarField> Foam::particleDragModels::Loth::CdKnRe
(const volScalarField& Re, const volScalarField& Kn) const
{
  return 24*(1 + 0.15*pow(Re, 0.687))*fKn(Kn);
}

Foam::tmp<Foam::volScalarField> Foam::particleDragModels::Loth::CdRare
(
  const volScalarField& Re,
  const volScalarField& Ma,
  const volScalarField& Kn,
  const volScalarField& S
) const
{
  volScalarField Ma4(pow(Ma, 4));

  return (CdKnRe(Re, Kn) + Re*Ma4*CdfmRe(Re, Ma, S))/(1 + Ma4);
}

//* * * * ** * * ** * * * * * * * Compressibility Regime * * * * * * * * * //
Foam::tmp<Foam::volScalarField> Foam::particleDragModels::Loth::CM
(const volScalarField& Ma) const
{
  volScalarField MaMax(max(Ma, SMALL));
  return
      neg0(MaMax - 1.5)*(1.65 + 0.65*tanh((4*MaMax - 3.4)()))
      + pos(MaMax - 1.5)*(2.18 - 0.13*tanh((0.9*MaMax - 2.7)()));
}

Foam::tmp<Foam::volScalarField> Foam::particleDragModels::Loth::GM
(const volScalarField& Ma) const
{
  volScalarField MaMax(max(Ma, SMALL));
  volScalarField Ma3(pow(MaMax, 3));
  return
      neg0(MaMax - 0.8)*(166*Ma3 + 3.29*sqr(MaMax) - 10.9*MaMax + 20)
      + pos(MaMax - 0.8)*(5 + 40/Ma3);
}

Foam::tmp<Foam::volScalarField> Foam::particleDragModels::Loth::HM
(const volScalarField& Ma) const
{
  volScalarField MaMax(max(Ma, SMALL));
  return
      neg0(MaMax - 1.0)*(0.0239*pow(MaMax, 3) + 0.212*pow(MaMax, 2) - 0.074*MaMax + 1)
      + pos(MaMax - 1.0)*(0.93 + 1/(3.5 + pow(MaMax, 5)));
}

Foam::tmp<Foam::volScalarField> Foam::particleDragModels::Loth::CdComp
(const volScalarField& Re, const volScalarField& Ma) const
{
  volScalarField CMMa(CM(Ma));
  CMMa.correctBoundaryConditions();
  return
      24*(1 + 0.15*pow(Re, 0.687))*HM(Ma)
      + Re*0.42*CMMa/(1 + 42500/max(pow(Re, 1.16*CMMa), SMALL)
          + GM(Ma)/max(sqrt(Re), SMALL));
}

Foam::tmp<Foam::volScalarField> Foam::particleDragModels::Loth::CdRe() const
{
    volScalarField Ma(pair_.magUr()/sqrt(gamma_*R_*pair_.continuous().thermo().T()));
    volScalarField Re(pair_.Re());
    volScalarField S(sqrt(gamma_/2)*Ma);
    volScalarField Kn(sqrt(constant::mathematical::pi)*S/max(Re, SMALL));
    Info << "min -> Ma: " << min(Ma).value() << " Re: " << min(Re).value() << " S: " << min(S).value() << " Kn: " << min(Kn).value() << endl;
    Info << "max -> Ma: " << max(Ma).value() << " Re: " << max(Re).value() << " S: " << max(S).value() << " Kn: " << max(Kn).value() << endl;

    return
        neg0(Re - 45)*CdRare(Re, Ma, Kn, S)
        + pos(Re - 45)*CdComp(Re, Ma);
}


// ************************************************************************* //
