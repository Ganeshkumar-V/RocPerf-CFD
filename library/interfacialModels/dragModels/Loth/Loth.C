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

// Foam::tmp<Foam::volScalarField> Foam::particleDragModels::Loth::erf
// (const volScalarField& gf) const
// {
//   volScalarField erfF
//   (
//     IOobject
//     (
//       "erf(" + gf.name() + ")",
//       gf.instance(),
//       gf.db(),
//       IOobject::NO_READ,
//       IOobject::NO_WRITE
//     ),
//     gf.mesh(),
//     gf.dimensions()
//   );
//
//   Internal::Field<scalar>& pF(erfF.primitiveFieldRef());
//   // volScalarField& pF(terf.ref());
//   forAll(pF, i)
//   {
//     pF[i] = std::erf(gf[i]);
//   }
//
//   Boundary& bFs(erfF.boundaryFieldRef());
//   const Boundary& gfbFs = gf.boundaryField();
//   forAll(bFs, i)
//   {
//     Field<scalar> bFpF = bFs[i];
//     const Field<scalar> gfbFpF = gfbFs[i];
//
//     forAll(bFpF, j)
//     {
//       bFpF[j] = std::erf(gfbFpF[j]);
//     }
//   }
//
//   return Foam::tmp<Foam::volScalarField>(new volScalarField ("terf", erfF));
// }
//
// Foam::tmp<Foam::volScalarField> Foam::particleDragModels::Loth::tanh
// (const volScalarField& gf) const
// {
//   volScalarField tanhF
//   (
//     IOobject
//     (
//         "tanh(" + gf.name() + ")",
//         gf.instance(),
//         gf.db(),
//         IOobject::NO_READ,
//         IOobject::NO_WRITE
//     ),
//     gf.mesh(),
//     dimensionedScalar("", dimless, 0)
//   );
//
//   Internal::Field<scalar>& pF = tanhF.primitiveFieldRef();
//   forAll(pF, i)
//   {
//     pF[i] = std::tanh(gf[i]);
//   }
//   // tanhF.correctBoundaryConditions();
//   Boundary& bFs = tanhF.boundaryFieldRef();
//   const Boundary& gfbFs = gf.boundaryField();
//   forAll(bFs, i)
//   {
//     Field<scalar> bFpF = bFs[i];
//     const Field<scalar> gfbFpF = gfbFs[i];
//     forAll(bFpF, j)
//     {
//       bFpF[j] = std::tanh(gfbFpF[j]);
//     }
//   }
//
//   return Foam::tmp<Foam::volScalarField>(new volScalarField ("ttanhF", tanhF));
// }

// *************************** Rarefied Flow ****************************//
Foam::scalar Foam::particleDragModels::Loth::fKn
(const Foam::scalar& Kn) const
{
  return 1/(1 + Kn*(2.514 + 0.8*exp(-0.55/max(Kn, SMALL))));
}

Foam::scalar Foam::particleDragModels::Loth::JM
(const scalar& Ma) const
{
  return
      neg0(Ma - 1)*(2.26 - 0.1/max(Ma, SMALL) + 0.14/max(pow(Ma, 3), SMALL))
      + pos(Ma - 1)*(1.6 + 0.25/max(Ma, SMALL) + 0.11/max(sqr(Ma), SMALL) + 0.44/max(pow(Ma, 3), SMALL));
}

Foam::scalar Foam::particleDragModels::Loth::Cdfm
(const scalar& S) const
{
  return
    (1 + 2*sqr(S))*exp(-sqr(S))/max((sqrt(constant::mathematical::pi)*pow(S, 3)), SMALL)
    + 2*sqrt(constant::mathematical::pi)/max((3*S), SMALL)
    + (4*pow(S, 4) + 4*sqr(S) - 1)*std::erf(S)/max((2*pow(S, 4)), SMALL);
}

Foam::scalar Foam::particleDragModels::Loth::CdfmRe
(const scalar& Re, const scalar& Ma, const scalar& S) const
{
  return Cdfm(S)/(1 + sqrt(Re/45)*(Cdfm(S)/JM(Ma) - 1));
}

Foam::scalar Foam::particleDragModels::Loth::CdKnRe
(const scalar& Re, const scalar& Kn) const
{
  return 24*(1 + 0.15*pow(Re, 0.687))*fKn(Kn);
}

Foam::scalar Foam::particleDragModels::Loth::CdRare
(
  const scalar& Re,
  const scalar& Ma,
  const scalar& Kn,
  const scalar& S
) const
{
  return (CdKnRe(Re, Kn) + Re*pow(Ma, 4)*CdfmRe(Re, Ma, S))/(1 + pow(Ma, 4));
}

//* * * * ** * * ** * * * * * * * Compressibility Regime * * * * * * * * //
Foam::scalar Foam::particleDragModels::Loth::CM
(const scalar& Ma) const
{
  return
      neg0(Ma - 1.5)*(1.65 + 0.65*std::tanh(4*Ma - 3.4))
      + pos(Ma - 1.5)*(2.18 - 0.13*std::tanh(0.9*Ma - 2.7));
}

Foam::scalar Foam::particleDragModels::Loth::GM
(const scalar& Ma) const
{
  return
      neg0(Ma - 0.8)*(166*pow(Ma, 3) + 3.29*sqr(Ma) - 10.9*Ma + 20)
      + pos(Ma - 0.8)*(5 + 40/max(pow(Ma, 3), SMALL));
}

Foam::scalar Foam::particleDragModels::Loth::HM
(const scalar& Ma) const
{
  return
      neg0(Ma - 1.0)*(0.0239*pow(Ma, 3) + 0.212*pow(Ma, 2) - 0.074*Ma + 1)
      + pos(Ma - 1.0)*(0.93 + 1/(3.5 + pow(Ma, 5)));
}

Foam::scalar Foam::particleDragModels::Loth::CdComp
(const scalar& Re, const scalar& Ma) const
{
  return
      24*(1 + 0.15*pow(Re, 0.687))*HM(Ma)
      + Re*0.42*CM(Ma)/(1 + 42500/max(pow(Re, 1.16*CM(Ma)), SMALL)
          + GM(Ma)/max(sqrt(Re), SMALL));
}

Foam::tmp<Foam::volScalarField> Foam::particleDragModels::Loth::CdRe() const
{
    const tmp<volScalarField> tT(pair_.continuous().thermo().T());
    const volScalarField& T(tT());

    // volScalarField M(max(pair_.magUr()/sqrt(gamma_*R_*T), SMALL));
    // volScalarField Re(max(pair_.Re(), SMALL));
    // volScalarField S(sqrt(gamma_/2)*M);
    // volScalarField Kn(sqrt(constant::mathematical::pi)*S/Re);

    // Implementation - 2
    const tmp<volScalarField> tM(max(pair_.magUr()/sqrt(gamma_*R_*T), SMALL));
    const volScalarField& M(tM());

    const tmp<volScalarField> tRe(max(pair_.Re(), SMALL));
    const volScalarField& Re(tRe());

    const tmp<volScalarField> tS(sqrt(gamma_/2)*M);
    const volScalarField& S(tS());

    const tmp<volScalarField> tKn(sqrt(constant::mathematical::pi)*S/Re);
    const volScalarField& Kn(tKn());

    volScalarField CdRe_
    (
      IOobject
      (
        "CdRe",
        pair_.phase1().mesh()
      ),
      pair_.phase1().mesh(),
      dimensionedScalar("", dimless, 0.0)
    );

    forAll(Re, i)
    {
        CdRe_[i] =
        (
          Re[i] <= 45 ? CdRare(Re[i], M[i], Kn[i], S[i])
          : CdComp(Re[i], M[i])
        );
    }
    CdRe_.correctBoundaryConditions();

    return Foam::tmp<Foam::volScalarField>(new volScalarField ("tCdRe", CdRe_));
}

// ************************************************************************* //
