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

#include "TurbulentLoth.H"
#include "mathematicalConstants.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace particleDragModels
{
    defineTypeNameAndDebug(TurbulentLoth, 0);
    addToRunTimeSelectionTable(particleDragModel, TurbulentLoth, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::particleDragModels::TurbulentLoth::TurbulentLoth
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

Foam::particleDragModels::TurbulentLoth::~TurbulentLoth()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> 
Foam::particleDragModels::TurbulentLoth::turbulentDispersion() const
{
    // Direction dependent parameter
    const volVectorField& Ur(pair_.Ur());
    const volScalarField& magUr(pair_.magUr());
    const volVectorField& Up(pair_.dispersed().U());
    const volScalarField& magUp(mag(Up));
    const volScalarField Cbeta(1.8 - 1.35*sqr((Ur&Up)/(magUr*magUp)));

    // Time scale of fluid turbulent motion
    const volScalarField Kpg
    (
        this->db().lookupObject<volScalarField>("Kpg." + pair_.dispersed().name())
    );

    const volScalarField Kg
    (
        this->db().lookupObject<volScalarField>("k." + pair_.continuous().name())
    );

    const volScalarField epsilonG
    (
        this->db().lookupObject<volScalarField>("epsilon." + pair_.continuous().name())
    );

    const volScalarField zetar(3.0*sqr(magUr)/(2.0*Kg));
    const volScalarField taupg((Kg/epsilonG)/sqrt(1.0 + Cbeta*zetar));

    return (1.0/3.0)*taupg*Kpg;
}

Foam::tmp<Foam::volVectorField> 
Foam::particleDragModels::TurbulentLoth::Ur() const
{
    //- Calculate Turbulent Dispersion Coefficient and drifting velocity
    const tmp<volScalarField> tDpg(turbulentDispersion());
    const volScalarField& Dpg(tDpg());

    const volScalarField& alphaP(pair_.dispersed());
    const volScalarField& alphaG(pair_.continuous());
    const volVectorField Ud(-Dpg*((1.0/alphaP)*fvc::grad(alphaP) - (1.0/alphaG)*fvc::grad(alphaG)));
    
    return pair_.Ur() - Ud;
}
// *************************** Rarefied Flow ****************************//
Foam::scalar Foam::particleDragModels::TurbulentLoth::fKn
(const Foam::scalar& Kn) const
{
  return 1/(1 + Kn*(2.514 + 0.8*exp(-0.55/max(Kn, SMALL))));
}

Foam::scalar Foam::particleDragModels::TurbulentLoth::JM
(const scalar& Ma) const
{
  return
      neg0(Ma - 1)*(2.26 - 0.1/max(Ma, SMALL) + 0.14/max(pow(Ma, 3), SMALL))
      + pos(Ma - 1)*(1.6 + 0.25/max(Ma, SMALL) + 0.11/max(sqr(Ma), SMALL) + 0.44/max(pow(Ma, 3), SMALL));
}

Foam::scalar Foam::particleDragModels::TurbulentLoth::Cdfm
(const scalar& S) const
{
  return
    (1 + 2*sqr(S))*exp(-sqr(S))/max((sqrt(constant::mathematical::pi)*pow(S, 3)), SMALL)
    + 2*sqrt(constant::mathematical::pi)/max((3*S), SMALL)
    + (4*pow(S, 4) + 4*sqr(S) - 1)*std::erf(S)/max((2*pow(S, 4)), SMALL);
}

Foam::scalar Foam::particleDragModels::TurbulentLoth::CdfmRe
(const scalar& Re, const scalar& Ma, const scalar& S) const
{
  return Cdfm(S)/(1 + sqrt(Re/45)*(Cdfm(S)/JM(Ma) - 1));
}

Foam::scalar Foam::particleDragModels::TurbulentLoth::CdKnRe
(const scalar& Re, const scalar& Kn) const
{
  return 24*(1 + 0.15*pow(Re, 0.687))*fKn(Kn);
}

Foam::scalar Foam::particleDragModels::TurbulentLoth::CdRare
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
Foam::scalar Foam::particleDragModels::TurbulentLoth::CM
(const scalar& Ma) const
{
  return
      neg0(Ma - 1.5)*(1.65 + 0.65*std::tanh(4*Ma - 3.4))
      + pos(Ma - 1.5)*(2.18 - 0.13*std::tanh(0.9*Ma - 2.7));
}

Foam::scalar Foam::particleDragModels::TurbulentLoth::GM
(const scalar& Ma) const
{
  return
      neg0(Ma - 0.8)*(166*pow(Ma, 3) + 3.29*sqr(Ma) - 10.9*Ma + 20)
      + pos(Ma - 0.8)*(5 + 40/max(pow(Ma, 3), SMALL));
}

Foam::scalar Foam::particleDragModels::TurbulentLoth::HM
(const scalar& Ma) const
{
  return
      neg0(Ma - 1.0)*(0.0239*pow(Ma, 3) + 0.212*pow(Ma, 2) - 0.074*Ma + 1)
      + pos(Ma - 1.0)*(0.93 + 1/(3.5 + pow(Ma, 5)));
}

Foam::scalar Foam::particleDragModels::TurbulentLoth::CdComp
(const scalar& Re, const scalar& Ma) const
{
  return
      24*(1 + 0.15*pow(Re, 0.687))*HM(Ma)
      + Re*0.42*CM(Ma)/(1 + 42500/max(pow(Re, 1.16*CM(Ma)), SMALL)
          + GM(Ma)/max(sqrt(Re), SMALL));
}

Foam::tmp<Foam::volScalarField> Foam::particleDragModels::TurbulentLoth::CdRe() const
{
    const tmp<volScalarField> tT(pair_.continuous().thermo().T());
    const volScalarField& T(tT());

    const tmp<volScalarField> tdp(pair_.dispersed().d());
    const volScalarField& dp(tdp());

    const tmp<volScalarField> tmu(pair_.continuous().nu()*pair_.continuous().rho());
    volScalarField mug(tmu());

    //- Calculate Turbulent Dispersion Coefficient and drifting velocity
    const tmp<volScalarField> tDpg(turbulentDispersion());
    const volScalarField& Dpg(tDpg());

    const volScalarField& alphaP(pair_.dispersed());
    const volScalarField& alphaG(pair_.continuous());
    const volVectorField Ud(-Dpg*((1.0/alphaP)*fvc::grad(alphaP) - (1.0/alphaG)*fvc::grad(alphaG)));
    
    const tmp<volVectorField> tUr(pair_.Ur() - Ud);

    const tmp<volScalarField> tMagUr(mag(tUr()));
    const volScalarField& magUr(tMagUr());

    //- Calculate Instantaneous relative velocity
    const volScalarField Kpg
    (
        this->db().lookupObject<volScalarField>("Kpg." + pair_.dispersed().name())
    );
    const volScalarField Kg
    (
        this->db().lookupObject<volScalarField>("k." + pair_.continuous().name())
    );
    const volScalarField Theta
    (
        this->db().lookupObject<volScalarField>("Theta." + pair_.dispersed().name())
    );

    const tmp<volScalarField> tmagUri(sqrt(sqr(magUr) + (2.0/3.0)*(Kg - Kpg) + Theta));
    const volScalarField& magUri(tmagUri());

    // Loth Drag Model
    const tmp<volScalarField> tM(max(magUri/sqrt(gamma_*R_*T), SMALL));
    const volScalarField& M(tM());

    const tmp<volScalarField> tRe(max(pair_.continuous().rho()*magUri*dp/mug, SMALL));
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

    // Update Boundary patches
    volScalarField::Boundary& CdReB(CdRe_.boundaryFieldRef());
    const volScalarField::Boundary& MB(M.boundaryField());
    const volScalarField::Boundary& ReB(Re.boundaryField());
    const volScalarField::Boundary& SB(S.boundaryField());
    const volScalarField::Boundary& KnB(Kn.boundaryField());

    forAll(CdReB, k)
    {
      Field<scalar>& pF(CdReB[k]);
      const Field<scalar>& MpF(MB[k]);
      const Field<scalar>& RepF(ReB[k]);
      const Field<scalar>& SpF(SB[k]);
      const Field<scalar>& KnpF(KnB[k]);

      forAll(pF, i)
      {
        pF[i] =
        (
          RepF[i] <= 45 ? CdRare(RepF[i], MpF[i], KnpF[i], SpF[i])
          : CdComp(RepF[i], MpF[i])
        );
      }
    }
    // CdRe_.correctBoundaryConditions();

    return Foam::tmp<Foam::volScalarField>(new volScalarField ("tCdRe", CdRe_));
}

// ************************************************************************* //
