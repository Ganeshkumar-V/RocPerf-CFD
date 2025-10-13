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

#include "KavanauRanzMarshall.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sharpInterfaceHeatTransferModels
{
    defineTypeNameAndDebug(KavanauRanzMarshall, 0);
    addToRunTimeSelectionTable(sharpInterfaceHeatTransferModel, KavanauRanzMarshall, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sharpInterfaceHeatTransferModels::KavanauRanzMarshall::KavanauRanzMarshall
(
    const dictionary& dict,
    const phasePair& pair
)
:
    sharpInterfaceHeatTransferModel(dict, pair),
    cutoff(dict.get<scalar>("cutoff")),
    gammaR
    (
      "gammaR", dimVelocity*dimVelocity/dimTemperature,
      dict.get<scalar>("gammaR")
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sharpInterfaceHeatTransferModels::KavanauRanzMarshall::~KavanauRanzMarshall()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::sharpInterfaceHeatTransferModels::KavanauRanzMarshall::Nu() const
{
    volScalarField Ma(pair_.magUr()/sqrt(gammaR*pair_.continuous().thermo().T()));
    volScalarField Re(pair_.Re());
    volScalarField Pr(pair_.Pr());
    volScalarField Nu0(scalar(2) + 0.6*sqrt(Re)*cbrt(Pr));

    return Nu0/(1 + 3.42*Nu0*Ma/max(Re*Pr, SMALL));
}

Foam::tmp<Foam::volScalarField>
Foam::sharpInterfaceHeatTransferModels::KavanauRanzMarshall::K(const scalar residualAlpha) const
{
    volScalarField Ma(pair_.magUr()/sqrt(gammaR*pair_.continuous().thermo().T()));
    volScalarField Re(pair_.Re());
    volScalarField Pr(pair_.Pr());
    volScalarField Nu0(scalar(2) + 0.6*sqrt(Re)*cbrt(Pr));
    volScalarField Nu(Nu0/(1 + 3.42*Nu0*Ma/max(Re*Pr, SMALL)));

    return
        6.0
       *pair_.dispersed()*pos(pair_.dispersed() - cutoff)
       *pair_.continuous().kappa()
       *Nu
       /sqr(pair_.dispersed().d());

    // Implementation - 2
    // volScalarField K
    // (
    //   IOobject("K.HeatTransfer", pair_.phase1().mesh()),
    //   pair_.phase1().mesh(),
    //   dimensionedScalar("", dimPower/(dimVolume*dimTemperature), 0.0)
    // );
    //
    // const tmp<volScalarField> tmagUr(pair_.magUr());
    // const volScalarField& magUr(tmagUr());
    //
    // const tmp<volScalarField> tmug(pair_.continuous().mu());
    // const volScalarField& mug(tmug());
    //
    // const tmp<volScalarField> trhog(pair_.continuous().rho());
    // const volScalarField& rhog(trhog());
    //
    // // const tmp<volScalarField> trhop(pair_.dispersed().rho());
    // // const volScalarField& rhop(trhop());
    //
    // const tmp<volScalarField> tkappac(pair_.continuous().kappa());
    // const volScalarField& kappac(tkappac());
    //
    // const tmp<volScalarField> tT(pair_.continuous().thermo().T());
    // const volScalarField& T(tT());
    //
    // const tmp<volScalarField> tCpg(pair_.continuous().thermo().Cpv());
    // const volScalarField& Cpg(tCpg());
    //
    // const tmp<volScalarField> tdp(pair_.dispersed().d());
    // const volScalarField& dp(tdp());
    //
    // const volScalarField& alphad(pair_.dispersed());
    // scalar Ma = 0.0;
    // scalar Re = 0.0;
    // scalar Pr = 0.0;
    // scalar Nu0 = 0.0;
    // scalar Nu = 0.0;
    //
    // forAll(K, i)
    // {
    //   Ma = magUr[i]/sqrt(gammaR.value()*T[i]);
    //   Re = max(rhog[i]*magUr[i]*dp[i]/mug[i], SMALL);
    //   Pr = mug[i]*Cpg[i]/kappac[i];
    //   Nu0 = 2.0 + 0.6*sqrt(Re)*cbrt(Pr);
    //   Nu = Nu0/(1 + 3.42*Nu0*Ma/(Re*Pr));
    //   K[i] = pos(alphad[i] - cutoff)*6.0*alphad[i]*kappac[i]*Nu/sqr(dp[i]);
    // }
    // K.correctBoundaryConditions();
    //
    // return Foam::tmp<Foam::volScalarField>(new volScalarField ("tK", K));
}

// ************************************************************************* //
