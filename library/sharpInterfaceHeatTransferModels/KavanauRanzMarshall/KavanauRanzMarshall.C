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
}

// ************************************************************************* //
