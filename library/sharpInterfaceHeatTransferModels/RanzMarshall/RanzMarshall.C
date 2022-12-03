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

#include "RanzMarshall.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sharpInterfaceHeatTransferModels
{
    defineTypeNameAndDebug(RanzMarshall, 0);
    addToRunTimeSelectionTable(sharpInterfaceHeatTransferModel, RanzMarshall, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sharpInterfaceHeatTransferModels::RanzMarshall::RanzMarshall
(
    const dictionary& dict,
    const phasePair& pair
)
:
    sharpInterfaceHeatTransferModel(dict, pair)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sharpInterfaceHeatTransferModels::RanzMarshall::~RanzMarshall()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::sharpInterfaceHeatTransferModels::RanzMarshall::K(const scalar residualAlpha) const
{
  // Info << "Reynolds Number: " << pair_.Re() << endl;
    // volScalarField Nu(scalar(2) + 0.6*sqrt(pair_.Re())*cbrt(pair_.Pr()));

    return
        6.0
       *pair_.dispersed()*pos(pair_.dispersed() - 1e-6)
       *pair_.continuous().kappa()
       *2.0
       /sqr(pair_.dispersed().d());
}


// ************************************************************************* //
