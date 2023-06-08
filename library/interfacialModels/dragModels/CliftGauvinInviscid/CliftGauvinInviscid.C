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

#include "CliftGauvinInviscid.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace particleDragModels
{
    defineTypeNameAndDebug(CliftGauvinInviscid, 0);
    addToRunTimeSelectionTable(particleDragModel, CliftGauvinInviscid, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::particleDragModels::CliftGauvinInviscid::CliftGauvinInviscid
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject
)
:
    particleDragModel(dict, pair, registerObject),
    residualRe_("residualRe", dimless, dict),
    As_(dimensionedScalar("As", dimPressure*dimTime/sqrt(dimTemperature), dict.get<scalar>("As"))),
    Ts_(dimensionedScalar("Ts", dimTemperature, dict.get<scalar>("Ts")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::particleDragModels::CliftGauvinInviscid::~CliftGauvinInviscid()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::particleDragModels::CliftGauvinInviscid::isInviscid() const
{
  return true;
}

Foam::tmp<Foam::volScalarField> Foam::particleDragModels::CliftGauvinInviscid::mu() const
{
    const tmp<volScalarField> tTc(pair_.continuous().thermo().T());
    const volScalarField& Tc(tTc());

    return As_*(sqrt(Tc)/(1.0 + Ts_/Tc));
}

Foam::tmp<Foam::volScalarField> Foam::particleDragModels::CliftGauvinInviscid::CdRe() const
{
    const tmp<volScalarField> tmagUr(pair_.magUr());
    const volScalarField& magUr(tmagUr());

    const tmp<volScalarField> trho(pair_.continuous().rho());
    const volScalarField& rho(trho());

    const tmp<volScalarField> td(pair_.dispersed().d());
    const volScalarField& d(td());

    volScalarField Re(rho*magUr*d/this->mu());

    return 24.0*(1.0 + 0.15*pow(Re, 0.687) + 0.0175*Re/(1.0 + 42500/max(pow(Re, 1.16), SMALL)));
}


// ************************************************************************* //
