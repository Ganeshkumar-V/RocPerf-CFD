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

#include "DrakeInviscid.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sharpInterfaceHeatTransferModels
{
    defineTypeNameAndDebug(DrakeInviscid, 0);
    addToRunTimeSelectionTable(sharpInterfaceHeatTransferModel, DrakeInviscid, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sharpInterfaceHeatTransferModels::DrakeInviscid::DrakeInviscid
(
    const dictionary& dict,
    const phasePair& pair
)
:
    sharpInterfaceHeatTransferModel(dict, pair),
    cutoff(dict.get<scalar>("cutoff")),
    As_(dimensionedScalar("As", dimPressure*dimTime/sqrt(dimTemperature), dict.get<scalar>("As"))),
    Ts_(dimensionedScalar("Ts", dimTemperature, dict.get<scalar>("Ts"))),
    Pr_
    (
      volScalarField
      (
          IOobject
          (
              "Pr",
              pair.phase1().mesh().time().timeName(),
              pair.phase1().mesh(),
              IOobject::NO_READ,
              IOobject::NO_WRITE
          ),
          pair.phase1().mesh(),
          dimensionedScalar(dimless, dict.get<scalar>("Pr"))
      )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sharpInterfaceHeatTransferModels::DrakeInviscid::~DrakeInviscid()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::sharpInterfaceHeatTransferModels::DrakeInviscid::K(const scalar residualAlpha) const
{
    const tmp<volScalarField> tmagUr(pair_.magUr());
    const volScalarField& magUr(tmagUr());

    const tmp<volScalarField> trho(pair_.continuous().rho());
    const volScalarField& rho(trho());

    const tmp<volScalarField> td(pair_.dispersed().d());
    const volScalarField& d(td());

    const tmp<volScalarField> tTc(pair_.continuous().thermo().T());
    const volScalarField& Tc(tTc());

    volScalarField mu(As_*(sqrt(Tc)/(1.0 + Ts_/Tc)));

    volScalarField Re(rho*magUr*d/mu);

    const tmp<volScalarField> tCp(pair_.continuous().thermo().Cpv());
    const volScalarField& Cp(tCp());

    volScalarField kappa(mu*Cp/Pr_);

    volScalarField Nu(scalar(2) + 0.459*pow(Re, 0.55)*pow(Pr_, 0.33));
    return
        6.0
       *pair_.dispersed()*pos(pair_.dispersed() - cutoff)
       *kappa
       *Nu
       /sqr(pair_.dispersed().d());
}

// ************************************************************************* //
