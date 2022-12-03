/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "interfaceTrackingModel.H"
#include "phasePair.H"
#include "processorFvPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(interfaceTrackingModel, 0);
    defineRunTimeSelectionTable(interfaceTrackingModel, dictionary);
}

const Foam::dimensionSet Foam::interfaceTrackingModel::dimDmdt =
    Foam::dimDensity/Foam::dimTime;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceTrackingModel::interfaceTrackingModel
(
    const dictionary& dict,
    const phasePair& pair
)
:
    pair_(pair),
    propellant_(dict.get<word>("propellant"))
{
  if (pair_.phase1().name() == propellant_ || pair_.phase2().name() == propellant_)
  {
    FatalErrorInFunction
     << "Propellant name should be different name other than the selected pair."
     << exit(FatalError);
  }
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
Foam::tmp<Foam::volScalarField>
Foam::interfaceTrackingModel::rhop() const
{
    return
    Foam::tmp<Foam::volScalarField>
    (
      new volScalarField("rhop", pair_.phase1().rho())
    );
}

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::interfaceTrackingModel>
Foam::interfaceTrackingModel::New
(
    const dictionary& dict,
    const phasePair& pair
)
{
    const word modelType(dict.get<word>("type"));

    Info<< "Selecting interfaceTrackingModel for "
        << pair << ": " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "interfaceTrackingModel",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << abort(FatalIOError);
    }

    return ctorPtr(dict, pair);
}

// ************************************************************************* //
