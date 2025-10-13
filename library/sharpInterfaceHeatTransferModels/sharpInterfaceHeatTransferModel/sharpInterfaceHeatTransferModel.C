/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "sharpInterfaceHeatTransferModel.H"
#include "phasePair.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sharpInterfaceHeatTransferModel, 0);
    defineRunTimeSelectionTable(sharpInterfaceHeatTransferModel, dictionary);
}

const Foam::dimensionSet Foam::sharpInterfaceHeatTransferModel::dimK(1, -1, -3, -1, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sharpInterfaceHeatTransferModel::sharpInterfaceHeatTransferModel
(
    const dictionary& dict,
    const phasePair& pair
)
:
    pair_(pair),
    residualAlpha_
    (
        "residualAlpha",
        dimless,
        dict.getOrDefault<scalar>
        (
            "residualAlpha",
            pair_.ordered()
          ? pair_.dispersed().residualAlpha().value()
          : pair_.phase1().residualAlpha().value()
        )
    ),
    kd_(pair_.phase1().thermo().kappa()),
    kc_(pair_.phase2().thermo().kappa())
{
  const fvMesh& mesh(pair_.phase1().mesh());

  dimensionedScalar tempd("Kd", kd_().dimensions(), dict.getOrDefault<scalar>("Kd", 0.0));
  dimensionedScalar tempc("Kc", kd_().dimensions(), dict.getOrDefault<scalar>("Kc", 0.0));

  kd_.ref() = volScalarField(IOobject("Kd", mesh), mesh, tempd);
  kc_.ref() = volScalarField(IOobject("Kc", mesh), mesh, tempc);
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::sharpInterfaceHeatTransferModel>
Foam::sharpInterfaceHeatTransferModel::New
(
    const dictionary& dict,
    const phasePair& pair
)
{
    const word modelType(dict.get<word>("type"));

    Info<< "Selecting sharpInterfaceHeatTransferModel for "
        << pair << ": " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "sharpInterfaceHeatTransferModel",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << abort(FatalIOError);
    }

    return ctorPtr(dict, pair);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::sharpInterfaceHeatTransferModel::K() const
{
    return K(residualAlpha_.value());
}

Foam::tmp<Foam::volScalarField>
Foam::sharpInterfaceHeatTransferModel::Nu() const
{
    return K(residualAlpha_.value());
}

const Foam::tmp<Foam::volScalarField>
Foam::sharpInterfaceHeatTransferModel::Kd() const
{
    return kd_;
}

const Foam::tmp<Foam::volScalarField>
Foam::sharpInterfaceHeatTransferModel::Kc() const
{
    return kc_;
}

// ************************************************************************* //
