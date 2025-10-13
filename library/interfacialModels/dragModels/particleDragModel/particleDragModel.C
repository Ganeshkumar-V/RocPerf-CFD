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

#include "particleDragModel.H"
#include "phasePair.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(particleDragModel, 0);
    defineRunTimeSelectionTable(particleDragModel, dictionary);
}

const Foam::dimensionSet Foam::particleDragModel::dimK(1, -3, -1, 0, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::particleDragModel::particleDragModel
(
    const phasePair& pair,
    const bool registerObject
)
:
    regIOobject
    (
        IOobject
        (
            IOobject::groupName(typeName, pair.name()),
            pair.phase1().mesh().time().timeName(),
            pair.phase1().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            registerObject
        )
    ),
    pair_(pair),
    factor_("", dimless, 1.0),
    Ur_
    (
      IOobject
      (
        IOobject::groupName("Ur", pair.name()),
        pair.phase1().mesh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
      ),
      pair.phase1().mesh(),
      dimensionedVector("Ur", dimVelocity, vector::zero)
    ),
    activeDrifting_(false)
{}


Foam::particleDragModel::particleDragModel
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject
)
:
    regIOobject
    (
        IOobject
        (
            IOobject::groupName(typeName, pair.name()),
            pair.phase1().mesh().time().timeName(),
            pair.phase1().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            registerObject
        )
    ),
    pair_(pair),
    factor_("", dimless, dict.getOrDefault<scalar>("factor", 1.0)),
    Ur_
    (
      IOobject
      (
        IOobject::groupName("Ur", pair.name()),
        pair.phase1().mesh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
      ),
      pair.Ur()
    ),
    activeDrifting_(dict.lookupOrDefault<bool>("drifting", false))
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::particleDragModel>
Foam::particleDragModel::New
(
    const dictionary& dict,
    const phasePair& pair
)
{
    const word modelType(dict.get<word>("type"));

    Info<< "Selecting particleDragModel for "
        << pair << ": " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "particleDragModel",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << abort(FatalIOError);
    }

    return ctorPtr(dict, pair, true);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::particleDragModel::~particleDragModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::particleDragModel::Ki() const
{
    const tmp<volScalarField> tmu(pair_.continuous().nu()*pair_.continuous().rho());
    volScalarField muc(tmu());

    if (isInviscid())
    {
      const tmp<volScalarField> tmc(mu());
      muc = tmc();
    }
    return
        0.75*factor_
       *CdRe()
       *muc
       /sqr(pair_.dispersed().d());
}


Foam::tmp<Foam::volScalarField> Foam::particleDragModel::K() const
{
    return max(pair_.dispersed(), pair_.dispersed().residualAlpha())*Ki();
}


Foam::tmp<Foam::surfaceScalarField> Foam::particleDragModel::Kf() const
{
    Info << "I am getting called!" << exit(FatalIOError);
    return
        max
        (
            fvc::interpolate(pair_.dispersed()),
            pair_.dispersed().residualAlpha()
        )*fvc::interpolate(Ki());
}

Foam::tmp<Foam::volScalarField> Foam::particleDragModel::mu() const
{
    return nullptr;
}

bool Foam::particleDragModel::isInviscid() const
{
    return false;
}

Foam::tmp<Foam::volVectorField> Foam::particleDragModel::Ur() const
{
    //- If drifting, return the relative velocity = Up - Ug - Ud
    //  - - - Appropriate turbulent particle drag model would replace this function
    
    //- If not drifting, return the default relative velocity
    return pair_.Ur();
}

bool Foam::particleDragModel::writeData(Ostream& os) const
{
    return os.good();
}


// ************************************************************************* //
