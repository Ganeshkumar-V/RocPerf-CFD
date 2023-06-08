/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2019 OpenFOAM Foundation
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

#include "JohnsonJacksonTheta.H"
#include "addToRunTimeSelectionTable.H"
#include "multiPhaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        JohnsonJacksonTheta
    );
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::JohnsonJacksonTheta::JohnsonJacksonTheta
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    restitutionCoefficient_("restitutionCoefficient", dimless, Zero),
    specularityCoefficient_("specularityCoefficient", dimless, Zero),
    phasename_("particles")
{}


Foam::JohnsonJacksonTheta::JohnsonJacksonTheta
(
    const JohnsonJacksonTheta& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    restitutionCoefficient_(ptf.restitutionCoefficient_),
    specularityCoefficient_(ptf.specularityCoefficient_),
    phasename_(ptf.phasename_)
{
}


Foam::JohnsonJacksonTheta::JohnsonJacksonTheta
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    restitutionCoefficient_("restitutionCoefficient", dimless, dict),
    specularityCoefficient_("specularityCoefficient", dimless, dict),
    phasename_(dict.get<word>("phase"))
{
    if
    (
        (restitutionCoefficient_.value() < 0)
     || (restitutionCoefficient_.value() > 1)
    )
    {
        FatalErrorInFunction
            << "The restitution coefficient has to be between 0 and 1"
            << abort(FatalError);
    }

    if
    (
        (specularityCoefficient_.value() < 0)
     || (specularityCoefficient_.value() > 1)
    )
    {
        FatalErrorInFunction
            << "The specularity coefficient has to be between 0 and 1"
            << abort(FatalError);
    }

    fvPatchScalarField::operator=
    (
        scalarField("value", dict, p.size())
    );
}


Foam::JohnsonJacksonTheta::JohnsonJacksonTheta
(
    const JohnsonJacksonTheta& ptf
)
:
    mixedFvPatchScalarField(ptf),
    restitutionCoefficient_(ptf.restitutionCoefficient_),
    specularityCoefficient_(ptf.specularityCoefficient_),
    phasename_(ptf.phasename_)
{}


Foam::JohnsonJacksonTheta::JohnsonJacksonTheta
(
    const JohnsonJacksonTheta& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    restitutionCoefficient_(ptf.restitutionCoefficient_),
    specularityCoefficient_(ptf.specularityCoefficient_),
    phasename_(ptf.phasename_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::JohnsonJacksonTheta::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);
}


void Foam::JohnsonJacksonTheta::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);
}


void Foam::JohnsonJacksonTheta::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // lookup all the fields on this patch
    const fvPatchScalarField& alpha
    (
        patch().lookupPatchField<volScalarField, scalar>
        (
            IOobject::groupName("alpha", phasename_)
        )
    );

    const fvPatchVectorField& U
    (
        patch().lookupPatchField<volVectorField, vector>
        (
            IOobject::groupName("U", phasename_)
        )
    );

    const fvPatchScalarField& gs0
    (
        patch().lookupPatchField<volScalarField, scalar>
        (
            IOobject::groupName("gs0", phasename_)
        )
    );

    const fvPatchScalarField& kappa
    (
        patch().lookupPatchField<volScalarField, scalar>
        (
            IOobject::groupName("kappa", phasename_)
        )
    );

    const scalarField Theta(patchInternalField());

    // lookup the packed volume fraction
    dimensionedScalar alphaMax
    (
        "alphaMax",
        dimless,
        db()
       .lookupObject<IOdictionary>
        (
            IOobject::groupName("turbulenceProperties", phasename_)
        )
       .subDict("RAS")
       .subDict("kineticTheoryCoeffs")
    );

    // calculate the reference value and the value fraction
    if (restitutionCoefficient_.value() != 1.0)
    {
        this->refValue() =
            (2.0/3.0)
           *specularityCoefficient_.value()
           *magSqr(U)
           /(scalar(1) - sqr(restitutionCoefficient_.value()));

        this->refGrad() = 0.0;

        scalarField c
        (
             constant::mathematical::pi*alpha*gs0
            *(scalar(1) - sqr(restitutionCoefficient_.value()))
            *sqrt(3*Theta)
            /max(4*kappa*alphaMax.value(), SMALL)
        );

        this->valueFraction() = c/(c + patch().deltaCoeffs());
    }

    // for a restitution coefficient of 1, the boundary degenerates to a fixed
    // gradient condition
    else
    {
        this->refValue() = 0.0;

        this->refGrad() =
            pos0(alpha - SMALL)
           *constant::mathematical::pi
           *specularityCoefficient_.value()
           *alpha
           *gs0
           *sqrt(3*Theta)
           *magSqr(U)
           /max(6*kappa*alphaMax.value(), SMALL);

        this->valueFraction() = 0;
    }

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::JohnsonJacksonTheta::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    os.writeEntry("restitutionCoefficient", restitutionCoefficient_);
    os.writeEntry("specularityCoefficient", specularityCoefficient_);
    os.writeEntry("phase", phasename_);
    writeEntry("value", os);
}


// ************************************************************************* //
