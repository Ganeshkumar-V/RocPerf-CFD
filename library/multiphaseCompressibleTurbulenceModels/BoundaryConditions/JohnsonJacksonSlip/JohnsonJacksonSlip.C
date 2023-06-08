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

#include "JohnsonJacksonSlip.H"
#include "addToRunTimeSelectionTable.H"
#include "multiPhaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        JohnsonJacksonSlip
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::JohnsonJacksonSlip::
JohnsonJacksonSlip
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    partialSlipFvPatchVectorField(p, iF),
    specularityCoefficient_("specularityCoefficient", dimless, Zero),
    phasename_("particles")
{}


Foam::JohnsonJacksonSlip::
JohnsonJacksonSlip
(
    const JohnsonJacksonSlip& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    partialSlipFvPatchVectorField(ptf, p, iF, mapper),
    specularityCoefficient_(ptf.specularityCoefficient_),
    phasename_(ptf.phasename_)
{}


Foam::JohnsonJacksonSlip::
JohnsonJacksonSlip
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    partialSlipFvPatchVectorField(p, iF),
    specularityCoefficient_("specularityCoefficient", dimless, dict),
    phasename_(dict.get<word>("phase"))
{
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

    fvPatchVectorField::operator=
    (
        vectorField("value", dict, p.size())
    );
}


Foam::JohnsonJacksonSlip::
JohnsonJacksonSlip
(
    const JohnsonJacksonSlip& ptf
)
:
    partialSlipFvPatchVectorField(ptf),
    specularityCoefficient_(ptf.specularityCoefficient_),
    phasename_(ptf.phasename_)
{}


Foam::JohnsonJacksonSlip::
JohnsonJacksonSlip
(
    const JohnsonJacksonSlip& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    partialSlipFvPatchVectorField(ptf, iF),
    specularityCoefficient_(ptf.specularityCoefficient_),
    phasename_(ptf.phasename_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::JohnsonJacksonSlip::autoMap
(
    const fvPatchFieldMapper& m
)
{
    partialSlipFvPatchVectorField::autoMap(m);
}


void Foam::JohnsonJacksonSlip::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    partialSlipFvPatchVectorField::rmap(ptf, addr);
}


void Foam::JohnsonJacksonSlip::updateCoeffs()
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

    const fvPatchScalarField& gs0
    (
        patch().lookupPatchField<volScalarField, scalar>
        (
            IOobject::groupName("gs0", phasename_)
        )
    );

    const scalarField nu
    (
        patch().lookupPatchField<volScalarField, scalar>
        (
            IOobject::groupName("nut", phasename_)
        )
    );

    const scalarField nuFric
    (
        patch().lookupPatchField<volScalarField, scalar>
        (
            IOobject::groupName("nuFric", phasename_)
        )
    );

    word ThetaName(IOobject::groupName("Theta", phasename_));

    const fvPatchScalarField& Theta
    (
        db().foundObject<volScalarField>(ThetaName)
      ? patch().lookupPatchField<volScalarField, scalar>(ThetaName)
      : alpha
    );

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

    // calculate the slip value fraction
    scalarField c
    (
        constant::mathematical::pi
       *alpha
       *gs0
       *specularityCoefficient_.value()
       *sqrt(3*Theta)
       /max(6*(nu - nuFric)*alphaMax.value(), SMALL)
    );

    this->valueFraction() = c/(c + patch().deltaCoeffs());

    partialSlipFvPatchVectorField::updateCoeffs();
}


void Foam::JohnsonJacksonSlip::write
(
    Ostream& os
) const
{
    fvPatchVectorField::write(os);
    os.writeEntry("specularityCoefficient", specularityCoefficient_);
    os.writeEntry("phase", phasename_);
    writeEntry("value", os);
}


// ************************************************************************* //
