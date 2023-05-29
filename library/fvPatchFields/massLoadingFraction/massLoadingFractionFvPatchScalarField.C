/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "massLoadingFractionFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::massLoadingFractionFvPatchScalarField::massLoadingFractionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    Name_("name"),
    gas_("gasPhase"),
    particle_("particlePhase"),
    phi_(0.0)
{}


Foam::massLoadingFractionFvPatchScalarField::massLoadingFractionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict, false),
    Name_(dict.get<word>("name")),
    gas_(dict.get<word>("gasPhase")),
    particle_(dict.get<word>("particlePhase")),
    phi_(dict.get<scalar>("massFraction"))
{
    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        Info << "value is missing in massLoadingFraction bounday condition" << exit(FatalError);
    }
}


Foam::massLoadingFractionFvPatchScalarField::massLoadingFractionFvPatchScalarField
(
    const massLoadingFractionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    Name_(ptf.Name_),
    gas_(ptf.gas_),
    particle_(ptf.particle_),
    phi_(ptf.phi_)
{}


Foam::massLoadingFractionFvPatchScalarField::massLoadingFractionFvPatchScalarField
(
    const massLoadingFractionFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    Name_(tppsf.Name_),
    gas_(tppsf.gas_),
    particle_(tppsf.particle_),
    phi_(tppsf.phi_)
{}


Foam::massLoadingFractionFvPatchScalarField::massLoadingFractionFvPatchScalarField
(
    const massLoadingFractionFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    Name_(tppsf.Name_),
    gas_(tppsf.gas_),
    particle_(tppsf.particle_),
    phi_(tppsf.phi_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::massLoadingFractionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
}


void Foam::massLoadingFractionFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);
}


void Foam::massLoadingFractionFvPatchScalarField::updateCoeffs
()
{
    if (updated())
    {
        return;
    }

    // Gas phase density
    const fvPatchScalarField& rhog
    (
        patch().lookupPatchField<volScalarField, scalar>
        (
            IOobject::groupName("thermo:rho", gas_)
        )
    );

    // particle phase density
    const fvPatchScalarField& rhop
    (
        patch().lookupPatchField<volScalarField, scalar>
        (
            IOobject::groupName("thermo:rho", particle_)
        )
    );

    if (Name_ == gas_)
    {
        operator==(1/(1 + (phi_/(1 - phi_))*rhog/rhop));
    }
    else
    {
        operator==(1/(1 + ((1 - phi_)/max(phi_, 1e-15))*rhop/rhog));
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}

void Foam::massLoadingFractionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeEntry("name", Name_);
    os.writeEntry("gasPhase", gas_);
    os.writeEntry("particlePhase", particle_);
    os.writeEntry("massFraction", phi_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        massLoadingFractionFvPatchScalarField
    );
}

// ************************************************************************* //
