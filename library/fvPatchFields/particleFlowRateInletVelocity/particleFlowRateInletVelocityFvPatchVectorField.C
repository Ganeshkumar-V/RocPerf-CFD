/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "particleFlowRateInletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "one.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::particleFlowRateInletVelocityFvPatchVectorField::
particleFlowRateInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    flowRate_(),
    rhoName_("rho"),
    rhoInlet_(0.0),
    volumetric_(false),
    phi_(0.0),
    continuous_("gas"),
    extrapolateProfile_(false)
{}


Foam::particleFlowRateInletVelocityFvPatchVectorField::
particleFlowRateInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict, false),
    rhoName_("rho"),
    rhoInlet_(dict.getOrDefault<scalar>("rhoInlet", -VGREAT)),
    volumetric_(false),
    phi_(dict.get<scalar>("massLoading")),
    continuous_(dict.get<word>("phase")),
    extrapolateProfile_
    (
        dict.getOrDefault<Switch>("extrapolateProfile", false)
    )
{
    if (dict.found("volumetricFlowRate"))
    {
        volumetric_ = true;
        flowRate_ =
            Function1<scalar>::New("volumetricFlowRate", dict, &db());
    }
    else if (dict.found("massFlowRate"))
    {
        volumetric_ = false;
        flowRate_ = Function1<scalar>::New("massFlowRate", dict, &db());
        rhoName_ = dict.getOrDefault<word>("rho", "rho");
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "Please supply either 'volumetricFlowRate' or"
            << " 'massFlowRate' and 'rho'" << nl
            << exit(FatalIOError);
    }

    // Value field require if mass based
    if (dict.found("value"))
    {
        fvPatchField<vector>::operator=
        (
            vectorField("value", dict, p.size())
        );
    }
    else
    {
        evaluate(Pstream::commsTypes::blocking);
    }
}


Foam::particleFlowRateInletVelocityFvPatchVectorField::
particleFlowRateInletVelocityFvPatchVectorField
(
    const particleFlowRateInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    flowRate_(ptf.flowRate_.clone()),
    rhoName_(ptf.rhoName_),
    rhoInlet_(ptf.rhoInlet_),
    volumetric_(ptf.volumetric_),
    phi_(ptf.phi_),
    continuous_(ptf.continuous_),
    extrapolateProfile_(ptf.extrapolateProfile_)
{}


Foam::particleFlowRateInletVelocityFvPatchVectorField::
particleFlowRateInletVelocityFvPatchVectorField
(
    const particleFlowRateInletVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    flowRate_(ptf.flowRate_.clone()),
    rhoName_(ptf.rhoName_),
    rhoInlet_(ptf.rhoInlet_),
    volumetric_(ptf.volumetric_),
    phi_(ptf.phi_),
    continuous_(ptf.continuous_),
    extrapolateProfile_(ptf.extrapolateProfile_)
{}


Foam::particleFlowRateInletVelocityFvPatchVectorField::
particleFlowRateInletVelocityFvPatchVectorField
(
    const particleFlowRateInletVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    flowRate_(ptf.flowRate_.clone()),
    rhoName_(ptf.rhoName_),
    rhoInlet_(ptf.rhoInlet_),
    volumetric_(ptf.volumetric_),
    phi_(ptf.phi_),
    continuous_(ptf.continuous_),
    extrapolateProfile_(ptf.extrapolateProfile_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class RhoType>
void Foam::particleFlowRateInletVelocityFvPatchVectorField::updateValues
(
    const RhoType& rho
)
{
    // const scalar t = db().time().timeOutputValue();
    //
    // const vectorField n(patch().nf());

    // lookup all the fields on this patch
    const fvPatchScalarField& alpha
    (
        patch().lookupPatchField<volScalarField, scalar>
        (
            IOobject::groupName("alpha", continuous_)
        )
    );

    const fvPatchScalarField& rhog
    (
        patch().lookupPatchField<volScalarField, scalar>
        (
            IOobject::groupName("rho", continuous_)
        )
    );

    const fvPatchVectorField& ug
    (
        patch().lookupPatchField<volVectorField, vector>
        (
            IOobject::groupName("U", continuous_)
        )
    );

    this->operator==((phi_/(1 - phi_))*(alpha/(1 - alpha))*rhog*ug/rho);
}


void Foam::particleFlowRateInletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (volumetric_ || rhoName_ == "none")
    {
        updateValues(one{});
    }
    else
    {
        // Mass flow-rate
        if (db().foundObject<volScalarField>(rhoName_))
        {
            const fvPatchField<scalar>& rhop =
                patch().lookupPatchField<volScalarField, scalar>(rhoName_);

            updateValues(rhop);
        }
        else
        {
            // Use constant density
            if (rhoInlet_ < 0)
            {
                FatalErrorInFunction
                    << "Did not find registered density field " << rhoName_
                    << " and no constant density 'rhoInlet' specified"
                    << exit(FatalError);
            }

            updateValues(rhoInlet_);
        }
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::particleFlowRateInletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    flowRate_->writeData(os);
    if (!volumetric_)
    {
        os.writeEntryIfDifferent<word>("rho", "rho", rhoName_);
        os.writeEntryIfDifferent<scalar>("rhoInlet", -VGREAT, rhoInlet_);
    }
    if (extrapolateProfile_)
    {
        os.writeEntry("extrapolateProfile", extrapolateProfile_);
    }
    writeEntry("value", os);
    os.writeEntry("massLoading", phi_);
    os.writeEntry("phase", continuous_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       particleFlowRateInletVelocityFvPatchVectorField
   );
}


// ************************************************************************* //
