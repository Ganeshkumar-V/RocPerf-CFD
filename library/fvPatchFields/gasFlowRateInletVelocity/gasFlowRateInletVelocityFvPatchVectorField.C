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

#include "gasFlowRateInletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "one.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gasFlowRateInletVelocityFvPatchVectorField::
gasFlowRateInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    flowRate_(),
    rhoName_("rho"),
    phaseName_("gas")
{}


Foam::gasFlowRateInletVelocityFvPatchVectorField::
gasFlowRateInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict, false),
    rhoName_("rho"),
    phaseName_("gas")
{
    if (dict.found("massFlowRate"))
    {
        flowRate_ = Function1<scalar>::New("massFlowRate", dict, &db());
        rhoName_ = dict.getOrDefault<word>("rho", "rho");
        phaseName_ = dict.getOrDefault<word>("phase", "gas");
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "Please supply 'massFlowRate' and 'rho'" << nl
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


Foam::gasFlowRateInletVelocityFvPatchVectorField::
gasFlowRateInletVelocityFvPatchVectorField
(
    const gasFlowRateInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    flowRate_(ptf.flowRate_.clone()),
    rhoName_(ptf.rhoName_),
    phaseName_(ptf.phaseName_)
{}


Foam::gasFlowRateInletVelocityFvPatchVectorField::
gasFlowRateInletVelocityFvPatchVectorField
(
    const gasFlowRateInletVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    flowRate_(ptf.flowRate_.clone()),
    rhoName_(ptf.rhoName_),
    phaseName_(ptf.phaseName_)
{}


Foam::gasFlowRateInletVelocityFvPatchVectorField::
gasFlowRateInletVelocityFvPatchVectorField
(
    const gasFlowRateInletVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    flowRate_(ptf.flowRate_.clone()),
    rhoName_(ptf.rhoName_),
    phaseName_(ptf.phaseName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class RhoType>
void Foam::gasFlowRateInletVelocityFvPatchVectorField::updateValues
(
    const RhoType& rho,
    const scalarField& alpha
)
{
    const scalar t = db().time().timeOutputValue();

    const vectorField n(patch().nf());

    const scalar avgU = -flowRate_->value(t)/gSum(alpha*rho*patch().magSf());
    operator==(avgU*n);

}


void Foam::gasFlowRateInletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Mass flow-rate
    if (db().foundObject<volScalarField>(rhoName_))
    {
        const fvPatchField<scalar>& rhop =
            patch().lookupPatchField<volScalarField, scalar>(rhoName_);

        const fvPatchScalarField& alpha
        (
            patch().lookupPatchField<volScalarField, scalar>
            (
                IOobject::groupName("alpha", phaseName_)
            )
        );

        updateValues(rhop, alpha);
    }
    else
    {
        FatalErrorInFunction
            << "Did not find registered density field " << rhoName_
            << exit(FatalError);
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::gasFlowRateInletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    flowRate_->writeData(os);
    os.writeEntryIfDifferent<word>("rho", "rho", rhoName_);
    os.writeEntryIfDifferent<word>("phase", "gas", phaseName_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       gasFlowRateInletVelocityFvPatchVectorField
   );
}


// ************************************************************************* //
