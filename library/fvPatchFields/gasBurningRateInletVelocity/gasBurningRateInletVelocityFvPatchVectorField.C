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

#include "gasBurningRateInletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "one.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gasBurningRateInletVelocityFvPatchVectorField::
gasBurningRateInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    phaseName_("gas"),
    rhoPropellant_(1000),
    a_(1.0),
    n_(1.0),
    X_(0.5)
{}


Foam::gasBurningRateInletVelocityFvPatchVectorField::
gasBurningRateInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict, false),
    phaseName_("gas")
{
    phaseName_ = dict.getOrDefault<word>("phase", "gas");
    rhoPropellant_ = dict.get<scalar>("rhoPropellant");
    a_ = dict.get<scalar>("a");
    n_ = dict.get<scalar>("n");
    X_ = dict.get<scalar>("X");

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


Foam::gasBurningRateInletVelocityFvPatchVectorField::
gasBurningRateInletVelocityFvPatchVectorField
(
    const gasBurningRateInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    phaseName_(ptf.phaseName_),
    rhoPropellant_(ptf.rhoPropellant_),
    a_(ptf.a_),
    n_(ptf.n_),
    X_(ptf.X_)
{}


Foam::gasBurningRateInletVelocityFvPatchVectorField::
gasBurningRateInletVelocityFvPatchVectorField
(
    const gasBurningRateInletVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    phaseName_(ptf.phaseName_),
    rhoPropellant_(ptf.rhoPropellant_),
    a_(ptf.a_),
    n_(ptf.n_),
    X_(ptf.X_)
{}


Foam::gasBurningRateInletVelocityFvPatchVectorField::
gasBurningRateInletVelocityFvPatchVectorField
(
    const gasBurningRateInletVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    phaseName_(ptf.phaseName_),
    rhoPropellant_(ptf.rhoPropellant_),
    a_(ptf.a_),
    n_(ptf.n_),
    X_(ptf.X_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class RhoType>
void Foam::gasBurningRateInletVelocityFvPatchVectorField::updateValues
(
    const RhoType& rho,
    const scalarField& alpha,
    const scalarField& P
)
{
    const vectorField n(patch().nf());

    const scalarField rb(a_*pow(P/10e6, n_));

    const scalarField avgU(-(X_*rhoPropellant_*rb)/gSum(alpha*rho));
    operator==(avgU*n);

}


void Foam::gasBurningRateInletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchField<scalar>& rho
    (
        patch().lookupPatchField<volScalarField, scalar>
        (
            IOobject::groupName("thermo:rho", phaseName_)
        )
    );

    const fvPatchScalarField& alpha
    (
        patch().lookupPatchField<volScalarField, scalar>
        (
            IOobject::groupName("alpha", phaseName_)
        )
    );

    const fvPatchScalarField& P
    (
        patch().lookupPatchField<volScalarField, scalar>("p")
    );

    updateValues(rho, alpha, P);

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::gasBurningRateInletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    os.writeEntryIfDifferent<word>("phase", "gas", phaseName_);
    os.writeEntry<scalar>("rhoPropellant", rhoPropellant_);
    os.writeEntry<scalar>("a", a_);
    os.writeEntry<scalar>("n", n_);
    os.writeEntry<scalar>("X", X_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       gasBurningRateInletVelocityFvPatchVectorField
   );
}


// ************************************************************************* //
