/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "subsonicNonReflectingNSCBCFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "EulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "backwardDdtScheme.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::subsonicNonReflectingNSCBCFvPatchField<Type>::subsonicNonReflectingNSCBCFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    advectiveFvPatchField<Type>(p, iF),
    psiName_("thermo:psi"),
    gamma_(0.0)
{}


template<class Type>
Foam::subsonicNonReflectingNSCBCFvPatchField<Type>::subsonicNonReflectingNSCBCFvPatchField
(
    const subsonicNonReflectingNSCBCFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    advectiveFvPatchField<Type>(ptf, p, iF, mapper),
    psiName_(ptf.psiName_),
    gamma_(ptf.gamma_)
{}


template<class Type>
Foam::subsonicNonReflectingNSCBCFvPatchField<Type>::subsonicNonReflectingNSCBCFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    advectiveFvPatchField<Type>(p, iF, dict),
    psiName_(dict.getOrDefault<word>("psi", "thermo:psi")),
    gamma_(dict.get<scalar>("gamma"))
{}


template<class Type>
Foam::subsonicNonReflectingNSCBCFvPatchField<Type>::subsonicNonReflectingNSCBCFvPatchField
(
    const subsonicNonReflectingNSCBCFvPatchField& ptpsf
)
:
    advectiveFvPatchField<Type>(ptpsf),
    psiName_(ptpsf.psiName_),
    gamma_(ptpsf.gamma_)
{}


template<class Type>
Foam::subsonicNonReflectingNSCBCFvPatchField<Type>::subsonicNonReflectingNSCBCFvPatchField
(
    const subsonicNonReflectingNSCBCFvPatchField& ptpsf,
    const DimensionedField<Type, volMesh>& iF
)
:
    advectiveFvPatchField<Type>(ptpsf, iF),
    psiName_(ptpsf.psiName_),
    gamma_(ptpsf.gamma_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::scalarField>
Foam::subsonicNonReflectingNSCBCFvPatchField<Type>::advectionSpeed() const
{
    // Lookup the velocity and compressibility of the patch
    const fvPatchField<scalar>& psip =
        this->patch().template
            lookupPatchField<volScalarField, scalar>(psiName_);

    const surfaceScalarField& phi =
        this->db().template lookupObject<surfaceScalarField>(this->phiName_);

    fvsPatchField<scalar> phip =
        this->patch().template
            lookupPatchField<surfaceScalarField, scalar>(this->phiName_);

    if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const fvPatchScalarField& rhop =
            this->patch().template
                lookupPatchField<volScalarField, scalar>(this->rhoName_);

        phip /= rhop;
    }

    // Calculate the speed of the field wave w
    // by summing the component of the velocity normal to the boundary
    // and the speed of sound (sqrt(gamma_/psi)).
    return phip/this->patch().magSf();
}

template<class Type>
Foam::tmp<Foam::scalarField>
Foam::subsonicNonReflectingNSCBCFvPatchField<Type>::soundSpeed() const
{
        const fvPatchField<scalar>& psip =
           this->patch().template
                lookupPatchField<volScalarField, scalar>(psiName_);

        return  sqrt(gamma_/psip);
}

template<class Type>
void Foam::subsonicNonReflectingNSCBCFvPatchField<Type>::updateCoeffs() const
{
    if (this->updated())
    {
        return;
    }

    const fvMesh& mesh = this->internalField().mesh();

    word ddtScheme
    (
        mesh.ddtScheme(this->internalField().name())
    );
    scalar deltaT = this->db().time().deltaTValue();

    const GeometricField<Type, fvPatchField, volMesh>& field =
        this->db().objectRegistry::template
        lookupObject<GeometricField<Type, fvPatchField, volMesh>>
        (
            this->internalField().name()
        );

    // Calculate the advection speed of the field wave
    // If the wave is incoming set the speed to 0.
    // const scalarField w(Foam::max(advectionSpeed(), scalar(0)));

    const  scalarField cP(soundSpeed());
    const  scalarField aP(advectionSpeed());

    const fvPatchScalarField& rhop =
            this->patch().template lookupPatchField<volScalarField, scalar>
            (
                rhoName_
            );

    const fvPatchVectorField& Up =
            this->patch().template lookupPatchField<volVectorField, vector>(UName_);



    label patchi = this->patch().index();

    // Non-reflecting outflow boundary
    // If lInf_ defined setup relaxation to the value fieldInf_.
    // if (lInf_ > 0)
    // {
    //
	  //   if
		//     (
		//      ddtScheme == fv::EulerDdtScheme<scalar>::typeName
		//      || ddtScheme == fv::CrankNicolsonDdtScheme<scalar>::typeName
		//     )
		//     {
		// 	    // Calculate the field relaxation coefficient k (See A2.2.2)
		// 	    const scalarField K(0.25*(1.0-sqr(aP/cP))*cP/lInf_);
    //
		// 	    this->valueFraction() = (1.0 + K*deltaT/2.0)/(1.0 + K*deltaT/2.0 + (aP+cP)/2.0*deltaT*this->patch().deltaCoeffs()) ;
		// 	    // ref-B.2.2
		// 	    this->refValue() =
		// 		    (
		// 		     field.oldTime().boundaryField()[patchi]
		// 		     + K * deltaT/2.0 * 101325 * fieldInf_
		// 		    )/( 1.0 + K * deltaT/2.0);
		// 	    // ref-B.2.3
		// 	    this->refGrad() = - 1.0 * rhop * cP * (this-> patch().nf() & Up.snGrad()) * fieldInf_;
		//     }
	  //   else if (ddtScheme == fv::backwardDdtScheme<scalar>::typeName)
	  //   {
		//     // Calculate the field relaxation coefficient k (See A2.2.2)
		//     const scalarField K(0.25*(1.0-sqr(aP/cP))*cP/lInf_);
    //
		//     this->valueFraction() = (1.5 + K*deltaT/2.0)/(1.5 + K*deltaT/2.0 + (aP+cP)/2.0*deltaT*this->patch().deltaCoeffs()) ;
		//     // ref-B.2.2
		//     this->refValue() =
		// 	    (
		// 	     2.0*field.oldTime().boundaryField()[patchi]
    //                        - 0.5*field.oldTime().oldTime().boundaryField()[patchi]
		// 	     + K * deltaT/2.0 * 101325 * fieldInf_
		// 	    )/( 1.5 + K * deltaT/2.0);
		//     // ref-B.2.3
		//     this->refGrad() = - 1.0 * rhop * cP * (this-> patch().nf() & Up.snGrad()) * fieldInf_;
	  //   }
    // }
    if (lInf_ > 0)
    {

	    if
		    (
		     ddtScheme == fv::EulerDdtScheme<scalar>::typeName
		     || ddtScheme == fv::CrankNicolsonDdtScheme<scalar>::typeName
		    )
		    {
			    // Calculate the field relaxation coefficient k (See A2.2.2)
			    const scalarField K(0.25*(1.0-sqr(aP/cP))*cP/lInf_);

			    this->valueFraction() = 1.0;
			    // ref-B.2.2
			    this->refValue() = field.oldTime().boundaryField()[patchi]*(K*deltaT + 1) - K*deltaT*fieldInf_;
				    // (
				    //  field.oldTime().boundaryField()[patchi]
				    //  + K * deltaT/2.0 * 101325 * fieldInf_
				    // )/( 1.0 + K * deltaT/2.0);
			    // ref-B.2.3
			    this->refGrad() = 0;
		    }
	    else if (ddtScheme == fv::backwardDdtScheme<scalar>::typeName)
	    {
		    // Calculate the field relaxation coefficient k (See A2.2.2)
		    const scalarField K(0.25*(1.0-sqr(aP/cP))*cP/lInf_);

		    this->valueFraction() = 1.0;
		    // ref-B.2.2
        this->refValue() = field.oldTime().boundaryField()[patchi]*(K*deltaT + 1) - K*deltaT*fieldInf_;

		    // ref-B.2.3
		    this->refGrad() = 0.0;
	    }
    }
    else
    {
        FatalErrorInFunction
                << "lInf_ must be above 0 "
                << exit(FatalError);
    }

    mixedFvPatchField<Type>::updateCoeffs();
}

template<class Type>
void Foam::subsonicNonReflectingNSCBCFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);

    os.writeEntryIfDifferent<word>("phi", "phi", this->phiName_);
    os.writeEntryIfDifferent<word>("rho", "rho", this->rhoName_);
    os.writeEntryIfDifferent<word>("psi", "thermo:psi", psiName_);

    os.writeEntry("gamma", gamma_);

    if (this->lInf_ > SMALL)
    {
        os.writeEntry("fieldInf", this->fieldInf_);
        os.writeEntry("lInf", this->lInf_);
    }

    this->writeEntry("value", os);
}


// ************************************************************************* //
