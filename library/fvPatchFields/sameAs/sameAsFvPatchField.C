/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "sameAsFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::sameAsFvPatchField<Type>::sameAsFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    Name_("field")
{}


template<class Type>
Foam::sameAsFvPatchField<Type>::sameAsFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict, false),
    Name_(dict.get<word>("field"))
{
    if (this->db().template foundObject<GeometricField<Type, fvPatchField, volMesh>>(Name_))
    {
        const fvPatchField<Type>& pF =
            this->patch().template lookupPatchField
            <
            GeometricField<Type, fvPatchField, volMesh>, Type
            >(Name_);

        fvPatchField<Type>::operator=(pF);
    }
    else
    {
        Warning << "Field " << Name_ << " does not exist. Skipping access." << endl;
    }
}


template<class Type>
Foam::sameAsFvPatchField<Type>::sameAsFvPatchField
(
    const sameAsFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    Name_(ptf.Name_)
{}


template<class Type>
Foam::sameAsFvPatchField<Type>::sameAsFvPatchField
(
    const sameAsFvPatchField& spf
)
:
    fixedValueFvPatchField<Type>(spf),
    Name_(spf.Name_)
{}


template<class Type>
Foam::sameAsFvPatchField<Type>::sameAsFvPatchField
(
    const sameAsFvPatchField& spf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(spf, iF),
    Name_(spf.Name_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::sameAsFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<Type>::autoMap(m);
}


template<class Type>
void Foam::sameAsFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<Type>::rmap(ptf, addr);
}


template<class Type>
void Foam::sameAsFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (this->db().template foundObject<GeometricField<Type, fvPatchField, volMesh>>(Name_))
    {
    const fvPatchField<Type>& pF =
        this->patch().template lookupPatchField
        <
        GeometricField<Type, fvPatchField, volMesh>, Type
        >(Name_);

        // Note: setting this field value using = operator (not ==)
        fvPatchField<Type>::operator=(pF);
    }
    else
    {
        Warning << "Field " << Name_ << " does not exist. Skipping access." << endl;
    }

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::sameAsFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    os.writeEntry("field", Name_);
    this->writeEntry("value", os);
}


// template<class Type>
// void Foam::sameAsFvPatchField<Type>::operator==
// (
//     const fvPatchField<Type>& ptf
// )
// {
//     const scalarField s(scalePtr_->value(this->db().time().timeOutputValue()));
//
//     forAll(s, facei)
//     {
//         const scalar si = s[facei];
//         if (mag(si) > ROOTVSMALL)
//         {
//             refValuePtr_->operator[](facei) = ptf[facei]/si;
//         }
//     }
//
//     Field<Type>::operator=(ptf);
// }
//
//
// template<class Type>
// void Foam::sameAsFvPatchField<Type>::operator==(const Field<Type>& tf)
// {
//     const scalarField s(scalePtr_->value(this->db().time().timeOutputValue()));
//
//     forAll(s, facei)
//     {
//         const scalar si = s[facei];
//         if (mag(si) > ROOTVSMALL)
//         {
//             refValuePtr_->operator[](facei) = tf[facei]/si;
//         }
//     }
//
//     Field<Type>::operator=(tf);
// }
//
//
// template<class Type>
// void Foam::sameAsFvPatchField<Type>::operator==(const Type& t)
// {
//     const scalarField s(scalePtr_->value(this->db().time().timeOutputValue()));
//
//     forAll(s, facei)
//     {
//         const scalar si = s[facei];
//         if (mag(si) > ROOTVSMALL)
//         {
//             refValuePtr_->operator[](facei) = t/si;
//         }
//     }
//
//     Field<Type>::operator=(t);
// }


// ************************************************************************* //
