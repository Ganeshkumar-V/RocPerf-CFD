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

Class
    Foam::sameAsFvPatchField

Group
    grpGenericBoundaryConditions

Description
    This condition applies a scalar multiplier to the value of another
    boundary condition.

Usage
    \table
        Property     | Description             | Required    | Default value
        scale        | Time varying scale      | yes         |
        patch        | patchField providing the raw patch value | yes |
    \endtable

    Example of the boundary condition specification to scale a reference
    velocity of (15 0 0)  supplied as a fixedValue by a table of values
    that ramps the scale from 0 to 1 over 1 second:
    \verbatim
    <patchName>
    {
        type            sameAs;

        scale table
        (
            (    0   0)
            (  1.0 1.0)
            (100.0 1.0)
        );

        refValue
        {
            type            fixedValue;
            value           uniform (15 0 0);
        }
    }
    \endverbatim

SourceFiles
    sameAsFvPatchField.C

SeeAlso
    Foam::PatchFunction1

\*---------------------------------------------------------------------------*/

#ifndef sameAsFvPatchField_H
#define sameAsFvPatchField_H

#include "fixedValueFvPatchField.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class sameAsFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class sameAsFvPatchField
:
    public fixedValueFvPatchField<Type>
{
protected:

    // Protected data

    //- Name of the field
    word Name_;


public:

    //- Runtime type information
    TypeName("sameAs");


    // Constructors

        //- Construct from patch and internal field
        sameAsFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        sameAsFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given a sameAsFvPatchField onto
        //- a new patch
        sameAsFvPatchField
        (
            const sameAsFvPatchField&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        sameAsFvPatchField(const sameAsFvPatchField&);

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new sameAsFvPatchField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        sameAsFvPatchField
        (
            const sameAsFvPatchField&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new sameAsFvPatchField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the name of the velocity field
            const word& Name() const
            {
                return Name_;
            }

            //- Return reference to the name of the velocity field
            //  to allow adjustment
            word& Name()
            {
                return Name_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const fvPatchFieldMapper& m);

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<Type>& ptf,
                const labelList& addr
            );


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


        // Force an assignment irrespective of form of patch

            // virtual void operator==(const fvPatchField<Type>& ptf);
            // virtual void operator==(const Field<Type>& tf);
            // virtual void operator==(const Type& t);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "sameAsFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
