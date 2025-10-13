/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::massLoadingFractionFvPatchScalarField

Group
    grpInletBoundaryConditions grpOutletBoundaryConditions

Description
    This boundary condition provides a condition of 
    volume fraction of phase based on the specified mass fraction:
        \f[
            \alpha_1 = 1/(1 + ((1 - \phi)/\phi)*(\rho_1 u_1 \dot S_f)/(\rho_2 u_2 \dot S_f))
        \f]
        where
        \vartable
            \alpha_1     | volume fraction of the phase 1
            \phi         | mass fraction of phase 1
            \gamma  | ratio of specific heats (Cp/Cv)
            \psi    | compressibility [m2/s2]
            G       | coefficient given by \f$\frac{\gamma}{1-\gamma}\f$
        \endvartable

    The modes of operation are set by the dimensions of the pressure field
    to which this boundary condition is applied, the \c psi entry and the value
    of \c gamma:
    \table
        Mode                    | dimensions | psi   | gamma
        incompressible subsonic | p/rho      |       |
        compressible subsonic   | p          | none  |
        compressible transonic  | p          | psi   | 1
        compressible supersonic | p          | psi   | > 1
    \endtable


Usage
    \table
        Property     | Description                | Required | Default value
        U            | Velocity field name        | no       | U
        phi          | Flux field name            | no       | phi
        rho          | Density field name         | no       | rho
        psi          | Compressibility field name | no       | none
        gamma        | (Cp/Cv)                    | no       | 1
        p0           | Total pressure             | yes      |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            totalPressure;
        p0              uniform 1e5;
    }
    \endverbatim

See also
    Foam::fixedValueFvPatchField

SourceFiles
    massLoadingFractionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef sameAsFvPatchScalarField_H
#define sameAsFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class massLoadingFractionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class massLoadingFractionFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Name of the field
        word Name_;

        //- Gas phase
        word gas_;

        //- Particle phase
        word particle_;

        //- mass fraction
        scalar phi_;


public:

    //- Runtime type information
    TypeName("massLoadingFraction");


    // Constructors

        //- Construct from patch and internal field
        massLoadingFractionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        massLoadingFractionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given massLoadingFractionFvPatchScalarField
        //  onto a new patch
        massLoadingFractionFvPatchScalarField
        (
            const massLoadingFractionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        massLoadingFractionFvPatchScalarField
        (
            const massLoadingFractionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new massLoadingFractionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        massLoadingFractionFvPatchScalarField
        (
            const massLoadingFractionFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new massLoadingFractionFvPatchScalarField(*this, iF)
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
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        // Evaluation functions

            //- Inherit updateCoeffs from fixedValueFvPatchScalarField
            using fixedValueFvPatchScalarField::updateCoeffs;

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
