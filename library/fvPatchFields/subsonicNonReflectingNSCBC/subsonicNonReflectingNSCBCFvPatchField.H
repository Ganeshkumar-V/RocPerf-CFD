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
    Foam::subsonicNonReflectingNSCBCFvPatchField

Group
    grpOutletBoundaryConditions

Description
    This boundary condition provides a wave transmissive outflow condition,
    based on solving DDt(W, field) = 0 at the boundary \c W is the wave velocity
    and \c field is the field to which this boundary condition is applied.

    The wave speed is calculated using:

        \f[
            w_p = \frac{\phi_p}{|Sf|} + \sqrt{\frac{\gamma}{\psi_p}}
        \f]

    where

    \vartable
        w_p     | patch wave speed
        \phi_p  | patch face flux
        \psi_p  | patch compressibility
        Sf      | patch face area vector
        \gamma  | ratio of specific heats
    \endvartable

Usage
    \table
        Property     | Description             | Required    | Default value
        phi          | flux field name         | no          | phi
        rho          | density field name      | no          | rho
        psi          | compressibility field name | no       | thermo:psi
        gamma        | ratio of specific heats (Cp/Cv) | yes |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            subsonicNonReflectingNSCBC;
        phi             phi;
        psi             psi;
        gamma           1.4;
    }
    \endverbatim

See also
    Foam::advectiveFvPatchField

SourceFiles
    subsonicNonReflectingNSCBCFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef subsonicNonReflectingNSCBCFvPatchField_H
#define subsonicNonReflectingNSCBCFvPatchField_H

#include "advectiveFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class subsonicNonReflectingNSCBCFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class subsonicNonReflectingNSCBCFvPatchField
:
    public advectiveFvPatchField<Type>
{

    // Private data

        //- Name of the compressibility field used to calculate the wave speed
        word psiName_;

        //- Heat capacity ratio
        scalar gamma_;


public:

    //- Runtime type information
    TypeName("subsonicNonReflectingNSCBC");


    // Constructors

        //- Construct from patch and internal field
        subsonicNonReflectingNSCBCFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        subsonicNonReflectingNSCBCFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given subsonicNonReflectingNSCBCFvPatchField
        //  onto a new patch
        subsonicNonReflectingNSCBCFvPatchField
        (
            const subsonicNonReflectingNSCBCFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        subsonicNonReflectingNSCBCFvPatchField
        (
            const subsonicNonReflectingNSCBCFvPatchField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new subsonicNonReflectingNSCBCFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        subsonicNonReflectingNSCBCFvPatchField
        (
            const subsonicNonReflectingNSCBCFvPatchField&,
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
                new subsonicNonReflectingNSCBCFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the heat capacity ratio
            scalar gamma() const
            {
                return gamma_;
            }

            //- Return reference to the heat capacity ratio to allow adjustment
            scalar& gamma()
            {
                return gamma_;
            }


        // Evaluation functions

            //- Calculate and return the advection speed at the boundary
            virtual tmp<scalarField> advectionSpeed() const;

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "subsonicNonReflectingNSCBCFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
