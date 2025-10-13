/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::burningRateInletVelocityFvPatchVectorField

Group
    grpInletBoundaryConditions

Description
    Velocity inlet boundary condition either correcting the extrapolated
    velocity or creating a uniform velocity field normal to the patch adjusted
    to match the specified flow rate

    For a mass-based flux:
    - the flow rate should be provided in kg/s
    - if \c rho is "none" the flow rate is in m3/s
    - otherwise \c rho should correspond to the name of the density field
    - if the density field cannot be found in the database, the user must
      specify the inlet density using the \c rhoInlet entry

Usage
    \table
        Property     | Description             | Required    | Default value
        massFlowRate | mass flow rate [kg/s]   | no          |
        rho          | density field name      | no          | rho
        rhoInlet     | inlet density           | no          |
    \endtable

    Example of the boundary condition specification for a mass flow rate:
    \verbatim
    <patchName>
    {
        type                burningRateInletVelocity;
        massFlowRate        0.2;
        rho                 rho;
        rhoInlet            1.0;
        value               uniform (0 0 0);
    }
    \endverbatim

    The \c flowRate entry is a \c Function1 of time, see Foam::Function1Types.

Note
    - \c rhoInlet is required for the case of a mass flow rate, where the
      density field is not available at start-up
    - The value is positive into the domain (as an inlet)
    - May not work correctly for transonic inlets
    - Strange behaviour with potentialFoam since the U equation is not solved

See also
    Foam::fixedValueFvPatchField
    Foam::Function1Types
    Foam::flowRateOutletVelocityFvPatchVectorField

SourceFiles
    burningRateInletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef burningRateInletVelocityFvPatchVectorField_H
#define burningRateInletVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class burningRateInletVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class burningRateInletVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private Data

        //- Burning rate pre-exponential factor
        scalar a_;

        //- Burning rate pressure exponent
        scalar n_;


    // Private Member Functions


public:

   //- Runtime type information
   TypeName("burningRateInletVelocity");


   // Constructors

        //- Construct from patch and internal field
        burningRateInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        burningRateInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  burningRateInletVelocityFvPatchVectorField
        //  onto a new patch
        burningRateInletVelocityFvPatchVectorField
        (
            const burningRateInletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        burningRateInletVelocityFvPatchVectorField
        (
            const burningRateInletVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new burningRateInletVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        burningRateInletVelocityFvPatchVectorField
        (
            const burningRateInletVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new burningRateInletVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member Functions

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
