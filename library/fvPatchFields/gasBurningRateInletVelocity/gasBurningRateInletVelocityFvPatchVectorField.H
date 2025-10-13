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
    Foam::gasBurningRateInletVelocityFvPatchVectorField

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
        type                gasBurningRateInletVelocity;
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
    gasBurningRateInletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef gasBurningRateInletVelocityFvPatchVectorField_H
#define gasBurningRateInletVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class gasBurningRateInletVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class gasBurningRateInletVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private Data

        //- Name of the phase
        word phaseName_;

        //- propellant density
        scalar rhoPropellant_;

        //- Burning rate pre-exponential factor
        scalar a_;

        //- Burning rate pressure exponent
        scalar n_;

        //- mass fraction
        scalar X_;


    // Private Member Functions

        //- Update the patch values given the appropriate density type and value
        template<class RhoType>
        void updateValues(const RhoType& rho, const scalarField& alpha, const scalarField& P);


public:

   //- Runtime type information
   TypeName("gasBurningRateInletVelocity");


   // Constructors

        //- Construct from patch and internal field
        gasBurningRateInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        gasBurningRateInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  gasBurningRateInletVelocityFvPatchVectorField
        //  onto a new patch
        gasBurningRateInletVelocityFvPatchVectorField
        (
            const gasBurningRateInletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        gasBurningRateInletVelocityFvPatchVectorField
        (
            const gasBurningRateInletVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new gasBurningRateInletVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        gasBurningRateInletVelocityFvPatchVectorField
        (
            const gasBurningRateInletVelocityFvPatchVectorField&,
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
                new gasBurningRateInletVelocityFvPatchVectorField(*this, iF)
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
