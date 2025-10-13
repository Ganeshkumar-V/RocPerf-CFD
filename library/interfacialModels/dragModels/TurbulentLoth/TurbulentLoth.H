/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::particleDragModels::TurbulentLoth

Description
    Schiller and Naumann drag model for dispersed bubbly flows.

SourceFiles
    TurbulentLoth.C

\*---------------------------------------------------------------------------*/

#ifndef TurbulentLoth_H
#define TurbulentLoth_H

#include "particleDragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace particleDragModels
{

/*---------------------------------------------------------------------------*\
                           Class TurbulentLoth Declaration
\*---------------------------------------------------------------------------*/

class TurbulentLoth
:
    public particleDragModel
{
    // Private data

        //- Residual Reynolds Number
        const dimensionedScalar residualRe_;

    // typedefs
    typedef DimensionedField<scalar, fvMesh> Internal;
    typedef FieldField<fvPatchField, scalar> Boundary;


public:

    //- Runtime type information
    TypeName("TurbulentLoth");


    // Constructors

        //- Construct from a dictionary and a phase pair
        TurbulentLoth
        (
            const dictionary& dict,
            const phasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~TurbulentLoth();

    dimensionedScalar R_;
    scalar gamma_;

    // Member Functions

    // Mathematical Functions
    // tmp<volScalarField> erf(const volScalarField& gf) const;
    // tmp<volScalarField> tanh(const volScalarField& gf) const;

    // Rarefied Flow Regime
    inline scalar CdKnRe(const scalar& Re, const scalar& Kn) const;
    inline scalar CdfmRe
    (
      const scalar& Re,
      const scalar& Ma,
      const scalar& S
    ) const;
    inline scalar Cdfm(const scalar& S) const;
    inline scalar JM(const scalar& Ma) const;
    inline scalar fKn(const scalar& Kn) const;
    inline scalar CdRare
    (
      const scalar& Re,
      const scalar& Ma,
      const scalar& Kn,
      const scalar& S
    ) const;

    // Compressibility Domainated Regime
    inline scalar CdComp
    (
      const scalar& Re,
      const scalar& Ma
    ) const;
    inline scalar CM(const scalar& Ma) const;
    inline scalar GM(const scalar& Ma) const;
    inline scalar HM(const scalar& Ma) const;

    // Turbulent Dispersion Coefficient
    tmp<volScalarField> turbulentDispersion() const;

    virtual tmp<volVectorField> Ur() const;

        //- Drag coefficient
        virtual tmp<volScalarField> CdRe() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace particleDragModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
