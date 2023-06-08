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
    Foam::particleDragModels::Loth

Description
    Schiller and Naumann drag model for dispersed bubbly flows.

SourceFiles
    Loth.C

\*---------------------------------------------------------------------------*/

#ifndef Loth_H
#define Loth_H

#include "particleDragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace particleDragModels
{

/*---------------------------------------------------------------------------*\
                           Class Loth Declaration
\*---------------------------------------------------------------------------*/

class Loth
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
    TypeName("Loth");


    // Constructors

        //- Construct from a dictionary and a phase pair
        Loth
        (
            const dictionary& dict,
            const phasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~Loth();

    dimensionedScalar R_;
    scalar gamma_;

    // Member Functions

    // Mathematical Functions
    tmp<volScalarField> erf(const volScalarField& gf) const;
    tmp<volScalarField> tanh(const volScalarField& gf) const;

    // Rarefied Flow Regime
    tmp<volScalarField> CdKnRe(const volScalarField& Re, const volScalarField& Kn) const;
    tmp<volScalarField> CdfmRe
    (
      const volScalarField& Re,
      const volScalarField& Ma,
      const volScalarField& S
    ) const;
    tmp<volScalarField> Cdfm(const volScalarField& S) const;
    tmp<volScalarField> JM(const volScalarField& Ma) const;
    tmp<volScalarField> fKn(const volScalarField& Kn) const;
    tmp<volScalarField> CdRare
    (
      const volScalarField& Re,
      const volScalarField& Ma,
      const volScalarField& Kn,
      const volScalarField& S
    ) const;

    // Compressibility Domainated Regime
    tmp<volScalarField> CdComp
    (
      const volScalarField& Re,
      const volScalarField& Ma
    ) const;
    tmp<volScalarField> CM(const volScalarField& Ma) const;
    tmp<volScalarField> GM(const volScalarField& Ma) const;
    tmp<volScalarField> HM(const volScalarField& Ma) const;


        //- Drag coefficient
        virtual tmp<volScalarField> CdRe() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace particleDragModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
