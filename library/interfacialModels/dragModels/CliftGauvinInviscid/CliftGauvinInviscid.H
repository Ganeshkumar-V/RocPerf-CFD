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
    Foam::particleDragModels::CliftGauvinInviscid

Description
    Schiller and Naumann drag model for dispersed bubbly flows.

SourceFiles
    CliftGauvinInviscid.C

\*---------------------------------------------------------------------------*/

#ifndef CliftGauvinInviscid_H
#define CliftGauvinInviscid_H

#include "particleDragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace particleDragModels
{

/*---------------------------------------------------------------------------*\
                           Class CliftGauvinInviscid Declaration
\*---------------------------------------------------------------------------*/

class CliftGauvinInviscid
:
    public particleDragModel
{
    // Private data

        //- Residual Reynolds Number
        const dimensionedScalar residualRe_;

        //- Viscosity coefficients from sutherland
        dimensionedScalar As_;
        dimensionedScalar Ts_;


public:

    //- Runtime type information
    TypeName("CliftGauvinInviscid");


    // Constructors

        //- Construct from a dictionary and a phase pair
        CliftGauvinInviscid
        (
            const dictionary& dict,
            const phasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~CliftGauvinInviscid();


    // Member Functions

        //- Drag coefficient
        virtual tmp<volScalarField> CdRe() const;

        //- Inviscid Functions
        virtual bool isInviscid() const;
        virtual tmp<volScalarField> mu() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace particleDragModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
