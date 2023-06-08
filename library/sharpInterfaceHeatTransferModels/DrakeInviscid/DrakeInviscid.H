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
    Foam::sharpInterfaceHeatTransferModels::DrakeInviscid

Description
    Ranz-Marshall correlation for turbulent heat transfer from the surface of a
    sphere to the surrounding fluid.

SourceFiles
    DrakeInviscid.C

\*---------------------------------------------------------------------------*/

#ifndef DrakeInviscid_H
#define DrakeInviscid_H

#include "sharpInterfaceHeatTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace sharpInterfaceHeatTransferModels
{

/*---------------------------------------------------------------------------*\
                          Class DrakeInviscid Declaration
\*---------------------------------------------------------------------------*/

class DrakeInviscid
:
    public sharpInterfaceHeatTransferModel
{
public:

    //- Runtime type information
    TypeName("DrakeInviscid");


    // Constructors

        //- Construct from components
        DrakeInviscid
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~DrakeInviscid();

    scalar cutoff;

    //- Thermal conductivity from sutherland
    dimensionedScalar As_;
    dimensionedScalar Ts_;
    volScalarField Pr_;
    // Member Functions

        //- The heat transfer function K used in the enthalpy equation
        tmp<volScalarField> K(const scalar residualAlpha) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace sharpInterfaceHeatTransferModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
