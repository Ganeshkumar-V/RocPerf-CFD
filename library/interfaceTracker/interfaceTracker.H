/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenFOAM Foundation
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
    Foam::interfaceTracker

Description

SourceFiles
    interfaceTracker.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceTracker_H
#define interfaceTracker_H

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "PtrListDictionary.H"
#include "phaseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;
/*---------------------------------------------------------------------------*\
                          Class interfaceTracker Declaration
\*---------------------------------------------------------------------------*/

class interfaceTracker
{
public:

    //- Runtime type information
    TypeName("interfaceTracker");

    //- Public Data
    volScalarField interface_;
    volScalarField Sf_;
    volScalarField Vf_;

    // Constructors

        //- Construct from components
        interfaceTracker(const volScalarField& alpha);


    //- Destructor
    virtual ~interfaceTracker();

    // Member Functions
    virtual void findInterface();

    virtual void interface();

    virtual void Sf();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
