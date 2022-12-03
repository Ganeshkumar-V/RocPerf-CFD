/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2018 OpenFOAM Foundation
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
    Foam::InterphaseHeatTransferPhaseSystem

Description
    Class which models interfacial heat transfer between a number of phases.
    Two heat transfer models are stored at each interface, one for each phase.
    This permits definition of an interface temperature with which heat transfer
    occurs. It also allows derived systems to define other thermodynamic
    properties at the interface and therefore represent phase changes.

See also
    OneResistanceHeatTransferPhaseSystem

SourceFiles
    InterphaseHeatTransferPhaseSystem.C

\*---------------------------------------------------------------------------*/

#ifndef InterphaseHeatTransferPhaseSystem_H
#define InterphaseHeatTransferPhaseSystem_H

#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class sharpInterfaceHeatTransferModel;

/*---------------------------------------------------------------------------*\
         Class InterphaseHeatTransferPhaseSystem Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseSystem>
class InterphaseHeatTransferPhaseSystem
:
    public BasePhaseSystem
{
protected:

    // Protected typedefs

        typedef HashTable
        <
            autoPtr<sharpInterfaceHeatTransferModel>,
            phasePairKey,
            phasePairKey::hash
        > sharpInterfaceHeatTransferModelTable;


    // Protected data

        // Sub Models

            //- Heat transfer models
            sharpInterfaceHeatTransferModelTable sharpInterfaceHeatTransferModels_;


public:

    // Constructors

        //- Construct from fvMesh
        InterphaseHeatTransferPhaseSystem(const fvMesh&);


    //- Destructor
    virtual ~InterphaseHeatTransferPhaseSystem();


    // Member Functions

        //- Return the heat transfer matrices
        virtual autoPtr<phaseSystem::heatTransferTable> heatTransfer() const;

        //- Correct the energy transport e.g. alphat and Tf
        virtual void correctEnergyTransport();

        //- Correct the interface thermodynamics
        virtual void correctInterfaceThermo();

        //- Read base phaseProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "InterphaseHeatTransferPhaseSystem.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
