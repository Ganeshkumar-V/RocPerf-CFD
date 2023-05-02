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
    Foam::PropellantCombustionPhaseSystem

Description
    Class which models non-thermally-coupled mass transfers; i.e.,
    representation changes, rather than phase changes.

SourceFiles
    PropellantCombustionPhaseSystem.C

\*---------------------------------------------------------------------------*/

#ifndef PropellantCombustionPhaseSystem_H
#define PropellantCombustionPhaseSystem_H

#include "phaseSystem.H"
#include "saturationModel.H"
#include "ChemicalEfficiency.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class interfaceTrackingModel;

/*---------------------------------------------------------------------------*\
                    Class PropellantCombustionPhaseSystem Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseSystem>
class PropellantCombustionPhaseSystem
:
    public BasePhaseSystem
{
protected:

    // Protected typedefs

        typedef HashTable
        <
            autoPtr<interfaceTrackingModel>,
            phasePairKey,
            phasePairKey::hash
        > interfaceTrackingModelTable;

        typedef HashPtrTable<volScalarField, phasePairKey, phasePairKey::hash>
            rDmdtTable;

        typedef HashPtrTable<volVectorField, phasePairKey, phasePairKey::hash>
            interfaceTable;

        typedef HashTable<scalar, phasePairKey, phasePairKey::hash>
            coeffTable;

    // Protected data

        dictionary phaseProperties;

        // Adiabatic Flame Temperature
        volScalarField Tad;

        // Burning rate
        volScalarField rb_;

        // propellant density
        dimensionedScalar rhoPropellant;

        // Regression Control Parameters
        bool regress_;
        word checkAndRegress_;
        scalar stopRegress_;
        scalar intialV_;
        dimensionedScalar retainer_;

        // Regresssion Parameters
        volScalarField alphaOld;
        volScalarField regressionAlpha;
        scalar epsilon;

        // Inlet velocity
        volVectorField Ug_;
        volVectorField Up_;

        // Sub Models

            //- Mass transfer models
            interfaceTrackingModelTable interfaceTrackingModels_;

            //- Mass transfer rates
            rDmdtTable rDmdt_;

            interfaceTable nHat_;

            //- The saturation model used to evaluate Tsat = Tf
            autoPtr<saturationModel> saturationModel_;

    // Protected member functions

        //- Return the representation mass transfer rate
        virtual tmp<volScalarField> rDmdt(const phasePairKey& key) const;


public:

    // Constructors

        //- Construct from fvMesh
        PropellantCombustionPhaseSystem(const fvMesh&);


    //- Destructor
    virtual ~PropellantCombustionPhaseSystem();

    typedef PtrListDictionary<phaseModel> phaseModelList;

    // Public Data
    chemicalEfficiency eta;

    // Member Functions

    //- Return the saturationModel
    const saturationModel& saturation() const;

        //- Return the mass transfer rate for a pair
        virtual tmp<volScalarField> dmdt(const phasePairKey& key) const;

        //- Return the mass transfer rates for each phase
        virtual PtrList<volScalarField> dmdts() const;

        //- Return the mass transfer matrices
        virtual autoPtr<phaseSystem::massTransferTable> massTransfer() const;

        //- Return the heat transfer matrices
        virtual autoPtr<phaseSystem::heatTransferTable> heatTransfer() const;

        //- Return the momentum transfer matrices
        virtual autoPtr<phaseSystem::momentumTransferTable> momentumTransfer();

        //- Regress Propellant Surface
        virtual void solve();

        //- Store old Times
        virtual void store();

        //- Correct the mass transfer rates
        virtual void correct();

        //- Calculate velocity
        void calculateVelocity();

        //- Read base phaseProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PropellantCombustionPhaseSystem.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
