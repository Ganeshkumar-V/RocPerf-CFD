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

\*---------------------------------------------------------------------------*/

#include "addToRunTimeSelectionTable.H"
#include "phaseSystem.H"
#include "multiPhaseSystem.H"
#include "multiphaseSystem.H"
#include "MomentumEnergyTransferPhaseSystem.H"
#include "InterphaseHeatTransferPhaseSystem.H"
#include "PropellantRegressionPhaseSystem.H"
#include "PropellantInterfacePhaseSystem.H"
#include "PropellantCombustionPhaseSystem.H"
#include "PropellantTransferPhaseSystem.H"
#include "EntrainedPropellantCombustionPhaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    typedef
        PropellantRegressionPhaseSystem
        <
            InterphaseHeatTransferPhaseSystem
            <
                MomentumEnergyTransferPhaseSystem<multiPhaseSystem>
            >
        >
        regressionMultiphaseSystem;

    addNamedToRunTimeSelectionTable
    (
        multiPhaseSystem,
        regressionMultiphaseSystem,
        dictionary,
        regressionMultiphaseSystem
    );

    typedef
        InterphaseHeatTransferPhaseSystem
        <
            MomentumEnergyTransferPhaseSystem<multiPhaseSystem>
        >
        gasParticleMultiphaseSystem;

    addNamedToRunTimeSelectionTable
    (
        multiPhaseSystem,
        gasParticleMultiphaseSystem,
        dictionary,
        gasParticleMultiphaseSystem
    );

    typedef
        PropellantInterfacePhaseSystem
        <
            InterphaseHeatTransferPhaseSystem
            <
                MomentumEnergyTransferPhaseSystem<multiPhaseSystem>
            >
        >
        interfaceRegressionMultiphaseSystem;

    addNamedToRunTimeSelectionTable
    (
        multiPhaseSystem,
        interfaceRegressionMultiphaseSystem,
        dictionary,
        interfaceRegressionMultiphaseSystem
    );

    typedef
        PropellantCombustionPhaseSystem
        <
            InterphaseHeatTransferPhaseSystem
            <
                MomentumEnergyTransferPhaseSystem<multiPhaseSystem>
            >
        >
        interfaceCombustionMultiphaseSystem;

    addNamedToRunTimeSelectionTable
    (
        multiPhaseSystem,
        interfaceCombustionMultiphaseSystem,
        dictionary,
        interfaceCombustionMultiphaseSystem
    );

    typedef
        PropellantTransferPhaseSystem
        <
            InterphaseHeatTransferPhaseSystem
            <
                MomentumEnergyTransferPhaseSystem<multiPhaseSystem>
            >
        >
        propellantTransferMultiphaseSystem;

    addNamedToRunTimeSelectionTable
    (
        multiPhaseSystem,
        propellantTransferMultiphaseSystem,
        dictionary,
        propellantTransferMultiphaseSystem
    );

    typedef
        EntrainedPropellantCombustionPhaseSystem
        <
            InterphaseHeatTransferPhaseSystem
            <
                MomentumEnergyTransferPhaseSystem<multiPhaseSystem>
            >
        >
        entrainedPropellantCombustionMultiphaseSystem;

    addNamedToRunTimeSelectionTable
    (
        multiPhaseSystem,
        entrainedPropellantCombustionMultiphaseSystem,
        dictionary,
        entrainedPropellantCombustionMultiphaseSystem
    );

}

// ************************************************************************* //
