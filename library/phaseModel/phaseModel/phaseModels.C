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

#include "rhoThermo.H"
#include "rhoReactionThermo.H"

#include "CombustionModel.H"

#include "phaseModel.H"
#include "ThermoPhaseModel.H"
#include "AnisothermalParticlePhaseModel.H"
#include "ThermalEnergyPhaseModel.H"
#include "PurePhaseModel.H"
#include "MultiComponentPhaseModel.H"
#include "InertPhaseModel.H"
#include "MovingPhaseModel.H"
#include "StaticPhaseModel.H"
#include "NoEnergyChangePhaseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef
        AnisothermalParticlePhaseModel
        <
            PurePhaseModel
            <
                InertPhaseModel
                <
                    MovingPhaseModel
                    <
                        ThermoPhaseModel<phaseModel, rhoThermo>
                    >
                >
            >
        >
        pureParticlePhaseModel;

    addNamedToRunTimeSelectionTable
    (
        phaseModel,
        pureParticlePhaseModel,
        phaseSystem,
        pureParticlePhaseModel
    );

    typedef
        NoEnergyChangePhaseModel
        <
            PurePhaseModel
            <
                InertPhaseModel
                <
                    MovingPhaseModel
                    <
                        ThermoPhaseModel<phaseModel, rhoThermo>
                    >
                >
            >
        >
        pureIsothermalParticlePhaseModel;

    addNamedToRunTimeSelectionTable
    (
        phaseModel,
        pureIsothermalParticlePhaseModel,
        phaseSystem,
        pureIsothermalParticlePhaseModel
    );

    typedef
        ThermalEnergyPhaseModel
        <
            PurePhaseModel
            <
                InertPhaseModel
                <
                    MovingPhaseModel
                    <
                        ThermoPhaseModel<phaseModel, rhoThermo>
                    >
                >
            >
        >
        pureThermalEnergyPhaseModel;

    addNamedToRunTimeSelectionTable
    (
        phaseModel,
        pureThermalEnergyPhaseModel,
        phaseSystem,
        pureThermalEnergyPhaseModel
    );

    typedef
        ThermalEnergyPhaseModel
        <
            MultiComponentPhaseModel
            <
                InertPhaseModel
                <
                    MovingPhaseModel
                    <
                        ThermoPhaseModel<phaseModel, rhoReactionThermo>
                    >
                >
            >
        >
        multiComponentThermalEnergyPhaseModel;

    addNamedToRunTimeSelectionTable
    (
        phaseModel,
        multiComponentThermalEnergyPhaseModel,
        phaseSystem,
        multiComponentThermalEnergyPhaseModel
    );

    typedef
        NoEnergyChangePhaseModel
        <
            PurePhaseModel
            <
                InertPhaseModel
                <
                    StaticPhaseModel
                    <
                        ThermoPhaseModel<phaseModel, rhoThermo>
                    >
                >
            >
        >
        pureIsothermalStaticPhaseModel;

    addNamedToRunTimeSelectionTable
    (
        phaseModel,
        pureIsothermalStaticPhaseModel,
        phaseSystem,
        pureIsothermalStaticPhaseModel
    );
}

// ************************************************************************* //
