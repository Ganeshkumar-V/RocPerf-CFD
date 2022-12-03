/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2019 OpenFOAM Foundation
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

#include "InterphaseHeatTransferPhaseSystem.H"
#include "sharpInterfaceHeatTransferModel.H"

#include "HashPtrTable.H"
#include "fvcDiv.H"
#include "fvmSup.H"
#include "fvMatrix.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::InterphaseHeatTransferPhaseSystem<BasePhaseSystem>::
InterphaseHeatTransferPhaseSystem
(
    const fvMesh& mesh
)
:
    BasePhaseSystem(mesh)
{
    this->generatePairsAndSubModels
    (
        "sharpInterfaceHeatTransfer",
        sharpInterfaceHeatTransferModels_
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::InterphaseHeatTransferPhaseSystem<BasePhaseSystem>::
~InterphaseHeatTransferPhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::heatTransferTable>
Foam::InterphaseHeatTransferPhaseSystem<BasePhaseSystem>::
heatTransfer() const
{
    autoPtr<phaseSystem::heatTransferTable> eqnsPtr
    (
        new phaseSystem::heatTransferTable()
    );

    phaseSystem::heatTransferTable& eqns = eqnsPtr();

    forAll(this->phaseModels_, phasei)
    {
        const phaseModel& phase = this->phaseModels_[phasei];

        eqns.set
        (
            phase.name(),
            new fvScalarMatrix(phase.thermo().he(), dimEnergy/dimTime)
        );
    }

    forAllConstIter
    (
        sharpInterfaceHeatTransferModelTable,
        sharpInterfaceHeatTransferModels_,
        sharpInterfaceHeatTransferModelIter
    )
    {
      scalar tol_ = 1e-6;
      const phasePair& pair
      (
          this->phasePairs_[sharpInterfaceHeatTransferModelIter.key()]
      );

      const phaseModel& phase1 = pair.dispersed();
      const phaseModel& phase2 = pair.continuous();

      const volScalarField& alpha2 = phase2;

      const volScalarField K(sharpInterfaceHeatTransferModelIter()->K());
      const volScalarField Kd(sharpInterfaceHeatTransferModelIter()->Kd());

      // Implicit  Source Treatment
      *eqns[phase1.name()] -=
            - pos(alpha2-tol_)*K/phase2.thermo().Cpv()*(phase2.thermo().he())
            + fvm::Sp(pos(alpha2-tol_)*K/phase1.thermo().Cpv(), phase1.thermo().he());
      *eqns[phase2.name()] +=
            - fvm::Sp(pos(alpha2-tol_)*K/phase2.thermo().Cpv(), phase2.thermo().he())
            + pos(alpha2-tol_)*K/phase1.thermo().Cpv()*(phase1.thermo().he());
    }

    return eqnsPtr;
}


template<class BasePhaseSystem>
void Foam::InterphaseHeatTransferPhaseSystem<BasePhaseSystem>::
correctEnergyTransport()
{
    BasePhaseSystem::correctEnergyTransport();
}


template<class BasePhaseSystem>
void Foam::InterphaseHeatTransferPhaseSystem<BasePhaseSystem>::
correctInterfaceThermo()
{}


template<class BasePhaseSystem>
bool Foam::InterphaseHeatTransferPhaseSystem<BasePhaseSystem>::read()
{
    if (BasePhaseSystem::read())
    {
        bool readOK = true;

        // Models ...

        return readOK;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
