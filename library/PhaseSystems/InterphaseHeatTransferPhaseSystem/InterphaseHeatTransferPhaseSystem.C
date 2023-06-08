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
    BasePhaseSystem(mesh),
    Tslip_
    (
      volScalarField
      (
          IOobject
          (
              "Tslip",
              this->mesh().time().timeName(),
              this->mesh(),
              IOobject::NO_READ,
              IOobject::AUTO_WRITE
          ),
          this->mesh(),
          dimensionedScalar(dimTemperature, 0.0)
      )
    )
{
    this->generatePairsAndSubModels
    (
        "sharpInterfaceHeatTransfer",
        sharpInterfaceHeatTransferModels_
    );

    forAllConstIter
    (
        sharpInterfaceHeatTransferModelTable,
        sharpInterfaceHeatTransferModels_,
        sharpInterfaceHeatTransferModelIter
    )
    {
      const phasePair& pair
      (
          this->phasePairs_[sharpInterfaceHeatTransferModelIter.key()]
      );

      const phaseModel& phase1 = pair.dispersed();
      const phaseModel& phase2 = pair.continuous();

      Tslip_ = phase1.thermo().T() - phase2.thermo().T();
    }
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
    autoPtr<phaseSystem::heatTransferTable> eqnsPtr =
            BasePhaseSystem::heatTransfer();

    phaseSystem::heatTransferTable& eqns = eqnsPtr();

    forAllConstIter
    (
        sharpInterfaceHeatTransferModelTable,
        sharpInterfaceHeatTransferModels_,
        sharpInterfaceHeatTransferModelIter
    )
    {
      const phasePair& pair
      (
          this->phasePairs_[sharpInterfaceHeatTransferModelIter.key()]
      );

      const phaseModel& phase1 = pair.dispersed();
      const phaseModel& phase2 = pair.continuous();

      const tmp<volScalarField> tK(sharpInterfaceHeatTransferModelIter()->K());
      const volScalarField& K(tK());

      const tmp<volScalarField> tCp1(phase1.thermo().Cpv());
      const volScalarField& Cp1(tCp1());

      const tmp<volScalarField> tCp2(phase2.thermo().Cpv());
      const volScalarField& Cp2(tCp2());

      const volScalarField& he1(phase1.thermo().he());
      const volScalarField& he2(phase2.thermo().he());

      // Explicit Source Treatment
      const tmp<volScalarField> tT1(phase1.thermo().T());
      const volScalarField& T1(tT1());

      const tmp<volScalarField> tT2(phase2.thermo().T());
      const volScalarField& T2(tT2());

      // *eqns[phase1.name()] -= K*(T1 - T2);
      // *eqns[phase2.name()] += K*(T1 - T2);

      *eqns[phase1.name()] -= K*(T1 - T2) - fvm::Sp(K/Cp1, he1) + K/Cp1*he1;
      *eqns[phase2.name()] += K*(T1 - T2) - fvm::Sp(K/Cp2, he2) + K/Cp2*he2;

      // // Implicit Source Treatment
      // *eqns[phase1.name()] -=
      //       - K/Cp2*he2 + fvm::Sp(K/Cp1, he1);
      // *eqns[phase2.name()] +=
      //       - fvm::Sp(K/Cp2, he2) + K/Cp1*he1;

      // Addition Heat diffusion in particle phase
      // *eqns[phase1.name()] -=
      //       - fvm::laplacian
      //         (
      //          fvc::interpolate(phase1*Kd/Cp1),
      //          he1
      //         );
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
void Foam::InterphaseHeatTransferPhaseSystem<BasePhaseSystem>::
store()
{
  forAllConstIter
  (
      sharpInterfaceHeatTransferModelTable,
      sharpInterfaceHeatTransferModels_,
      sharpInterfaceHeatTransferModelIter
  )
  {
    const phasePair& pair
    (
        this->phasePairs_[sharpInterfaceHeatTransferModelIter.key()]
    );

    const phaseModel& phase1 = pair.dispersed();
    const phaseModel& phase2 = pair.continuous();

    Tslip_ = phase1.thermo().T() - phase2.thermo().T();
  }
}

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
