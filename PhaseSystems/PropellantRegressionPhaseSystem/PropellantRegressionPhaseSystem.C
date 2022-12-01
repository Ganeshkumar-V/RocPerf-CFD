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

\*---------------------------------------------------------------------------*/

#include "PropellantRegressionPhaseSystem.H"
#include "interfaceTrackingModel.H"
#include "fvmSup.H"
#include "phaseSystem.H"
#include "fvmLaplacian.H"

// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::PropellantRegressionPhaseSystem<BasePhaseSystem>::rDmdt
(
    const phasePairKey& key
) const
{
    if (!rDmdt_.found(key))
    {
        return phaseSystem::dmdt(key);
    }

    const scalar rDmdtSign(Pair<word>::compare(rDmdt_.find(key).key(), key));
    return rDmdtSign**rDmdt_[key];
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::PropellantRegressionPhaseSystem<BasePhaseSystem>::PropellantRegressionPhaseSystem
(
    const fvMesh& mesh
)
:
    BasePhaseSystem(mesh),
    molecularWeights_(this->subDict("molecularWeight")),
    eqR_(this->template get<scalar>("equivalenceRatio")),
    Tad
    (
      volScalarField
      (
        IOobject("Tadiabatic", mesh), mesh,
        dimensionedScalar("", dimTemperature, this->template get<scalar>("Tad"))
      )
    ),
    Ug
    (
      volVectorField
      (
        IOobject("Ugas", mesh), mesh,
        dimensionedVector("", dimVelocity, this->template get<vector>("Ugas"))
      )
    ),
    Up
    (
      volVectorField
      (
        IOobject("Uparticle", mesh), mesh,
        dimensionedVector("", dimVelocity, this->template get<vector>("Uparticle"))
      )
    ),
    saturationModel_
    (
        saturationModel::New(this->subDict("saturationModel"), mesh)
    )
{
    this->generatePairsAndSubModels
    (
        "interfaceTracking",
        interfaceTrackingModels_
    );

    forAllConstIter
    (
        interfaceTrackingModelTable,
        interfaceTrackingModels_,
        interfaceTrackingModelIter
    )
    {
        const phasePair& pair = this->phasePairs_[interfaceTrackingModelIter.key()];

        // Initially assume no mass transfer
        rDmdt_.set
        (
            pair,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("rDmdt", pair.name()),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                this->mesh(),
                dimensionedScalar(dimDensity/dimTime)
            )
        );
    }

    // Coefficient of mass transfer
    forAllConstIter
    (
      interfaceTrackingModelTable,
      interfaceTrackingModels_,
      interfaceTrackingModelIter
    )
    {
      // Molecular weights of Species
      scalar Al = molecularWeights_.get<scalar>("Al");
      scalar Al2O3 = molecularWeights_.get<scalar>("Al2O3");
      scalar H2O = molecularWeights_.get<scalar>("H2O");
      scalar H2 = molecularWeights_.get<scalar>("H2");

      // number of moles of fuel
      scalar zeta = Al/(eqR_*H2O);
      scalar coeff = 1.0;

      if (eqR_ <= 1.0) // Lean or Stoichiometric Mixture
      {
        coeff = Al2O3/(2*Al*(1 + 1/eqR_));
      }
      else  // Rich mixture
      {
        coeff = 1.0 - H2*zeta/(Al*(1 + 1/eqR_));
      }

      this->coeff_.set
      (
        interfaceTrackingModelIter.key(),
        coeff
      );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::PropellantRegressionPhaseSystem<BasePhaseSystem>::
~PropellantRegressionPhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
template<class BasePhaseSystem>
const Foam::saturationModel&
Foam::PropellantRegressionPhaseSystem<BasePhaseSystem>::saturation() const
{
    return saturationModel_();
}

template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::PropellantRegressionPhaseSystem<BasePhaseSystem>::dmdt
(
    const phasePairKey& key
) const
{
    // return BasePhaseSystem::dmdt(key) + this->rDmdt(key);
    return BasePhaseSystem::dmdt(key);
}


template<class BasePhaseSystem>
Foam::PtrList<Foam::volScalarField>
Foam::PropellantRegressionPhaseSystem<BasePhaseSystem>::dmdts() const
{
    PtrList<volScalarField> dmdts(BasePhaseSystem::dmdts());

    // Fill the mass transfer rates with zero
    this->fillFields("dmdt", dimDensity/dimTime, dmdts);

    forAllConstIter(rDmdtTable, rDmdt_, rDmdtIter)
    {
        const phasePair& pair = this->phasePairs_[rDmdtIter.key()];
        const volScalarField& rDmdt = *rDmdtIter();

        const scalar coeff = coeff_[rDmdtIter.key()];
        this->addField(pair.phase1(), "dmdt", coeff*rDmdt, dmdts);
        this->addField(pair.phase2(), "dmdt", (1.0 - coeff)*rDmdt, dmdts);
    }
    return dmdts;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::massTransferTable>
Foam::PropellantRegressionPhaseSystem<BasePhaseSystem>::massTransfer() const
{
    // Create a mass transfer matrix for each species of each phase
    autoPtr<phaseSystem::massTransferTable> eqnsPtr
    (
        new phaseSystem::massTransferTable()
    );
    // (No Species Present)
    return eqnsPtr;
}

template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::heatTransferTable>
Foam::PropellantRegressionPhaseSystem<BasePhaseSystem>::heatTransfer() const
{
  autoPtr<phaseSystem::heatTransferTable> eqnsPtr =
          BasePhaseSystem::heatTransfer();

  phaseSystem::heatTransferTable& eqns = eqnsPtr();

  forAllConstIter(rDmdtTable, rDmdt_, rDmdtIter)
  {
    const phasePair& pair = this->phasePairs_[rDmdtIter.key()];
    const volScalarField& rDmdt = *rDmdtIter();

    const scalar coeff = coeff_[rDmdtIter.key()];

    const phaseModel& phase1 = pair.phase1();
    const phaseModel& phase2 = pair.phase2();

    // specific Heat Capacity
    const volScalarField Cp1(phase1.thermo().Cp());
    const volScalarField Cp2(phase2.thermo().Cp());

    // Equations
    fvScalarMatrix& eqn1 = *eqns[phase1.name()];
    fvScalarMatrix& eqn2 = *eqns[phase2.name()];

    // Enthalpy Source
    eqn1 += - fvm::Sp(coeff*rDmdt, eqn1.psi())
            + coeff*rDmdt*Cp1*Tad;
    eqn2 += - fvm::Sp((1.0 - coeff)*rDmdt, eqn2.psi())
            + (1.0 - coeff)*rDmdt*Cp2*Tad;

    // // Stabilization Terms
    // eqn1 += coeff*rDmdt*phase1.thermo().he()
    //         - fvm::Sp(coeff*rDmdt, eqn1.psi());
    // eqn2 += (1.0 - coeff)*rDmdt*phase2.thermo().he()
    //         - fvm::Sp((1.0 - coeff)*rDmdt, eqn2.psi());


    // Kinetic Energy Source
    eqn1 += - coeff*rDmdt*phase1.K()
            + coeff*rDmdt*(0.5*magSqr(Up));
    eqn2 += - (1.0 - coeff)*rDmdt*phase2.K()
            + (1.0 - coeff)*rDmdt*(0.5*magSqr(Ug));
  }

  return eqnsPtr;
}

template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::momentumTransferTable>
Foam::PropellantRegressionPhaseSystem<BasePhaseSystem>::momentumTransfer()
{
  autoPtr<phaseSystem::momentumTransferTable> eqnsPtr =
                BasePhaseSystem::momentumTransfer();

  phaseSystem::momentumTransferTable& eqns = eqnsPtr();

  forAllConstIter(rDmdtTable, rDmdt_, rDmdtIter)
  {
    const phasePair& pair = this->phasePairs_[rDmdtIter.key()];
    const volScalarField& rDmdt = *rDmdtIter();

    const scalar coeff = coeff_[rDmdtIter.key()];

    const phaseModel& phase1 = pair.phase1();
    const phaseModel& phase2 = pair.phase2();

    // Equations
    fvVectorMatrix& eqn1 = *eqns[phase1.name()];
    fvVectorMatrix& eqn2 = *eqns[phase2.name()];

    // Momentum Source
    eqn1 += - fvm::Sp(coeff*rDmdt, eqn1.psi())
            + coeff*rDmdt*Up;
    eqn2 += - fvm::Sp((1.0 - coeff)*rDmdt, eqn2.psi())
            + (1.0 - coeff)*rDmdt*Ug;
  }

  return eqnsPtr;
}

template<class BasePhaseSystem>
void Foam::PropellantRegressionPhaseSystem<BasePhaseSystem>::solve()
{
  // Regress Propellant surface (Manipulate propellant volume fraction)
  forAllIter
  (
    interfaceTrackingModelTable,
    interfaceTrackingModels_,
    interfaceTrackingModelIter
  )
  {
    word propellant = "alpha." + interfaceTrackingModelIter()->propellant_;
    volScalarField& alpha = this->db().template lookupObjectRef<volScalarField>(propellant);
    interfaceTrackingModelIter()->regress(alpha);
  }

  // Solve other phase volume fraction equations
  BasePhaseSystem::solve();
}

template<class BasePhaseSystem>
void Foam::PropellantRegressionPhaseSystem<BasePhaseSystem>::correct()
{
    BasePhaseSystem::correct();

    //- Finds burning rate (rb = aP^n)
    forAllIter
    (
        interfaceTrackingModelTable,
        interfaceTrackingModels_,
        interfaceTrackingModelIter
    )
    {
        interfaceTrackingModelIter()->correct();
        *rDmdt_[interfaceTrackingModelIter.key()]
                  = dimensionedScalar(dimDensity/dimTime);
    }

    //- return burning Rate, As and propellant density -> find mdot of propellant
    forAllConstIter
    (
        interfaceTrackingModelTable,
        interfaceTrackingModels_,
        interfaceTrackingModelIter
    )
    {
        *rDmdt_[interfaceTrackingModelIter.key()]
              = interfaceTrackingModelIter()->rb()
                *interfaceTrackingModelIter()->As()
                *interfaceTrackingModelIter()->rhop();
    }

}

template<class BasePhaseSystem>
bool Foam::PropellantRegressionPhaseSystem<BasePhaseSystem>::read()
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
