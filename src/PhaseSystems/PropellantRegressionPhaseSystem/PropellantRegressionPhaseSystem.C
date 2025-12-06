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
    solveParticle(true),
    Tad
    (
      volScalarField
      (
        IOobject("Tadiabatic", mesh), mesh,
        dimensionedScalar("", dimTemperature, this->template get<scalar>("Tad"))
      )
    ),
    Hs1
    (
      volScalarField
      (
        IOobject("Hs1", mesh), mesh,
        dimensionedScalar("", dimVelocity*dimVelocity, 0)
      )
    ),
    Hs2
    (
      volScalarField
      (
        IOobject("Hs2", mesh), mesh,
        dimensionedScalar("", dimVelocity*dimVelocity, 0)
      )
    ),
    RTf
    (
      volScalarField
      (
        IOobject("RTf", mesh), mesh,
        dimensionedScalar("", dimEnergy/dimMass, 0)
      )
    ),
    rb_
    (
      volScalarField
      (
        IOobject("burningRate", mesh), mesh,
        dimensionedScalar("", dimVelocity,0)
      )
    ),
    alphaOld
    (
      volScalarField
      (
        IOobject
        (
          "alphaOld",
          mesh.time().timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless,0)
      )
    ),
    rhoPropellant("rhoprop", dimDensity, this->template get<scalar>("propellantRho")),
    Ug_
    (
      volVectorField
      (
        IOobject("Ugas", mesh), mesh,
        dimensionedVector("", dimVelocity, vector(0, 0, 0))
      )
    ),
    Up_
    (
      volVectorField
      (
        IOobject("Uparticle", mesh), mesh,
        dimensionedVector("", dimVelocity, vector(0, 0, 0))
      )
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

        // Set Source terms
        const phaseModel& phase1 = pair.phase1();
        const phaseModel& phase2 = pair.phase2();
        Hs1 = phase1.thermo().he(phase1.thermo().p(), Tad);
        Hs2 = phase2.thermo().he(phase2.thermo().p(), Tad);

        // Set Propellant volume fraction field
        word propellant = "alpha." + interfaceTrackingModelIter()->propellant_;
        alphaOld = this->db().template lookupObject<volScalarField>(propellant);
    }

    // Coefficient of mass transfer
    forAllConstIter
    (
      interfaceTrackingModelTable,
      interfaceTrackingModels_,
      interfaceTrackingModelIter
    )
    {
        // Mass fraction of particles in combustion products
        scalar Xp = this->template get<scalar>("Xp");
        if(Xp == 1.0)
        {
            FatalErrorInFunction 
              << "Mass fraction of particles in combustion products cannot be 1.0" 
              << exit(FatalError);
        }
        if(Xp == 0.0)
        {
            solveParticle = false;
        }

        this->coeff_.set
        (
            interfaceTrackingModelIter.key(),
            Xp
        );

        // Calculate RT of the gas phase at the propellant surface
        const phasePair& pair = this->phasePairs_[interfaceTrackingModelIter.key()];
        const phaseModel& phase2 = pair.phase2();
        RTf = dimensionedScalar("R", dimEnergy/dimMoles/dimTemperature, 8314.5)*Tad/phase2.thermo().W();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::PropellantRegressionPhaseSystem<BasePhaseSystem>::
~PropellantRegressionPhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

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

    phaseSystem::massTransferTable& eqns = eqnsPtr();

    forAll(this->phaseModels_, phasei)
    {
        const phaseModel& phase = this->phaseModels_[phasei];

        const PtrList<volScalarField>& Yi = phase.Y();

        forAll(Yi, i)
        {
            eqns.set
            (
                Yi[i].name(),
                new fvScalarMatrix(Yi[i], dimMass/dimTime)
            );
        }
    }

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

    // Equations
    fvScalarMatrix& eqn1 = *eqns[phase1.name()];
    fvScalarMatrix& eqn2 = *eqns[phase2.name()];

    // Energy Source 1
    eqn1 += - fvm::Sp(coeff*rDmdt, eqn1.psi())
            + coeff*rDmdt*Hs1;
    eqn2 += - fvm::Sp((1.0 - coeff)*rDmdt, eqn2.psi())
            + (1.0 - coeff)*rDmdt*Hs2;

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
            + coeff*rDmdt*Up_;
    eqn2 += - fvm::Sp((1.0 - coeff)*rDmdt, eqn2.psi())
            + (1.0 - coeff)*rDmdt*Ug_;
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

        interfaceTrackingModelIter()->regress(alpha, alphaOld);
    }
  
    // Solve other phase volume fraction equations if required
    if (solveParticle)
    {
        BasePhaseSystem::solve();
    }
    else
    {
        forAllIter
        (
            interfaceTrackingModelTable,
            interfaceTrackingModels_,
            interfaceTrackingModelIter
        )
        {
            const phasePair& pair = this->phasePairs_[interfaceTrackingModelIter.key()];

            word propellant = "alpha." + interfaceTrackingModelIter()->propellant_;
            const volScalarField& alphaPropellant = this->db().template lookupObject<volScalarField>(propellant);

            forAll(this->phases(), phasei)
            {
                phaseModel& phase = this->phases()[phasei];
                volScalarField& alpha = phase;
                
                if (phase.stationary())
                {
                    Info << phase.name() << " fraction, min, max = "
                         << phase.weightedAverage(phase.mesh().V()).value()
                         << ' ' << min(phase).value()
                         << ' ' << max(phase).value()
                         << endl;

                    continue;
                };
                
                if (phase.name() == pair.phase2().name()) // Gas phase
                {
                    alpha = 1.0 - alphaPropellant;
                    phase.alphaPhiRef() = fvc::interpolate(alpha)*phase.phi();
                    phase.alphaRhoPhiRef() =
                    fvc::interpolate(phase.rho())*phase.alphaPhi();
                }
                else // particle phase
                {
                    alpha = 0*alphaPropellant;
                    phase.alphaPhiRef() = 0*phase.phi();
                    phase.alphaRhoPhiRef() =
                    fvc::interpolate(phase.rho())*phase.alphaPhi();
                }
                phase.clip(SMALL, 1 - SMALL);

                Info << phase.name() << " fraction, min, max = "
                     << phase.weightedAverage(phase.mesh().V()).value()
                     << ' ' << min(phase).value()
                     << ' ' << max(phase).value()
                     << endl;
            }
            
        }
    }
    
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
        rb_ = interfaceTrackingModelIter()->rb();
        *rDmdt_[interfaceTrackingModelIter.key()]
              = interfaceTrackingModelIter()->dmdt()()*rhoPropellant;
    }

    // calculate velocity of the gas and particle source
    calculateVelocity();

}

template<class BasePhaseSystem>
void Foam::PropellantRegressionPhaseSystem<BasePhaseSystem>::store()
{
    BasePhaseSystem::store();

    forAllIter
    (
        interfaceTrackingModelTable,
        interfaceTrackingModels_,
        interfaceTrackingModelIter
    )
    {
        word propellant = "alpha." + interfaceTrackingModelIter()->propellant_;
        alphaOld = this->db().template lookupObject<volScalarField>(propellant);
    }
}

template<class BasePhaseSystem>
void Foam::PropellantRegressionPhaseSystem<BasePhaseSystem>::calculateVelocity()
{
    //- Calculate velocity of the gas and particles entering
    //                the combustion chamber

    forAllIter
    (
        interfaceTrackingModelTable,
        interfaceTrackingModels_,
        interfaceTrackingModelIter
    )
    {
        const phasePair& pair = this->phasePairs_[interfaceTrackingModelIter.key()];
        const phaseModel& phase1 = pair.phase1();
        const phaseModel& phase2 = pair.phase2();
        const scalar Xp = coeff_[interfaceTrackingModelIter.key()];

        const tmp<volScalarField> tAs = interfaceTrackingModelIter()->As();
        const tmp<volScalarField> trb = interfaceTrackingModelIter()->rb();
        const volScalarField& dmdt = *rDmdt_[interfaceTrackingModelIter.key()];

        const volScalarField& As(tAs());
        const volScalarField& rb(trb());

        const tmp<volScalarField> trhogf(phase2.thermo().p()/RTf);
        const volScalarField& rhogf(trhogf());
        
        const tmp<volScalarField> talphagf(1.0/(1.0 + (rhogf/phase1.rho())*(Xp/(1 - Xp))));
        const volScalarField& alphagf(talphagf());

        forAll(Ug_, i)
        {
            if(As[i] != 0)
            {
                Ug_[i] = 
                (
                    (rb[i] - (1 - Xp)*dmdt[i]/
                        (alphagf[i]*As[i]*rhogf[i]))*vector(1, 0, 0)
                );
            }
        }
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
