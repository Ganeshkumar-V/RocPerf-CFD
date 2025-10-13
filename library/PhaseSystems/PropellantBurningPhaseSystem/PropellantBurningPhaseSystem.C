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

#include "PropellantBurningPhaseSystem.H"
#include "burningRateModel.H"
#include "fvmSup.H"
#include "phaseSystem.H"
#include "fvmLaplacian.H"

// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::PropellantBurningPhaseSystem<BasePhaseSystem>::rDmdt
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
Foam::PropellantBurningPhaseSystem<BasePhaseSystem>::PropellantBurningPhaseSystem
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
    zeta(this->template get<scalar>("Efficiency")),
    rb_
    (
      volScalarField
      (
        IOobject
        (
          "rb",
          mesh.time().timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimVelocity,0)
      )
    ),
    R_("R", dimEnergy/dimMass/dimTemperature, 0),
    rhoPropellant("rhoprop", dimDensity, this->template get<scalar>("propellantRho")),
    rhoParticle("rhopar", dimDensity, this->template get<scalar>("particleRho")),
    alphaRhoAl("alphaRhoAl", dimDensity, 0),
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
          IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless,0)
      )
    ),
    regressionAlpha
    (
      volScalarField
      (
        IOobject
        (
          "regressionAlpha",
          mesh.time().timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless,0)
      )
    ),
    epsilon(this->template get<scalar>("trapingFactor")),
    Ug_
    (
      volVectorField
      (
        IOobject
        (
          "Ug",
          mesh.time().timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("", dimVelocity, vector(0, 0, 0))
      )
    ),
    Up_
    (
      volVectorField
      (
        IOobject
        (
          "Up",
          mesh.time().timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("", dimVelocity, vector(0, 0, 0))
      )
    ),
    saturationModel_
    (
        saturationModel::New(this->subDict("saturationModel"), mesh)
    )
{
    this->generatePairsAndSubModels
    (
        "burningRate",
        burningRateModels_
    );

    forAllConstIter
    (
        burningRateModelTable,
        burningRateModels_,
        burningRateModelIter
    )
    {
        const phasePair& pair = this->phasePairs_[burningRateModelIter.key()];

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

        nHat_.set
        (
          pair,
          new volVectorField
          (
            IOobject
            (
              IOobject::groupName("nHat", pair.name()),
              this->mesh().time().timeName(),
              this->mesh(),
              IOobject::READ_IF_PRESENT,
              IOobject::AUTO_WRITE
            ),
            this->mesh(),
            dimensionedVector(dimless)
          )
        );
    }

    // molecularWeights
    MW_.Al = molecularWeights_.get<scalar>("Al");
    MW_.Al2O3 = molecularWeights_.get<scalar>("Al2O3");
    MW_.H2O = molecularWeights_.get<scalar>("H2O");
    MW_.H2 = molecularWeights_.get<scalar>("H2");
    AAlC_ = this->template get<scalar>("activeAlContent");
    K_ = (MW_.Al/MW_.Al2O3)*(1/AAlC_ - 1);
    MW_.Prop = MW_.Al + K_*MW_.Al2O3 + (1.5/eqR_)*MW_.H2O;
    regress_ = this->template get<bool>("regression");

    R_.value() = 8314.5/MW_.H2;
    alphaRhoAl.value() = rhoPropellant.value()*eqR_/(1 + eqR_);

    // Old Alpha
    forAllIter
    (
      burningRateModelTable,
      burningRateModels_,
      burningRateModelIter
    )
    {
      word propellant = "alpha." + burningRateModelIter()->propellant_;
      alphaOld = this->db().template lookupObject<volScalarField>(propellant);
      regressionAlpha = alphaOld;
    }

    // Coefficient of mass transfer
    forAllConstIter
    (
      burningRateModelTable,
      burningRateModels_,
      burningRateModelIter
    )
    {
      scalar gcoeff = 1.0;
      scalar pcoeff = 1.0;

      if (eqR_ <= 1.0) // Lean or Stoichiometric Mixture
      {
        // pcoeff = AAlC_*MW_.Al2O3/(2*MW_.Al*(1 + 1/eqR_)) + (1 - AAlC_);
        // gcoeff = AAlC_*(1.0 - MW_.Al2O3/(2*MW_.Al*(1 + 1/eqR_)));
        pcoeff = epsilon*(K_ + 0.5)*MW_.Al2O3/MW_.Prop;
        gcoeff = (1.5*MW_.H2 + ((1.5/eqR_) - 1.5)*MW_.H2O)/MW_.Prop;
      }
      else  // Rich mixture
      {
        FatalErrorInFunction << "Equivalence Ratio is greater than 1.0 : "
        << " Fuel Rich Mixture Coefficients is not implemented"
        << exit(FatalError);
      }

      this->pCoeff_.set
      (
        burningRateModelIter.key(),
        pcoeff
      );
      this->gCoeff_.set
      (
        burningRateModelIter.key(),
        gcoeff
      );
    }

    // Exit if efficiency is 0%
    if (zeta == 0)
    {
        FatalErrorInFunction
            << "Efficiency cannot be 0%\n"
            << "Choose any number -> (0, 1.0]" << exit(FatalError);
    }
    if (zeta > 1.0)
    {
        FatalErrorInFunction
            << "Efficiency cannot be more than 100%\n"
            << "Choose any number -> (0, 1.0]" << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::PropellantBurningPhaseSystem<BasePhaseSystem>::
~PropellantBurningPhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
template<class BasePhaseSystem>
const Foam::saturationModel&
Foam::PropellantBurningPhaseSystem<BasePhaseSystem>::saturation() const
{
    return saturationModel_();
}

template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::PropellantBurningPhaseSystem<BasePhaseSystem>::dmdt
(
    const phasePairKey& key
) const
{
    // return BasePhaseSystem::dmdt(key) + this->rDmdt(key);
    return BasePhaseSystem::dmdt(key);
}


template<class BasePhaseSystem>
Foam::PtrList<Foam::volScalarField>
Foam::PropellantBurningPhaseSystem<BasePhaseSystem>::dmdts() const
{
    PtrList<volScalarField> dmdts(BasePhaseSystem::dmdts());

    // Fill the mass transfer rates with zero
    this->fillFields("dmdt", dimDensity/dimTime, dmdts);

    forAllConstIter(rDmdtTable, rDmdt_, rDmdtIter)
    {
        const phasePair& pair = this->phasePairs_[rDmdtIter.key()];
        const volScalarField& rDmdt = *rDmdtIter();

        const scalar pcoeff = pCoeff_[rDmdtIter.key()];
        const scalar gcoeff = gCoeff_[rDmdtIter.key()];

        this->addField(pair.phase1(), "dmdt", pcoeff*rDmdt, dmdts);
        this->addField(pair.phase2(), "dmdt", gcoeff*rDmdt, dmdts);

        // Subtract for Propellant Phase
        this->addField(this->phases()[2], "dmdt", -rDmdt, dmdts);
    }
    return dmdts;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::massTransferTable>
Foam::PropellantBurningPhaseSystem<BasePhaseSystem>::massTransfer() const
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

    forAllConstIter
    (
      burningRateModelTable,
      burningRateModels_,
      burningRateModelIter
    )
    {
        const phasePair& pair(this->phasePairs_[burningRateModelIter.key()]);

        const phaseModel& phase = pair.continuous();
        const volScalarField dmdt(this->rDmdt(pair));
        const dimensionedScalar coeff(dimless, gCoeff_[burningRateModelIter.key()]);

        const PtrList<volScalarField>& Yi = phase.Y();

        const dimensionedScalar fH2(dimless, 1.5*MW_.H2/MW_.Prop);
        const dimensionedScalar fH2O(dimless, (1.5/eqR_ - 1.5)*MW_.H2O/MW_.Prop);

        if (min(dmdt).value() < 0)
        {
          Info << "min(dmdt): " << min(dmdt).value() << endl;
          FatalErrorInFunction
              << "Mass Transfer(): dmdt or one of the factors are negative"
              << exit(FatalError);
        }

        if (eqns.size() != 0)
        {
          *eqns[Yi[0].name()] += dmdt*fH2;
          *eqns[Yi[1].name()] += dmdt*fH2O;
        }

    }
    return eqnsPtr;
}

template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::heatTransferTable>
Foam::PropellantBurningPhaseSystem<BasePhaseSystem>::heatTransfer() const
{
  autoPtr<phaseSystem::heatTransferTable> eqnsPtr =
          BasePhaseSystem::heatTransfer();

  phaseSystem::heatTransferTable& eqns = eqnsPtr();

  forAllConstIter(rDmdtTable, rDmdt_, rDmdtIter)
  {
    const phasePair& pair = this->phasePairs_[rDmdtIter.key()];
    const volScalarField& rDmdt = *rDmdtIter();

    const scalar pcoeff = pCoeff_[rDmdtIter.key()];
    const scalar gcoeff = gCoeff_[rDmdtIter.key()];

    const phaseModel& phase1 = pair.phase1();
    const phaseModel& phase2 = pair.phase2();

    // Enthalpy source
    const tmp<volScalarField> ths1(phase1.thermo().he(phase1.thermo().p(), Tad));
    const volScalarField& hs1(ths1());

    const tmp<volScalarField> ths2(phase2.thermo().he(phase2.thermo().p(), Tad));
    const volScalarField& hs2(ths2());

    // Equations
    fvScalarMatrix& eqn1 = *eqns[phase1.name()];
    fvScalarMatrix& eqn2 = *eqns[phase2.name()];

    eqn1 += - fvm::Sp(pcoeff*rDmdt, eqn1.psi())
            + pcoeff*rDmdt*hs1;
    eqn2 += - fvm::Sp(gcoeff*rDmdt, eqn2.psi())
            + gcoeff*rDmdt*hs2;
  }

  return eqnsPtr;
}

template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::momentumTransferTable>
Foam::PropellantBurningPhaseSystem<BasePhaseSystem>::momentumTransfer()
{
  autoPtr<phaseSystem::momentumTransferTable> eqnsPtr =
                BasePhaseSystem::momentumTransfer();

  phaseSystem::momentumTransferTable& eqns = eqnsPtr();

  forAllConstIter(rDmdtTable, rDmdt_, rDmdtIter)
  {
    const phasePair& pair = this->phasePairs_[rDmdtIter.key()];
    const volScalarField& rDmdt = *rDmdtIter();

    const scalar pcoeff = pCoeff_[rDmdtIter.key()];
    const scalar gcoeff = gCoeff_[rDmdtIter.key()];

    const phaseModel& phase1 = pair.phase1();
    const phaseModel& phase2 = pair.phase2();

    // Equations
    fvVectorMatrix& eqn1 = *eqns[phase1.name()];
    fvVectorMatrix& eqn2 = *eqns[phase2.name()];

    // Momentum Source
    eqn1 += - fvm::Sp(pcoeff*rDmdt, eqn1.psi())
            + pcoeff*rDmdt*Up_;
    eqn2 += - fvm::Sp(gcoeff*rDmdt, eqn2.psi())
            + gcoeff*rDmdt*Ug_;
  }

  return eqnsPtr;
}

template<class BasePhaseSystem>
void Foam::PropellantBurningPhaseSystem<BasePhaseSystem>::solve()
{
  // Regress Propellant surface (Manipulate propellant volume fraction)
  if (regress_)
  {
    forAllIter
    (
      burningRateModelTable,
      burningRateModels_,
      burningRateModelIter
    )
    {
      word propellant = "alpha." + burningRateModelIter()->propellant_;
      volScalarField& alpha = this->db().template lookupObjectRef<volScalarField>(propellant);
      burningRateModelIter()->regress(alpha, alphaOld);
      alpha.clip(SMALL, 1-SMALL);
    }
  }

  // Solve other phase volume fraction equations
  BasePhaseSystem::solve();
}

template<class BasePhaseSystem>
void Foam::PropellantBurningPhaseSystem<BasePhaseSystem>::correct()
{
    BasePhaseSystem::correct();

    //- Finds burning rate (rb = aP^n)
    forAllIter
    (
        burningRateModelTable,
        burningRateModels_,
        burningRateModelIter
    )
    {
        // burningRateModelIter()->correct();
        *rDmdt_[burningRateModelIter.key()]
                  = dimensionedScalar(dimDensity/dimTime);
    }

    //- return burning Rate, As and propellant density -> find mdot of propellant
    forAllConstIter
    (
        burningRateModelTable,
        burningRateModels_,
        burningRateModelIter
    )
    {
        *rDmdt_[burningRateModelIter.key()]
              = burningRateModelIter()->dmdt()*rhoPropellant;
        rb_ = burningRateModelIter()->rb();
        *nHat_[burningRateModelIter.key()]
              = burningRateModelIter()->nHat();
    }

    // calculate velocity of the gas and particle source
    calculateVelocity();
}

template<class BasePhaseSystem>
void Foam::PropellantBurningPhaseSystem<BasePhaseSystem>::store()
{
  forAllIter
  (
    burningRateModelTable,
    burningRateModels_,
    burningRateModelIter
  )
  {
    word propellant = "alpha." + burningRateModelIter()->propellant_;
    alphaOld = this->db().template lookupObject<volScalarField>(propellant);
  }
}

template<class BasePhaseSystem>
void Foam::PropellantBurningPhaseSystem<BasePhaseSystem>::calculateVelocity()
{
    //- Calculate velocity of the gas and particles comes
    //                into the combustion chamber

    // density of gas phase entering,
    const tmp<volScalarField> trhog(this->phases()[0].thermo().p()/(Tad*R_));
    const volScalarField& rhog(trhog());

    // volume fraction of particle phase
    const tmp<volScalarField> talphap(1/(1 + (rhoParticle/rhog)*(zeta/(8 + 9*eqR_))));
    const volScalarField& alphap(talphap());

    // Velocity of gas and particle phase
    forAllConstIter(interfaceTable, nHat_, nHatIter)
    {
      const volVectorField& nHat = *nHatIter();

      Ug_ = -alphaRhoAl*rb_*nHat/(9.0*eqR_*(1 - alphap)*rhog);
      Up_ = Ug_*zeta;

      Ug_ += rb_*nHat;
      Up_ += rb_*nHat;
    }
}

template<class BasePhaseSystem>
bool Foam::PropellantBurningPhaseSystem<BasePhaseSystem>::read()
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
