/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "multiphaseKineticTheoryModel.H"
#include "mathematicalConstants.H"
#include "fvOptions.H"
#include "processorFvPatch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RASModels::multiphaseKineticTheoryModel::multiphaseKineticTheoryModel
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& phase,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity
    <
        RASModel<EddyDiffusivity<phaseCompressibleTurbulenceModel>>
    >
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        phase,
        propertiesName
    ),

    phase_(phase),
    kineticCoeffDict_(coeffDict_.subDict("kineticTheoryCoeffs")),
    otherPhase_
    (
      this->db().lookupObject<phaseModel>
      (
        coeffDict_.get<word>("otherPhase")
      )
    ),
    viscosityModel_
    (
        kineticTheoryModels::viscosityModel::New
        (
            kineticCoeffDict_
        )
    ),
    conductivityModel_
    (
        kineticTheoryModels::conductivityModel::New
        (
            kineticCoeffDict_
        )
    ),
    radialModel_
    (
        kineticTheoryModels::radialModel::New
        (
            kineticCoeffDict_
        )
    ),
    granularPressureModel_
    (
        kineticTheoryModels::granularPressureModel::New
        (
            kineticCoeffDict_
        )
    ),
    frictionalStressModel_
    (
        kineticTheoryModels::frictionalStressModel::New
        (
            kineticCoeffDict_
        )
    ),

    equilibrium_(kineticCoeffDict_.get<bool>("equilibrium")),
    e_("e", dimless, kineticCoeffDict_),
    alphaMax_("alphaMax", dimless, kineticCoeffDict_),
    alphaMinFriction_("alphaMinFriction", dimless, kineticCoeffDict_),
    residualAlpha_("residualAlpha", dimless, kineticCoeffDict_),
    maxNut_("maxNut", dimViscosity, 1000, kineticCoeffDict_),

    ThetaSmall("", dimVelocity*dimVelocity, kineticCoeffDict_.get<scalar>("minTheta")),

    Theta_
    (
        IOobject
        (
            IOobject::groupName("Theta", phase.name()),
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),

    lambda_
    (
        IOobject
        (
            IOobject::groupName("lambda", phase.name()),
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedScalar(dimViscosity, Zero)
    ),

    gs0_
    (
        IOobject
        (
            IOobject::groupName("gs0", phase.name()),
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedScalar(dimless, Zero)
    ),

    kappa_
    (
        IOobject
        (
            IOobject::groupName("kappa", phase.name()),
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedScalar(dimDynamicViscosity, Zero)
    ),

    nuFric_
    (
        IOobject
        (
            IOobject::groupName("nuFric", phase.name()),
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar(dimViscosity, Zero)
    ),

    cutoff_(kineticCoeffDict_.get<scalar>("cutoff"))
{
    if (type == typeName)
    {
        printCoeffs(type);
    }

    volScalarField Pp(
        max(Theta_, ThetaSmall)*
        granularPressureModel_->granularPressureCoeff
        (
            alpha,
            radialModel_->g0(alpha, alphaMinFriction_, alphaMax_),
            rho,
            e_
        )
    );
    Info << "min/max Pp: " << min(Pp).value() << " - " << max(Pp).value() << endl;
    Info << "min/max P: " << min(otherPhase_.thermo().p()).value() << " - " << max(otherPhase_.thermo().p()).value() << endl;
    Info << "min/max Theta: " << min(Theta_).value() << " - " << max(Theta_).value() << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RASModels::multiphaseKineticTheoryModel::~multiphaseKineticTheoryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
namespace Foam
{
void ImposeWall(volScalarField& psi, const volScalarField& alpha)
{
  const labelList& Own(psi.mesh().owner());
  const labelList& Nei(psi.mesh().neighbour());

  scalar One(0.999); //1.0 - SMALL);

  forAll(Own, celli)
  {
    if ((alpha[Own[celli]] < One) && (alpha[Nei[celli]] >= One))
    {
      psi[Nei[celli]] = psi[Own[celli]];
    }
    else if ((alpha[Nei[celli]] < One) && (alpha[Own[celli]] >= One))
    {
      psi[Own[celli]] = psi[Nei[celli]];
    }
  }

  forAll(psi.mesh().boundary(),patchi)
  {
    if (isType<processorFvPatch>(psi.mesh().boundary()[patchi]))
    {
      const processorPolyPatch& pp
          = refCast<const processorPolyPatch>(psi.mesh().boundaryMesh()[patchi]);
      if (pp.owner())
      {
        const scalarField& psibF(psi.boundaryField()[patchi].patchNeighbourField());
        const scalarField& alphaF(alpha.boundaryField()[patchi].patchInternalField());
        const scalarField& alphabF(alpha.boundaryField()[patchi].patchNeighbourField());
        const labelList& fC(psi.mesh().boundary()[patchi].faceCells());
        forAll(fC, i)
        {
          if ((alphabF[i] < One) && (alphaF[i] >= One))
          {
            psi[fC[i]] = psibF[i];
          }
        }
      }
    }
  }

}

void detectZero(const volScalarField& psi)
{
    const volScalarField::Internal& psiIF(psi.internalField());
    const word name(psi.name());

    forAll(psiIF, i)
    {
        if (psiIF[i] == 0)
        {
            Info << name << " at " << i << " -> " << psiIF[i] << endl;
        }
    }

    const volScalarField::Boundary& psibF(psi.boundaryField());
    forAll(psibF, bFi)
    {
        const scalarField& sF(psi.boundaryField()[bFi]);
        const word patchName(psi.mesh().boundary()[bFi].name());

        forAll(sF, sFi)
        {
            if (sF[sFi] == 0)
            {
                Info << name << " -> " << patchName << " at " << sFi << " -> " << 0 << endl;
            }
        }
    }
}

}

bool Foam::RASModels::multiphaseKineticTheoryModel::read()
{
    if
    (
        eddyViscosity
        <
            RASModel<EddyDiffusivity<phaseCompressibleTurbulenceModel>>
        >::read()
    )
    {
        coeffDict().readEntry("equilibrium", equilibrium_);
        e_.readIfPresent(coeffDict());
        alphaMax_.readIfPresent(coeffDict());
        alphaMinFriction_.readIfPresent(coeffDict());

        viscosityModel_->read();
        conductivityModel_->read();
        radialModel_->read();
        granularPressureModel_->read();
        frictionalStressModel_->read();

        return true;
    }

    return false;
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::multiphaseKineticTheoryModel::k() const
{
    NotImplemented;
    return nut_;
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::multiphaseKineticTheoryModel::epsilon() const
{
    NotImplemented;
    return nut_;
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::multiphaseKineticTheoryModel::omega() const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::volSymmTensorField>
Foam::RASModels::multiphaseKineticTheoryModel::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                IOobject::groupName("R", U_.group()),
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
          - (nut_)*dev(twoSymm(fvc::grad(U_)))
          - (lambda_*fvc::div(phi_))*symmTensor::I
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::multiphaseKineticTheoryModel::pPrime() const
{
    const volScalarField& rho = phase_.rho();

    tmp<volScalarField> tpPrime
    (
        Theta_
       *granularPressureModel_->granularPressureCoeffPrime
        (
            alpha_,
            radialModel_->g0(alpha_, alphaMinFriction_, alphaMax_),
            radialModel_->g0prime(alpha_, alphaMinFriction_, alphaMax_),
            rho,
            e_
        )
     // +  frictionalStressModel_->frictionalPressurePrime
     //    (
     //        phase_,
     //        alphaMinFriction_,
     //        alphaMax_
     //    )
    );

    volScalarField::Boundary& bpPrime =
        tpPrime.ref().boundaryFieldRef();

    forAll(bpPrime, patchi)
    {
        if (!bpPrime[patchi].coupled())
        {
            bpPrime[patchi] == 0;
        }
    }

    return tpPrime;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::RASModels::multiphaseKineticTheoryModel::pPrimef() const
{
    return fvc::interpolate(pPrime());
}


Foam::tmp<Foam::volSymmTensorField>
Foam::RASModels::multiphaseKineticTheoryModel::devRhoReff() const
{
    return devRhoReff(U_);
}


Foam::tmp<Foam::volSymmTensorField>
Foam::RASModels::multiphaseKineticTheoryModel::devRhoReff
(
    const volVectorField& U
) const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                IOobject::groupName("devRhoReff", U.group()),
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
          - (rho_*nut_)
           *dev(twoSymm(fvc::grad(U)))
          - ((rho_*lambda_)*fvc::div(phi_))*symmTensor::I
        )
    );
}


Foam::tmp<Foam::fvVectorMatrix>
Foam::RASModels::multiphaseKineticTheoryModel::divDevRhoReff
(
    volVectorField& U
) const
{
    return
    (
      - fvm::laplacian(rho_*nut_, U)
      - fvc::div
        (
            (rho_*nut_)*dev2(T(fvc::grad(U)))
          + ((rho_*lambda_)*fvc::div(phi_))
           *dimensioned<symmTensor>("I", dimless, symmTensor::I)
        )
    );
}


void Foam::RASModels::multiphaseKineticTheoryModel::correct()
{
    // Local references
    volScalarField alpha(max(alpha_, scalar(0)));
    const volScalarField& rho = phase_.rho();

    Theta_.clip(ThetaSmall, max(Theta_));

    const scalar sqrtPi = sqrt(constant::mathematical::pi);
    // dimensionedScalar ThetaSmall("ThetaSmall", Theta_.dimensions(), 1e-6);
    // dimensionedScalar ThetaSmallSqrt(sqrt(ThetaSmall));

    tmp<volScalarField> tda(phase_.d());
    const volScalarField& da = tda();

    tmp<volTensorField> tgradU(fvc::grad(U_));
    const volTensorField& gradU(tgradU());
    volSymmTensorField D(symm(gradU));

    // Calculating the radial distribution function
    gs0_ = radialModel_->g0(alpha, alphaMinFriction_, alphaMax_);

    if (!equilibrium_)
    {
        const surfaceScalarField& alphaRhoPhi = alphaRhoPhi_;
        const volVectorField& Uc_ = otherPhase_.U();
        const volVectorField& U = U_;
        // refCast<const twoPhaseSystem>(phase_.fluid()).otherPhase(phase_).U();

        // Particle viscosity (Table 3.2, p.47)
        nut_ = viscosityModel_->nu(alpha, Theta_, gs0_, rho, da, e_);

        volScalarField ThetaSqrt("sqrtTheta", sqrt(Theta_));

        // Bulk viscosity  p. 45 (Lun et al. 1984).
        lambda_ = (4.0/3.0)*sqr(alpha)*da*gs0_*(1 + e_)*ThetaSqrt/sqrtPi;

        // Stress tensor, Definitions, Table 3.1, p. 43
        volSymmTensorField tau
        (
            rho*(2*nut_*D + (lambda_ - (2.0/3.0)*nut_)*tr(D)*I)
        );

        // Dissipation (Eq. 3.24, p.50)
        volScalarField gammaCoeff
        (
            "gammaCoeff",
            12*(1 - sqr(e_))
           *max(sqr(alpha), residualAlpha_)
           *rho*gs0_*(1.0/da)*ThetaSqrt/sqrtPi
        );

        // Drag
        volScalarField beta
        (
            this->db().lookupObject<volScalarField>("Kd.particlesInGas")
        );

        // Eq. 3.25, p. 50 Js = J1 - J2
        volScalarField J1("J1", 3*beta);
        volScalarField J2
        (
            "J2",
            0.25*sqr(beta)*da*magSqr(U - Uc_)
           /(
               max(alpha, residualAlpha_)*rho
              *sqrtPi*ThetaSqrt
            )
        );

        // particle pressure - coefficient in front of Theta (Eq. 3.22, p. 45)
        volScalarField PsCoeff
        (
            granularPressureModel_->granularPressureCoeff
            (
                alpha,
                gs0_,
                rho,
                e_
            )
        );

        // 'thermal' conductivity (Table 3.3, p. 49)
        kappa_ = conductivityModel_->kappa(alpha, Theta_, gs0_, rho, da, e_);

        fv::options& fvOptions(fv::options::New(mesh_));

        // Construct the granular temperature equation (Eq. 3.20, p. 44)
        // NB. note that there are two typos in Eq. 3.20:
        //     Ps should be without grad
        //     the laplacian has the wrong sign
        fvScalarMatrix ThetaEqn
        (
            1.5*
            (
                fvm::ddt(alpha, rho, Theta_)
              + fvm::div(alphaRhoPhi, Theta_)
              - fvc::Sp(fvc::ddt(alpha, rho) + fvc::div(alphaRhoPhi), Theta_)
            )
          - fvm::laplacian(kappa_, Theta_, "laplacian(kappa,Theta)")
         ==
          - fvm::SuSp((PsCoeff*I) && gradU, Theta_)
          + (tau && gradU)
          + fvm::Sp(-gammaCoeff, Theta_)
          + fvm::Sp(-J1, Theta_)
          + fvm::Sp(J2/Theta_, Theta_)
          + fvOptions(alpha, rho, Theta_)
        );

        // Fix empty cell values
        label emptyCellSize = 0;
        forAll(alpha, i)
        {
          emptyCellSize = alpha[i] <= cutoff_ ? emptyCellSize + 1 : emptyCellSize;
        }
        Info << "EmptyCellSize: " << emptyCellSize << endl;
        labelList emptyCells(emptyCellSize);
        label j = 0;
        forAll(alpha, i)
        {
          if(alpha[i] <= cutoff_)
          {
            emptyCells[j] = i;
            j++;
          }
        }
        scalarField Theta0(emptyCellSize, 1e-15);
        ThetaEqn.setValues(emptyCells, Theta0);

        ThetaEqn.relax();
        // fvOptions.constrain(ThetaEqn);

        ThetaEqn.solve();

        // Impose wall on propellant surface
        // const volScalarField& alphaProp(this->db().lookupObject<volScalarField>("alpha.propellant"));
        // ImposeWall(Theta_, alphaProp);
        // Theta_ = pos0(alpha - minAlpha_)*Theta_;
        fvOptions.correct(Theta_);
    }
    else
    {
        // Equilibrium => dissipation == production
        // Eq. 4.14, p.82
        volScalarField K1("K1", 2*(1 + e_)*rho*gs0_);
        volScalarField K3
        (
            "K3",
            0.5*da*rho*
            (
                (sqrtPi/(3*(3.0 - e_)))
               *(1 + 0.4*(1 + e_)*(3*e_ - 1)*alpha*gs0_)
               +1.6*alpha*gs0_*(1 + e_)/sqrtPi
            )
        );

        volScalarField K2
        (
            "K2",
            4*da*rho*(1 + e_)*alpha*gs0_/(3*sqrtPi) - 2*K3/3.0
        );

        volScalarField K4("K4", 12*(1 - sqr(e_))*rho*gs0_/(da*sqrtPi));

        volScalarField trD
        (
            "trD",
            alpha/(alpha + residualAlpha_)
           *fvc::div(phi_)
        );
        volScalarField tr2D("tr2D", sqr(trD));
        volScalarField trD2("trD2", tr(D & D));

        volScalarField t1("t1", K1*alpha + rho);
        volScalarField l1("l1", -t1*trD);
        volScalarField l2("l2", sqr(t1)*tr2D);
        volScalarField l3
        (
            "l3",
            4.0
           *K4
           *alpha
           *(2*K3*trD2 + K2*tr2D)
        );

        Theta_ = sqr
        (
            (l1 + sqrt(l2 + l3))
           /(2*max(alpha*K4, residualAlpha_*dimensionedScalar("", K4.dimensions(), 1.0)))
        );

        kappa_ = conductivityModel_->kappa(alpha, Theta_, gs0_, rho, da, e_);
    }

    Theta_.max(ThetaSmall.value());
    Theta_.min(100);
    {
        // particle viscosity (Table 3.2, p.47)
        nut_ = viscosityModel_->nu(alpha, Theta_, gs0_, rho, da, e_);

        volScalarField ThetaSqrt("sqrtTheta", sqrt(Theta_));

        // Bulk viscosity  p. 45 (Lun et al. 1984).
        lambda_ = (4.0/3.0)*sqr(alpha)*da*gs0_*(1 + e_)*ThetaSqrt/sqrtPi;

        // Frictional pressure
        volScalarField pf
        (
            frictionalStressModel_->frictionalPressure
            (
                phase_,
                alphaMinFriction_,
                alphaMax_
            )
        );

        nuFric_ = frictionalStressModel_->nu
        (
            phase_,
            alphaMinFriction_,
            alphaMax_,
            pf/rho,
            D
        );

        // Limit viscosity and add frictional viscosity
        nut_.min(maxNut_);
        nuFric_ = min(nuFric_, maxNut_ - nut_);
        nut_ += nuFric_;
    }

    // if (debug)
    // {
        Info<< typeName << ':' << nl
            << "    max(Theta) = " << max(Theta_).value() << nl
            << "    max(nut) = " << max(nut_).value() << endl;
    // }
}


// ************************************************************************* //
