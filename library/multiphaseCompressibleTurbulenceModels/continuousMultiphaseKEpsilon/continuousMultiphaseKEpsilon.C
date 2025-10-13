/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2017 OpenFOAM Foundation
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

#include "continuousMultiphaseKEpsilon.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
continuousMultiphaseKEpsilon<BasicTurbulenceModel>::continuousMultiphaseKEpsilon
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    kEpsilon<BasicTurbulenceModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName,
        type
    ),
    Kpg_
    (
        IOobject
        (
            IOobject::groupName("Kpg", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("Kpg", dimVelocity*dimVelocity, 0.0)
    ),
    particlePhaseName_
    (
      this->coeffDict_.template get<word>("particlePhase")
    ),
    CEpsilon3_("CEpsilon3", dimless, 1.2)
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<fvScalarMatrix>
continuousMultiphaseKEpsilon<BasicTurbulenceModel>::kSource() const
{
    const phaseModel& particlePhase_(this->db().template lookupObject<phaseModel>("alpha." + particlePhaseName_));

    const volScalarField& alpha = particlePhase_;

    // Drag Force
    volScalarField rhoPbeta 
    (
        this->db().template lookupObject<volScalarField>("Kd.particlesInGas")
    );

    return - fvm::Sp(2.0*alpha*rhoPbeta, this->k_) + alpha*rhoPbeta*Kpg_;
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix>
continuousMultiphaseKEpsilon<BasicTurbulenceModel>::epsilonSource() const
{
    const phaseModel& particlePhase_(this->db().template lookupObject<phaseModel>("alpha."+particlePhaseName_));

    const volScalarField& alpha = particlePhase_;

    // Drag Force
    volScalarField rhoPbeta
    (
        this->db().template lookupObject<volScalarField>("Kd.particlesInGas")
    );

    return fvm::Sp(CEpsilon3_*alpha*rhoPbeta*(Kpg_/max(this->k_, dimensionedScalar("", dimVelocity*dimVelocity, SMALL)) - 2.0), this->epsilon_);
}

template<class BasicTurbulenceModel>
void continuousMultiphaseKEpsilon<BasicTurbulenceModel>::correctKpg()
{
    const phaseModel& particlePhase_ = this->db().template lookupObject<phaseModel>("alpha."+particlePhaseName_);

    const tmp<volScalarField> tdp(particlePhase_.d());
    const volScalarField& dp(tdp());

    const tmp<volScalarField> trhoP(particlePhase_.rho());
    const volScalarField& rhoP(trhoP());

    // Drag Force
    volScalarField rhoPbeta
    (
        this->db().template lookupObject<volScalarField>("Kd.particlesInGas")
    );

    // Relative velocity
    volVectorField Ur
    (
        this->db().template lookupObject<volVectorField>("Ur.particlesInGas")
    );

    // Fluctuating kinetic energy of particle phase
    volScalarField Theta
    (
        this->db().template lookupObject<volScalarField>("Theta." + particlePhase_.name())
    );

    Kpg_ = (rhoPbeta*dp*sqr(mag(Ur))/rhoP)/(4.0*sqrt(constant::mathematical::pi*Theta)); // Koch Relation
    // Kpg_ = Csf/sqrt(this->k_*Theta); // Isotropic turbulence
}

template<class BasicTurbulenceModel>
void continuousMultiphaseKEpsilon<BasicTurbulenceModel>::correct()
{
    // Correct gas-particle velocity correlation
    correctKpg();

    // Solve the turbulence equations
    kEpsilon<BasicTurbulenceModel>::correct();
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
