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

#include "StaticPhaseModel.H"

// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

template<class BasePhaseModel>
template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::StaticPhaseModel<BasePhaseModel>::zeroField
(
    const word& name,
    const dimensionSet& dims,
    const bool cache
) const
{
    return tmp<GeometricField<Type, PatchField, GeoMesh>>
    (
        new GeometricField<Type, PatchField, GeoMesh>
        (
            IOobject
            (
                IOobject::groupName(name, this->name()),
                this->mesh().time().timeName(),
                this->mesh()
            ),
            this->mesh(),
            dimensioned<Type>(dims, Zero)
        )
    );
}


template<class BasePhaseModel>
template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::StaticPhaseModel<BasePhaseModel>::zeroVolField
(
    const word& name,
    const dimensionSet& dims,
    const bool cache
) const
{
    return zeroField<Type, fvPatchField, volMesh>(name, dims, cache);
}


template<class BasePhaseModel>
template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::StaticPhaseModel<BasePhaseModel>::zeroSurfaceField
(
    const word& name,
    const dimensionSet& dims,
    const bool cache
) const
{
    return zeroField<Type, fvsPatchField, surfaceMesh>(name, dims, cache);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::StaticPhaseModel<BasePhaseModel>::StaticPhaseModel
(
    const phaseSystem& fluid,
    const word& phaseName,
    const label index
)
:
    BasePhaseModel(fluid, phaseName, index),
    U_
    (
        IOobject
        (
            IOobject::groupName("U", this->name()),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        fluid.mesh(),
        dimensionedVector(dimVelocity, vector(0, 0, 0))
    ),
    phi_
    (
        IOobject
        (
            IOobject::groupName("phi", phaseModel::name()),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        fluid.mesh(),
        dimensionedScalar(dimensionSet(0, 3, -1, 0, 0), Zero)
    ),
    alphaPhi_
    (
        IOobject
        (
            IOobject::groupName("alphaPhi", this->name()),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        fluid.mesh(),
        dimensionedScalar(dimensionSet(0, 3, -1, 0, 0), Zero)
    ),
    alphaRhoPhi_
    (
        IOobject
        (
            IOobject::groupName("alphaRhoPhi", this->name()),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        fluid.mesh(),
        dimensionedScalar(dimensionSet(1, 0, -1, 0, 0), Zero)
    ),
    divU_(nullptr)
{
  phi_.setOriented();
  alphaPhi_.setOriented();
  alphaRhoPhi_.setOriented();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::StaticPhaseModel<BasePhaseModel>::~StaticPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasePhaseModel>
bool Foam::StaticPhaseModel<BasePhaseModel>::stationary() const
{
    return true;
}

template<class BasePhaseModel>
bool Foam::StaticPhaseModel<BasePhaseModel>::moving() const
{
    return false;
}


template<class BasePhaseModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::StaticPhaseModel<BasePhaseModel>::UEqn()
{
    FatalErrorInFunction
        << "Cannot construct a momentum equation for a static phase"
        << exit(FatalError);

    return tmp<fvVectorMatrix>();
}


template<class BasePhaseModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::StaticPhaseModel<BasePhaseModel>::UfEqn()
{
    FatalErrorInFunction
        << "Cannot construct a momentum equation for a static phase"
        << exit(FatalError);

    return tmp<fvVectorMatrix>();
}


template<class BasePhaseModel>
Foam::tmp<Foam::volVectorField>
Foam::StaticPhaseModel<BasePhaseModel>::U() const
{
    return tmp<volVectorField>(U_);
}


template<class BasePhaseModel>
Foam::volVectorField&
Foam::StaticPhaseModel<BasePhaseModel>::URef()
{
    FatalErrorInFunction
        << "Cannot access the velocity of a static phase"
        << exit(FatalError);

    return const_cast<volVectorField&>(volVectorField::null());
}


template<class BasePhaseModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::StaticPhaseModel<BasePhaseModel>::phi() const
{
    return tmp<surfaceScalarField>(phi_);
}


template<class BasePhaseModel>
Foam::surfaceScalarField&
Foam::StaticPhaseModel<BasePhaseModel>::phiRef()
{
    FatalErrorInFunction
        << "Cannot access the flux of a static phase"
        << exit(FatalError);

    return const_cast<surfaceScalarField&>(surfaceScalarField::null());
}


template<class BasePhaseModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::StaticPhaseModel<BasePhaseModel>::alphaPhi() const
{
    return tmp<surfaceScalarField>(alphaPhi_);
}


template<class BasePhaseModel>
Foam::surfaceScalarField&
Foam::StaticPhaseModel<BasePhaseModel>::alphaPhiRef()
{
    FatalErrorInFunction
        << "Cannot access the volumetric flux of a static phase"
        << exit(FatalError);

    return const_cast<surfaceScalarField&>(surfaceScalarField::null());
}


template<class BasePhaseModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::StaticPhaseModel<BasePhaseModel>::alphaRhoPhi() const
{
    return tmp<surfaceScalarField>(alphaRhoPhi_);
}


template<class BasePhaseModel>
Foam::surfaceScalarField&
Foam::StaticPhaseModel<BasePhaseModel>::alphaRhoPhiRef()
{
    FatalErrorInFunction
        << "Cannot access the mass flux of a static phase"
        << exit(FatalError);

    return const_cast<surfaceScalarField&>(surfaceScalarField::null());
}


template<class BasePhaseModel>
Foam::tmp<Foam::volVectorField>
Foam::StaticPhaseModel<BasePhaseModel>::DUDt() const
{
    return zeroVolField<vector>("DUDt", dimVelocity/dimTime);
}


template<class BasePhaseModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::StaticPhaseModel<BasePhaseModel>::DUDtf() const
{
    return zeroSurfaceField<scalar>("DUDtf", dimVelocity*dimArea/dimTime);
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::StaticPhaseModel<BasePhaseModel>::continuityError() const
{
    return zeroVolField<scalar>("continuityError", dimDensity/dimTime);
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::StaticPhaseModel<BasePhaseModel>::continuityErrorFlow() const
{
    return zeroVolField<scalar>("continuityErrorFlow", dimDensity/dimTime);
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::StaticPhaseModel<BasePhaseModel>::continuityErrorSources() const
{
    return zeroVolField<scalar>("continuityErrorSources", dimDensity/dimTime);
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::StaticPhaseModel<BasePhaseModel>::K() const
{
    return zeroVolField<scalar>("K", sqr(dimVelocity));
}

template<class BasePhaseModel>
Foam::tmp<Foam::volSymmTensorField>
Foam::StaticPhaseModel<BasePhaseModel>::devRhoReff() const
{
    Info << "devRhoReff() NotImplemented for stationary phase !" << exit(FatalError);
    return tmp<volSymmTensorField>();
}

template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::StaticPhaseModel<BasePhaseModel>::divU() const
{
    return divU_.valid() ? tmp<volScalarField>(divU_()) : tmp<volScalarField>();
}


template<class BasePhaseModel>
void Foam::StaticPhaseModel<BasePhaseModel>::divU
(
    tmp<volScalarField> divU
)
{
    divU_ = divU;
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::StaticPhaseModel<BasePhaseModel>::mut() const
{
    return zeroVolField<scalar>("continuityError", dimDynamicViscosity);
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::StaticPhaseModel<BasePhaseModel>::muEff() const
{
    return this->thermo().mu();
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::StaticPhaseModel<BasePhaseModel>::nut() const
{
    return zeroVolField<scalar>("continuityError", dimViscosity);
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::StaticPhaseModel<BasePhaseModel>::nuEff() const
{
    return this->thermo().nu();
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::StaticPhaseModel<BasePhaseModel>::kappaEff() const
{
    return this->thermo().kappa();
}


template<class BasePhaseModel>
Foam::tmp<Foam::scalarField>
Foam::StaticPhaseModel<BasePhaseModel>::kappaEff(const label patchi) const
{
    return this->thermo().kappa(patchi);
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::StaticPhaseModel<BasePhaseModel>::alphaEff() const
{
    return this->thermo().alpha();
}


template<class BasePhaseModel>
Foam::tmp<Foam::scalarField>
Foam::StaticPhaseModel<BasePhaseModel>::alphaEff(const label patchi) const
{
    return this->thermo().alpha(patchi);
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::StaticPhaseModel<BasePhaseModel>::k() const
{
    return zeroVolField<scalar>("k", sqr(dimVelocity));
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::StaticPhaseModel<BasePhaseModel>::pPrime() const
{
    return zeroVolField<scalar>("pPrime", dimPressure);
}


// ************************************************************************* //
