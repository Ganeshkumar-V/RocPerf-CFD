/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "burningRateDy.H"
#include "phasePair.H"
#include "phaseSystem.H"
#include "addToRunTimeSelectionTable.H"
#include "processorFvPatch.H"
#include "cutCellIso.H"
#include "volPointInterpolation.H"
#include "wallPolyPatch.H"
#include "fvcSmooth.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace interfaceTrackingModels
{
    defineTypeNameAndDebug(burningRateDy, 0);
    addToRunTimeSelectionTable(interfaceTrackingModel, burningRateDy, dictionary);
}
}
// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

void Foam::interfaceTrackingModels::burningRateDy::updateInterface()
{
    const volScalarField& alpha = propellantPhase_;
    interfaceArea_ = dimensionedScalar(dimless/dimLength, 0)*alpha;
    interfaceAreaVector_ = vector(0, 0, 0)*alpha;

    scalarField ap
    (
        volPointInterpolation::New(mesh_).interpolate(alpha)
    );

    cutCellIso cutCell(mesh_, ap);

    forAll(interfaceArea_, celli)
    {
        label status = cutCell.calcSubCell(celli, isoAlpha_);
        if (status == 0) // cell is cut
        {
            interfaceArea_[celli] = pos(cutCell.faceArea().x())
                          *cutCell.faceArea().x()/mesh_.V()[celli];
            interfaceAreaVector_[celli].x() = 1.0;
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceTrackingModels::burningRateDy::burningRateDy
(
    const dictionary& dict,
    const phasePair& pair
)
:
    interfaceTrackingModel(dict, pair),
    mesh_(pair.phase1().mesh()),
    n(dict.get<scalar>("n")),
    a
    (
      "",
      pow(dimLength*dimTime*dimTime/dimMass, n)*dimVelocity,
      dict.get<scalar>("a")
    ),
    propellantPhase_
    (
      pair.phase1().db().lookupObject<phaseModel>("alpha." + this->propellant_)
    ),
    interfaceArea_
    (
      IOobject
      (
        "interfaceArea",
        this->mesh_.time().timeName(),
        this->mesh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
      ),
      this->mesh_,
      dimensionedScalar(dimless/dimLength, Zero)
    ),
    interfaceAreaVector_
    (
      IOobject
      (
        "interfaceArea",
        this->mesh_.time().timeName(),
        this->mesh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
      ),
      this->mesh_,
      dimensionedVector(dimless, vector(0, 0, 0))
    ),
    isoAlpha_(dict.getOrDefault<scalar>("isoAlpha", 0.5)),
    rb_
    (
      volScalarField
      (
        IOobject("rb", pair_.phase1().mesh()),
        pair_.phase1().mesh(),
        dimensionedScalar("", dimVelocity, 0.0)
      )
    ),
    crb_("", dimVelocity, dict.getOrDefault<scalar>("rb", -1))
{
  Info << "Propellant Phase is " << propellantPhase_.name()
            << " index - " << propellantPhase_.index() << endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interfaceTrackingModels::burningRateDy::~burningRateDy()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::interfaceTrackingModels::burningRateDy::correct()
{

  // rb = constant
  if (crb_.value() != -1)
  {
    rb_.ref() = crb_;
  }
  else
  {
    // rb = aP^n
    const phaseModel& phase = pair_.phase1();
    const volScalarField& p = phase.db().lookupObject<volScalarField>("p");
    rb_.ref() = (a*pow(p()/1e6, n))*1e-2;
  }

  // Find Interface
  updateInterface();
}

// find interface and stores interface in owner cell
void Foam::interfaceTrackingModels::burningRateDy::regress
(
  volScalarField& alpha
)
{}

void Foam::interfaceTrackingModels::burningRateDy::findInterface
(
  const volScalarField& alpha
)
{}

Foam::tmp<Foam::volScalarField>
Foam::interfaceTrackingModels::burningRateDy::rb() const
{
    return Foam::tmp<Foam::volScalarField>(new volScalarField("trb", rb_));
}

Foam::tmp<Foam::volScalarField>
Foam::interfaceTrackingModels::burningRateDy::As() const
{
    return Foam::tmp<Foam::volScalarField>(new volScalarField("tAs", interfaceArea_));
}

Foam::tmp<Foam::volVectorField>
Foam::interfaceTrackingModels::burningRateDy::nHat() const
{
    return Foam::tmp<Foam::volVectorField>(new volVectorField("nHat", interfaceAreaVector_));
}
// ************************************************************************* //
