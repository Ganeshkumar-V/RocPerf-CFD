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

#include "subCellularInterfaceMotion.H"
#include "phasePair.H"
#include "phaseSystem.H"
#include "addToRunTimeSelectionTable.H"
#include "processorFvPatch.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace interfaceTrackingModels
{
    defineTypeNameAndDebug(subCellularInterfaceMotion, 0);
    addToRunTimeSelectionTable(interfaceTrackingModel, subCellularInterfaceMotion, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceTrackingModels::subCellularInterfaceMotion::subCellularInterfaceMotion
(
    const dictionary& dict,
    const phasePair& pair
)
:
    interfaceTrackingModel(dict, pair),
    n(dict.get<scalar>("n")),
    a("", pow(dimLength*dimTime*dimTime/dimMass, n)*dimVelocity, dict.get<scalar>("a")),
    interface_
    (
      volScalarField
      (
        IOobject("interface", pair_.phase1().mesh()),
        pair_.phase1().mesh(),
        dimensionedScalar("", dimless, 0.0)
      )
    ),
    rb_
    (
      volScalarField
      (
        IOobject("rb", pair_.phase1().mesh()),
        pair_.phase1().mesh(),
        dimensionedScalar("", dimVelocity, 0.0)
      )
    ),
    As_
    (
      volScalarField
      (
        IOobject("As", pair_.phase1().mesh()),
        pair_.phase1().mesh(),
        dimensionedScalar("", dimArea/dimVolume, 0.0)
      )
    ),
    crb_("", dimVelocity, dict.getOrDefault<scalar>("rb", -1))
{
  const phaseModel& phase = pair_.phase1();
  const volScalarField& alpha
        = phase.db().lookupObject<volScalarField>("alpha." + propellant_);
  this->findInterface(alpha);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interfaceTrackingModels::subCellularInterfaceMotion::~subCellularInterfaceMotion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::interfaceTrackingModels::subCellularInterfaceMotion::correct()
{

  // rb = constant
  if (crb_.value() != -1)
  {
    rb_.ref() = interface_.ref()*crb_;
  }
  else
  {
    // rb = aP^n
    const phaseModel& phase = pair_.phase1();
    const volScalarField& p = phase.db().lookupObject<volScalarField>("p");
    rb_.ref() = interface_.ref()*(a*pow(p()/1e6, n))*1e-2;
  }

}

void Foam::interfaceTrackingModels::subCellularInterfaceMotion::regress
(
  volScalarField& alpha
)
{
  const volScalarField& alpha0 = alpha.oldTime();

  const fvMesh& mesh = alpha.mesh();
  const labelList& Own = mesh.owner();
  const labelList& Nei = mesh.neighbour();
  const surfaceScalarField& Sf = mesh.magSf();
  const scalar dt = mesh.time().deltaTValue();
  const scalarField& V = mesh.V();

  forAll(Own, i)
  {
    // case:1 Interface is present at the center of the owner cell
    if (alpha0[Own[i]] == 0.5 && alpha0[Nei[i]] == 1)
    {
      interface_[Own[i]] = 1;
      As_[Own[i]] = Sf[i]/V[Own[i]];  // Area of face between owner and neighbour
      alpha[Own[i]] = alpha0[Own[i]] - rb_[Own[i]]*Sf[i]*dt/V[Own[i]];
      if (alpha[Own[i]] < 0)
      {
        scalar Vr = -alpha[Own[i]]*V[Own[i]];
        alpha[Nei[i]] = alpha0[Nei[i]] - Vr/V[Nei[i]];
        alpha[Own[i]] = 0;
        if (alpha[Nei[i]] < 0)
        {
          FatalErrorInFunction
            << "Regression is very fast!\n"
            << "Hint: Reduce time step."
            << exit(FatalError);
        }
      }
    }
    // case:2 Interface is present in between owner center and face
    else if ((alpha0[Own[i]] < 0.5) && (alpha0[Nei[i]] == 1))
    {
      interface_[Own[i]] = 1;
      As_[Own[i]] = Sf[i]/V[Own[i]];  // Area of face between owner and neighbour
      alpha[Own[i]] = alpha0[Own[i]] - rb_[Own[i]]*Sf[i]*dt/V[Own[i]];
      if (alpha[Own[i]] < 0)
      {
        scalar Vr = -alpha[Own[i]]*V[Own[i]];
        alpha[Nei[i]] = alpha0[Nei[i]] - Vr/V[Nei[i]];
        alpha[Own[i]] = 0;
        if (alpha[Nei[i]] < 0)
        {
          FatalErrorInFunction
            << "Regression is very fast!\n"
            << "Hint: Reduce time step."
            << exit(FatalError);
        }
      }
    }
    // case:3 Interface is present exactly at the center
    else if ((alpha0[Own[i]] == 0) && (alpha0[Nei[i]] == 1))
    {
      interface_[Own[i]] = 1;
      As_[Own[i]] = Sf[i]/V[Nei[i]];  // Area of face between owner and neighbour
      alpha[Nei[i]] = alpha0[Nei[i]] - rb_[Own[i]]*Sf[i]*dt/V[Nei[i]];
      if (alpha[Nei[i]] < 0)
      {
        FatalErrorInFunction
          << "Regression is very fast!\n"
          << "Hint: Reduce time step."
          << exit(FatalError);
      }
    }
    // case:4 Interface is present in between face and neighbour center
    else if ((alpha0[Own[i]] == 0) && (alpha0[Nei[i]] > 0.5))
    {
      interface_[Own[i]] = 1;
      As_[Own[i]] = Sf[i]/V[Nei[i]];  // Area of face between owner and neighbour
      alpha[Nei[i]] = alpha0[Nei[i]] - rb_[Own[i]]*Sf[i]*dt/V[Nei[i]];
      if (alpha[Nei[i]] < 0)
      {
        FatalErrorInFunction
          << "Regression is very fast!\n"
          << "Hint: Reduce time step."
          << exit(FatalError);
      }
    }
    // case:5 Interface is not present (Ignore)
    else continue;
  }
}

void Foam::interfaceTrackingModels::subCellularInterfaceMotion::findInterface
(
  const volScalarField& alpha
)
{
  // -If found interface -> interface_ = 1
  // use owner neighbour approach

  const fvMesh& mesh = alpha.mesh();
  const labelList& Own = mesh.owner();
  const labelList& Nei = mesh.neighbour();

  forAll(Own, i)
  {
    // case:1 Interface is present at the center of the owner cell
    if (alpha[Own[i]] == 0.5 && alpha[Nei[i]] == 1)
    {
      interface_[Own[i]] = 1;
    }
    // case:2 Interface is present in between owner center and face
    else if ((alpha[Own[i]] < 0.5) && (alpha[Nei[i]] == 1))
    {
      interface_[Own[i]] = 1;
    }
    // case:3 Interface is present exactly at the center
    else if ((alpha[Own[i]] == 0) && (alpha[Nei[i]] == 1))
    {
      interface_[Own[i]] = 1;
    }
    // case:4 Interface is present in between face and neighbour center
    else if ((alpha[Own[i]] == 0) && (alpha[Nei[i]] > 0.5))
    {
      interface_[Own[i]] = 1;
    }
    // case:5 Interface is not present (Ignore)
    else continue;
  }

}

Foam::tmp<Foam::volScalarField>
Foam::interfaceTrackingModels::subCellularInterfaceMotion::rb() const
{
    return Foam::tmp<Foam::volScalarField>(new volScalarField("trb", rb_));
}

Foam::tmp<Foam::volScalarField>
Foam::interfaceTrackingModels::subCellularInterfaceMotion::As() const
{
    return Foam::tmp<Foam::volScalarField>(new volScalarField("tAs", As_));
}
// ************************************************************************* //
