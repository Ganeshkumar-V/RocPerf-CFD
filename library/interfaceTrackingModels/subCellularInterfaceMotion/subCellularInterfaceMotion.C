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
    f(dict.get<scalar>("f")),
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
    dmdt_
    (
      volScalarField
      (
        IOobject("dmdt", pair_.phase1().mesh()),
        pair_.phase1().mesh(),
        dimensionedScalar("", dimVelocity/dimLength, 0.0)
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

Foam::scalar Foam::interfaceTrackingModels::subCellularInterfaceMotion::rb(scalar P)
{
  if (P > 1e6)
    return (a*(pow(P/1e6, n))).value()*1e-2;
  else
    return (a*(pow(P/1e6, n*f))).value()*1e-2;
}

// void Foam::interfaceTrackingModels::subCellularInterfaceMotion::regress
// (
//   volScalarField& alpha,
//   const volScalarField& alphaOld
// )
// {
//   const phaseModel& phase = pair_.phase1();
//   const volScalarField& p = phase.db().lookupObject<volScalarField>("p");
//
//   const volScalarField& alpha0 = alphaOld;
//
//   const fvMesh& mesh = alpha.mesh();
//   const labelList& Own = mesh.owner();
//   const labelList& Nei = mesh.neighbour();
//   const surfaceScalarField& Sf = mesh.magSf();
//   const scalar dt = mesh.time().deltaTValue();
//   const scalarField& V = mesh.V();
//   const scalar One(1 - SMALL);
//   const scalar Zero(SMALL);
//   interface_ = dimensionedScalar(dimless, 0.0);
//   As_ = dimensionedScalar(As_.dimensions(), 0.0);
//   dmdt_ = dimensionedScalar(dmdt_.dimensions(), 0.0);
//   rb_ = dimensionedScalar(rb_.dimensions(), 0.0);
//
//   // Internal Cells
//   forAll(Own, i)
//   {
//     // case:1 Interface is present at the center of the owner cell (or)
//     // case:2 Interface is present in between owner center and face
//     if (alpha0[Own[i]] <= 0.5 && alpha0[Nei[i]] == One)
//     {
//       interface_[Own[i]] = 1;
//       As_[Own[i]] = Sf[i]/V[Own[i]];  // Area of face between owner and neighbour
//       // rb_[Own[i]] = (a*pow(p[Own[i]]/1e6, n)).value()*1e-2;  // burning Rate
//       rb_[Own[i]] = rb(p[Own[i]]);  // burning Rate
//       dmdt_[Own[i]] = rb_[Own[i]]*As_[Own[i]];
//       alpha[Own[i]] = alpha0[Own[i]] - rb_[Own[i]]*As_[Own[i]]*dt;
//       if (alpha[Own[i]] < 0)
//       {
//           scalar Vr = -alpha[Own[i]]*V[Own[i]];
//           alpha[Nei[i]] = alpha0[Nei[i]] - Vr/V[Nei[i]];
//           alpha[Own[i]] = 0;
//           if (alpha[Nei[i]] < 0)
//           {
//               FatalErrorInFunction
//                 << "Regression is very fast!\n"
//                 << "Hint: Reduce time step."
//                 << exit(FatalError);
//           }
//       }
//     }
//     // case:3 Interface is present exactly at the face
//     // case:4 Interface is present in between face and neighbour center
//     else if ((alpha0[Own[i]] == Zero) && (alpha0[Nei[i]] >= 0.5))
//     {
//       interface_[Nei[i]] = 1;
//       As_[Nei[i]] = Sf[i]/V[Nei[i]];  // Area of face between owner and neighbour
//       // rb_[Nei[i]] = (a*pow(p[Nei[i]]/1e6, n)).value()*1e-2;  // burning Rate
//       rb_[Nei[i]] = rb(p[Nei[i]]);  // burning Rate
//
//       // Source Term Distribution - Technique 1
//       dmdt_[Own[i]] = 2*(alpha0[Nei[i]] - 0.5)*(rb_[Nei[i]]*As_[Nei[i]]);
//       dmdt_[Nei[i]] = (1.0 - 2*(alpha0[Nei[i]] - 0.5))*rb_[Nei[i]]*As_[Nei[i]];
//       As_[Own[i]] = As_[Nei[i]];
//       rb_[Own[i]] = rb_[Nei[i]];
//
//       // Source Term Distribution - Technique 2
//       // dmdt_[Own[i]] = alpha0[Nei[i]]*(rb_[Nei[i]]*As_[Nei[i]]);
//       // dmdt_[Nei[i]] = (1.0 -alpha0[Nei[i]])*rb_[Nei[i]]*As_[Nei[i]];
//       // As_[Own[i]] = As_[Nei[i]];
//       // rb_[Own[i]] = rb_[Nei[i]];
//
//       // Source Term Distribution - Technique 3
//       // dmdt_[Own[i]] = 0.2*(rb_[Nei[i]]*As_[Nei[i]]);
//       // dmdt_[Nei[i]] = 0.8*rb_[Nei[i]]*As_[Nei[i]];
//       // As_[Own[i]] = As_[Nei[i]];
//       // rb_[Own[i]] = rb_[Nei[i]];
//
//       // No source Term Distribution (put Everything in neighbour)
//       // dmdt_[Nei[i]] = rb_[Nei[i]]*As_[Nei[i]];
//
//       alpha[Nei[i]] = alpha0[Nei[i]] - rb_[Nei[i]]*As_[Nei[i]]*dt;
//       if (alpha[Nei[i]] < 0)
//       {
//         FatalErrorInFunction
//           << "Regression is very fast!\n"
//           << "Hint: Reduce time step."
//           << exit(FatalError);
//       }
//     }
//     // case:5 Interface is not present (Ignore)
//     else continue;
//   }
//
  // Boundary Patches
  // forAll(mesh.boundary(), patchi)
  // {
  //   const fvPatch& patch = mesh.boundary()[patchi];
  //   const labelList& fC = patch.faceCells();
  //   const scalarField pSf(patch.magSf());
  //   forAll(fC, celli)
  //   {
  //     if
  //     (
  //         (alpha0[fC[celli]] <= 0.5) &&
  //         (alpha0[fC[celli]] > Zero) &&
  //         (interface_[fC[celli]] == 0)
  //     )
  //     {
  //       interface_[fC[celli]] = 1.0;
  //       As_[fC[celli]] = pSf[celli]/V[fC[celli]];
  //       rb_[fC[celli]] = (a*pow(p[fC[celli]]/1e6, n)).value()*1e-2;  // burning Rate
  //       dmdt_[fC[celli]] = rb_[fC[celli]]*As_[fC[celli]];
  //
  //       alpha[fC[celli]] = alpha0[fC[celli]] - rb_[fC[celli]]*As_[fC[celli]]*dt;
  //       if (alpha[fC[celli]] < 0)
  //       {
  //         alpha[fC[celli]] = Zero;
  //         As_[fC[celli]] = (alpha0[fC[celli]] - alpha[fC[celli]])
  //                           /(rb_[fC[celli]]*dt);
  //       }
  //     }
  //     else continue;
  //   }
  // }
// }

Foam::label Foam::interfaceTrackingModels::subCellularInterfaceMotion::findNeighbour
(
  const volScalarField& alpha,
  scalar NEI
)
{
  const fvMesh& mesh = alpha.mesh();
  const labelList& Own = mesh.owner();
  const labelList& Nei = mesh.neighbour();
  const scalar One(1 - SMALL);

  forAll(Own, i)
  {
    if (Own[i] == NEI)
    {
      Info << "Nei: " << Nei[i] << " alpha: " << alpha[Nei[i]] << endl;
      if (alpha[Nei[i]] == One)
      {
        return Nei[i];
      }
    }
  }
  return -1;
}


void Foam::interfaceTrackingModels::subCellularInterfaceMotion::regress
(
  volScalarField& alpha,
  const volScalarField& alphaOld
)
{
  const phaseModel& phase = pair_.phase1();
  const volScalarField& p = phase.db().lookupObject<volScalarField>("p");

  const volScalarField& alpha0 = alphaOld;

  const fvMesh& mesh = alpha.mesh();
  const labelList& Own = mesh.owner();
  const labelList& Nei = mesh.neighbour();
  const surfaceScalarField& Sf = mesh.magSf();
  const scalar dt = mesh.time().deltaTValue();
  const scalarField& V = mesh.V();
  // const scalar One(1 - SMALL);
  const scalar Zero(SMALL);

  interface_ = dimensionedScalar(dimless, 0.0);
  As_ = dimensionedScalar(As_.dimensions(), 0.0);
  dmdt_ = dimensionedScalar(dmdt_.dimensions(), 0.0);
  rb_ = dimensionedScalar(rb_.dimensions(), 0.0);

  // Internal Cells
  forAll(Own, i)
  {
    // case:1 Interface is present in the Neighbour Cell
    if
    (
      (alpha0[Own[i]] == Zero && alpha0[Nei[i]] > Zero) && 
      (Nei[i] == Own[i] + 1)
    )
    {
      interface_[Nei[i]] = 1;

      As_[Nei[i]] = Sf[i]/V[Nei[i]];  // Area of face between owner and neighbour
      rb_[Nei[i]] = rb(p[Nei[i]]);  // burning Rate
      dmdt_[Nei[i]] = (1 - alpha[Nei[i]])*rb_[Nei[i]]*As_[Nei[i]];

      As_[Own[i]] = As_[Nei[i]];
      rb_[Own[i]] = rb_[Nei[i]];
      dmdt_[Own[i]] = alpha[Nei[i]]*rb_[Nei[i]]*As_[Nei[i]];

      alpha[Nei[i]] = alpha0[Nei[i]] - rb_[Nei[i]]*As_[Nei[i]]*dt;
      if (alpha[Nei[i]] < 0)
      {
          scalar Vr = -alpha[Nei[i]]*V[Nei[i]];
          alpha[Nei[i]] = SMALL;

          // Find Neighbour cell
          label NNei = findNeighbour(alpha0, Nei[i]);
          alpha[NNei] = alpha0[NNei] - Vr/V[NNei];
          if (alpha[NNei] < 0)
          {
              FatalErrorInFunction
                << "Regression is very fast!\n"
                << "Hint: Reduce time step."
                << exit(FatalError);
          }
      }
    }
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
  const scalar One(1 - SMALL);
  const scalar Zero(SMALL);
  interface_ = dimensionedScalar(dimless, 0.0);

  forAll(Own, i)
  {
    // case:1 Interface is present at the center of the owner cell
    if (alpha[Own[i]] == 0.5 && alpha[Nei[i]] == One)
    {
      interface_[Own[i]] = 1;
    }
    // case:2 Interface is present in between owner center and face
    else if ((alpha[Own[i]] < 0.5) && (alpha[Nei[i]] == One))
    {
      interface_[Own[i]] = 1;
    }
    // case:3 Interface is present exactly at the center
    else if ((alpha[Own[i]] == Zero) && (alpha[Nei[i]] == One))
    {
      interface_[Own[i]] = 1;
    }
    // case:4 Interface is present in between face and neighbour center
    else if ((alpha[Own[i]] == Zero) && (alpha[Nei[i]] > 0.5))
    {
      interface_[Own[i]] = 1;
    }
    // case:5 Interface is not present (Ignore)
    else continue;
  }

}

Foam::tmp<Foam::volScalarField>
Foam::interfaceTrackingModels::subCellularInterfaceMotion::interface() const
{
    return Foam::tmp<Foam::volScalarField>(new volScalarField("tinterface", interface_));
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

Foam::tmp<Foam::volVectorField>
Foam::interfaceTrackingModels::subCellularInterfaceMotion::nHat() const
{
    return As_*dimensionedScalar(dimLength, 1.0)*vector(1, 0, 0)/max(Foam::mag(As_*dimensionedScalar(dimLength, 1.0)), SMALL);
}

Foam::tmp<Foam::volScalarField>
Foam::interfaceTrackingModels::subCellularInterfaceMotion::dmdt() const
{
    return Foam::tmp<Foam::volScalarField>(new volScalarField("tdmdt", dmdt_));
}
// ************************************************************************* //
