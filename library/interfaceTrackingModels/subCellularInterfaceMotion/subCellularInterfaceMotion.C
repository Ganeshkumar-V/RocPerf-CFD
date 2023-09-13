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
    crb_("", dimVelocity, dict.getOrDefault<scalar>("rb", -1)),
    pBufnB(Pstream::commsTypes::nonBlocking),
    transfer_(0),
    transferAlpha_(0, 0),
    communicate_(0)
{
  const phaseModel& phase = pair_.phase1();
  const volScalarField& alpha
        = phase.db().lookupObject<volScalarField>("alpha." + propellant_);
  this->findInterface(alpha);

  forAll(alpha.mesh().boundary(), patchi)
  {
    if (isType<processorFvPatch>(alpha.mesh().boundary()[patchi]))
    {
      // Check for shared cell
      const processorPolyPatch& pp
          = refCast<const processorPolyPatch>(alpha.mesh().boundaryMesh()[patchi]);
      if (pp.owner())
      {
          transferAlpha_ = alpha.boundaryField()[patchi];
      }
    }
  }

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

  // Initialize to allow transfer/communicate_ if waranted
  transfer_ = 0;
  communicate_ = 0;
  forAll(mesh.boundary(),patchi)
  {
    if (isType<processorFvPatch>(mesh.boundary()[patchi]))
    {
      const processorPolyPatch& pp
          = refCast<const processorPolyPatch>(mesh.boundaryMesh()[patchi]);
      if (pp.owner())
      {
        transferAlpha_ = alpha0.boundaryField()[patchi].patchNeighbourField();
      }
    }
  }

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
      dmdt_[Nei[i]] = (1 - alpha0[Nei[i]])*rb_[Nei[i]]*As_[Nei[i]];

      As_[Own[i]] = As_[Nei[i]];
      rb_[Own[i]] = rb_[Nei[i]];
      dmdt_[Own[i]] = alpha0[Nei[i]]*rb_[Nei[i]]*As_[Nei[i]];

      scalar newalpha = alpha0[Nei[i]] - rb_[Nei[i]]*As_[Nei[i]]*dt;
      if (newalpha < 0)
      {
          // Find Neighbour of Neighbour cell
          label NNei = findNeighbour(alpha0, Nei[i]);

          if (NNei != -1)
          {
                scalar Vr = -newalpha*V[Nei[i]];
                alpha[Nei[i]] = SMALL;

                alpha[NNei] = alpha0[NNei] - Vr/V[NNei];
                if (alpha[NNei] < 0)
                {
                    FatalErrorInFunction
                      << "Regression is very fast!\n"
                      << "Hint: Reduce time step."
                      << exit(FatalError);
                }
          }
          else
          {
               // check for processor shared neighbour cells and regress
               forAll(mesh.boundary(),patchi)
               {
                 if (isType<processorFvPatch>(mesh.boundary()[patchi]))
                 {
                   const processorPolyPatch& pp
                       = refCast<const processorPolyPatch>(mesh.boundaryMesh()[patchi]);
                   const scalarField& nf(alpha0.boundaryField()[patchi].patchNeighbourField());

                   transferAlpha_ = alpha.boundaryField()[patchi].patchNeighbourField();

                   if (pp.owner())
                   {
                     const labelList& fC(mesh.boundary()[patchi].faceCells());
                     forAll(fC, celli)
                     {
                       if ((fC[celli] == Nei[i]) && (nf[celli] == 1 - SMALL))
                       {
                         transfer_ = 1.0; // To allow processor exchange information
                         // Interface Transfer
                         scalar Vr = -newalpha*V[Nei[i]];
                         alpha[Nei[i]] = SMALL;

                         transferAlpha_[celli] = nf[celli] - Vr/V[Nei[i]];
                         Pout << "transferCell: " << transferAlpha_ << endl;
                         // This should be V[Nf[celli]]! But how to get neighbouring cell's volume?
                       }
                     }
                   }
                 }
               }
          }
      }
      else
      {
          alpha[Nei[i]] = newalpha;
      }
    }
  }

  // Check the need for transfer
  reduce(transfer_, sumOp<scalar>());

  if (transfer_ > 0)
  {
    // Send and Receive Information about interface transfer
    forAll(mesh.boundary(), patchi)
    {
      // check for processor patch
      if (isType<processorFvPatch>(mesh.boundary()[patchi]))
      {
        const processorPolyPatch& patch
            = refCast<const processorPolyPatch>(mesh.boundaryMesh()[patchi]);

        if (patch.owner())
        {
          // Send alpha Information
          UOPstream sendToNeighbour(patch.neighbProcNo(), pBufnB);
          sendToNeighbour << transferAlpha_;
        }
      }
    }
    pBufnB.finishedSends();

    forAll(mesh.boundary(), patchi)
    {
      // check for processor patch
      if (isType<processorFvPatch>(mesh.boundary()[patchi]))
      {
        const processorPolyPatch& patch
            = refCast<const processorPolyPatch>(mesh.boundaryMesh()[patchi]);

        if (patch.neighbour())
        {
          // Receive alpha Information
          scalarField transferCellReceive(transferAlpha_.size(), 0);
          UIPstream recvFromOwner(patch.neighbProcNo(), pBufnB);
          recvFromOwner >> transferCellReceive;
          const labelList& faceCells(patch.faceCells());
          const scalarField alpha0NF(alpha0.boundaryField()[patchi].patchNeighbourField());
          forAll(faceCells, i)
          {
              if (alpha0NF[i] != transferCellReceive[i])
              alpha[faceCells[i]] = transferCellReceive[i];
          }
        }
      }
    }
  }


  if ((alpha0[0] != 1 - SMALL) && (alpha0[0] != SMALL)) communicate_ = 1.0;
  // check the need for communiation of sources
  reduce(communicate_, sumOp<scalar>());

  // Perform calculation to communicate Sources back to owner
  if (communicate_ > 0)
  {
    scalar dmdtComm = 0;
    scalar AsComm = 0;
    scalar rbComm = 0;

    if ((alpha0[0] != 1 - SMALL) && (alpha0[0] != SMALL))
    {
      interface_[0] = 1;
      As_[0] = Sf[0]/V[0];
      rb_[0] = rb(p[0]);
      dmdt_[0] = (1 - alpha0[0])*rb_[0]*As_[0];

      dmdtComm = alpha0[0]*rb_[0]*As_[0];
      AsComm = As_[0];
      rbComm = rb_[0];

      alpha[0] = alpha0[0] - rb_[0]*As_[0]*dt;
      if (alpha[0] < 0)
      {
        alpha[1] = alpha0[1] + alpha[0]*V[0]/V[1];
        alpha[0] = SMALL;
      }
    }

    // Send Data
    forAll(mesh.boundary(), patchi)
    {
      // check for processor patch
      if (isType<processorFvPatch>(mesh.boundary()[patchi]))
      {
        const processorPolyPatch& patch
            = refCast<const processorPolyPatch>(mesh.boundaryMesh()[patchi]);

        if (patch.neighbour())
        {
          // Send alpha Information
          UOPstream sendToOwner(patch.neighbProcNo(), pBufnB);
          sendToOwner << dmdtComm;
          sendToOwner << AsComm;
          sendToOwner << rbComm;
        }
      }
    }
    pBufnB.finishedSends();

    // Receive Data
    forAll(mesh.boundary(), patchi)
    {
      if (isType<processorFvPatch>(mesh.boundary()[patchi]))
      {
        const processorPolyPatch& patch
            = refCast<const processorPolyPatch>(mesh.boundaryMesh()[patchi]);

        if (patch.owner())
        {
          // Receive from Neigbour
          UIPstream recvFromNeigh(patch.neighbProcNo(), pBufnB);
          recvFromNeigh >> dmdt_[dmdt_.size() - 1];
          recvFromNeigh >> As_[dmdt_.size() - 1];
          recvFromNeigh >> rb_[dmdt_.size() - 1];
        }
      }
    }
  }

  dmdt_.correctBoundaryConditions();
  As_.correctBoundaryConditions();
  rb_.correctBoundaryConditions();
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
