/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenFOAM Foundation and Ganeshkumar V, IIT Gandhinagar
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

#include "Surface.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::Surface::Surface
(
    const fvMesh& mesh,
    const word name,
    const dictionary& dict
)
:
    alpha_
    (
        volScalarField
        (
            IOobject
            (
                IOobject::groupName("alpha", name),
                mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    ),
    alphaOld_(alpha_),
    interface_
    (
        volScalarField
        (
          IOobject
          (
            IOobject::groupName("interface", name),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
          ),
          mesh,
          dimensionedScalar("", dimless, 0)
        )
    ),
    interfaceOwners_(nullptr),
    interfaceNeighbours_(nullptr),
    As_
    (
        volScalarField
        (
          IOobject
          (
            IOobject::groupName("As", name),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
          ),
          mesh,
          dimensionedScalar("", dimArea/dimVolume, 0)
        )
    ),
    nHat_
    (
        volVectorField
        (
          IOobject
          (
            IOobject::groupName("nHat", name),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
          ),
          mesh,
          dimensionedVector("", dimless, vector(0, 0, 0))
        )
    ),
    rb_
    (
        volScalarField
        (
          IOobject
          (
            IOobject::groupName("rb", name),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
          ),
          mesh,
          dimensionedScalar("", dimVelocity, 0)
        )
    ),
    n_(dict.get<scalar>("n")),
    f_(dict.get<scalar>("f")),
    Pmin_(dict.get<scalar>("Pmin")),
    a_(dict.get<scalar>("a"))
{
    // Find interface
    alpha_.clip(SMALL, 1 - SMALL);
    alphaOld_.clip(SMALL, 1 - SMALL);

    this->findInterface();

    label interfaceCells(sum(interface_).value());
    interfaceOwners_.reset(new labelList(interfaceCells, -1));
    interfaceNeighbours_.reset(new labelList(interfaceCells, -1));

    findInterfaceCells();
    Info << "Interface: " << sum(interface_) << endl;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::Surface::rb(const scalar P)
{
    if (P < Pmin_)
        return (a_*(pow(Pmin_/1e6, n_))).value()*1e-2;
    else
        return (a_*(pow(P/1e6, n_))).value()*1e-2;
    // if (P > 1e6)
    //     return (a_*(pow(P/1e6, n_))).value()*1e-2;
    // else
    //     return (a_*(pow(P/1e6, n_*f_))).value()*1e-2;
}

Foam::label Foam::Surface::findNeighbour
(
    const volScalarField& alpha,
    scalar NEI
)
{
    // Finds neighbour cell
    const fvMesh& mesh = alpha.mesh();
    const labelList& Own = mesh.owner();
    const labelList& Nei = mesh.neighbour();
    const scalar One(1 - SMALL);

    forAll(Own, i)
    {
        if (Own[i] == NEI)
        {
            //Info << "Nei: " << Nei[i] << " alpha: " << alpha[Nei[i]] << endl;
            if (alpha[Nei[i]] == One)
            {
                return Nei[i];
            }
        }
    }
    return -1;
}

Foam::scalar Foam::Surface::findNeighbourSurfaceArea
(
    const volScalarField& alpha,
    scalar Owner, 
    scalar Neighbour
)
{
    // Finds neighbour cell
    const fvMesh& mesh = alpha.mesh();
    const surfaceScalarField& Sf = mesh.magSf();
    const labelList& Own = mesh.owner();
    const labelList& Nei = mesh.neighbour();
    const scalar One(1 - SMALL);
    
    if (Neighbour != -1)
    {    
        forAll(Own, i)
        {
            if ((Own[i] == Owner) && (Nei[i] == Neighbour))
            {
                return Sf[i];
            }
        }
    }
    else
    {
        // Neighbour of interface cell is not found! 
        // Owner of the inerface is expected to be the boundary cell --- Searching for boundary patches
        const fvPatchList& patches = mesh.boundary();
        forAll(patches, i)
        {
            const fvPatch& patch = patches[i];
            const labelList& fC = patch.faceCells();
            const scalarField& Sf = patch.magSf();
            forAll(fC, j)
            {
                if (fC[j] == Owner)
                {
                    return Sf[j];
                }
            }
        }
    }
    Info << "--------- Neighbour Surface Area not found! ------------- " << exit(FatalError);
    return 0;
}

void Foam::Surface::findInterface()
{
    // -If found interface -> interface_ = 1
    // use owner neighbour approach
    const fvMesh& mesh = alpha_.mesh();
    const volScalarField& alpha = alpha_;
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
        else if ((alpha[Own[i]] == Zero) && (alpha[Nei[i]] >= 0.5))
        {
            interface_[Own[i]] = 1;
        }
        // case:5 Interface is present in the boundary cell
        else if ((alpha[Own[i]] == Zero) && (alpha[Nei[i]] < 0.5))
        {
            // Additional check for boundary cells
            const fvPatchList& patches = mesh.boundary();
            forAll(patches, j)
            {
                const fvPatch& patch = patches[j];
                const labelList& fC = patch.faceCells();
                forAll(fC, k)
                {
                    if (fC[k] == Nei[i])
                    {
                        interface_[Nei[i]] = 1;
                    }
                }
            }
        }
        // case:6 Interface is not present (Ignore)
        else  continue;
    }
}

void Foam::Surface::findInterfaceCells()
{

    const volScalarField& alpha0 = alphaOld_;
    const fvMesh& mesh = alpha_.mesh();
    const labelList& Own = mesh.owner();
    const labelList& Nei = mesh.neighbour();
    const scalar Zero(SMALL);
    interface_ = dimensionedScalar(dimless, 0.0);

    // interface owners and neighbours
    labelList& iOwners(interfaceOwners_());
    labelList& iNeighbours(interfaceNeighbours_());
    iOwners = -1;
    iNeighbours = -1;

    // Internal Cells
    label j = 0;
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
          iOwners[j] = Own[i];
          iNeighbours[j] = Nei[i];
          j++;
      }
    }

}

Foam::tmp<Foam::volScalarField> Foam::Surface::regressInterface
(
    const volScalarField& p, const labelList& bed
)
{
    // dmdt Field
    tmp<volScalarField> tdmdt
    (
        new volScalarField
        (
            IOobject("dmdt", p.mesh()),
            p.mesh(),
            dimensionedScalar("", dimless/dimTime, 0.0)
        )
    );
    volScalarField& dmdt_(tdmdt.ref());

    const volScalarField& alpha0 = alphaOld_;
    volScalarField& alpha = alpha_;
    const fvMesh& mesh = alpha_.mesh();
    const labelList& Own = mesh.owner();
    const labelList& Nei = mesh.neighbour();
    const surfaceScalarField& Sf = mesh.magSf();
    const scalar dt = mesh.time().deltaTValue();
    const scalarField& V = mesh.V();
    const scalar Zero(SMALL);

    interface_ = dimensionedScalar(dimless, 0.0);
    As_ = dimensionedScalar(As_.dimensions(), 0.0);
    rb_ = dimensionedScalar(rb_.dimensions(), 0.0);
    nHat_ = dimensionedVector(dimless, vector(0, 0, 0));

    // interface owners and neighbours
    labelList& iOwners(interfaceOwners_());
    labelList& iNeighbours(interfaceNeighbours_());
    iOwners = -1;
    iNeighbours = -1;

    // Internal Cells
    label k = 0;
    // Info << "Name: " << alpha.name() << " ";
    forAll(Own, i)
    {
      // case:1 Interface is present in the Neighbour Cell
      if
      (
          (alpha0[Own[i]] == Zero && alpha0[Nei[i]] > Zero) &&
          (Nei[i] == Own[i] + 1)
      )
      {
//	Info << "Flame Regress: Cell: " << Nei[i]  << " - " << alpha0[Nei[i]]; 
          interface_[Nei[i]] = 1;
          iOwners[k] = Own[i];
          iNeighbours[k] = Nei[i];

          // Constant Area ------------------
          // As_[Nei[i]] = Sf[i]/V[Nei[i]];  // Area of face between owner and neighbour
          
          // Linear Interpolation of Area -----------------
          scalar Asi = Sf[i];
          label NNei = findNeighbour(alpha0, Nei[i]);
          scalar Asip1 = findNeighbourSurfaceArea(alpha0, Nei[i], NNei);
          scalar Sfj = alpha0[Nei[i]]*Asi + (1 - alpha0[Nei[i]])*Asip1;
          As_[Nei[i]] = Sfj/V[Nei[i]];  // Area of face between owner and neighbour
          // ----------------------------------------------
          
          rb_[Nei[i]] = rb(p[bed[k]]);  // burning Rate
          dmdt_[Nei[i]] = rb_[Nei[i]]*As_[Nei[i]];
          nHat_[Nei[i]] = vector(1, 0, 0);
          //Info << "Bed Regress: " << Nei[i] << " ( " << dmdt_[Nei[i]] << " ) ";
          scalar newalpha = alpha0[Nei[i]] - rb_[Nei[i]]*As_[Nei[i]]*dt;
          if (newalpha < 0)
          {
              scalar Vr = -newalpha*V[Nei[i]];
              alpha[Nei[i]] = SMALL;

              // Find Neighbour of Neighbour cell
              bool isFound = false;
            //   label NNei = findNeighbour(alpha0, Nei[i]);

              if (NNei != -1)
              {
                    alpha[NNei] = alpha0[NNei] - Vr/V[NNei];
                    isFound = true;
                    if (alpha[NNei] < 0)
                    {
                        FatalErrorInFunction
                          << "Regression is very fast!\n"
                          << "Hint: Reduce time step."
                          << exit(FatalError);
                    }
              }
              if (isFound == false) // No adjacent cells have been found and hence stoping the regression here.
              {
                   // Correcting source terms for termination
                   dmdt_[Nei[i]] = alpha0[Nei[i]]*V[Nei[i]]/dt;
                   dmdt_[Own[i]] = 0.0;
              }
          }
          else
          {
              alpha[Nei[i]] = newalpha;
          }
	//Info << " New Alpha: " << alpha[Nei[i]] << endl;
          k++;
      }
    }
    Info << endl;
    return tdmdt;
}

Foam::tmp<Foam::volScalarField> Foam::Surface::regressInterface
(
    const volScalarField& p,
    const volScalarField& dmdt,
    const labelList& flame,
    const scalar fp,
    const scalar MR
)
{

    // dmdt Field
    tmp<volScalarField> tdmdt
    (
        new volScalarField
        (
            IOobject("dmdt", p.mesh()),
            p.mesh(),
            dimensionedScalar("", dimless/dimTime, 0.0)
        )
    );
    volScalarField& dmdt_(tdmdt.ref());

    // Regress interface based on dmdt
    const volScalarField& alpha0 = alphaOld_;
    volScalarField& alpha = alpha_;
    const fvMesh& mesh = alpha.mesh();
    const labelList& Own = mesh.owner();
    const labelList& Nei = mesh.neighbour();
    const surfaceScalarField& Sf = mesh.magSf();
    const scalar dt = mesh.time().deltaTValue();
    const scalarField& V = mesh.V();
    const scalar Zero(SMALL);

    interface_ = dimensionedScalar(dimless, 0.0);
    As_ = dimensionedScalar(As_.dimensions(), 0.0);
    rb_ = dimensionedScalar(rb_.dimensions(), 0.0);
    nHat_ = dimensionedVector(dimless, vector(0, 0, 0));

    // interface owners and neighbours
    labelList& iOwners(interfaceOwners_());
    labelList& iNeighbours(interfaceNeighbours_());
    iOwners = -1;
    iNeighbours = -1;
    // Info << "Name: " << alpha.name() << " ";

    // Internal Cells
    label k = 0;
    forAll(Own, i)
    {
        // case:1 Interface is present in the Neighbour Cell
        if
        (
            (alpha0[Own[i]] == Zero && alpha0[Nei[i]] > Zero) &&
            (Nei[i] == Own[i] + 1)
        )
        {
//		Info << "Bed Regress -> Cell: " << Nei[i] << " - " << alpha0[Nei[i]];
            scalar dmdtflame = 0;
            scalar Vflame = 0;

            if (flame.size() > 0)
            {
              if (flame[k] != -1)
              {
                  dmdtflame = dmdt[flame[k]];
                  Vflame = V[flame[k]];
              }
            }

            interface_[Nei[i]] = 1;
            iOwners[k] = Own[i];
            iNeighbours[k] = Nei[i];

            As_[Nei[i]] = Sf[i]/V[Nei[i]];  // Area of face between owner and neighbour
            
            // Linear Interpolation of Area -----------------
            scalar Asi = Sf[i];
            label NNei = findNeighbour(alpha0, Nei[i]);
            scalar Asip1 = findNeighbourSurfaceArea(alpha0, Nei[i], NNei);
            scalar Sfj = alpha0[Nei[i]]*Asi + (1 - alpha0[Nei[i]])*Asip1;
            As_[Nei[i]] = Sfj/V[Nei[i]];  // Area of face between owner and neighbour
            // ----------------------------------------------

            rb_[Nei[i]] = rb(p[Nei[i]]);  // burning Rate
            dmdt_[Nei[i]] = (1 - alpha0[Nei[i]])*dmdtflame*Vflame/V[Nei[i]];
            nHat_[Nei[i]] = vector(1, 0, 0);

            As_[Own[i]] = As_[Nei[i]];
            rb_[Own[i]] = rb_[Nei[i]];
            dmdt_[Own[i]] = alpha0[Nei[i]]*dmdtflame*Vflame/V[Own[i]];
            nHat_[Own[i]] = vector(1, 0, 0);
            Info << " ( " << Nei[i] << " " << Own[i] << " ) -> "
              << " ( " << dmdt_[Nei[i]] << " " << dmdt_[Own[i]] << " ) "
              << " Flame: dmdt: " << dmdt[flame[k]] << " vF: " << V[flame[k]]
              << " fp: " << fp << endl;

            scalar newalpha = alpha0[Nei[i]] - (1.0 - MR/fp)*dmdt[flame[k]]*V[flame[k]]*dt/V[Nei[i]];
            if (newalpha < 0)
            {
                scalar Vr = -newalpha*V[Nei[i]];
                alpha[Nei[i]] = SMALL;

                // Find Neighbour of Neighbour cell
                bool isFound = false;
                // label NNei = findNeighbour(alpha0, Nei[i]);

                if (NNei != -1)
                {
                     alpha[NNei] = alpha0[NNei] - Vr/V[NNei];
                     isFound = true;
                     if (alpha[NNei] < 0)
                     {
                          FatalErrorInFunction
                            << "Regression is very fast!\n"
                            << "Hint: Reduce time step."
                            << exit(FatalError);
                     }
                }
                if (isFound == false) // No adjacent cells have been found and hence stoping the regression here.
                {
                     // Correcting source terms for termination
                     dmdt_[Nei[i]] = alpha0[Nei[i]]*V[Nei[i]]/dt;
                     dmdt_[Own[i]] = 0.0;
                }
            }
            else
            {
                alpha[Nei[i]] = newalpha;
            }
//		Info << " New Alpha: " << alpha[Nei[i]] << endl;
            k++;
        }
    }

    return tdmdt;
}
