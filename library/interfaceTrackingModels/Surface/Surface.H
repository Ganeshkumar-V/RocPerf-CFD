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
#ifndef Surface_H
#define Surface_H

#include "phaseSystem.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class Surface
\*---------------------------------------------------------------------------*/

class Surface
{
protected:

  // Protected data

      // Volume fraction used for identifying interface
      volScalarField alpha_;

      // Volume fraction at old time step
      volScalarField alphaOld_;

      // cells containing surface
      volScalarField interface_;
      autoPtr<labelList> interfaceOwners_;
      autoPtr<labelList> interfaceNeighbours_;

      // surface area
      volScalarField As_;

      // surface normal
      volVectorField nHat_;

      // burning rate
      volScalarField rb_;

      // Burning rate Parameters
      scalar n_;
      scalar f_;
      scalar Pmin_;
      dimensionedScalar a_;

public:

  // Constructors

    //- Construct from dictionary
    Surface
    (
        const fvMesh& mesh,
        const word name,
        const dictionary& dict
    );

  // Methods

    // Public functions

    // Burning Rate
    scalar rb(const scalar P);

    // Find neighbour cell for cell transition
    label findNeighbour(const volScalarField& alpha, scalar NEI);
    scalar findNeighbourSurfaceArea(const volScalarField& alpha, scalar Owner, scalar Neighbour);

    // Find interface and correct interface, As and nHat Fields
    void findInterface();
    void findInterfaceCells();

    // Find and regress
    tmp<volScalarField> regressInterface
    (
        const volScalarField& p,
        const labelList& bed
    );

    // Regress based on mdot
    tmp<volScalarField> regressInterface
    (
        const volScalarField& p,
        const volScalarField& dmdt,
        const labelList& flame,
        const scalar fp,
        const scalar MR
    );


    // Public Access Functions

    tmp<volScalarField> interface() const
    {
        return tmp<volScalarField>(new volScalarField("interface", interface_));
    }

    tmp<volScalarField> As() const
    {
        return tmp<volScalarField>(new volScalarField("As", As_));
    }

    tmp<volVectorField> nHat() const
    {
        return tmp<volVectorField>(new volVectorField("nHat", nHat_));
    }

    tmp<volScalarField> rb() const
    {
        return tmp<volScalarField>(new volScalarField("rb", rb_));
    }

    const labelList& iOwners() const
    {
        return interfaceOwners_();
    }

    const labelList& iNeighbours() const
    {
        return interfaceNeighbours_();
    }

    const volScalarField& alpha()
    {
        return alpha_;
    }

    void store()
    {
        alphaOld_ = alpha_;
    }

};

}

#endif
