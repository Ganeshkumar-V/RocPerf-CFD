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

Class
    Foam::interfaceTrackingModels::burningRate

Description

SourceFiles
    burningRate.C

\*---------------------------------------------------------------------------*/

#ifndef burningRate_H
#define burningRate_H

#include "interfaceTrackingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace interfaceTrackingModels
{

/*---------------------------------------------------------------------------*\
                          Class burningRate Declaration
\*---------------------------------------------------------------------------*/

class burningRate
:
    public interfaceTrackingModel
{
private:

  // Private Member Functions

      //- Update interface
      void updateInterface();

      //- mesh
      const fvMesh& mesh_;

public:

    //- Runtime type information
    TypeName("burningRate");

    scalar n;
    dimensionedScalar a;

    const phaseModel& propellantPhase_;

    //- Interface area and Vector
    volScalarField interfaceArea_;
    volVectorField interfaceAreaVector_;

    //- Interface Iso-value
    scalar isoAlpha_;
    scalar ztol_;

    volScalarField rb_;

    dimensionedScalar crb_;
    // Constructors

        //- Construct from components
        burningRate
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~burningRate();

    // Member Functions
    virtual void correct();

    virtual void regress(volScalarField& alpha, const volScalarField& alphaOld);

    virtual void findInterface(const volScalarField& alpha);

    //- The burning rate
    virtual tmp<volScalarField> rb() const;

    //- Interface Area / Volume
    virtual tmp<volScalarField> As() const;

    //- Interface Area / Volume
    virtual tmp<volVectorField> nHat() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace interfaceTrackingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
