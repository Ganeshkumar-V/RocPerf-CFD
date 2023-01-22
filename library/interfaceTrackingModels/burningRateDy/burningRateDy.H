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
    Foam::interfaceTrackingModels::burningRateDy

Description

SourceFiles
    burningRateDy.C

\*---------------------------------------------------------------------------*/

#ifndef burningRateDy_H
#define burningRateDy_H

#include "interfaceTrackingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace interfaceTrackingModels
{

/*---------------------------------------------------------------------------*\
                          Class burningRateDy Declaration
\*---------------------------------------------------------------------------*/

class burningRateDy
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
    TypeName("burningRateDy");

    scalar n;
    dimensionedScalar a;

    const phaseModel& propellantPhase_;

    //- Interface area and Vector
    volScalarField interfaceArea_;
    volVectorField interfaceAreaVector_;

    //- Interface Iso-value
    scalar isoAlpha_;

    volScalarField rb_;

    dimensionedScalar crb_;
    // Constructors

        //- Construct from components
        burningRateDy
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~burningRateDy();

    // Member Functions
    virtual void correct();

    virtual void regress(volScalarField& alpha);

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
