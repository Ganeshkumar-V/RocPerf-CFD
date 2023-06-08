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

Class
    Foam::interfaceTrackingModel

Description

SourceFiles
    interfaceTrackingModel.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceTrackingModel_H
#define interfaceTrackingModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "PtrListDictionary.H"
#include "phaseModel.H"

namespace Foam
{

// Forward Declarations
class phasePair;

/*---------------------------------------------------------------------------*\
                        Class interfaceTrackingModel Declaration
\*---------------------------------------------------------------------------*/

class interfaceTrackingModel
{
protected:

    // Protected Data

        //- Phase pair
        const phasePair& pair_;



public:

    //- Runtime type information
    TypeName("interfaceTrackingModel");

    typedef PtrListDictionary<phaseModel> phaseModelList;

    //- Propellant phase
    word propellant_;

    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            interfaceTrackingModel,
            dictionary,
            (
                const dictionary& dict,
                const phasePair& pair
            ),
            (dict, pair)
        );


    // Static Data Members

        //- Mass transfer rate dimensions
        static const dimensionSet dimDmdt;


    // Constructors

        //- Construct from a dictionary and a phase pair
        interfaceTrackingModel
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~interfaceTrackingModel() = default;


    // Selectors

        static autoPtr<interfaceTrackingModel> New
        (
            const dictionary& dict,
            const phasePair& pair
        );


    // Member Functions
    virtual void correct() = 0;
    virtual void regress(volScalarField& alpha, const volScalarField& alphaOld) = 0;

    virtual void findInterface(const volScalarField& alpha) = 0;

    //- The burning rate
    virtual tmp<volScalarField> rb() const = 0;

    //- Interface Area/Volume
    virtual tmp<volScalarField> As() const = 0;

    //- Interface Area vector
    virtual tmp<volVectorField> nHat() const;

    //- Propellant Density
    virtual tmp<volScalarField> rhop() const;

    //- Mass Transfer Rate
    virtual tmp<volScalarField> dmdt() const;

    //- Interface Cell Field
    virtual tmp<volScalarField> interface() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
