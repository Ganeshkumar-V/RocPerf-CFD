/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "calculateThrust.H"
#include "fluidThermo.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(calculateThrust, 0);
    addToRunTimeSelectionTable(functionObject, calculateThrust, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::calculateThrust::calc()
{
    const surfaceScalarField& phiGas =
                lookupObject<surfaceScalarField>("phi.gas");
    const surfaceScalarField& phiParticle =
                lookupObject<surfaceScalarField>("phi.particles");

    const volScalarField& alphaGas =
                lookupObject<volScalarField>("alpha.gas");
    const volScalarField& alphaParticle =
                lookupObject<volScalarField>("alpha.particles");

    const volVectorField& UGas =
                lookupObject<volVectorField>("U.gas");
    const volVectorField& UParticle =
                lookupObject<volVectorField>("U.particles");

    const fvMesh& mesh(alphaGas.mesh());

    const volScalarField rhoG = alphaGas;
    const volScalarField rhoP = alphaParticle;

    forAll(mesh.boundary(), patchi)
    {
      if (mesh.boundary()[patchi].name() == "outlet")
      {
        scalar mDotG = sum
        (
          rhoG.boundaryField()[patchi]*alphaGas.boundaryField()[patchi]
          *phiGas.boundaryField()[patchi]
        );

        scalar mDotP = sum
        (
          rhoP.boundaryField()[patchi]*alphaParticle.boundaryField()[patchi]
          *phiParticle.boundaryField()[patchi]
        );

        vectorField momG
          (rhoG.boundaryField()[patchi]*alphaGas.boundaryField()[patchi]
          *phiGas.boundaryField()[patchi]*UGas.boundaryField()[patchi]);

        vectorField momP
          (rhoP.boundaryField()[patchi]*alphaParticle.boundaryField()[patchi]
          *phiParticle.boundaryField()[patchi]*UParticle.boundaryField()[patchi]);

        fileName outputFile("output.txt");
        OFstream os(outputFile);
        os  << "This is first Line \n" << "This is second line." << endl;
      }
    }
    return true;
}

// * * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::calculateThrust::write()
{}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::calculateThrust::calculateThrust
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict)
{}


// ************************************************************************* //
