/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
    Copyright (C) 2019-2020 DLR
    Copyright (C) 2025 Ganeshkumar V
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

#include "generalizedCylinderImplicitFunction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace implicitFunctions
    {
        defineTypeNameAndDebug(generalizedCylinderImplicitFunction, 0);
        addToRunTimeSelectionTable
        (
            implicitFunction,
            generalizedCylinderImplicitFunction,
            dict
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::implicitFunctions::generalizedCylinderImplicitFunction::generalizedCylinderImplicitFunction
(
    const point& origin,
    const vector& direction
)
:
    origin_(origin),
    radius_(),
    direction_(normalised(direction)),
    project_(tensor::I - direction_*direction_) // outer product
{}


Foam::implicitFunctions::generalizedCylinderImplicitFunction::generalizedCylinderImplicitFunction
(
    const dictionary& dict
)
:
    // __INTEL_COMPILER bug with inheriting constructors?? (issue #1821)
    origin_(dict.get<point>("origin")),
    direction_(dict.get<vector>("direction").normalise()),
    project_(tensor::I - direction_*direction_) // outer product
{
    if (dict.found("radius"))
    {
        radius_ = Function1<scalar>::New("radius", dict);
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "Please supply 'radius' as a scalar function (Function1)" << nl
            << exit(FatalIOError);
    }
}


// ************************************************************************* //
