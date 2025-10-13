/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

Application
    postProcess

Description
    Execute the set of functionObjects specified in the selected dictionary
    (which defaults to system/controlDict) or on the command-line for the
    selected set of times on the selected set of fields.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "profiling.H"
#include "timeSelector.H"
#include "ReadFields.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "uniformDimensionedFields.H"
#include "fileFieldSelection.H"
#include "mapPolyMesh.H"
#include "fvCFD.H"
#include "multiPhaseSystem.H"
#include "pimpleControl.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "phasePair.H"
#include "GeometricField.H"
#include "scalarIOField.H"
#include "processorFvPatch.H"
#include <stdio.h>
#include <sstream>

using namespace Foam;

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Execute the set of functionObjects specified in the selected"
        " dictionary or on the command-line for the"
        " selected set of times on the selected set of fields"
    );

    timeSelector::addOptions();
    #include "addProfilingOption.H"
    #include "addRegionOption.H"
    #include "addFunctionObjectOptions.H"

    // Set functionObject post-processing mode
    functionObject::postProcess = true;

    #include "setRootCase.H"

    if (args.found("list"))
    {
        functionObjectList::list();
        return 0;
    }

    #include "createTime.H"
    // instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createNamedMesh.H"
    #include "readGravitationalAcceleration.H"
    #include "readhRef.H"

    Info<< "Creating phaseSystem\n" << endl;

    autoPtr<multiPhaseSystem> fluidPtr
    (
        multiPhaseSystem::New(mesh)
    );
    multiPhaseSystem& fluid = fluidPtr();
    multiPhaseSystem::phaseModelList& phases = fluid.phases();

    #include "gh.H"

    // gas phase
    const phaseModel& gasPhase(phases[0]);
    const tmp<volScalarField>& tmu(gasPhase.thermo().mu());
    const volScalarField& mu(tmu());
    dimensionedScalar dRef("", dimLength, fluid.get<scalar>("dRef"));
    volScalarField Re
    (
        IOobject
        (
          "Re",
          runTime.timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
        ),
        (gasPhase.rho()*mag(gasPhase.U())*dRef/mu)()
    );
    Re.write();

    // particle phase
    const phaseModel& particlePhase(phases[1]);

    Info<< "End\n" << endl;

    return 0;
}
