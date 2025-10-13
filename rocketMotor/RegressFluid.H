{
Info << "RegressFluid: Solving for propellant volume fraction" << endl;
// Get alpha source
volScalarField alphaSource(phaseChange.alphaSource());
volScalarField distributedAlphaSource(0*alphaSource);

const scalarField& V(mesh.V());
const scalar dt(runTime.deltaT().value());

forAll(alphaSource, celli)
{
    if (alphaSource[celli]*dt > alpha1OldTime[celli])
    {
        label adjacent = findAdjacent(celli, beta1OldTime);
        if (adjacent == -1)
        {
            if (meshObject::debug)
            {
                Pout << " Adjacent cell not found in present processor to pass on alpha source. \n "
                 << " Celli: " << celli << " \n "
                 << " alpha1: " << alpha1OldTime[celli] << "\n "
                 << " alphaSource*dt: " << alphaSource[celli]*dt << endl;
            }
            
            distributedAlphaSource[celli] += alphaSource[celli];
        }
        else
        {
            distributedAlphaSource[adjacent] = (alphaSource[celli]*dt - alpha1OldTime[celli])*V[celli]/(V[adjacent]*dt);
            distributedAlphaSource[celli] = alpha1OldTime[celli]/dt;
            if (meshObject::debug)
            {
                Pout << "Moving the source from " << celli << " to " << adjacent << endl;
                Pout << " Celli: " << celli << " \n "
                << " alpha1: " << alpha1OldTime[celli] << "\n "
                << " alphaSource*dt: " << alphaSource[celli]*dt << endl;
                Pout << "alphaSource[celli]: " << distributedAlphaSource[celli]*dt << " alphaSource[adjacent]: " << distributedAlphaSource[adjacent]*dt << endl;
            }
            Pout << "~~~~~~~~~~~~ CELL TRANSITION ~~~~~~~~~~~~~~ " << endl;
            
            if (distributedAlphaSource[adjacent]*dt > alpha1OldTime[adjacent])
            {
                if (meshObject::debug)
                {
                    Warning << " cell: " << adjacent << " is going to become negative. "
                        << " alphaSource*dt: " << distributedAlphaSource[adjacent]*dt
                        << " alpha1: " << alpha1[adjacent]
                        << " ----- Please Reduce the time step!" << exit(FatalError);
                }
                Pout << "Time step is too large! " << " Please reduce the time step. " << exit(FatalError);
            }
        }

     }
    else
    {
        distributedAlphaSource[celli] += alphaSource[celli];

        // check for negative alpha in the adjacent cell
        if (distributedAlphaSource[celli]*dt > alpha1OldTime[celli])
        {
            
            if (meshObject::debug)
            {
                Warning << " cell: " << celli << " is going to become negative. "
                    << " alphaSource: " << distributedAlphaSource[celli]
                    << " alpha1: " << alpha1[celli] 
                    << " ----- Please Reduce the time step!" << exit(FatalError);
            }
            Pout << "Time step is too large! " << " Please reduce the time step. " << exit(FatalError);
        }
    }
}

// For processor boundary patches
distributedAlphaSource.correctBoundaryConditions(); // Send/Recv processor boundary information

{
    const volScalarField::Boundary& alpha1OldTimebF(alpha1OldTime.boundaryField());
    const volScalarField::Boundary& beta1OldTimebF(beta1OldTime.boundaryField());
    const volScalarField::Boundary& VbF(cellV.boundaryField());
    volScalarField::Boundary& distributedAlphaSourcebF(distributedAlphaSource.boundaryFieldRef());
    forAll(distributedAlphaSourcebF, bFi)
    {
        word BCtype = mesh.boundaryMesh().types()[bFi];
        if ( BCtype == "processor" )    // For Processor Boundary
        {
            const processorPolyPatch& pp = refCast<const processorPolyPatch>( mesh.boundaryMesh()[bFi]);
            scalarField VNF(VbF[bFi].patchNeighbourField());
            const scalarField alpha1OldTimeNF(alpha1OldTimebF[bFi].patchNeighbourField());
            const scalarField beta1OldTimeNF(beta1OldTimebF[bFi].patchNeighbourField());
            const UList<label>& fC(mesh.boundaryMesh()[bFi].faceCells());

            scalarField distributedAlphaSourceNF(distributedAlphaSourcebF[bFi].patchNeighbourField());
            forAll(fC, i)
            {
                // Receiving patch
                if (distributedAlphaSourceNF[i]*dt > alpha1OldTimeNF[i])
                {
                    if ((beta1OldTimeNF[i] < 1.0) && (beta1OldTime[fC[i]] == 1.0))
                    {
                        distributedAlphaSource[fC[i]] = (distributedAlphaSourceNF[i]*dt - alpha1OldTimeNF[i])*VNF[i]/(V[fC[i]]*dt);
                        if (meshObject::debug)
                        {
                            Pout << "Receiving source from: " << pp.neighbProcNo() << endl;
                            Pout << "AlphaSource-Neighbour[celli]: " << fC[i] << " -> "<< distributedAlphaSourceNF[i] << endl;
                            Pout << "New Alphasource: " << distributedAlphaSource[fC[i]] << " added to " << fC[i] << " alpha here is: " << alpha1OldTime[fC[i]] << endl;
                        }
                        Pout << "\n~~~~~~~~~~~~ INTER PROCESSOR -> CELL TRANSITION ~~~~~~~~~~~~~~ " << endl;

                        if (distributedAlphaSource[fC[i]]*dt > alpha1OldTime[fC[i]])
                        {
                            Pout << "Time step is too large! " << " Please reduce the time step. " << exit(FatalError);
                        }
                    }
                }

                // Sending Patch
                if (distributedAlphaSource[fC[i]]*dt > alpha1OldTime[fC[i]])
                {
                    if (beta1OldTimeNF[i] == 1.0)
                    {
                        if (meshObject::debug)
                        {
                            Pout << "Sending source from: " << pp.myProcNo() << endl;
                            Pout << "AlphaSource[i]: " << i << " -> "<< distributedAlphaSource[fC[i]] << endl;
                        }
                        distributedAlphaSource[fC[i]] = alpha1OldTime[fC[i]]/dt;
                    }

                }
            }
        }
    }
}
distributedAlphaSource.correctBoundaryConditions();

// For doamin boundary patches
{
    const volScalarField::Boundary& alpha1OldTimebF(alpha1OldTime.boundaryField());
    const volScalarField::Boundary& beta1OldTimebF(beta1OldTime.boundaryField());
    const volScalarField::Boundary& VbF(cellV.boundaryField());
    volScalarField::Boundary& distributedAlphaSourcebF(distributedAlphaSource.boundaryFieldRef());
    forAll(distributedAlphaSourcebF, bFi)
    {
        word BCtype = mesh.boundaryMesh().names()[bFi];
        if ( BCtype == bottomPatchName )    // For Non-processor Boundary
        {
            const UList<label>& fC(mesh.boundaryMesh()[bFi].faceCells());
            forAll(fC, i)
            {
                if (distributedAlphaSource[fC[i]]*dt > alpha1OldTime[fC[i]])
                {
                    // Info Statements
                    Info << " Boundary patch name: " << BCtype << endl;
                    Pout << "Interface is at the doamin boundary! So emptying the cell! " << endl;
                    Pout << "Celli: " << fC[i] << endl << "alphaSource*dt: " << distributedAlphaSource[fC[i]]*dt << endl;
                    distributedAlphaSource[fC[i]] = alpha1OldTime[fC[i]]/dt;
                }
            }
        }
    }
    
}
distributedAlphaSource.correctBoundaryConditions();

forAll(alpha1, celli)
{
    alpha1[celli] = alpha1OldTime[celli] - distributedAlphaSource[celli]*dt;
    if (mag(alpha1[celli]) < SMALL)
    {
        alpha1[celli] = 0.0;
    }
}

alpha1.correctBoundaryConditions();
distributedAlphaSource.correctBoundaryConditions();

// Concentrated mass source
phaseChange.psi0Ref() = distributedAlphaSource*phase1.thermo().rho();
normalizeVolumeFraction(alpha1, alpha2, beta1, initialVOF);

Info<< alpha1.name() << " Liquid volume fraction = "
                << alpha1.weightedAverage(mesh.V()).value()
                << "  Min(alpha1) = " << min(alpha1).value()
                << "  Max(alpha1) = " << max(alpha1).value()
                << endl;
}
