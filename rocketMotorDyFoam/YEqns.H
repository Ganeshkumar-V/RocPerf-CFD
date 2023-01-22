{
    autoPtr<phaseSystem::massTransferTable>
        massTransferPtr(fluid.massTransfer());

    phaseSystem::massTransferTable&
        massTransfer(massTransferPtr());

    forAll(fluid.multiComponentPhases(), multiComponentPhasei)
    {
        phaseModel& phase = fluid.multiComponentPhases()[multiComponentPhasei];

        UPtrList<volScalarField>& Y = phase.YActiveRef();
        const volScalarField& alpha = phase;
        const volScalarField& rho = phase.rho();

        forAll(Y, i)
        {
            fvScalarMatrix YiEqn
            (
              fvm::ddt(alpha, rho, Y[i])
              + fvm::div
              (
                phase.alphaRhoPhi(), Y[i],
                "div(" + phase.alphaRhoPhi()().name() + ",Yi)"
              )
             ==
               *massTransfer[Y[i].name()]
              + fvOptions(alpha, rho, Y[i])
            );

            YiEqn.relax();
            YiEqn.solve(mesh.solver("Yi"));

            // Clip Yi's
            Y[i].clip(SMALL, 1 - SMALL);
        }
    }
}
