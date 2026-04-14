{
    /*
    TSolidus = alpha1*Tsolidus1 + (1.0 - alpha1)*Tsolidus2;
    TLiquidus = alpha1*Tliquidus1 + (1.0 - alpha1)*Tliquidus2;
    //beta = alpha1*beta1 + (1.0 - alpha1)*beta2;
    //rhok = (1.0 - max(epsilon1*(beta)*(T - TSolidus), 0.0));

    scalarField& alpha1I = alpha1;
    scalarField& TSolidusI = TSolidus;
    scalarField& TLiquidusI = TLiquidus;
    const scalar Tsolidus1Value = Tsolidus1.value();
    const scalar Tliquidus1Value = Tliquidus1.value();
    forAll(alpha1I, celli)
    {
        if (alpha1I[celli] > 0.05)
        {
            TSolidusI[celli] = Tsolidus1Value;
            TLiquidusI[celli] = Tliquidus1Value;
        }
    }

    TSolidus.correctBoundaryConditions();
    TLiquidus.correctBoundaryConditions();
    //beta.correctBoundaryConditions();  // ADDED
    //rhok.correctBoundaryConditions();  // ADDED

    LatentHeat = alpha1*LatentHeat1 + (1.0 - alpha1)*LatentHeat2;
    LatentHeat.correctBoundaryConditions();
    */

    #include "updateKappaCp.H"

    gh = (g & mesh.C());
    ghf = (g & mesh.Cf());
}
