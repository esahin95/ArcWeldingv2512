{
    int iter = 0;
    scalarField residual(alphaLiquid.primitiveField().size(), 1.0);
    scalar maxResidual = 1.0;
    scalar meanResidual = 1.0;

    TRHS *= 0.0;

    rhoCp = rho * fluid.Cp();
    rhoCp.correctBoundaryConditions();  
    
    rhoPhiCp = fvc::interpolate(fluid.Cp()) * rhoPhi;

    do
    {
        alphaLiquid.storePrevIter();
        T.storePrevIter();
        TRHS =
            LatentHeat*((fvc::ddt(rho, alphaLiquid) + fvc::div(rhoPhi, alphaLiquid)));
        
        TRHS.correctBoundaryConditions();  

        /*
        Qv =
            0.82 * Lvap * Mm * p0
           * Foam::exp(Lvap * Mm * ((T - Tvap) / (R * T * Tvap)))
           / Foam::pow(2 * M_PI * Mm * R * T, 0.5);
        
        Qv.correctBoundaryConditions();  
        */
        
        fvScalarMatrix TEqn
        (
            fvm::ddt(rhoCp,T)
          + fvm::div(rhoPhiCp, T)
          - fvm::Sp(fvc::ddt(rhoCp) + fvc::div(rhoPhiCp), T)
          - fvm::laplacian(kappa, T)
         ==
          - TRHS
        );

        TEqn.relax();

        TEqn.solve();
        
        T.correctBoundaryConditions();  

        Tcorr = (Tliq - Tsol) * alphaLiquid + Tsol;
        Tcorr.correctBoundaryConditions();  

        alphaLiquid =
            max
            (
                min
                (
                    alphaLiquid + (epsilonRel * fluid.Cp() / Lmelt) * (T - Tcorr),
                    scalar(1)
                ),
                scalar(0)
            );
        
        alphaLiquid.correctBoundaryConditions();  

        residual =
            mag
            (
                alphaLiquid.primitiveField()
              - alphaLiquid.prevIter().primitiveField()
            );

        maxResidual = gMax(residual);

        meanResidual = gSum(residual*mesh.V())/gSum(mesh.V().field());

        Info<< "Correcting alphaLiquid, mean residual = " << meanResidual
            << ", max residual = " << maxResidual
            << endl;
    }
    while
    (
        (
            ++iter < minTCorr
         || (maxResidual > epsilonTol.value())
        )
     && iter <= maxTCorr
    );
}
