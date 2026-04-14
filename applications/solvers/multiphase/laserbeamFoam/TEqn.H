{
    int iter = 0;
    scalarField residual(alphaSolid.primitiveField().size(), 1.0);
    scalar maxResidual = 1.0;
    scalar meanResidual = 1.0;

    rhoCp = rho*cp;
    rhoCp.correctBoundaryConditions();  
    
    rhophicp = fvc::interpolate(cp)*rhoPhi;

    do
    {
        //epsilon1.storePrevIter();
        alphaSolid.storePrevIter();
        //T.storePrevIter();

        /*
        fvScalarMatrix TEqn
        (
            fvm::ddt(rhoCp,T)
          + fvm::div(rhophicp, T)
          - fvm::Sp(fvc::ddt(rhoCp) + fvc::div(rhophicp), T)
          - fvm::laplacian(kappa, T)
          - fvm::Sp(Qrad_inpl*gradAlphaDamper, T)
          //- (tau && gradU) // Viscous dissipation
         ==
            laser.deposition()
          - Qv*mag(gradAlpha)*thermalDamper
          + (Qrad_expl - Qv)*gradAlphaDamper
          - TRHS
        );
        */

        // Reduced temperature equation
        fvScalarMatrix TEqn
        (
            fvm::ddt(rhoCp,T)
          + fvm::div(rhophicp, T)
          - fvm::Sp(fvc::ddt(rhoCp) + fvc::div(rhophicp), T)
          - fvm::laplacian(kappa, T)
         ==
            LatentHeat1*(fvc::ddt(rho, alphaSolid) + fvc::div(rhoPhi, alphaSolid)) 
         //- LatentHeat*((fvc::ddt(rho, epsilon1) + fvc::div(rhoPhi, epsilon1)))
        );

        TEqn.relax();

        TEqn.solve();
        
        T.correctBoundaryConditions();
        
        /*
        Tcorr = alpha1 * Tliquidus1 - alphaSolid * (Tliquidus1 - Tsolidus1);
        Tcorr.correctBoundaryConditions();  

        alphaSolid = max
        (
            min
            (
                alphaSolid + epsilonRel*cp/LatentHeat1*(Tcorr - T * alpha1),
                alpha1
            ),
            0.0
        );
        alphaSolid.correctBoundaryConditions(); 
        */

        Tcorr = Tliquidus1 - solidFraction1 * (Tliquidus1 - Tsolidus1);
        //Tcorr.correctBoundaryConditions(); 
        solidFraction1 = 
            max
            (
                min
                (
                    solidFraction1 + epsilonRel * cp / LatentHeat1 * (Tcorr - T),
                    1.0
                ),
                0.0
            );
        alphaSolid = clamp(alpha1, zero_one{}) * solidFraction1;
        alphaSolid.correctBoundaryConditions(); 

        /*
        Tcorr = (TLiquidus - TSolidus)*epsilon1 + TSolidus;
        Tcorr.correctBoundaryConditions();  
        epsilon1 =
            max
            (
                min
                (
                    epsilon1 + (epsilonRel*cp/LatentHeat)*(T - Tcorr),
                    scalar(1)
                ),
                scalar(0)
            );        
        epsilon1.correctBoundaryConditions();  

        residual =
            mag
            (
                epsilon1.primitiveField()
              - epsilon1.prevIter().primitiveField()
            );

        maxResidual = gMax(residual);

        meanResidual = gSum(residual*mesh.V())/gSum(mesh.V().field());

        Info<< "Correcting epsilon1, mean residual = " << meanResidual
            << ", max residual = " << maxResidual
            << endl;
        */

        residual =
            mag
            (
                alphaSolid.primitiveField()
              - alphaSolid.prevIter().primitiveField()
            );

        maxResidual = gMax(residual);

        meanResidual = gSum(residual*mesh.V())/gSum(mesh.V().field());

        Info<< "Correcting alphaSolid, mean residual = " << meanResidual
            << ", max residual = " << maxResidual
            << endl;
    }
    while
    (
        (
            ++iter < minTCorr
         // || (meanResidual > epsilonTol.value())
         || (maxResidual > epsilonTol.value())
        )
     && iter <= maxTCorr
    );

    gradT = fvc::grad(T);
    gradT.correctBoundaryConditions();  // ADDED - CRITICAL!
}
