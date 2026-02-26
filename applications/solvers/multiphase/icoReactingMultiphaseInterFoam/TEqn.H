{
    radiation->correct();
    rhoCp =  rho*fluid.Cp();

    const surfaceScalarField rhoCpPhi(fvc::interpolate(fluid.Cp())*rhoPhi);

    const volScalarField kappaEff
    (
        "kappaEff",
        fluid.kappa() + fluid.Cp()*turbulence->mut()/fluid.Prt()
    );

    fvScalarMatrix TSource
    (
        fluid.heatTransfer(T)
        + radiation->ST(T)
        + fvOptions(rhoCp, T)
    );

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix TEqn
        (
            fvm::ddt(rhoCp, T)
          + fvm::div(rhoCpPhi, T, "div(phi,T)")
          - fvm::Sp(fvc::ddt(rhoCp) + fvc::div(rhoCpPhi), T)
          - fvm::laplacian(kappaEff, T,  "laplacian(kappa,T)")
          == 
            TSource
        );

        TEqn.relax();

        fvOptions.constrain(TEqn);

        TEqn.solve(T.select(pimple.finalNonOrthogonalIter()));

        if (pimple.finalNonOrthogonalIter())
        {            
            T.relax();
            
            fvOptions.correct(T);

            fluid.correct();

            Info<< "min/max(T) = " << min(T).value() << ", " << max(T).value() << endl;
        }
    }
}
