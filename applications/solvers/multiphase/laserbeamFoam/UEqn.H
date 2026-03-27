MRF.correctBoundaryVelocity(U);

DC = DClarge*Foam::pow((1.0 - alphaLiquid), 2) /
     (Foam::pow(alphaLiquid, 3) + DCsmall);
DC.correctBoundaryConditions();

if (damperSwitch)
{
    // To account for high density difference
    damper = fvc::interpolate(rho * (2.0 / (rho1 + rho2)));
}

pVap = 0.54 * p0 * Foam::exp(Lvap * Mm * ((T - Tvap) / (R * T * Tvap)));
pVap.correctBoundaryConditions();  // ADDED

fvVectorMatrix UEqn
(
    fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
    + MRF.DDt(rho, U)
    + turbulence->divDevRhoReff(rho, U)
    + fvm::Sp(DC, U)
    ==
    fvOptions(rho, U)
);

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve
    (
        UEqn
        ==
        fvc::reconstruct
        (
            (
                fluid.surfaceTensionForce() * damper
                + fvc::interpolate(pVap) * fvc::snGrad(alpha1) * damper
                - fvc::snGrad(p_rgh)
            ) * mesh.magSf()
        )
    );

    fvOptions.correct(U);
}
