{
    volScalarField& he = thermo.he();

    fvScalarMatrix EEqn
    (
        fvm::div(phi, he)
        + (
            he.name() == "e"
            ? fvc::div(phi, volScalarField("Ekp", 0.5 * magSqr(U) + p / rho))
            : fvc::div(phi, volScalarField("K", 0.5 * magSqr(U)))
        )
        - fvm::laplacian(turbulence->alphaEff(), he)
        ==
        fvOptions(rho, he)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    eresidual = EEqn.solve().initialResidual();

    fvOptions.correct(he);

    thermo.correct(); // Here are calculated both temperature and density based on P,U and he.
}
