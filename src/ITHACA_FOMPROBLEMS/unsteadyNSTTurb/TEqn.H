// Solve temperature equation
{
    alphat = turbulence->nut() / Prt;
    alphat.correctBoundaryConditions();
    volScalarField alphaEff("alphaEff",  nu / Pr + alphat);

    fvScalarMatrix TEqn
    (
        fvm::ddt(T)
        + fvm::div(phi, T)
        - fvm::laplacian(alphaEff, T)
    );
    TEqn.relax();
    TEqn.solve();

}
