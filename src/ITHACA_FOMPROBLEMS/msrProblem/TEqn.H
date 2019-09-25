
alphat.correctBoundaryConditions();


tmp<fvScalarMatrix> tTEqn
(

    fvm::div(phi, T)
    ==
    fvm::laplacian(turbulence->nu() / Pr + alphat, T)
    +(v / CpRef) * (decLam1 * dec1 + decLam2 * dec2 + decLam3* dec3)
    + (SP * (1 - decbetaTot)*v / CpRef)*flux

    );
fvScalarMatrix& TEqn = tTEqn.ref();
//TEqn.relax();
//TEqn.solve();
//Tresidual=TEqn.solve().initialResidual();





