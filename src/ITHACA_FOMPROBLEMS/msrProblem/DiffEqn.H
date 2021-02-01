
tmp<fvScalarMatrix> tNEqn
(

    - fvm::laplacian(D, flux)
    ==
    fvm::Sp(NSF * (1 - betaTot), flux)
    - fvm::Sp(A, flux)
    + lam1* prec1
    + lam2* prec2
    + lam3* prec3
    + lam4* prec4
    + lam5* prec5
    + lam6* prec6
    + lam7* prec7
    + lam8* prec8

    );

fvScalarMatrix& NEqn = tNEqn.ref();
NEqn.relax();
//NEqn.solve();
fluxresidual = NEqn.solve().initialResidual();




