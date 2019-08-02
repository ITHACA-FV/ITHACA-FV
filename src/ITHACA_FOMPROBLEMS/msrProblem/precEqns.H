
difft.correctBoundaryConditions();

tmp<fvScalarMatrix> tprecEq1
(
    fvm::div(phi, prec1)
    ==
    fvm::laplacian(turbulence->nu() + difft, prec1)
    - fvm::Sp(lam1, prec1)
    + beta1* NSF* flux
    );

fvScalarMatrix& precEq1 = tprecEq1.ref();
precEq1.relax();
//precEq1.solve();
prec1residual = precEq1.solve().initialResidual();

tmp<fvScalarMatrix> tprecEq2
(
    fvm::div(phi, prec2)
    ==
    fvm::laplacian(turbulence->nu() + difft, prec2)
    - fvm::Sp(lam2, prec2)
    + beta2* NSF* flux
    );

fvScalarMatrix& precEq2 = tprecEq2.ref();
precEq2.relax();
//precEq2.solve();
prec2residual = precEq2.solve().initialResidual();


tmp<fvScalarMatrix> tprecEq3
(
    fvm::div(phi, prec3)
    ==
    fvm::laplacian(turbulence->nu() + difft, prec3)
    - fvm::Sp(lam3, prec3)
    + beta3* NSF* flux
    );

fvScalarMatrix& precEq3 = tprecEq3.ref();
precEq3.relax();
// precEq3.solve();
prec3residual = precEq3.solve().initialResidual();


tmp<fvScalarMatrix> tprecEq4
(
    fvm::div(phi, prec4)
    ==
    fvm::laplacian(turbulence->nu() + difft, prec4)
    - fvm::Sp(lam4, prec4)
    + beta4* NSF* flux
    );

fvScalarMatrix& precEq4 = tprecEq4.ref();
precEq4.relax();
//precEq4.solve();
prec4residual = precEq4.solve().initialResidual();


tmp<fvScalarMatrix> tprecEq5
(

    fvm::div(phi, prec5)
    ==
    fvm::laplacian(turbulence->nu() + difft, prec5)
    - fvm::Sp(lam5, prec5)
    + beta5* NSF* flux
    );

fvScalarMatrix& precEq5 = tprecEq5.ref();
precEq5.relax();
//precEq5.solve();
prec5residual = precEq5.solve().initialResidual();


tmp<fvScalarMatrix> tprecEq6
(
    fvm::div(phi, prec6)
    ==
    fvm::laplacian(turbulence->nu() + difft, prec6)
    - fvm::Sp(lam6, prec6)
    + beta6* NSF* flux
    );

fvScalarMatrix& precEq6 = tprecEq6.ref();
precEq6.relax();
//precEq6.solve();
prec6residual = precEq6.solve().initialResidual();


tmp<fvScalarMatrix> tprecEq7
(
    fvm::div(phi, prec7)
    ==
    fvm::laplacian(turbulence->nu() + difft, prec7)
    - fvm::Sp(lam7, prec7)
    + beta7* NSF* flux
    );

fvScalarMatrix& precEq7 = tprecEq7.ref();
precEq7.relax();
//precEq7.solve();
prec7residual = precEq7.solve().initialResidual();


tmp<fvScalarMatrix> tprecEq8
(
    fvm::div(phi, prec8)
    ==
    fvm::laplacian(turbulence->nu() + difft, prec8)
    - fvm::Sp(lam8, prec8)
    + beta8* NSF* flux
    );

fvScalarMatrix& precEq8 = tprecEq8.ref();
precEq8.relax();
//precEq8.solve();
prec8residual = precEq8.solve().initialResidual();


