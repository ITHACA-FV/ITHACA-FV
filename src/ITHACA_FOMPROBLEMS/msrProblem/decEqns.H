
difft.correctBoundaryConditions();

tmp<fvScalarMatrix> tdecEq1
(
    fvm::div(phi, dec1)
    ==
    fvm::laplacian(turbulence->nu() / Sc + difft, dec1)
    - fvm::Sp(decLam1, dec1)
    + decBeta1* SP* flux
    );

fvScalarMatrix& decEq1 = tdecEq1.ref();
decEq1.relax();
dec1residual = decEq1.solve().initialResidual();



tmp<fvScalarMatrix> tdecEq2
(
    fvm::div(phi, dec2)
    ==
    fvm::laplacian(turbulence->nu() / Sc + difft, dec2)
    - fvm::Sp(decLam2, dec2)
    + decBeta2* SP* flux
    );

fvScalarMatrix& decEq2 = tdecEq2.ref();
decEq2.relax();
dec2residual = decEq2.solve().initialResidual();



tmp<fvScalarMatrix> tdecEq3
(
    fvm::div(phi, dec3)
    ==
    fvm::laplacian(turbulence->nu() / Sc + difft, dec3)
    - fvm::Sp(decLam3, dec3)
    + decBeta3* SP* flux
    );

fvScalarMatrix& decEq3 = tdecEq3.ref();
decEq3.relax();
dec3residual = decEq3.solve().initialResidual();

