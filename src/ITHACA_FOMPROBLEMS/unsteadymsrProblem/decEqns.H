/*---------------------------------------------------------------------------*\
     ██╗████████╗██╗  ██╗ █████╗  ██████╗ █████╗       ███████╗██╗   ██╗
     ██║╚══██╔══╝██║  ██║██╔══██╗██╔════╝██╔══██╗      ██╔════╝██║   ██║
     ██║   ██║   ███████║███████║██║     ███████║█████╗█████╗  ██║   ██║
     ██║   ██║   ██╔══██║██╔══██║██║     ██╔══██║╚════╝██╔══╝  ╚██╗ ██╔╝
     ██║   ██║   ██║  ██║██║  ██║╚██████╗██║  ██║      ██║      ╚████╔╝
     ╚═╝   ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝      ╚═╝       ╚═══╝

 * In real Time Highly Advanced Computational Applications for Finite Volumes
 * Copyright (C) 2017 by the ITHACA-FV authors
-------------------------------------------------------------------------------

  License
  This file is part of ITHACA-FV

  ITHACA-FV is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ITHACA-FV is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with ITHACA-FV. If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/
difft.correctBoundaryConditions();

fvScalarMatrix decEq1
(
    fvm::ddt(dec1)
    ==
    - fvm::div(phi, dec1)
    + fvm::laplacian(turbulence->nu() / Sc + difft, dec1)
    - fvm::Sp(decLam1, dec1)
    + decBeta1* SP* flux
    );

decEq1.relax();
decEq1.solve();




fvScalarMatrix decEq2
(
    fvm::ddt(dec2)
    ==
    - fvm::div(phi, dec2)
    + fvm::laplacian(turbulence->nu() / Sc + difft, dec2)
    - fvm::Sp(decLam2, dec2)
    + decBeta2* SP* flux
    );

decEq2.relax();
decEq2.solve();




fvScalarMatrix decEq3
(
    fvm::ddt(dec3)
    ==
    - fvm::div(phi, dec3)
    + fvm::laplacian(turbulence->nu() / Sc + difft, dec3)
    - fvm::Sp(decLam3, dec3)
    + decBeta3* SP* flux
    );

decEq3.relax();
decEq3.solve();


