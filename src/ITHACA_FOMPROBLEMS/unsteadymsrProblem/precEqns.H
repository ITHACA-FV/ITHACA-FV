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

fvScalarMatrix precEq1
(
    fvm::ddt(prec1)
    ==
    - fvm::div(phi, prec1)
    + fvm::laplacian(turbulence->nu() / Sc + difft, prec1)
    - fvm::Sp(lam1, prec1)
    + beta1* NSF* flux
    );

precEq1.relax();
precEq1.solve();


fvScalarMatrix precEq2
(
    fvm::ddt(prec2)
    ==
    - fvm::div(phi, prec2)
    + fvm::laplacian(turbulence->nu() / Sc + difft, prec2)
    - fvm::Sp(lam2, prec2)
    + beta2* NSF* flux
    );

precEq2.relax();
precEq2.solve();


fvScalarMatrix precEq3
(
    fvm::ddt(prec3)
    ==
    - fvm::div(phi, prec3)
    + fvm::laplacian(turbulence->nu() / Sc + difft, prec3)
    - fvm::Sp(lam3, prec3)
    + beta3* NSF* flux
    );

precEq3.relax();
precEq3.solve();


fvScalarMatrix precEq4
(
    fvm::ddt(prec4)
    ==
    - fvm::div(phi, prec4)
    + fvm::laplacian(turbulence->nu() / Sc + difft, prec4)
    - fvm::Sp(lam4, prec4)
    + beta4* NSF* flux
    );

precEq4.relax();
precEq4.solve();


fvScalarMatrix precEq5
(
    fvm::ddt(prec5)
    ==
    - fvm::div(phi, prec5)
    + fvm::laplacian(turbulence->nu() / Sc + difft, prec5)
    - fvm::Sp(lam5, prec5)
    + beta5* NSF* flux
    );

precEq5.relax();
precEq5.solve();


fvScalarMatrix precEq6
(
    fvm::ddt(prec6)
    ==
    - fvm::div(phi, prec6)
    + fvm::laplacian(turbulence->nu() / Sc + difft, prec6)
    - fvm::Sp(lam6, prec6)
    + beta6* NSF* flux
    );

precEq6.relax();
precEq6.solve();


fvScalarMatrix precEq7
(
    fvm::ddt(prec7)
    ==
    - fvm::div(phi, prec7)
    + fvm::laplacian(turbulence->nu() / Sc + difft, prec7)
    - fvm::Sp(lam7, prec7)
    + beta7* NSF* flux
    );

precEq7.relax();
precEq7.solve();


fvScalarMatrix precEq8
(
    fvm::ddt(prec8)
    ==
    - fvm::div(phi, prec8)
    + fvm::laplacian(turbulence->nu() / Sc + difft, prec8)
    - fvm::Sp(lam8, prec8)
    + beta8* NSF* flux
    );

precEq8.relax();
precEq8.solve();


