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

/// \file
/// Source file of the laplacianProblem class.


#include "laplacianProblem.H"

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructors
laplacianProblem::laplacianProblem() {}
laplacianProblem::laplacianProblem(int argc, char *argv[])
{
#include "setRootCase.H"
#include "createTime.H"
#include "createMesh.H"
#include "createFields.H"
}


// * * * * * * * * * * * * * * Full Order Methods * * * * * * * * * * * * * * //
// Method to performa a truthSolve
void laplacianProblem::truthSolve()
{
  volScalarField& T = _T();
  volScalarField& S = _S();
  fvScalarMatrix lhs(theta[0]*operator_list[0]);
  for (int i = 1; i < operator_list.size(); i++)
  {
    lhs += theta[i] * operator_list[i];
  }
  solve(lhs == -S);
  exportSolution(T, name(counter), "./ITHACAoutput/Offline/");
  Tfield.append(T);
  counter++;
}
// Perform the projection onto the POD modes
void laplacianProblem::project(label Nmodes)
{
  NTmodes = Nmodes;
  A_matrices.resize(operator_list.size());
  source.resize(Nmodes, 1);
  volScalarField& S = _S();
  for (int i = 0; i < operator_list.size(); i++)
  {
    A_matrices[i].resize(Nmodes, Nmodes);
    for (int j = 0; j < Nmodes; j++)
    {
      if (i == 0)
      {
        source(j, 0) = fvc::domainIntegrate( Tmodes[j] * S).value();
      }
      for (int k = 0; k < Nmodes; k++)
      {
        A_matrices[i](j, k) = fvc::domainIntegrate( Tmodes[j] * fvc::laplacian(nu_list[i], Tmodes[k])).value();
      }
    }
  }
  /// Export the A matrices
  ITHACAstream::exportMatrix(A_matrices, "A", "python", "./ITHACAoutput/Matrices/");
  ITHACAstream::exportMatrix(A_matrices, "A", "matlab", "./ITHACAoutput/Matrices/");
  ITHACAstream::exportMatrix(A_matrices, "A", "eigen", "./ITHACAoutput/Matrices/A_matrices");
  /// Export the source term
  ITHACAstream::exportMatrix(source, "S", "python", "./ITHACAoutput/Matrices/");
  ITHACAstream::exportMatrix(source, "S", "matlab", "./ITHACAoutput/Matrices/");
  ITHACAstream::exportMatrix(source, "S", "eigen", "./ITHACAoutput/Matrices/");

}




