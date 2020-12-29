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
laplacianProblem::laplacianProblem(int argc, char* argv[])
{
    _args = autoPtr<argList>
            (
                new argList(argc, argv)
            );

    if (!_args->checkRootCase())
    {
        Foam::FatalError.exit();
    }

    argList& args = _args();
#include "createTime.H"
#include "createMesh.H"
#include "createFields.H"
    offline = ITHACAutilities::check_off();
    podex = ITHACAutilities::check_pod();
}


// * * * * * * * * * * * * * * Full Order Methods * * * * * * * * * * * * * * //
// Method to performa a truthSolve
void laplacianProblem::truthSolve(List<scalar> mu_now, word folder)
{
    volScalarField& T = _T();
    volScalarField& S = _S();
    fvScalarMatrix lhs(theta[0]*operator_list[0]);

    for (label i = 1; i < operator_list.size(); i++)
    {
        lhs += theta[i] * operator_list[i];
    }

    solve(lhs == -S);
    ITHACAstream::exportSolution(T, name(counter), folder);
    Tfield.append(tmp<volScalarField>(T));
    counter++;
    writeMu(mu_now);
    // --- Fill in the mu_samples with parameters (mu) to be used for the PODI sample points
    mu_samples.conservativeResize(mu_samples.rows() + 1, mu_now.size());

    for (label i = 0; i < mu_now.size(); i++)
    {
        mu_samples(mu_samples.rows() - 1, i) = mu_now[i];
    }

    // Resize to Unitary if not initialized by user (i.e. non-parametric problem)
    if (mu.cols() == 0)
    {
        mu.resize(1, 1);
    }

    if (mu_samples.rows() == mu.cols())
    {
        ITHACAstream::exportMatrix(mu_samples, "mu_samples", "eigen",
                                   folder);
    }
}
// Perform the projection onto the POD modes
void laplacianProblem::project(label Nmodes)
{
    NTmodes = Nmodes;
    A_matrices.resize(operator_list.size());
    source.resize(Nmodes, 1);
    volScalarField& S = _S();

    for (label i = 0; i < operator_list.size(); i++)
    {
        A_matrices[i].resize(Nmodes, Nmodes);

        for (label j = 0; j < Nmodes; j++)
        {
            if (i == 0)
            {
                source(j, 0) = fvc::domainIntegrate( Tmodes[j] * S).value();
            }

            for (label k = 0; k < Nmodes; k++)
            {
                A_matrices[i](j, k) = fvc::domainIntegrate( Tmodes[j] * fvc::laplacian(
                                          nu_list[i], Tmodes[k])).value();
            }
        }
    }

    /// Export the A matrices
    ITHACAstream::exportMatrix(A_matrices, "A", "python",
                               "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(A_matrices, "A", "matlab",
                               "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(A_matrices, "A", "eigen",
                               "./ITHACAoutput/Matrices/A_matrices");
    /// Export the source term
    ITHACAstream::exportMatrix(source, "S", "python", "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(source, "S", "matlab", "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(source, "S", "eigen", "./ITHACAoutput/Matrices/");
}




