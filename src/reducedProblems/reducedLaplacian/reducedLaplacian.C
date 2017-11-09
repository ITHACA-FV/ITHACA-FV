/*---------------------------------------------------------------------------*\
Copyright (C) 2017 by the ITHACA-FV authors

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
/// Source file of the reducedLaplacian class

#include "reducedLaplacian.H"


// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructor
reducedLaplacian::reducedLaplacian()
{

}

reducedLaplacian::reducedLaplacian(laplacianProblem& problem)
{
    source = problem.source;
    NTmodes = problem.NTmodes;
    A_matrices = problem.A_matrices;
    Tmodes = problem.Tmodes;
}

void reducedLaplacian::solveOnline(Eigen::MatrixXd mu)
{
    if (mu.cols() != A_matrices.size())
    {
        Info << "wrong dimension of online parameters" << endl;
        exit(0);
    }
    Eigen::MatrixXd A;
    A.setZero(NTmodes, NTmodes);
    for (int i = 0; i < A_matrices.size() ; i++)
    {
        A += A_matrices[i] * mu(0, i);
    }
    Eigen::MatrixXd x;
    x = A.colPivHouseholderQr().solve(-source);
    online_solution.conservativeResize(count_online_solve, NTmodes + 1);
    online_solution(count_online_solve - 1, 0) = count_online_solve;
    online_solution.row(count_online_solve - 1).tail(NTmodes) = x.transpose();
    count_online_solve += 1;
}

void reducedLaplacian::reconstruct(laplacianProblem& problem, fileName folder, int printevery)
{
    mkDir(folder);
    system("ln -s ../../constant " + folder + "/constant");
    system("ln -s ../../0 " + folder + "/0");
    system("ln -s ../../system " + folder + "/system");

    int counter = 0;
    int nextwrite = 0;
    int counter2 = 1;

    for (label i = 0; i < online_solution.rows(); i++)
    {
        if (counter == nextwrite)
        {
            volScalarField T_rec("T_rec", Tmodes[0] * 0);
            for (label j = 0; j < NTmodes; j++)
            {
                T_rec += Tmodes[j] * online_solution(i, j + 1);
            }
            problem.exportSolution(T_rec, name(online_solution(i, 0)), folder);

            nextwrite += printevery;
            counter2 ++;

        }
        counter++;
    }
}




// ************************************************************************* //

