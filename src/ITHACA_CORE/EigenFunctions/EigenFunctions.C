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

#include "EigenFunctions.H"

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void EigenFunctions::sortEigenvalues(Eigen::VectorXd& eigenvalues,
                                     Eigen::MatrixXd& eigenvectors)
{
    labelList order;
    scalarField eigenValues(eigenvalues.size());

    for (label i = 0; i < eigenvalues.size(); i++)
    {
        eigenValues[i] = eigenvalues(i);
    }

    sortedOrder(eigenValues, order);
    scalarField eigenValues2(eigenValues);

    for (label i = 0; i < order.size(); i++)
    {
        eigenvalues(i) = eigenValues[order[order.size() - i - 1]];
    }

    Eigen::MatrixXd eigenvectors2 = eigenvectors;

    for (label i = 0; i < eigenvalues.size(); i++)
    {
        for (label k = 0; k < eigenvalues.size(); k++)
        {
            eigenvectors2(i, k) = eigenvectors(k, order[order.size() - i - 1]);
        }
    }

    eigenvectors = eigenvectors2;
}
