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

namespace EigenFunctions
{
void sortEigenvalues(Eigen::VectorXd& eigenvalues,
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


Eigen::VectorXd ExpSpaced(double first, double last, int n)
{
    Eigen::VectorXd vector(n);
    double m = (double) 1 / (n * 1.0 - 1);
    double quotient = std::pow(last / first, m);
    vector(0) = first;

    for (int i = 1; i < n; i++)
    {
        vector(i) = vector(i - 1) * quotient;
    }

    return vector;
}

template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
vectorTensorProduct(const Eigen::Matrix<T, Eigen::Dynamic, 1>&
                    g,
                    const Eigen::Tensor<T, 3 >& c,
                    const Eigen::Matrix<T, Eigen::Dynamic, 1>& a)
{
    int prodDim = c.dimension(0);
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> prod;
    prod.resize(prodDim, 1);

    for (int i = 0; i < prodDim; i++)
    {
        prod(i, 0) = g.transpose() *
                     SliceFromTensor(c, 0, i) * a;
    }

    return prod;
}

template Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
vectorTensorProduct<>(
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& g,
    const Eigen::Tensor<double, 3 >& c,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& a);

template Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>
vectorTensorProduct(
    const Eigen::Matrix<int, Eigen::Dynamic, 1>& g, const Eigen::Tensor<int, 3 >& c,
    const Eigen::Matrix<int, Eigen::Dynamic, 1>& a);

template Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>
vectorTensorProduct(
    const Eigen::Matrix<float, Eigen::Dynamic, 1>& g,
    const Eigen::Tensor<float, 3 >& c,
    const Eigen::Matrix<float, Eigen::Dynamic, 1>& a);
}
