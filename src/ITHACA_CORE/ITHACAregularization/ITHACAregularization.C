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
/// Source file of the ITHACAregularization class, it contains the implementation of
/// several methods for regularization.

#include "ITHACAregularization.H"


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
namespace ITHACAregularization
{

Eigen::VectorXd  TSVD(Eigen::MatrixXd A,
                      Eigen::MatrixXd b, int filter)
{
    M_Assert(b.cols() == 1, "The b input in TSVD must have only one column");
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A,
                                          Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::MatrixXd U = svd.matrixU();
    Eigen::MatrixXd V = svd.matrixV();
    Eigen::VectorXd x = Eigen::VectorXd::Zero(b.size());

    for (label i = 0; i < filter; i++)
    {
        x += (U.col(i).dot(b.col(0)) / svd.singularValues()(i)) * (V.col(i));
    }

    return x;
}

Eigen::VectorXd  TSVD(Eigen::MatrixXd A,
                      Eigen::MatrixXd b, double noiseVariance, word parameterMethod)
{
    int filter;

    if (parameterMethod == "DP")
    {
        Info << "\nRegularization parameter selected by Discrepancy principle" << endl;
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(A,
                                              Eigen::ComputeThinU | Eigen::ComputeThinV);
        Eigen::MatrixXd U = svd.matrixU();
        Eigen::VectorXd bVect = b.col(0);
        double min;

        for (int col = 0; col < U.cols() - 1; col++)
        {
            double f = 0;

            for (int i = col + 1; i < U.cols(); i++)
            {
                Eigen::VectorXd tempU = U.col(i);
                f += tempU.col(i).dot(bVect) * tempU.col(i).dot(bVect);
            }

            f += 2 * noiseVariance * col;

            if (col == 0)
            {
                min = f;
                filter = col + 1;
            }
            else if (min > f)
            {
                min = f;
                filter = col + 1;
            }

            Info << "debug : f = " << f << endl;
            Info << "debug : min = " << min << endl;
            Info << "debug : k = " << filter << endl;
        }
    }
    else if (parameterMethod == "UPRE")
    {
        Info << "\nRegularization parameter selected by Discrepancy principle" << endl;
    }
    else
    {
        Info << "Regularization parameter selection methods available are:" << endl
             << "DP, UPRE" << endl;
        exit(1);
    }

    return ITHACAregularization::TSVD(A, b, filter);
}

Eigen::VectorXd  Tikhonov(Eigen::MatrixXd A,
                          Eigen::MatrixXd b, double regularizationParameter)
{
    M_Assert(b.cols() == 1, "The b input in TSVD must have only one column");
    Eigen::MatrixXd Anew = A.transpose() * A + regularizationParameter *
                           Eigen::MatrixXd::Identity(A.rows(), A.cols());
    Eigen::MatrixXd bNew = A.transpose() * b;
    Eigen::VectorXd x = A.inverse() * b;
    std::cout << "x = \n" << x.transpose() << std::endl;
    x = Anew.inverse() * bNew;
    std::cout << "xNew = \n" << x.transpose() << std::endl;
    return x;
}
}
