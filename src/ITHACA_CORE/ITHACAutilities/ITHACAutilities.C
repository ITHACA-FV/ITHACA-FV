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

#include "ITHACAutilities.H"
#include "ITHACAstream.H"
#include "ITHACAparameters.H"
#include "turbulentTransportModel.H"

/// \file
/// Source file of the ITHACAutilities namespace.

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

namespace ITHACAutilities
{

Eigen::MatrixXd rand(label rows, label cols, double min,
                     double max)
{
    std::srand(static_cast<long unsigned int>
               (std::chrono::high_resolution_clock::now().time_since_epoch().count()));
    Eigen::MatrixXd matr = Eigen::MatrixXd::Random(rows, cols);
    matr = (matr.array() + 1) / 2;
    matr = matr.array() * (max - min);
    matr = matr.array() + min;
    return matr;
}

Eigen::MatrixXd rand(label rows, Eigen::MatrixXd minMax)
{
    std::srand(static_cast<long unsigned int>
               (std::chrono::high_resolution_clock::now().time_since_epoch().count()));
    label cols = minMax.rows();
    Eigen::MatrixXd matr = Eigen::MatrixXd::Random(rows, cols);
    matr = (matr.array() + 1) / 2;

    for (label i = 0; i < cols; i++)
    {
        matr.col(i) = matr.col(i).array() * (minMax(i, 1) - minMax(i, 0));
        matr.col(i) = matr.col(i).array() + (minMax(i, 0));
    }

    return matr;
}


bool isInteger(double ratio)
{
    bool checkResult = 0;

    if (abs(round(ratio) - ratio) < std::sqrt(SMALL))
    {
        checkResult = true;
    }
    else
    {
        checkResult = false;
    }

    return checkResult;
}

bool isTurbulent()
{
    bool checkTurb;
    ITHACAparameters* para = ITHACAparameters::getInstance();
    auto& tur =
        para->mesh.lookupObject<incompressible::turbulenceModel>("turbulenceProperties");

    if (tur.type() == "Stokes" || tur.type() == "Maxwell"
            || tur.type() == "laminarModel")
    {
        checkTurb = false;
    }
    else
    {
        checkTurb = true;
    }

    return checkTurb;
}

template<typename T>
List<T> combineList(List<List<T>> & doubleList)
{
    List<T> a = ListListOps::combine<List<T>>(doubleList,
                accessOp<List<T>>());
#if OPENFOAM >= 1812
    inplaceUniqueSort(a);
#else
    labelList order;
    uniqueOrder(a, order);
    List<T> b(order.size());

    for (label i = 0; i < order.size(); ++i)
    {
        b[i] = a[order[i]];
    }

    a.resize(order.size());
    a = b;
#endif
    return a;
}

template List<label> combineList(List<List<label>> & doubleList);


// Using the Eigen library, using the SVD decomposition method to solve the
// matrix pseudo-inverse, the default error er is 0
Eigen::MatrixXd pinv_eigen_based(Eigen::MatrixXd& origin, const float er)
{
    // perform svd decomposition
    Eigen::JacobiSVD<Eigen::MatrixXd> svd_holder(origin,
            Eigen::ComputeThinU | Eigen::ComputeThinV);
    // Build SVD decomposition results
    Eigen::MatrixXd U = svd_holder.matrixU();
    Eigen::MatrixXd V = svd_holder.matrixV();
    Eigen::MatrixXd D = svd_holder.singularValues();
    // Build the S matrix
    Eigen::MatrixXd S(V.cols(), U.cols());
    S.setZero();

    for (unsigned int i = 0; i < D.size(); ++i)
    {
        if (D(i, 0) > er)
        {
            S(i, i) = 1 / D(i, 0);
        }
        else
        {
            S(i, i) = 0;
        }
    }

    return V* S* U.transpose();
}


Eigen::MatrixXd invertMatrix(Eigen::MatrixXd& matrixToInvert,
                             const word inversionMethod)
{
    Info << "Inversion method : " << inversionMethod << endl;

    if (inversionMethod == "pinv_eigen_based")
    {
        return pinv_eigen_based(matrixToInvert);
    }
    else if (inversionMethod == "direct")
    {
        return matrixToInvert.inverse();
    }
    else if (inversionMethod == "fullPivLu")
    {
        return matrixToInvert.fullPivLu().inverse();
    }
    else if (inversionMethod == "partialPivLu")
    {
        return matrixToInvert.partialPivLu().inverse();
    }
    else if (inversionMethod == "householderQr")
    {
        return matrixToInvert.householderQr().solve(Eigen::MatrixXd::Identity(
                    matrixToInvert.rows(), matrixToInvert.cols()));
    }
    else if (inversionMethod == "colPivHouseholderQr")
    {
        return matrixToInvert.colPivHouseholderQr().inverse();
    }
    else if (inversionMethod == "fullPivHouseholderQr")
    {
        return matrixToInvert.fullPivHouseholderQr().inverse();
    }
    else if (inversionMethod == "completeOrthogonalDecomposition")
    {
        return matrixToInvert.completeOrthogonalDecomposition().pseudoInverse();
    }
    else if (inversionMethod == "jacobiSvd")
    {
        return matrixToInvert.jacobiSvd(Eigen::ComputeThinU |
                                        Eigen::ComputeThinV).solve(Eigen::MatrixXd::Identity(matrixToInvert.rows(),
                                                matrixToInvert.cols()));
    }
    else if (inversionMethod == "llt")
    {
        return matrixToInvert.llt().solve(Eigen::MatrixXd::Identity(
                                              matrixToInvert.rows(), matrixToInvert.cols()));
    }
    else if (inversionMethod == "ldlt")
    {
        Eigen::LLT<Eigen::MatrixXd> lltOfA(matrixToInvert);
        return lltOfA.solve(Eigen::MatrixXd::Identity(matrixToInvert.rows(),
                            matrixToInvert.cols()));
    }
    else if (inversionMethod == "bdcSvd")
    {
        return matrixToInvert.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(
                   Eigen::MatrixXd::Identity(matrixToInvert.rows(), matrixToInvert.cols()));
    }
    else
    {
        Info << "Unkwown inversion method, solving with : completeOrthogonalDecomposition"
             << endl;
        return matrixToInvert.completeOrthogonalDecomposition().pseudoInverse();
    }
}


}
