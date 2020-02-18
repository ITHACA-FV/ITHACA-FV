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
/// Source file of the reducedProblem class.

#include "ReducedProblem.H"

// ******************** //
// class reducedProblem //
// ******************** //

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructor

reducedProblem::reducedProblem()
{
}

reducedProblem::reducedProblem(reductionProblem& problem)
    :
    problem(&problem)
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void reducedProblem::solveOnline()
{
    Info << "The method reducedProblem::solveOnline in reducedProblem.C is a virtual method"
         << endl;
    Info << "It must be overridden, exiting the code" << endl;
    exit(0);
}

Eigen::MatrixXd reducedProblem::solveLinearSys(List<Eigen::MatrixXd> LinSys,
        Eigen::MatrixXd x, Eigen::VectorXd& residual, const Eigen::MatrixXd& bc,
        const std::string solverType)

{
    M_Assert(solverType == "fullPivLu" || solverType == "partialPivLu"
             || solverType == "householderQR" || solverType == "colPivHouseholderQR"
             || solverType == "fullPivHouseholderQR"
             || solverType == "CompleteOrthogonalDecomposition" || solverType == "llt"
             || solverType == "ldlt" || solverType == "bdcSvd"
             || solverType == "jacobiSvd", "solver not defined");
    Eigen::MatrixXd y;

    for (label i = 0; i < bc.size(); i++)
    {
        LinSys[0].row(i) *= 0;
        LinSys[0](i, i) = 1;
        LinSys[1](i, 0) = bc(i);
    }

    residual = LinSys[0] * x - LinSys[1];

    if (solverType == "fullPivLu")
    {
        y = LinSys[0].fullPivLu().solve(LinSys[1]);
    }
    else if (solverType == "partialPivLu")
    {
        y = LinSys[0].partialPivLu().solve(LinSys[1]);
    }
    else if (solverType == "householderQr")
    {
        y = LinSys[0].householderQr().solve(LinSys[1]);
    }
    else if (solverType == "colPivHouseholderQr")
    {
        y = LinSys[0].colPivHouseholderQr().solve(LinSys[1]);
    }
    else if (solverType == "fullPivHouseholderQr")
    {
        y = LinSys[0].fullPivHouseholderQr().solve(LinSys[1]);
    }
    else if (solverType == "completeOrthogonalDecomposition")
    {
        y = LinSys[0].completeOrthogonalDecomposition().solve(LinSys[1]);
    }
    else if (solverType == "llt")
    {
        y = LinSys[0].llt().solve(LinSys[1]);
    }
    else if (solverType == "ldlt")
    {
        y = LinSys[0].ldlt().solve(LinSys[1]);
    }
    else if (solverType == "bdcSvd")
    {
        y = LinSys[0].bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(
                LinSys[1]);
    }
    else if (solverType == "jacobiSvd")
    {
        y = LinSys[0].jacobiSvd().solve(LinSys[1]);
    }

    return y;
}

Eigen::MatrixXd reducedProblem::solveLinearSys(List<Eigen::MatrixXd> LinSys,
        Eigen::MatrixXd x, Eigen::VectorXd& residual, const std::string solverType)
{
    const Eigen::MatrixXd& bc = Eigen::MatrixXd::Zero(0, 0);
    Eigen::MatrixXd y = reducedProblem::solveLinearSys(LinSys, x, residual, bc,
                        solverType);
    return y;
}

// ****************** //
// class onlineInterp //
// ****************** //

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructor

onlineInterp::onlineInterp()
{
    Nmu_samples =
        ITHACAstream::readMatrix("./ITHACAoutput/Offline/mu_samples_mat.txt").cols();
}

Eigen::MatrixXd onlineInterp::getInterpCoeffRBF(std::vector<SPLINTER::RBFSpline>
        rbfVec, Eigen::MatrixXd mu_interp)
{
    M_Assert(mu_interp.cols() == Nmu_samples,
             "Matrix 'mu_interp' must have same number and order of columns (i.e. parameters) as the matrix 'mu_samples'.");
    label Nmodes = rbfVec.size();
    // Nsamples (snapshots) to evaluate the interpolator at
    label Nsamples = mu_interp.rows();
    Eigen::MatrixXd muInterpRefined(mu_interp.rows(), 0);

    // Discard columns of mu_interp that have constant value for all elements, and use a refined matrix muInterpRefined
    for (label i = 0; i < mu_interp.cols(); i++)
    {
        double firstVal = mu_interp(0, i);
        bool isConst = (mu_interp.col(i).array() == firstVal).all();

        // Consider only columns with non-constant values
        if (isConst == 0)
        {
            muInterpRefined.conservativeResize(muInterpRefined.rows(),
                                               muInterpRefined.cols() + 1);
            muInterpRefined.col(muInterpRefined.cols() - 1) = mu_interp.col(i);
        }
    }

    label Nmu = muInterpRefined.cols();
    Eigen::MatrixXd coeff_interp(Nmodes, Nsamples);

    for (label i = 0; i < Nmodes; i++)
    {
        for (label j = 0; j < Nsamples; j++)
        {
            SPLINTER::DenseVector x(Nmu);

            for (label k = 0; k < Nmu; k++)
            {
                x(k) = muInterpRefined(j, k);
            }

            coeff_interp(i, j) = rbfVec[i].eval(x);
        }
    }

    return coeff_interp;
}


Eigen::MatrixXd onlineInterp::getInterpCoeffSPL(std::vector<SPLINTER::BSpline>
        splVec, Eigen::MatrixXd mu_interp)
{
    M_Assert(mu_interp.cols() == Nmu_samples,
             "Matrix 'mu_interp' must have same number and order of columns (i.e. parametes) as the matrix 'mu_samples'.");
    label Nmodes = splVec.size();
    // Nsamples (snapshots) to evaluate the interpolator at
    label Nsamples = mu_interp.rows();
    Eigen::MatrixXd muInterpRefined(mu_interp.rows(), 0);

    // Discard columns of mu_interp that have constant value for all elements, and use a refined matrix muInterpRefined
    for (label i = 0; i < mu_interp.cols(); i++)
    {
        double firstVal = mu_interp(0, i);
        bool isConst = (mu_interp.col(i).array() == firstVal).all();

        // Consider only columns with non-constant values
        if (isConst == 0)
        {
            muInterpRefined.conservativeResize(muInterpRefined.rows(),
                                               muInterpRefined.cols() + 1);
            muInterpRefined.col(muInterpRefined.cols() - 1) = mu_interp.col(i);
        }
    }

    label Nmu = muInterpRefined.cols();
    Eigen::MatrixXd coeff_interp(Nmodes, Nsamples);

    for (label i = 0; i < Nmodes; i++)
    {
        for (label j = 0; j < Nsamples; j++)
        {
            SPLINTER::DenseVector x(Nmu);

            for (label k = 0; k < Nmu; k++)
            {
                x(k) = muInterpRefined(j, k);
            }

            coeff_interp(i, j) = splVec[i].eval(x);
        }
    }

    return coeff_interp;
}

// ************************************************************************* //


