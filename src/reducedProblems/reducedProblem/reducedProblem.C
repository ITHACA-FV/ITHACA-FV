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

#include "reducedProblem.H"

// ******************** //
// class reducedProblem //
// ******************** //

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructor

reducedProblem::reducedProblem()
{
}

reducedProblem::reducedProblem(reductionProblem& problem)
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
             "Matrix 'mu_interp' must have same number and order of columns (i.e. parametes) as the matrix 'mu_samples'.");
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


