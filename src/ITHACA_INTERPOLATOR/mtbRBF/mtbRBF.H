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

Class
    mtbRBF
Description
    RBF interpolation class using the Mathtoolbox library.
SourceFiles
    mtbRBF.C
\*---------------------------------------------------------------------------*/

#ifndef mtbRBF_H
#define mtbRBF_H

#include <Eigen/Core>
#include <memory>
#include <string>
#include <mathtoolbox/rbf-interpolation.hpp>
#include "dictionary.H"
#include "Ostream.H"
#include <cmath>

class mtbRBF
{
public:
    // Constructor
    // kernelType can be "gaussian", "linear", "thin_plate_spline", "cubic"
    mtbRBF(const Foam::word& kernelType = "gaussian",
        bool usePolynomialTerm = true,
        bool useRegularization = true,
        Foam::scalar lambda = 0.001,
        bool normalize = true);

    // Constructor from OpenFOAM dictionary
    mtbRBF(const Foam::dictionary& dict);

    // Destructor
    ~mtbRBF();

    // Fit function
    void fit(const Eigen::MatrixXd& X, const Eigen::VectorXd& y);

    // Predict function
    Foam::scalar predict(const Eigen::VectorXd& x);

    // Predict for multiple points
    Eigen::VectorXd predict(const Eigen::MatrixXd& X);

    // Print model info
    void printInfo();

private:
    std::unique_ptr<mathtoolbox::RbfInterpolator> impl_;
    Foam::word kernelType_;
    bool useRegularization_;
    bool usePolynomialTerm_;
    Foam::scalar lambda_;
    bool normalize_;
    Foam::scalar epsilon_;

    Eigen::VectorXd xMean_;
    Eigen::VectorXd xStd_;
    Foam::scalar yMean_;
    Foam::scalar yStd_;
};

#endif
