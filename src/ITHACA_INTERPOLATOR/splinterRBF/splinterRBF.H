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
    splinterRBF
Description
    RBF interpolation class using the SPLINTER library.
SourceFiles
    splinterRBF.C
\*---------------------------------------------------------------------------*/

#ifndef splinterRBF_H
#define splinterRBF_H

#include <Eigen/Core>
#include <memory>
#include <string>
#include <rbfspline.h>
#include <datatable.h>
#include "dictionary.H"
#include <iostream>
#include <Ostream.H>
#include <cmath>

class splinterRBF
{
public:
    // Constructor
    // kernelType can be "gaussian", "thin_plate_spline", "multiquadric", "inverse_multiquadric", "inverse_quadratic"
    splinterRBF(const Foam::word& kernelType = "gaussian",
        bool normalize = true,
        Foam::scalar epsilon = 1.0);

    // Constructor from OpenFOAM dictionary
    splinterRBF(const Foam::dictionary& dict);

    // Destructor
    ~splinterRBF();

    // Fit function
    void fit(const Eigen::MatrixXd& X, const Eigen::VectorXd& y);

    // Predict function
    Foam::scalar predict(const Eigen::VectorXd& x);

    // Predict for multiple points
    Eigen::VectorXd predict(const Eigen::MatrixXd& X);

    // Print model info
    void printInfo();

private:
    std::unique_ptr<SPLINTER::RBFSpline> impl_;
    Foam::word kernelType_;
    bool normalize_;
    Foam::scalar epsilon_;

    Eigen::VectorXd xMean_;
    Eigen::VectorXd xStd_;
    Foam::scalar yMean_;
    Foam::scalar yStd_;

    SPLINTER::RadialBasisFunctionType getKernelType(const Foam::word& kernelType);
};

#endif
