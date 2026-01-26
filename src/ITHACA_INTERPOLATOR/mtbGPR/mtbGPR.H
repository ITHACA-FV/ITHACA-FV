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
    mtbGPR
Description
    Gaussian Process Regression wrapper around the Mathtoolbox implementation.
SourceFiles
    mtbGPR.C
\*---------------------------------------------------------------------------*/

#ifndef mtbGPR_H
#define mtbGPR_H

#include <Eigen/Core>
#include <memory>
#include <mathtoolbox/gaussian-process-regression.hpp>
#include "dictionary.H"
#include "Ostream.H"

class mtbGPR
{
public:
    // Constructor from OpenFOAM dictionary
    explicit mtbGPR(const Foam::dictionary& dict);

    // Destructor
    ~mtbGPR();

    // Fit function
    void fit(const Eigen::MatrixXd& X, const Eigen::VectorXd& y);

    // Predict function (mean only)
    Foam::scalar predict(const Eigen::VectorXd& x);

    // Predict for multiple points
    Eigen::VectorXd predict(const Eigen::MatrixXd& X);

    // Print model info
    void printInfo() const;

private:
    mathtoolbox::GaussianProcessRegressor::KernelType parseKernelType(const Foam::word& kernelWord) const;

    std::unique_ptr<mathtoolbox::GaussianProcessRegressor> impl_;
    Foam::word kernelTypeWord_;
    bool useDataNormalization_;
    bool optimizeHyperparams_;
    Foam::scalar kernelScale_;
    Foam::scalar lengthScale_;
    Foam::scalar noise_; // observation noise hyperparameter
};

#endif
