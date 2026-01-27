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
    ithacaInterpolator Declaration
Description
    Radial Basis Function Interpolator interface class that wraps different RBF packages.
Source files
    ithacaInterpolator.C
\*---------------------------------------------------------------------------*/

#ifndef ithacaInterpolator_H
#define ithacaInterpolator_H

#include <Eigen/Core>
#include <memory>
#include <string>
#include "dictionary.H"

// Forward declarations
class mtbGPR;
class mtbRBF;
class splinterRBF;

class ithacaInterpolator
{
public:
    // Constructor from OpenFOAM dictionary
    ithacaInterpolator(const Foam::dictionary& dict);

    // Destructor
    ~ithacaInterpolator();

    // Fit function
    void fit(const Eigen::MatrixXd& X, const Eigen::VectorXd& y);

    // Predict function
    double predict(const Eigen::VectorXd& x);

    // Predict for multiple points
    Eigen::VectorXd predict(const Eigen::MatrixXd& X);

    // Print model information
    void printInfo() const;

private:
    std::unique_ptr<mtbGPR> mtbGPR_;
    std::unique_ptr<mtbRBF> mtb_;
    std::unique_ptr<splinterRBF> splinter_;
    std::string algorithm_;
    std::string package_;
};

#endif
