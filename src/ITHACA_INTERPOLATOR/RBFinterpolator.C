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

#include "RBFinterpolator.H"
#include "mtbRBF.H"
#include "splinterRBF.H"
#include "error.H"

RBFinterpolator::RBFinterpolator(const Foam::dictionary& dict)
{
    package_ = dict.lookupOrDefault<Foam::word>("package", "mathtoolbox");

    if (package_ == "mathtoolbox")
    {
        mtb_ = std::make_unique<mtbRBF>(dict);
    }
    else if (package_ == "splinter")
    {
        splinter_ = std::make_unique<splinterRBF>(dict);
    }
    else
    {
        FatalErrorInFunction
            << "Unknown package: " << package_
            << ". Valid options are: mathtoolbox, splinter"
            << Foam::exit(Foam::FatalError);
    }
}

RBFinterpolator::~RBFinterpolator() = default;

void RBFinterpolator::fit(const Eigen::MatrixXd& X, const Eigen::VectorXd& y)
{
    if (package_ == "mathtoolbox")
    {
        mtb_->fit(X, y);
    }
    else if (package_ == "splinter")
    {
        splinter_->fit(X, y);
    }
    else
    {
        FatalErrorInFunction
            << "Unknown package: " << package_
            << Foam::exit(Foam::FatalError);
    }
}

double RBFinterpolator::predict(const Eigen::VectorXd& x)
{
    if (package_ == "mathtoolbox")
    {
        return mtb_->predict(x);
    }
    else if (package_ == "splinter")
    {
        return splinter_->predict(x);
    }    
    else
    {
        FatalErrorInFunction
            << "Unknown package: " << package_
            << Foam::exit(Foam::FatalError);
    }
    return 0.0; // unreachable, keeps compiler happy
}

Eigen::VectorXd RBFinterpolator::predict(const Eigen::MatrixXd& X)
{
    if (package_ == "mathtoolbox")
    {
        return mtb_->predict(X);
    }
    else if (package_ == "splinter")
    {
        return splinter_->predict(X);
    }
    else
    {
        FatalErrorInFunction
            << "Unknown package: " << package_
            << Foam::exit(Foam::FatalError);
    }
    return Eigen::VectorXd(); // unreachable, keeps compiler happy
}

void RBFinterpolator::printInfo() const
{
    if (package_ == "mathtoolbox")
    {
        mtb_->printInfo();
    }
    else if (package_ == "splinter")
    {
        splinter_->printInfo();
    }
    else
    {
        FatalErrorInFunction
            << "Unknown package: " << package_
            << Foam::exit(Foam::FatalError);
    }
}
