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

#include "ithacaInterpolator.H"
#include "mtbGPR.H"
#include "mtbRBF.H"
#include "splinterRBF.H"
#include "error.H"

ithacaInterpolator::ithacaInterpolator(const Foam::dictionary& dict)
{
    package_ = dict.lookupOrDefault<Foam::word>("package", "mathtoolbox");
    algorithm_ = dict.lookupOrDefault<Foam::word>("algorithm", "RBF");

    if (package_ == "mathtoolbox")
    {
        if (algorithm_ == "RBF" || algorithm_ == "rbf")
        {
            mtb_ = std::make_unique<mtbRBF>(dict);
        }
        else if (algorithm_ == "GPR" || algorithm_ == "gpr")
        {
            mtbGPR_ = std::make_unique<mtbGPR>(dict);
        }
        else
        {
            FatalErrorInFunction
                << "Unknown algorithm for mathtoolbox package: " << algorithm_
                << ". Valid options are: RBF, GPR"
                << Foam::exit(Foam::FatalError);
        }
    }
    else if (package_ == "splinter")
    {
        if (algorithm_ != "RBF" && algorithm_ != "rbf")
        {
            FatalErrorInFunction
                << "Unknown algorithm for splinter package: " << algorithm_
                << ". Valid option is: RBF"
                << Foam::exit(Foam::FatalError);
        }
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

ithacaInterpolator::~ithacaInterpolator() = default;

void ithacaInterpolator::fit(const Eigen::MatrixXd& X, const Eigen::VectorXd& y)
{
    if (package_ == "mathtoolbox")
    {
        if (mtb_)
        {
            mtb_->fit(X, y);
        }
        else if (mtbGPR_)
        {
            mtbGPR_->fit(X, y);
        }
        else
        {
            FatalErrorInFunction
                << "mathtoolbox algorithm not initialized"
                << Foam::exit(Foam::FatalError);
        }
    }
    else if (package_ == "splinter")
    {
        if (algorithm_ != "RBF" && algorithm_ != "rbf")
        {
            FatalErrorInFunction
                << "Unknown algorithm for splinter package: " << algorithm_
                << ". Valid option is: RBF"
                << Foam::exit(Foam::FatalError);
        }
        splinter_->fit(X, y);
    }
    else
    {
        FatalErrorInFunction
            << "Unknown package: " << package_
            << Foam::exit(Foam::FatalError);
    }
}

double ithacaInterpolator::predict(const Eigen::VectorXd& x)
{
    if (package_ == "mathtoolbox")
    {
        if (mtb_)
        {
            return mtb_->predict(x);
        }
        else if (mtbGPR_)
        {
            return mtbGPR_->predict(x);
        }
        else
        {
            FatalErrorInFunction
                << "mathtoolbox algorithm not initialized"
                << Foam::exit(Foam::FatalError);
        }
    }
    else if (package_ == "splinter")
    {
        if (algorithm_ != "RBF" && algorithm_ != "rbf")
        {
            FatalErrorInFunction
                << "Unknown algorithm for splinter package: " << algorithm_
                << ". Valid option is: RBF"
                << Foam::exit(Foam::FatalError);
        }
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

Eigen::VectorXd ithacaInterpolator::predict(const Eigen::MatrixXd& X)
{
    if (package_ == "mathtoolbox")
    {
        if (mtb_)
        {
            return mtb_->predict(X);
        }
        else if (mtbGPR_)
        {
            return mtbGPR_->predict(X);
        }
        else
        {
            FatalErrorInFunction
                << "mathtoolbox algorithm not initialized"
                << Foam::exit(Foam::FatalError);
        }
    }
    else if (package_ == "splinter")
    {
        if (algorithm_ != "RBF" && algorithm_ != "rbf")
        {
            FatalErrorInFunction
                << "Unknown algorithm for splinter package: " << algorithm_
                << ". Valid option is: RBF"
                << Foam::exit(Foam::FatalError);
        }
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

void ithacaInterpolator::printInfo() const
{
    if (package_ == "mathtoolbox")
    {
        if (mtb_)
        {
            mtb_->printInfo();
        }
        else if (mtbGPR_)
        {
            mtbGPR_->printInfo();
        }
        else
        {
            FatalErrorInFunction
                << "mathtoolbox algorithm not initialized"
                << Foam::exit(Foam::FatalError);
        }
    }
    else if (package_ == "splinter")
    {
        if (algorithm_ != "RBF" && algorithm_ != "rbf")
        {
            FatalErrorInFunction
                << "Unknown algorithm for splinter package: " << algorithm_
                << ". Valid option is: RBF"
                << Foam::exit(Foam::FatalError);
        }
        splinter_->printInfo();
    }
    else
    {
        FatalErrorInFunction
            << "Unknown package: " << package_
            << Foam::exit(Foam::FatalError);
    }
}
