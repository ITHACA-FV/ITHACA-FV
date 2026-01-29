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

#include "splinterRBF.H"

splinterRBF::splinterRBF(const Foam::word& kernelType, bool normalize, Foam::scalar epsilon)
:   kernelType_(kernelType),
    normalize_(normalize),
    epsilon_(epsilon)
{
}

splinterRBF::splinterRBF(const Foam::dictionary& dict)
{
    kernelType_ = dict.lookupOrDefault<Foam::word>("kernel", "gaussian");
    normalize_ = dict.lookupOrDefault<bool>("normalize", true);
    epsilon_ = dict.lookupOrDefault<Foam::scalar>("epsilon", 1.0);
}

splinterRBF::~splinterRBF() = default;

SPLINTER::RadialBasisFunctionType splinterRBF::getKernelType(const Foam::word& kernelType)
{
    if (kernelType == "gaussian")
    {
        return SPLINTER::RadialBasisFunctionType::GAUSSIAN;
    }
    else if (kernelType == "linear")
    {
        return SPLINTER::RadialBasisFunctionType::LINEAR;
    }
    else if (kernelType == "cubic")
    {
        return SPLINTER::RadialBasisFunctionType::CUBIC;
    }
    else if (kernelType == "quintic")
    {
        return SPLINTER::RadialBasisFunctionType::QUINTIC;
    }
    else if (kernelType == "multiquadric")
    {
        return SPLINTER::RadialBasisFunctionType::MULTIQUADRIC;
    }
    else if (kernelType == "inverse_multiquadric")
    {
        return SPLINTER::RadialBasisFunctionType::INVERSE_MULTIQUADRIC;
    }
    else if (kernelType == "inverse_quadratic")
    {
        return SPLINTER::RadialBasisFunctionType::INVERSE_QUADRIC;
    }
    else if (kernelType == "thin_plate_spline")
    {
        return SPLINTER::RadialBasisFunctionType::THIN_PLATE_SPLINE;
    }
    else
    {
        FatalErrorInFunction
            << "Unknown kernel type: " << kernelType
            << ". Valid options are: gaussian, linear, cubic, quintic, multiquadric, "
            << "inverse_multiquadric, inverse_quadratic, thin_plate_spline"
            << Foam::exit(Foam::FatalError);
        return SPLINTER::RadialBasisFunctionType::GAUSSIAN; // Unreachable, but suppresses warning
    }
}

void splinterRBF::fit(const Eigen::MatrixXd& X, const Eigen::VectorXd& y)
{
    SPLINTER::DataTable data;

    if (normalize_)
    {
        // Input normalization
        xMean_ = X.rowwise().mean();
        if (X.cols() > 1)
        {
            Eigen::MatrixXd centered = X.colwise() - xMean_;
            xStd_ = (centered.array().square().rowwise().sum() / (X.cols() - 1)).sqrt();
        }
        else
        {
            xStd_ = Eigen::VectorXd::Ones(X.rows());
        }

        for(int i=0; i<xStd_.size(); ++i) {
            if(xStd_(i) < 1e-16) xStd_(i) = 1.0; 
        }

        Eigen::MatrixXd X_norm = (X.colwise() - xMean_).array().colwise() / xStd_.array();

        // Output normalization
        yMean_ = y.mean();
        if (y.size() > 1)
        {
            Eigen::VectorXd yCentered = y.array() - yMean_;
            yStd_ = std::sqrt(yCentered.array().square().sum() / (y.size() - 1));
        }
        else
        {
            yStd_ = 1.0;
        }
        
        if(yStd_ < 1e-16) yStd_ = 1.0;

        Eigen::VectorXd y_norm = (y.array() - yMean_) / yStd_;

        for (int i = 0; i < X_norm.cols(); ++i)
        {
            data.addSample(X_norm.col(i), y_norm(i));
        }
    }
    else
    {
        for (int i = 0; i < X.cols(); ++i)
        {
            data.addSample(X.col(i), y(i));
        }
    }
    
    impl_ = std::make_unique<SPLINTER::RBFSpline>(data, getKernelType(kernelType_), epsilon_);
}

Foam::scalar splinterRBF::predict(const Eigen::VectorXd& x)
{
    if (normalize_)
    {
        Eigen::VectorXd x_norm = (x - xMean_).array() / xStd_.array();
        Foam::scalar y_pred_norm = impl_->eval(x_norm);
        return y_pred_norm * yStd_ + yMean_;
    }
    return impl_->eval(x);
}

Eigen::VectorXd splinterRBF::predict(const Eigen::MatrixXd& X)
{
    Eigen::VectorXd result(X.cols());

    if (normalize_)
    {
        Eigen::MatrixXd X_norm = (X.colwise() - xMean_).array().colwise() / xStd_.array();
        
        for (int i = 0; i < X.cols(); ++i)
        {
            result(i) = impl_->eval(X_norm.col(i));
        }
        return result.array() * yStd_ + yMean_;
    }

    for (int i = 0; i < X.cols(); ++i)
    {
        result(i) = impl_->eval(X.col(i));
    }
    return result;
}

void splinterRBF::printInfo()
{
    Foam::Info << "splinterRBF Model Info:" << Foam::endl;
    Foam::Info << "\t Kernel Type: " << kernelType_ << Foam::endl;
    Foam::Info << "\t Normalize: " << (normalize_ ? "true" : "false") << Foam::endl;
    Foam::Info << "\t Epsilon: " << epsilon_ << Foam::endl;
}
