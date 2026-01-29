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

#include "mtbRBF.H"

namespace {
    struct QuinticRbfKernel {
        double operator()(double r) const {
            return std::pow(r, 5);
        }
    };

    struct MultiquadricRbfKernel {
        double epsilon;
        MultiquadricRbfKernel(double eps = 1.0) : epsilon(eps) {}
        double operator()(double r) const {
            return std::sqrt(1.0 + (epsilon * r) * (epsilon * r));
        }
    };

    struct InverseMultiquadricRbfKernel {
        double epsilon;
        InverseMultiquadricRbfKernel(double eps = 1.0) : epsilon(eps) {}
        double operator()(double r) const {
            return 1.0 / std::sqrt(1.0 + (epsilon * r) * (epsilon * r));
        }
    };

    struct InverseQuadraticRbfKernel {
        double epsilon;
        InverseQuadraticRbfKernel(double eps = 1.0) : epsilon(eps) {}
        double operator()(double r) const {
            return 1.0 / (1.0 + (epsilon * r) * (epsilon * r));
        }
    };

    std::function<double(const double)> getKernelFunction(const std::string& kernelType, double epsilon)
    {
        if (kernelType == "gaussian")
        {
            return mathtoolbox::GaussianRbfKernel(epsilon);
        }
        else if (kernelType == "linear")
        {
            return mathtoolbox::LinearRbfKernel();
        }
        else if (kernelType == "cubic")
        {
            return mathtoolbox::CubicRbfKernel();
        }
        else if (kernelType == "quintic")
        {
            return QuinticRbfKernel();
        }
        else if (kernelType == "multiquadric")
        {
            return MultiquadricRbfKernel(epsilon);
        }
        else if (kernelType == "inverse_multiquadric")
        {
            return InverseMultiquadricRbfKernel(epsilon);
        }
        else if (kernelType == "inverse_quadratic")
        {
            return InverseQuadraticRbfKernel(epsilon);
        }
        else if (kernelType == "thin_plate_spline")
        {
            return mathtoolbox::ThinPlateSplineRbfKernel();
        }
        else
        {
            FatalErrorInFunction
                << "Unknown kernel type: " << kernelType
                << ". Valid options are: gaussian, linear, cubic, quintic, multiquadric, "
                << "inverse_multiquadric, inverse_quadratic, thin_plate_spline"
                << Foam::exit(Foam::FatalError);
            return mathtoolbox::GaussianRbfKernel(1.0); // Unreachable, but suppresses warning
        }
    }
}

mtbRBF::mtbRBF(const Foam::word& kernelType, bool usePolynomialTerm, bool useRegularization, Foam::scalar lambda, bool normalize)
:   kernelType_(kernelType),
    useRegularization_(useRegularization),
    usePolynomialTerm_(usePolynomialTerm),
    lambda_(lambda),
    normalize_(normalize)
{
    auto kernel = getKernelFunction(kernelType_, 1.0);
    impl_ = std::make_unique<mathtoolbox::RbfInterpolator>(kernel, usePolynomialTerm_);
}

mtbRBF::mtbRBF(const Foam::dictionary& dict)
{
    kernelType_ = dict.lookupOrDefault<Foam::word>("kernel", "gaussian");    
    usePolynomialTerm_ = dict.lookupOrDefault<bool>("polynomial", true);
    useRegularization_ = dict.lookupOrDefault<bool>("regularization", true);
    lambda_ = dict.lookupOrDefault<Foam::scalar>("lambda", 0.001);
    normalize_ = dict.lookupOrDefault<bool>("normalize", true);
    epsilon_ = dict.lookupOrDefault<Foam::scalar>("epsilon", 1.0);

    auto kernel = getKernelFunction(kernelType_, epsilon_);
    impl_ = std::make_unique<mathtoolbox::RbfInterpolator>(kernel, usePolynomialTerm_);
}

mtbRBF::~mtbRBF() = default;

void mtbRBF::fit(const Eigen::MatrixXd& X, const Eigen::VectorXd& y)
{
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

        impl_->SetData(X_norm, y_norm);
    }
    else
    {
        impl_->SetData(X, y);
    }
    impl_->CalcWeights(useRegularization_, lambda_);
}

Foam::scalar mtbRBF::predict(const Eigen::VectorXd& x)
{
    if (normalize_)
    {
        Eigen::VectorXd x_norm = (x - xMean_).array() / xStd_.array();
        Foam::scalar y_pred_norm = impl_->CalcValue(x_norm);
        return y_pred_norm * yStd_ + yMean_;
    }
    return impl_->CalcValue(x);
}

Eigen::VectorXd mtbRBF::predict(const Eigen::MatrixXd& X)
{
    if (normalize_)
    {
        Eigen::MatrixXd X_norm = (X.colwise() - xMean_).array().colwise() / xStd_.array();
        
        Eigen::VectorXd result(X.cols());
        for (int i = 0; i < X.cols(); ++i)
        {
            result(i) = impl_->CalcValue(X_norm.col(i));
        }
        return result.array() * yStd_ + yMean_;
    }

    Eigen::VectorXd result(X.cols());
    for (int i = 0; i < X.cols(); ++i)
    {
        result(i) = impl_->CalcValue(X.col(i));
    }
    return result;
}

void mtbRBF::printInfo()
{
    Foam::Info << "mtbRBF Model Info:" << Foam::endl;
    Foam::Info << "\t useKernel: " << kernelType_ << Foam::endl;
    Foam::Info << "\t usePolynomialTerm: " << (usePolynomialTerm_ ? "true" : "false") << Foam::endl;
    Foam::Info << "\t useRegularization: " << (useRegularization_ ? "true" : "false") << Foam::endl;
    Foam::Info << "\t lambda: " << lambda_ << Foam::endl;
    Foam::Info << "\t normalize: " << (normalize_ ? "true" : "false") << Foam::endl;
}