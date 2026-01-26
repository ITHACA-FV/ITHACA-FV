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

#include "mtbGPR.H"
#include "error.H"

mtbGPR::mtbGPR(const Foam::dictionary& dict)
{
    kernelTypeWord_ = dict.lookupOrDefault<Foam::word>("kernel", "matern");
    useDataNormalization_ = dict.lookupOrDefault<bool>("normalize", true);
    optimizeHyperparams_ = dict.lookupOrDefault<bool>("optimizeHyperparams", true);
    kernelScale_ = dict.lookupOrDefault<Foam::scalar>("kernelScale", 0.10);
    lengthScale_ = dict.lookupOrDefault<Foam::scalar>("lengthScale", 0.10);
    noise_ = dict.lookupOrDefault<Foam::scalar>("noise", 1e-5);
}

mtbGPR::~mtbGPR() = default;

mathtoolbox::GaussianProcessRegressor::KernelType mtbGPR::parseKernelType(const Foam::word& kernelWord) const
{
    const Foam::word lower = kernelWord;
    if (lower == "matern")
    {
        return mathtoolbox::GaussianProcessRegressor::KernelType::ArdMatern52;
    }
    else if (lower == "squared_exp")
    {
        return mathtoolbox::GaussianProcessRegressor::KernelType::ArdSquaredExp;
    }

    FatalErrorInFunction
        << "Unknown GPR kernel: " << kernelWord
        << ". Valid options are: matern, squared_exp"
        << Foam::exit(Foam::FatalError);
    return mathtoolbox::GaussianProcessRegressor::KernelType::ArdMatern52; // unreachable
}

void mtbGPR::fit(const Eigen::MatrixXd& X, const Eigen::VectorXd& y)
{
    if (X.rows() == 0)
    {
        FatalErrorInFunction << "Input matrix has zero rows" << Foam::exit(Foam::FatalError);
    }
    if (X.cols() != y.size())
    {
        FatalErrorInFunction
            << "Input size mismatch: cols(X) = " << X.cols() << ", size(y) = " << y.size()
            << Foam::exit(Foam::FatalError);
    }

    const auto kernelType = parseKernelType(kernelTypeWord_);

    Eigen::VectorXd kernelHyperparams = Eigen::VectorXd::Constant(X.rows() + 1, lengthScale_);
    kernelHyperparams[0] = kernelScale_;

    impl_ = std::make_unique<mathtoolbox::GaussianProcessRegressor>(X, y, kernelType, useDataNormalization_);

    if (optimizeHyperparams_)
    {
        impl_->PerformMaximumLikelihood(kernelHyperparams, noise_);
    }
    else
    {
        impl_->SetHyperparams(kernelHyperparams, noise_);
    }
}

Foam::scalar mtbGPR::predict(const Eigen::VectorXd& x)
{
    if (!impl_)
    {
        FatalErrorInFunction << "mtbGPR used before calling fit()" << Foam::exit(Foam::FatalError);
    }
    return impl_->PredictMean(x);
}

Eigen::VectorXd mtbGPR::predict(const Eigen::MatrixXd& X)
{
    if (!impl_)
    {
        FatalErrorInFunction << "mtbGPR used before calling fit()" << Foam::exit(Foam::FatalError);
    }

    Eigen::VectorXd result(X.cols());
    for (int i = 0; i < X.cols(); ++i)
    {
        result(i) = impl_->PredictMean(X.col(i));
    }
    return result;
}

void mtbGPR::printInfo() const
{
    Foam::Info << "mtbGPR Model Info:" << Foam::endl;
    Foam::Info << "\t kernel: " << kernelTypeWord_ << Foam::endl;
    Foam::Info << "\t normalize: " << (useDataNormalization_ ? "true" : "false") << Foam::endl;
    Foam::Info << "\t optimizeHyperparams: " << (optimizeHyperparams_ ? "true" : "false") << Foam::endl;
    Foam::Info << "\t kernelScale: " << kernelScale_ << Foam::endl;
    Foam::Info << "\t lengthScale: " << lengthScale_ << Foam::endl;
    Foam::Info << "\t noise: " << noise_ << Foam::endl;
}
