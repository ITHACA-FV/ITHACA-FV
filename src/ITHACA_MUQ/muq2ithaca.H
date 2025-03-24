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
    muq2ithaca
Description
    Set of function to convert ITHACA objects into MUQ2 objects and vice versa
SourceFiles
    muq2ithaca.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the muq2ithaca namespace. It contains functions to
/// \dir
/// Directory containing the header and source files for the muq2ithaca namespace.

#ifndef muq2ithaca_H
#define muq2ithaca_H

#include <iostream>
#include <Eigen/Eigen>
#include "ITHACAassert.H"
#include "Foam2Eigen.H"

#include "MUQ/Modeling/Distributions/Gaussian.h"
#include "MUQ/Modeling/Distributions/Density.h"

#include "MUQ/SamplingAlgorithms/SamplingProblem.h"
#include "MUQ/SamplingAlgorithms/SingleChainMCMC.h"
#include "MUQ/SamplingAlgorithms/MCMCFactory.h"

#include <boost/property_tree/ptree.hpp>


namespace ITHACAmuq
{
namespace muq2ithaca
{
//--------------------------------------------------------------------------
/// @brief      Ensemble Kalman Filter
///
/// @param[in]  prior           Samples of the prior
/// @param[in]  measurements    Measured data
/// @param[in]  measurementsCov Covariance matrix for the measurements, gaussian noise with zero mean is assumed
/// @param[in]  observedState   Ensemble of the observed state
///
/// @return     Ensamble of the posterior
///
Eigen::MatrixXd EnsembleKalmanFilter(Eigen::MatrixXd prior,
                                     Eigen::VectorXd measurements,
                                     Eigen::MatrixXd measurementsCov,
                                     Eigen::MatrixXd observedState);

//--------------------------------------------------------------------------
/// @brief      Ensemble Kalman Filter
///
/// @param[in]  prior           Samples of the prior
/// @param[in]  measurements    Measured data
/// @param[in]  measurementsCov Covariance matrix for the measurements, gaussian noise with zero mean is assumed
/// @param[in]  observedState   Ensemble of the observed state
///
/// @return     Ensamble of the posterior
///
Eigen::MatrixXd EnsembleKalmanFilter(PtrList<volScalarField>& prior,
                                     Eigen::VectorXd measurements,
                                     Eigen::MatrixXd measurementsCov,
                                     Eigen::MatrixXd observedState);

//--------------------------------------------------------------------------
/// @brief      Returns quantile for a vector of samples
///
/// @param[in]  samps   Samples
/// @param[in]  p       Quantile parameter
/// @param[in]  method  Method used to compute quantile according to Hyndman and Fan (1986)
///
/// @return     Quantile
///
double quantile(Eigen::VectorXd samps, double p, int method = 1);

//--------------------------------------------------------------------------
/// @brief      Returns quantile for a matrix of samples
///
/// @param[in]  samps   Samples, each variable on a row
/// @param[in]  p       Quantile parameter
/// @param[in]  method  Method used to compute quantile according to Hyndman and Fan (1986)
///
/// @return     Quantile
///
Eigen::VectorXd quantile(Eigen::MatrixXd samps, double p, int method = 1);
}
}

#endif
