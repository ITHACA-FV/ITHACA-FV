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
Description
    Example of a state reconstruction in 1D heat transfer problem using EnKF
SourceFiles
    02enKF_1DinverseHeatTransfer.C
\*---------------------------------------------------------------------------*/

#include <iostream>
#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include "IOmanip.H"
#include "Time.H"
#include "laplacianProblem.H"
#include "ITHACAutilities.H"
#include <Eigen/Dense>
#define _USE_MATH_DEFINES
#include <cmath>
#include "Foam2Eigen.H"

#include "MUQ/Modeling/Distributions/Gaussian.h"

#include "muq2ithaca.H"

#include "02enKF_1DinverseHeatTransfer.H"

using namespace SPLINTER;


int main(int argc, char* argv[])
{
    solverPerformance::debug = 1; //No verbose output
    // Number of elements in the ensemble
    int Nseeds = 50;
    // Create the example object of the 02enKF_1DinverseHeatTransfer type
    EnKF_1DinverseHeatTransfer example(argc, argv, Nseeds);
    // Reading parameters from file
    ITHACAparameters* para = ITHACAparameters::getInstance(example._mesh(),
                             example._runTime());
    example.k = para->ITHACAdict->lookupOrDefault<double>("thermalConductivity", 0);
    M_Assert(example.k > 0, "thermalConductivity, k, not specified");
    example.rho = para->ITHACAdict->lookupOrDefault<double>("density", 0);
    M_Assert(example.rho > 0, "Density, rho, not specified");
    example.Cp = para->ITHACAdict->lookupOrDefault<double>("heatCapacity", 0);
    M_Assert(example.Cp > 0, "heatCapacity, Cp, not specified");
    // Mesh and temperature field setup
    fvMesh& mesh = example._mesh();
    volScalarField& T = example._T();
    int stateSize = T.size();
    scalar initialField = T.internalField()[0];
    // Performes true solution
    example.solveDirect();
    int Ntimes = example.Ntimes;
    // Setting up the densities
    // Gaussian densities are assumed
    example.priorSetup(initialField, 0.7);
    example.modelErrorSetup(0.0, 0.7);
    example.measNoiseSetup(0.0, 0.05);
    Eigen::MatrixXd posteriorSamples(stateSize, Nseeds);
    Eigen::MatrixXd priorSamples(stateSize, Nseeds);
    // Sampling prior density
    example.priorSampling();
    Eigen::MatrixXd posteriorMean(stateSize, Ntimes);
    Eigen::MatrixXd minConfidence = posteriorMean;
    Eigen::MatrixXd maxConfidence = minConfidence;
    posteriorMean.col(0) = posteriorSamples.rowwise().mean();
    // Performs reconstruction
    example.reconstruct();
    return 0;
}
//--------
/// \dir  02enKF_1DinverseHeatTransfer Folder of the UQ tutorial 2
/// \file
/// \brief Implementation of state reconstruction in a 1D heat transfer problem

/// \example 02enKF_1DinverseHeatTransfer.C
/// \section intro_invProb 1D heat transfer state reconstruction
/// The tutorial consists in the state reconstruction of a 1D heat transfer problem where we have some measurement points inside the domain but partially known boundary conditions.
///
/// The equations of the true problem are
///   \f{eqnarray*}{
///        \rho C_p\frac{\partial T}{\partial t} + k \Delta T  &=& 0 &\text{ in } (0,1) \times (0, t_f]\\
///        \nabla T \cdot \mathbf{n}  &=& g_{true}(t) &\text{ on } x=0 \times (0, t_f]\\
///        T  &=& 0.5 &\text{ on } x=1 \times (0, t_f]\\
///        T  &=& T_0 &\text{ on } (0,1) \times t=0
///   \f}
/// where \f$\rho\f$, \f$C_p\f$, \f$k\f$ and \f$T_0\f$ are defined by the user in the ITHACAdict and 0/T file while \f$g_{true}(t) = 5 t + 2\f$ and \f$T_0 = 1\f$.
///
/// The objective of the tutorial is to reconstruct the true solution having some noisy temperature measurements defined in the measurementsDict and the model
///   \f{eqnarray*}{
///        \rho C_p\frac{\partial T}{\partial t} + k \Delta T  &=& 0 &\text{ in } (0,1) \times (0, t_f]\\
///        \nabla T \cdot \mathbf{n}  &=& g(t) &\text{ on } x=0 \times (0, t_f]\\
///        T  &=& 0.5 &\text{ on } x=1 \times (0, t_f]\\
///        T  &=& T_0 &\text{ on } (0,1) \times t=0
///   \f}
/// where the boundary condition \f$g(t)=5t\f$ is different from the right one \f$g_{true}\f$.
///
/// To solve this problem, we use the Ensemble Kalman Filter implemented in the ITHACAmuq::muq2ithaca::EnsembleKalmanFilter method.
/// We start by assuming a Gaussian prior for the state (the temperature in this case)
/// \f[
/// \omega_{prior} \sim \mathcal{N}(T_0, \sigma_{prior}I).
/// \f]
/// Then, we sample it creating the prior ensemble
/// \f[
/// E_{prior}^0.
/// \f]
///
/// Now, at each timestep, we perform a forecasting step.
/// In the forecasting step, we solve the direct problem with the "wrong" heat flux for each of the member of the ensemble adding to the solution the model error
/// \f[
/// \omega_{model} \sim \mathcal{N}(0, \sigma_{model}I).
/// \f]
/// This way, we produced at each timestep a posterior ensamble
/// \f[
/// E_{post}^n.
/// \f]
///
/// If there are no measurements available at this timestep
/// \f[
/// E_{post}^n = E_{prior}^{n-1}.
/// \f]
/// However, when we have available measurements, we perform a Kalman filter correction step.
/// The Kalman filter updates the posterior ensamble such that
/// \f[
/// \hat{E}_{post}^n = E_{post}^{n} + Z M,
/// \f]
/// where
/// \f[
/// Z = \frac{1}{Nseeds - 1} A  (HA)^T,
/// \f]
/// \f[
/// A = E_{prior}(i) - \mu_{prior},
/// \f]
/// \mu_{prior} being a matrix having the mean value on the ensemble repeated on each column.
/// \f$HA\f$ is similar to \f$A\f$ but with the value of the state in the observation point,
/// \f[
/// M = P Y,
/// \f]
/// \f[
/// P = \frac{1}{Nseeds - 1} (HA)  (HA)^T  + \sigma_{meas} I,
/// \f]
/// and \f$Y\f$ is the matrix containing the difference between the measured and computed state at the measurement points.
///
/// The image below shows the true end reconstructed temperature at \f$x=0.5\f$ together with the 95\% quantile
/// \image html 02enKF_1DinverseHeatTransfer.png
///
///
/// \section plaincode The plain program
/// Here there's the plain code
///
