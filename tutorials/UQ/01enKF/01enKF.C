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
    Test of the Ensemble Kalman filter implementation
SourceFiles
    01enKF.C
\*---------------------------------------------------------------------------*/


#include "MUQ/Modeling/Distributions/Gaussian.h"
#include "MUQ/Modeling/Distributions/Density.h"

#include <iostream>
#include "fvCFD.H"
#include "IOmanip.H"
#include "ITHACAutilities.H"
#include <Eigen/Dense>
#define _USE_MATH_DEFINES
#include <cmath>
#include "Foam2Eigen.H"
#include "muq2ithaca.H"

int main(int argc, char* argv[])
{
    word outputFolder = "./ITHACAoutput/";
    int Nseeds = 1000;
    Eigen::MatrixXd A = ITHACAstream::readMatrix("A_mat.txt");
    Eigen::MatrixXd Aw = ITHACAstream::readMatrix("Awrong_mat.txt");
    Eigen::MatrixXd H = ITHACAstream::readMatrix("observation_mat.txt");
    M_Assert(ITHACAstream::readMatrix("initialState_mat.txt").cols() == 1,
             "Wrong initialState input");
    Eigen::VectorXd x0 = ITHACAstream::readMatrix("initialState_mat.txt").col(0);
    int stateSize = A.rows();
    int obsSize = H.rows();
    std::cout <<
              "In this tutorial we have a dynamical system in the form:\ndx/dt = A * x" <<
              std::endl;
    std::cout << "with A = \n" << A << std::endl;
    std::cout << "We observe the state x by mean of the observation matrix \nH = \n"
              << H << std::endl;
    std::cout <<
              "The objective is to reconstruct the vector state knowing H and x0 = \n" <<
              x0.transpose() <<
              "\nbut having a wrong A" << std::endl;
    std::cout << "A_wrong =\n" << Aw << std::endl;
    int Ntimes = 201;
    int sampleDeltaStep = 10;
    double endTime = 10;
    double deltaTime = endTime / Ntimes;
    Eigen::VectorXd time = Eigen::VectorXd::LinSpaced(Ntimes, 0, endTime);
    Eigen::VectorXd xOld = x0;
    Eigen::MatrixXd X(stateSize, Ntimes);
    X.col(0) = x0;
    int sampleFlag = sampleDeltaStep;
    int Nsamples = (Ntimes - 1) / sampleDeltaStep;
    int sampleI = 0;
    Eigen::MatrixXd obs(obsSize, Nsamples);

    for (int timeI = 0; timeI < Ntimes - 1; timeI++)
    {
        Eigen::VectorXd xNew = (A * deltaTime + Eigen::MatrixXd::Identity(A.rows(),
                                A.cols()))  * xOld;
        xOld = xNew;
        Eigen::VectorXd dNew = H * xNew;
        X.col(timeI + 1) = xNew;
        sampleFlag--;

        if (sampleFlag == 0)
        {
            sampleFlag = sampleDeltaStep;
            obs.col(sampleI) = dNew;
            sampleI++;
        }
    }

    M_Assert(Nsamples == sampleI, "Something went wrong in the sampling");
    ITHACAstream::exportMatrix(time, "time", "eigen", outputFolder);
    ITHACAstream::exportMatrix(X, "X", "eigen", outputFolder);
    Eigen::VectorXd x = x0;
    Eigen::VectorXd prior_mu = x * 0.0;
    Eigen::MatrixXd prior_cov = Eigen::MatrixXd::Identity(stateSize,
                                stateSize) * 0.5;
    auto priorDensity = std::make_shared<muq::Modeling::Gaussian>(prior_mu,
                        prior_cov);
    Eigen::VectorXd modelError_mu = x * 0.0;
    Eigen::MatrixXd modelError_cov = Eigen::MatrixXd::Identity(stateSize,
                                     stateSize) * 0.7;
    auto modelErrorDensity = std::make_shared<muq::Modeling::Gaussian>
                             (modelError_mu, modelError_cov);
    Eigen::MatrixXd posteriorSamples(stateSize, Nseeds);
    Eigen::MatrixXd priorSamples(stateSize, Nseeds);

    for (int i = 0; i < Nseeds; i++)
    {
        priorSamples.col(i) = priorDensity->Sample();
    }

    posteriorSamples = priorSamples;
    Eigen::MatrixXd meas_cov = Eigen::MatrixXd::Identity(obsSize, obsSize) * 0.3;
    auto measNoise = std::make_shared<muq::Modeling::Gaussian>
                     (Eigen::VectorXd::Zero(obsSize), meas_cov);
    Eigen::MatrixXd posteriorMean(stateSize, Ntimes);
    Eigen::MatrixXd minConfidence = posteriorMean;
    Eigen::MatrixXd maxConfidence = minConfidence;
    posteriorMean.col(0) = posteriorSamples.rowwise().mean();
    sampleFlag = sampleDeltaStep;
    sampleI = 0;
    Eigen::MatrixXd forwardSamples(stateSize, Nseeds);

    for (int timeI = 0; timeI < Ntimes - 1; timeI++)
    {
        priorSamples = posteriorSamples;
        Eigen::MatrixXd forwardSamplesOld = forwardSamples;

        //Forecast step
        for (int i = 0; i < Nseeds; i++)
        {
            forwardSamples.col(i) = (A * deltaTime + Eigen::MatrixXd::Identity(A.rows(),
                                     A.cols())) * priorSamples.col(i) + modelErrorDensity->Sample();
        }

        sampleFlag--;

        if (sampleFlag == 0)
        {
            sampleFlag = sampleDeltaStep;
            Eigen::VectorXd meas = obs.col(sampleI);
            //Kalman filter
            posteriorSamples = ITHACAmuq::muq2ithaca::EnsembleKalmanFilter(forwardSamples,
                               meas, meas_cov, H * forwardSamples);
            sampleI++;
        }
        else
        {
            posteriorSamples = forwardSamples;
        }

        posteriorMean.col(timeI + 1) = posteriorSamples.rowwise().mean();
        minConfidence.col(timeI + 1) = ITHACAmuq::muq2ithaca::quantile(posteriorSamples,
                                       0.05);
        maxConfidence.col(timeI + 1) = ITHACAmuq::muq2ithaca::quantile(posteriorSamples,
                                       0.95);
    }

    ITHACAstream::exportMatrix(posteriorMean, "posteriorMean", "eigen",
                               outputFolder);
    ITHACAstream::exportMatrix(minConfidence, "minConfidence", "eigen",
                               outputFolder);
    ITHACAstream::exportMatrix(maxConfidence, "maxConfidence", "eigen",
                               outputFolder);
    return 0;
}
