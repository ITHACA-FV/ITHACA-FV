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
    Fang2017filter_wDF
Description
    A class that implements the Ensamble Kalman Filter in the version by Fang et al. 2017
SourceFiles
    Fang2017filter_wDF.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the Fang2017filter_wDF class. It contains functions to
/// \dir
/// Directory containing the header and source files for the Fang2017filter_wDF class.

#ifndef fang2017filter_wDF_H
#define fang2017filter_wDF_H

#include "ITHACAutilities.H"
#include "muq2ithaca.H"
#include "ensembleClass.H"

namespace ITHACAmuq
{
//--------------------------------------------------------------------------
/// @brief      Class for Fang2017filter_wDFs
///
class Fang2017filter_wDF
{
    int Nsamples = 0;

    int stateSize = 0;
    int parameterSize = 0;
    int observationSize = 0;

    int timeStepI = 0;
    double startTime;
    double deltaTime;
    double endTime;
    int Ntimes;

    bool initialStateFlag = 0;
    bool parameterPriorFlag = 0;
    bool modelErrorFlag = 0;
    bool measurementNoiseFlag = 0;

    /// Mean of the initial state
    Eigen::VectorXd initialStateMean;

    /// Covariance of the initial state
    Eigen::MatrixXd initialStateCov;

    /// Mean of the prior parameter
    Eigen::VectorXd parameterPriorMean;

    /// Covariance of the prior parameter
    Eigen::MatrixXd parameterPriorCov;

    /// Parameters mean
    Eigen::MatrixXd parameterMean;

    /// State mean
    Eigen::MatrixXd stateMean;

    /// Timesteps vector
    Eigen::VectorXd timeVector;

    /// Timestep at which observations start
    int observationStart;

    /// Number of timesteps between observations
    int observationDelta;

    /// Vector that has 1 for the timestesp for which I have samples
    Eigen::VectorXi observationBoolVec;

    /// Matrix containing on each column the observations
    Eigen::MatrixXd observations;

    /// Initial state density
    std::shared_ptr<muq::Modeling::Gaussian> initialStateDensity;

    /// Parameter prior density
    std::shared_ptr<muq::Modeling::Gaussian> parameterPriorDensity;


    int innerLoopI = 0;

    /// Set to 1 to use a univariate sampling for the initial state density
    /// (use with big states)
    bool univariateInitStateDensFlag = 0;


    public:
        // Constructors
        Fang2017filter_wDF();

        Fang2017filter_wDF(int _Nsamples);

        virtual ~Fang2017filter_wDF() = default;
            
        /// Model noise density
        std::shared_ptr<muq::Modeling::Gaussian> modelErrorDensity;

        /// Measurement noise density
        std::shared_ptr<muq::Modeling::Gaussian> measNoiseDensity;

        /// State ensemble
        ensemble stateEns;

        /// Parameter ensemble
        ensemble parameterEns;

        /// Observation ensemble
        ensemble observationEns;

        /// Joint state and parameter ensemble
        ensemble jointEns;

        /// State maximum confidence level
        Eigen::MatrixXd state_maxConf;

        /// State minimum confidence level
        Eigen::MatrixXd state_minConf;

        /// parameter maximum confidence level
        Eigen::MatrixXd parameter_maxConf;

        /// parameter minimum confidence level
        Eigen::MatrixXd parameter_minConf;

        // Functions

        //--------------------------------------------------------------------------
        /// Return number of samples per ensamble
        int getNumberOfSamples();

        //--------------------------------------------------------------------------
        /// Return time
        double getTime();

        //--------------------------------------------------------------------------
        /// Return time for input timestep
        double getTime(int _timeStepI);

        //--------------------------------------------------------------------------
        /// Return timestep
        int getTimeStep();

        //--------------------------------------------------------------------------
        /// Return time vector
        Eigen::VectorXd getTimeVector();

        //--------------------------------------------------------------------------
        /// Return state mean
        Eigen::MatrixXd getStateMean();

        //--------------------------------------------------------------------------
        /// Return parameter mean
        Eigen::MatrixXd getParameterMean();

        //--------------------------------------------------------------------------
        /// Return parameter Max Confidence
        Eigen::MatrixXd getParameterMaxConf();   // Kabir

        //--------------------------------------------------------------------------
        /// Return parameter Min Confidence
        Eigen::MatrixXd getParameterMinConf();   // Kabir

        //--------------------------------------------------------------------------
        /// Set the observations matrix
        void setObservations(Eigen::MatrixXd _observations);

        //--------------------------------------------------------------------------
        /// Setup the time vector 
        void setTime(double _startTime, double _deltaTime, double _endTime);

        //--------------------------------------------------------------------------
        /// Setup the observation vector 
        void setObservationTime(int _observationStart, int _observationDelta);

        //--------------------------------------------------------------------------
        /// Setup of the model error distribution
        void setModelError(double cov, bool univariate = 0);

        //--------------------------------------------------------------------------
        /// Setup of the measurement noise distribution
        void setMeasNoise(double cov);

        //--------------------------------------------------------------------------
        /// Create initial state ensemble
        void setInitialStateDensity(Eigen::VectorXd _mean, Eigen::MatrixXd _cov,
                bool _univariateFlag = 0);

        //--------------------------------------------------------------------------
        /// Create initial state ensemble 
        void sampleInitialState();

        //--------------------------------------------------------------------------
        /// Create parameter ensemble 
        void setParameterPriorDensity(Eigen::VectorXd _mean, Eigen::MatrixXd _cov);

        //--------------------------------------------------------------------------
        /// Create parameter ensemble 
        void sampleParameterDist();

        //--------------------------------------------------------------------------
        /// General class to sample from an input density 
        Eigen::MatrixXd ensembleFromDensity(
                std::shared_ptr<muq::Modeling::Gaussian> _density);

        //--------------------------------------------------------------------------
        /// Concatenate state and parameter ensambles to create the joint ensamble 
        void buildJointEns();

        //--------------------------------------------------------------------------
        /// 
        virtual void stateProjection() = 0;

        //--------------------------------------------------------------------------
        /// 
        virtual void observeState() = 0;

        //--------------------------------------------------------------------------
        /// 
        void setObservationSize(int _size);

        //--------------------------------------------------------------------------
        /// 
        void setStateSize(int _size);

        //--------------------------------------------------------------------------
        /// 
        int getStateSize();

        //--------------------------------------------------------------------------
        /// 
        void setParameterSize(int _size);

        //--------------------------------------------------------------------------
        /// 
        int getParameterSize();

        //--------------------------------------------------------------------------
        /// 
        void updateJointEns(Eigen::VectorXd _observation);

        //--------------------------------------------------------------------------
        /// Run the filtering
        void run(int innerLoopMax, word outputFolder);
};
}
#endif
