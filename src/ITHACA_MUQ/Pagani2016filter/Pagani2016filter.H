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
    Pagani2016filter
Description
    A class that implements the Ensamble Kalman Filter in the version by Pagani et al. 2016
SourceFiles
    Pagani2016filter.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the Pagani2016filter class. It contains functions to
/// \dir
/// Directory containing the header and source files for the Pagani2016filter class.

#ifndef pagani2016filter_H
#define pagani2016filter_H

#include "ITHACAutilities.H"
#include "muq2ithaca.H"
#include "ensembleClass.H"




namespace ITHACAmuq
{
//--------------------------------------------------------------------------
/// @brief      Class for Pagani2016filters
///
class Pagani2016filter
{
    int Nseeds = 0;

    int timeSampI = 0;
    int stateSize = 0;
    int parameterSize = 0;
    int observationSize = 0;

    int timeStepI = 0;
    double startTime;
    double deltaTime;
    double endTime;
    int Ntimes = 0;

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
    Eigen::MatrixXd trueObservations;

    /// Initial state density
    std::shared_ptr<muq::Modeling::Gaussian> initialStateDensity;

    /// Parameter prior density
    std::shared_ptr<muq::Modeling::Gaussian> parameterPriorDensity;

   /// Parameter-observation cross covariance
   Eigen::MatrixXd parameterObservation_crossCov;

   /// State-observation cross covariance
   Eigen::MatrixXd stateObservation_crossCov;

   /// Observation covariance
   Eigen::MatrixXd observation_Cov;


    int innerLoopI = 0;


    public:
        // Constructors
        Pagani2016filter();

        Pagani2016filter(int _Nseeds);

        virtual ~Pagani2016filter() = default;

        /// Measurement noise density
        std::shared_ptr<muq::Modeling::Gaussian> measNoiseDensity;

        /// Parameter error density
        std::shared_ptr<muq::Modeling::Gaussian> parameterErrorDensity;

        /// State ensemble
        ensemble stateEns;

        /// Old state ensemble
        ensemble oldStateEns;

        /// Parameter ensemble
        ensemble parameterEns;

        /// Observation ensemble
        ensemble observationEns;

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
        /// Return the size of the observation vector
        int getObservationSize();

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
        /// Compute state mean
        void computeStateMean(int _index);

        //--------------------------------------------------------------------------
        /// Compute parameter mean
        void computeParameterMean(int _index);

        //--------------------------------------------------------------------------
        /// Return 
        int getObservationCounter();

        //--------------------------------------------------------------------------
        /// Return timestep at corresponding to the input observation step
        ///
        int getObservationTimestep(int _timeSampleI);

        //--------------------------------------------------------------------------
        /// Set the observations matrix
        void setTrueObservations(Eigen::MatrixXd _observations);

        //--------------------------------------------------------------------------
        /// Setup the time vector 
        void setTime(double _startTime, double _deltaTime, double _endTime);

        //--------------------------------------------------------------------------
        /// Setup the observation vector 
        void setObservationTime(int _observationStart, int _observationDelta);

        //--------------------------------------------------------------------------
        /// Setup of the measurement noise distribution
        void setMeasNoise(double cov);

        //--------------------------------------------------------------------------
        /// Setup of the parameter error distribution
        void setParameterError(double cov);

        //--------------------------------------------------------------------------
        /// Create initial state ensemble 
        void setInitialStateDensity(Eigen::VectorXd _mean, Eigen::MatrixXd _cov);

        //--------------------------------------------------------------------------
        /// Create initial state ensemble 
        void sampleInitialState();

        //--------------------------------------------------------------------------
        /// Create initial parameter ensemble 
        void sampleInitialParameter();

        //--------------------------------------------------------------------------
        /// Create parameter ensemble 
        void setParameterPriorDensity(Eigen::VectorXd _mean, Eigen::MatrixXd _cov);

        //--------------------------------------------------------------------------
        /// General class to sample from an input density 
        Eigen::MatrixXd ensembleFromDensity(std::shared_ptr<muq::Modeling::Gaussian> 
                _density);

        //--------------------------------------------------------------------------
        /// Project the state in between two sampling steps and fills the state and 
        /// parameter mean and confidence level, their ensembles and the observationEns
        /// It also must make timeStepI evolving 
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
        /// Perform Kalman filter update 
        void update();

        //--------------------------------------------------------------------------
        /// Run the filtering
        void run(word outputFolder);
};
}
#endif
