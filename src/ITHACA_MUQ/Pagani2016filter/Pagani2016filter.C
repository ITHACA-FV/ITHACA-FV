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


/// \file
/// Source file of the Pagani2016filter class.

#include "Pagani2016filter.H"

namespace ITHACAmuq
{
// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //
// Constructor
Pagani2016filter::Pagani2016filter() {}
Pagani2016filter::Pagani2016filter(int _Nseeds)
{
    Nseeds = _Nseeds;
}

// * * * * * * * * * * * * * * Filter Methods * * * * * * * * * * * * * * //

// Return number of samples per ensamble
int Pagani2016filter::getNumberOfSamples()
{
    return Nseeds;
}

int Pagani2016filter::getObservationSize()
{
    return trueObservations.rows();
}

// Return time
double Pagani2016filter::getTime()
{
    return timeVector(timeStepI);
}

//--------------------------------------------------------------------------
/// Return time for input timestep
double Pagani2016filter::getTime(int _timeStepI)
{
    return timeVector(_timeStepI);
}

//--------------------------------------------------------------------------
/// Return timestep
int Pagani2016filter::getTimeStep()
{
    return timeStepI;
}

//--------------------------------------------------------------------------
/// Return time vector
Eigen::VectorXd Pagani2016filter::getTimeVector()
{
    return timeVector;
}

//--------------------------------------------------------------------------
/// Return state mean
Eigen::MatrixXd Pagani2016filter::getStateMean()
{
    return stateMean;
}

//--------------------------------------------------------------------------
/// Return parameter mean
Eigen::MatrixXd Pagani2016filter::getParameterMean()
{
    return parameterMean;
}

//--------------------------------------------------------------------------
/// Compute state mean
void Pagani2016filter::computeStateMean(int _index)
{
    Info << "debug : stateMean.cols() = " << stateMean.cols() << endl;
    Info << "debug : _index = " << _index << endl;
    stateMean.col(_index) = stateEns.mean();
}

//--------------------------------------------------------------------------
/// Compute parameter mean
void Pagani2016filter::computeParameterMean(int _index)
{
    parameterMean.col(_index) = parameterEns.mean();
}

//--------------------------------------------------------------------------
/// Return observation caounter timeSampleI
int Pagani2016filter::getObservationCounter()
{
    return timeSampI;
}

//--------------------------------------------------------------------------
/// Return state mean
int Pagani2016filter::getObservationTimestep(int _timeSampleI)
{
    return observationStart + observationDelta * _timeSampleI;
}

//--------------------------------------------------------------------------
/// Set the observations matrix
void Pagani2016filter::setTrueObservations(Eigen::MatrixXd _observations)
{
    trueObservations = _observations;
    std::cout << "trueObservations = \n" << trueObservations << std::endl;
    for(int colI = 0; colI < trueObservations.cols(); colI++)
    {
        trueObservations.col(colI) += measNoiseDensity->Sample();
    }
    std::cout << "trueObservations w. error = \n" << trueObservations << std::endl;
}

//--------------------------------------------------------------------------
/// Setup the time vector 
void Pagani2016filter::setTime(double _startTime, double _deltaTime, double _endTime)
{
    M_Assert(_endTime > _startTime, "endTime must be bigger than startTime");
    startTime = _startTime;
    deltaTime = _deltaTime;
    endTime = _endTime;
    Ntimes = round((endTime - startTime) / deltaTime);
    Info << "startTime = " << startTime << endl;
    Info << "endTime = " << endTime << endl;
    Info << "deltaTime = " << deltaTime << endl;
    Info << "Ntimes = " << Ntimes << endl;
    timeVector = Eigen::VectorXd::LinSpaced(Ntimes + 1, startTime, endTime);
}

//--------------------------------------------------------------------------
/// Setup the observation vector 
void Pagani2016filter::setObservationTime(int _observationStart, int _observationDelta)
{
    M_Assert(timeVector.size() > 0, "Setup the timeVector before setting up the observations vector");
    M_Assert(_observationStart > 0, "First observation timestep can't be 0");
    observationStart = _observationStart;
    observationDelta = _observationDelta;
    Info << "First observation at time = " << timeVector(observationStart) << " s" << endl;
    Info << "Observations taken every " << observationDelta << " timesteps" << endl;
    observationBoolVec = Eigen::VectorXi::Zero(timeVector.size() - 1);
    for(int i = observationStart - 1; i < Ntimes; i += observationDelta)
    {
        observationBoolVec(i) = 1;
    }
}

//--------------------------------------------------------------------------
/// Setup of the measurement noise distribution
void Pagani2016filter::setMeasNoise(double cov)
{
    M_Assert(observationSize > 0, "Read measurements before setting up the measurement noise");
    Eigen::VectorXd measNoise_mu = Eigen::VectorXd::Zero(observationSize);
    Eigen::MatrixXd measNoise_cov = Eigen::MatrixXd::Identity(observationSize,
                                    observationSize) * cov;
    measNoiseDensity = std::make_shared<muq::Modeling::Gaussian>(measNoise_mu,
                       measNoise_cov);
    measurementNoiseFlag = 1;
}

//--------------------------------------------------------------------------
/// Setup of the measurement noise distribution
void Pagani2016filter::setParameterError(double cov)
{
    M_Assert(parameterSize > 0, 
            "Set parameter size before setting up the parameter error");
    Eigen::VectorXd parameterError_mu = Eigen::VectorXd::Zero(parameterSize);
    Eigen::MatrixXd parameterError_cov = Eigen::MatrixXd::Identity(parameterSize,
                                    parameterSize) * cov;
    parameterErrorDensity = std::make_shared<muq::Modeling::Gaussian>(parameterError_mu,
                       parameterError_cov);
}

//--------------------------------------------------------------------------
/// Create initial state ensemble 
void Pagani2016filter::setInitialStateDensity(Eigen::VectorXd _mean, Eigen::MatrixXd _cov)
{
    if(stateSize == 0)
    {
        stateSize = _mean.size();
    }
    else
    {
        std::string message = "State has size = " + std::to_string(stateSize)
            + " while input mean vector has size = "
            + std::to_string(_mean.size());

        M_Assert(stateSize == _mean.size(), message.c_str());
    }
    M_Assert(_cov.rows() == stateSize && _cov.cols() == stateSize, "To initialize the state ensemble use mean and cov with same dimentions");
    initialStateDensity = std::make_shared<muq::Modeling::Gaussian>(_mean, _cov);
    initialStateFlag = 1;
}

//--------------------------------------------------------------------------
/// Create initial state ensemble 
void Pagani2016filter::sampleInitialState()
{
    M_Assert(initialStateFlag == 1, "Initialize the initial state density before sampling it");
    stateEns.assignSamples(ensembleFromDensity(initialStateDensity));
}

//--------------------------------------------------------------------------
/// Create initial parameter ensemble 
void Pagani2016filter::sampleInitialParameter()
{
    M_Assert(parameterPriorFlag == 1, "Initialize the initial parameter density before sampling it");
    parameterEns.assignSamples(ensembleFromDensity(parameterPriorDensity));
}

//--------------------------------------------------------------------------
/// Create parameter ensemble 
void Pagani2016filter::setParameterPriorDensity(Eigen::VectorXd _mean, Eigen::MatrixXd _cov)
{
    if(parameterSize == 0)
    {
        parameterSize = _mean.size();
    }
    else
    {
        std::string message = "The input mean has size = " + std::to_string(_mean.size())
            + " but the parameterSize is "
            + std::to_string(parameterSize);

        M_Assert(parameterSize == _mean.size(), message.c_str());
    }
    M_Assert(_cov.rows() == parameterSize && _cov.cols() == parameterSize, "Use mean and cov with same dimentions");
    parameterPriorMean = _mean;
    parameterPriorCov = _cov;
    parameterPriorDensity = std::make_shared<muq::Modeling::Gaussian>(_mean, _cov);
    parameterPriorFlag = 1;
}

//--------------------------------------------------------------------------
/// General class to sample from an input density 
Eigen::MatrixXd Pagani2016filter::ensembleFromDensity(std::shared_ptr<muq::Modeling::Gaussian> _density)
{
    M_Assert(Nseeds > 0, "Number of samples not set up correctly");
    Eigen::MatrixXd output(_density->Sample().size(), Nseeds);
    for (int i = 0; i < Nseeds; i++)
    {
        output.col(i) = _density->Sample();
    }

    return output;
}

//--------------------------------------------------------------------------
/// 
void Pagani2016filter::setObservationSize(int _size)
{
    observationSize = _size;
    Info << "Observation size = " << observationSize << endl;
}

//--------------------------------------------------------------------------
/// 
void Pagani2016filter::setStateSize(int _size)
{
    stateSize = _size;
    Info << "State size = " << stateSize << endl;
}

//--------------------------------------------------------------------------
/// 
int Pagani2016filter::getStateSize()
{
    return stateSize;
}

//--------------------------------------------------------------------------
/// 
void Pagani2016filter::setParameterSize(int _size)
{
    parameterSize = _size;
    Info << "Parameter size = " << parameterSize << endl;
}

//--------------------------------------------------------------------------
/// 
int Pagani2016filter::getParameterSize()
{
    return parameterSize;
}

//--------------------------------------------------------------------------
/// Perform Kalman filter update 
void Pagani2016filter::update()
{
    Eigen::MatrixXd parameterObservation_crossCov = 
        parameterEns.crossCov(observationEns.getSamples());
    Eigen::MatrixXd stateObservation_crossCov = 
        stateEns.crossCov(observationEns.getSamples());
    Eigen::MatrixXd observation_Cov = observationEns.cov();

    Eigen::MatrixXd P = measNoiseDensity->ApplyCovariance(
            Eigen::MatrixXd::Identity(observationSize, observationSize));
    P += observation_Cov;
    std::cout << "debug : P = \n" << P << std::endl;
    P = P.inverse(); //TODO check if invertible
    std::cout << "debug : P inverse = \n" << P << std::endl;
    
    Eigen::MatrixXd parameterUpdateMat = parameterObservation_crossCov * P;
    Eigen::MatrixXd stateUpdateMat = stateObservation_crossCov * P;
    for(int seedI = 0; seedI < Nseeds; seedI++)
    {
        Eigen::VectorXd observationDelta = 
            trueObservations.col(timeSampI) - observationEns.getSample(seedI);
        parameterEns.assignSample(seedI, 
            parameterEns.getSample(seedI) + parameterUpdateMat * observationDelta);
        stateEns.assignSample(seedI, 
            stateEns.getSample(seedI) + stateUpdateMat * observationDelta);
    }
    std::cout << "debug : parameterEns.mean() =\n" << parameterEns.mean() << std::endl;
    std::cout << "debug : stateEns.mean() =\n" << stateEns.mean() << std::endl;
}

//--------------------------------------------------------------------------
/// Run the filtering
void Pagani2016filter::run(word outputFolder)
{
    M_Assert(initialStateFlag == 1, "Set up the initial state density");
    M_Assert(parameterPriorFlag == 1, "Set up the parameter prior density");
    M_Assert(measurementNoiseFlag == 1, "Set up the measurement noise");
    M_Assert(stateSize > 0, "Set state size");
    M_Assert(observationSize > 0, "Set observation size");
    M_Assert(parameterSize > 0, "Set parameter size");

    // Initialization
    sampleInitialState();
    sampleInitialParameter();
    int NtimeObservations = observationBoolVec.sum();
    stateMean.resize(stateSize, Ntimes);
    state_maxConf.resize(stateSize, Ntimes);
    state_minConf.resize(stateSize, Ntimes);
    parameter_maxConf.resize(parameterSize, Ntimes);
    parameter_minConf.resize(parameterSize, Ntimes);
    parameterMean.resize(parameterSize, Ntimes);


    for(timeSampI = 0; timeSampI < NtimeObservations; timeSampI++)
    {
        Info << "timeSamp " << timeSampI << endl;
        oldStateEns.assignSamples(stateEns.getSamples());

        stateProjection();

        update();
    }
    ITHACAstream::exportMatrix(stateMean, "stateMean", "eigen", outputFolder);
    ITHACAstream::exportMatrix(parameterMean, "parameterMean", "eigen", outputFolder);
    ITHACAstream::exportMatrix(parameter_maxConf, "parameter_maxConf", "eigen", 
            outputFolder);
    ITHACAstream::exportMatrix(parameter_minConf, "parameter_minConf", "eigen", 
            outputFolder);
    Info << "\n*****************************************************" << endl;
    Info << "Kalman filter run ENDED" << endl;
    Info << "\n*****************************************************" << endl;

}
}
