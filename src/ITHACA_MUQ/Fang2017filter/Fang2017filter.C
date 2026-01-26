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
/// Source file of the Fang2017filter class.

#include "Fang2017filter.H"

namespace ITHACAmuq
{
// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //
// Constructor
Fang2017filter::Fang2017filter() {}
Fang2017filter::Fang2017filter(int _Nsamples)
{
    Nsamples = _Nsamples;
}

// * * * * * * * * * * * * * * Filter Methods * * * * * * * * * * * * * * //

// Return number of samples per ensamble
int Fang2017filter::getNumberOfSamples()
{
    return Nsamples;
}

// Return time
double Fang2017filter::getTime()
{
    return timeVector(timeStepI);
}

//--------------------------------------------------------------------------
/// Return time for input timestep
double Fang2017filter::getTime(int _timeStepI)
{
    return timeVector(_timeStepI);
}

//--------------------------------------------------------------------------
/// Return timestep
int Fang2017filter::getTimeStep()
{
    return timeStepI;
}

//--------------------------------------------------------------------------
/// Return time vector
Eigen::VectorXd Fang2017filter::getTimeVector()
{
    return timeVector;
}

//--------------------------------------------------------------------------
/// Return state mean
Eigen::MatrixXd Fang2017filter::getStateMean()
{
    return stateMean;
}

//--------------------------------------------------------------------------
/// Set the observations matrix
void Fang2017filter::setObservations(Eigen::MatrixXd _observations)
{
    observations = _observations;
}

//--------------------------------------------------------------------------
/// Setup the time vector
void Fang2017filter::setTime(double _startTime, double _deltaTime,
                             double _endTime)
{
    M_Assert(_endTime > _startTime, "endTime must be bigger than startTime");
    startTime = _startTime;
    deltaTime = _deltaTime;
    endTime = _endTime;
    Ntimes = (endTime - startTime) / deltaTime;
    Info << "startTime = " << startTime << endl;
    Info << "endTime = " << endTime << endl;
    Info << "deltaTime = " << deltaTime << endl;
    Info << "Ntimes = " << Ntimes << endl;
    timeVector = Eigen::VectorXd::LinSpaced(Ntimes + 1, startTime, endTime);
}

//--------------------------------------------------------------------------
/// Setup the observation vector
void Fang2017filter::setObservationTime(int _observationStart,
                                        int _observationDelta)
{
    M_Assert(timeVector.size() > 0,
             "Setup the timeVector before setting up the observations vector");
    M_Assert(_observationStart > 0, "First observation timestep can't be 0");
    observationStart = _observationStart;
    observationDelta = _observationDelta;
    Info << "First observation at time = " << timeVector(observationStart) << " s"
         << endl;
    Info << "Observations taken every " << observationDelta << " timesteps" << endl;
    observationBoolVec = Eigen::VectorXi::Zero(timeVector.size() - 1);

    for (int i = observationStart - 1; i < Ntimes; i += observationDelta)
    {
        observationBoolVec(i) = 1;
    }
}

//--------------------------------------------------------------------------
/// Setup of the model error distribution
void Fang2017filter::setModelError(double cov, bool univariate)
{
    M_Assert(stateSize > 0, "Set the stateSize before setting up the model error");
    Eigen::VectorXd modelError_mu;
    Eigen::MatrixXd modelError_cov;

    if (univariate)
    {
        modelError_mu = Eigen::VectorXd::Zero(1);
        modelError_cov = Eigen::MatrixXd::Identity(1,
                         1) * cov;
    }
    else
    {
        modelError_mu = Eigen::VectorXd::Zero(stateSize);
        modelError_cov = Eigen::MatrixXd::Identity(stateSize,
                         stateSize) * cov;
    }

    modelErrorDensity = std::make_shared<muq::Modeling::Gaussian>(modelError_mu,
                        modelError_cov);
    modelErrorFlag = 1;
}

//--------------------------------------------------------------------------
/// Setup of the measurement noise distribution
void Fang2017filter::setMeasNoise(double cov)
{
    M_Assert(observationSize > 0,
             "Read measurements before setting up the measurement noise");
    Eigen::VectorXd measNoise_mu = Eigen::VectorXd::Zero(observationSize);
    Eigen::MatrixXd measNoise_cov = Eigen::MatrixXd::Identity(observationSize,
                                    observationSize) * cov;
    measNoiseDensity = std::make_shared<muq::Modeling::Gaussian>(measNoise_mu,
                       measNoise_cov);
    measurementNoiseFlag = 1;
}

//--------------------------------------------------------------------------
/// Create initial state ensemble
void Fang2017filter::setInitialStateDensity(Eigen::VectorXd _mean,
        Eigen::MatrixXd _cov)
{
    if (stateSize == 0)
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

    M_Assert(_cov.rows() == stateSize
             && _cov.cols() == stateSize,
             "To initialize the state ensemble use mean and cov with same dimentions");
    initialStateDensity = std::make_shared<muq::Modeling::Gaussian>(_mean, _cov);
    initialStateFlag = 1;
}

//--------------------------------------------------------------------------
/// Create initial state ensemble
void Fang2017filter::sampleInitialState()
{
    M_Assert(initialStateFlag == 1,
             "Initialize the initial state density before sampling it");
    stateEns.assignSamples(ensembleFromDensity(initialStateDensity));
}

//--------------------------------------------------------------------------
/// Create parameter ensemble
void Fang2017filter::setParameterPriorDensity(Eigen::VectorXd _mean,
        Eigen::MatrixXd _cov)
{
    if (parameterSize == 0)
    {
        parameterSize = _mean.size();
    }
    else
    {
        std::string message = "The input mean has size = " + std::to_string(
                                  _mean.size())
                              + " but the parameterSize is "
                              + std::to_string(parameterSize);
        M_Assert(parameterSize == _mean.size(), message.c_str());
    }

    M_Assert(_cov.rows() == parameterSize
             && _cov.cols() == parameterSize, "Use mean and cov with same dimentions");
    parameterPriorMean = _mean;
    parameterPriorCov = _cov;
    parameterPriorDensity = std::make_shared<muq::Modeling::Gaussian>(_mean, _cov);
    parameterPriorFlag = 1;
}

//--------------------------------------------------------------------------
/// Create parameter ensemble
void Fang2017filter::sampleParameterDist()
{
    M_Assert(parameterPriorFlag == 1, "Set up the parameter prior density");
    parameterEns.assignSamples(ensembleFromDensity(parameterPriorDensity));
    parameterMean.col(timeStepI) = parameterEns.mean();
}

//--------------------------------------------------------------------------
/// General class to sample from an input density
Eigen::MatrixXd Fang2017filter::ensembleFromDensity(
    std::shared_ptr<muq::Modeling::Gaussian> _density)
{
    M_Assert(Nsamples > 0, "Number of samples not set up correctly");
    Eigen::MatrixXd output(_density->Sample().size(), Nsamples);

    for (int i = 0; i < Nsamples; i++)
    {
        output.col(i) = _density->Sample();
    }

    return output;
}

//--------------------------------------------------------------------------
/// Concatenate state and parameter ensambles to create the joint ensamble
void Fang2017filter::buildJointEns()
{
    Eigen::MatrixXd stateSamps = stateEns.getSamples();
    Eigen::MatrixXd paramSamps = parameterEns.getSamples();
    M_Assert(stateSamps.cols() == paramSamps.cols(),
             "State and parameter ensambles must have same number of samples");
    Eigen::MatrixXd temp(stateSamps.rows() + paramSamps.rows(), stateSamps.cols());
    temp << stateSamps,
         paramSamps;
    jointEns.assignSamples(temp);
}

//--------------------------------------------------------------------------
///
void Fang2017filter::setObservationSize(int _size)
{
    observationSize = _size;
    Info << "Observation size = " << observationSize << endl;
}

//--------------------------------------------------------------------------
///
void Fang2017filter::setStateSize(int _size)
{
    stateSize = _size;
    Info << "State size = " << stateSize << endl;
}

//--------------------------------------------------------------------------
///
int Fang2017filter::getStateSize()
{
    return stateSize;
}

//--------------------------------------------------------------------------
///
void Fang2017filter::setParameterSize(int _size)
{
    parameterSize = _size;
    Info << "Parameter size = " << parameterSize << endl;
}

//--------------------------------------------------------------------------
///
int Fang2017filter::getParameterSize()
{
    return parameterSize;
}

//--------------------------------------------------------------------------
///
void Fang2017filter::updateJointEns(Eigen::VectorXd _observation)
{
    M_Assert(_observation.size() == observationSize,
             "Observation has wrong dimentions");
    int ensSize = jointEns.getSize();
    //Eigen::MatrixXd temp = _observation - observationEns.getSamples().colwise();
    //TODO deal with invertibility of observationEns.cov()
    Eigen::MatrixXd autoCovInverse = observationEns.cov().inverse();
    //Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(observationEns.cov());
    //auto P_rank = lu_decomp.rank();
    //if (P_rank < P.cols() || P_rank < P.rows())
    //{
    //    std::cout << "Pseudo inverse of P should be implemented in the EnKF, exiting" <<
    //              std::endl;
    //    std::cout << P << std::endl;
    //    exit(10);
    //}
    //else
    //{
    //    P = P.inverse();
    //}
    Eigen::MatrixXd crossCov = jointEns.crossCov(observationEns.getSamples());

    for (int i = 0; i < ensSize; i++)
    {
        Eigen::VectorXd newSamp = jointEns.getSample(i) + crossCov * autoCovInverse *
                                  (_observation - observationEns.getSample(i));
        jointEns.assignSample(i, newSamp);
    }
}

//--------------------------------------------------------------------------
/// Run the filtering
void Fang2017filter::run(int innerLoopMax, word outputFolder)
{
    M_Assert(initialStateFlag == 1, "Set up the initial state density");
    M_Assert(parameterPriorFlag == 1, "Set up the parameter prior density");
    M_Assert(modelErrorFlag == 1, "Set up the model error");
    M_Assert(measurementNoiseFlag == 1, "Set up the measurement noise");
    M_Assert(stateSize > 0, "Set state size");
    M_Assert(observationSize > 0, "Set observation size");
    M_Assert(parameterSize > 0, "Set parameter size");
    timeStepI = 0;
    sampleInitialState();
    stateMean.resize(stateSize, Ntimes);
    state_maxConf.resize(stateSize, Ntimes);
    state_minConf.resize(stateSize, Ntimes);
    parameter_maxConf.resize(parameterSize, Ntimes);
    parameter_minConf.resize(parameterSize, Ntimes);
    parameterMean.resize(parameterSize, Ntimes);

    while (timeStepI < Ntimes)
    {
        Info << "timeStep " << timeStepI << endl;
        innerLoopI = 0;
        oldStateEns.assignSamples(stateEns.getSamples());

        while (innerLoopI < innerLoopMax)
        {
            Info << "Inner loop " << innerLoopI << endl;

            // Projection
            if (innerLoopI == 0)
            {
                if (timeStepI == 0)
                {
                    setParameterPriorDensity(parameterPriorMean, parameterPriorCov);
                    sampleParameterDist();
                }
                else
                {
                    setParameterPriorDensity(
                        parameterMean.col(timeStepI - 1), parameterPriorCov);
                    sampleParameterDist();
                }

                std::cout << "\ndebug : parameterPriorMean = " <<
                          parameterPriorMean << std::endl;
                std::cout << "\ndebug : parameterMean.col(" << timeStepI << ") =\n" <<
                          parameterMean.col(timeStepI) << std::endl;
            }
            else
            {
                std::cout << "\ndebug : parameterMean before loop =\n" <<
                          parameterMean.col(timeStepI) << std::endl;
                setParameterPriorDensity(parameterMean.col(timeStepI), parameterPriorCov);
                sampleParameterDist();
                std::cout << "\ndebug : parameterMean after loop =\n" <<
                          parameterMean.col(timeStepI) << std::endl;
            }

            stateProjection();
            buildJointEns();

            if (observationBoolVec(timeStepI) == 1)
            {
                Eigen::MatrixXd measNoiseSamps = ensembleFromDensity(measNoiseDensity);
                observeState();
                std::cout << "\ndebug : observation =\n" <<
                          observations.col(observationBoolVec.head(timeStepI + 1).sum() - 1) <<
                          std::endl;
                updateJointEns(
                    observations.col(
                        observationBoolVec.head(timeStepI + 1).sum() - 1));
            }

            parameterMean.col(timeStepI) = jointEns.mean().tail(parameterSize);
            innerLoopI++;
        }

        stateEns.assignSamples(jointEns.getSamples().topRows(stateSize));
        parameterEns.assignSamples(jointEns.getSamples().bottomRows(parameterSize));
        stateMean.col(timeStepI) = jointEns.mean().head(stateSize);
        state_maxConf.col(timeStepI) = muq2ithaca::quantile(stateEns.getSamples(),
                                       0.95);
        state_minConf.col(timeStepI) = muq2ithaca::quantile(stateEns.getSamples(),
                                       0.05);
        parameter_maxConf.col(timeStepI) = muq2ithaca::quantile(
                                               parameterEns.getSamples(),
                                               0.95);
        parameter_minConf.col(timeStepI) = muq2ithaca::quantile(
                                               parameterEns.getSamples(),
                                               0.05);
        timeStepI++;
    }

    ITHACAstream::exportMatrix(stateMean, "stateMean", "eigen", outputFolder);
    ITHACAstream::exportMatrix(parameterMean, "parameterMean", "eigen",
                               outputFolder);
    ITHACAstream::exportMatrix(parameter_maxConf, "parameter_maxConf", "eigen",
                               outputFolder);
    ITHACAstream::exportMatrix(parameter_minConf, "parameter_minConf", "eigen",
                               outputFolder);
}
}
