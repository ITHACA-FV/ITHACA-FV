# Tutorial UQ-01

## Introduction
This small example tests the Ensemble Kalman Filter (EnKF) wrapper included in ITHACA‑FV/MUQ. A simple linear dynamical system $\dot{x} = A x$ is considered and only partial observations are available
through a matrix $H$. The goal is to reconstruct the full state starting from an incorrect model matrix $A_{wrong}$.

## Problem setup
The state evolves according to

$$ \frac{dx}{dt} = A x $$

with initial condition $x(0)=x_0$. Observations at discrete times are obtained as

$$ y = H x. $$

Matrices required by the tutorial are stored in text files under the tutorial root:

* `A_mat.txt` – true system matrix
* `Awrong_mat.txt` – incorrect matrix used for forecasting
* `observation_mat.txt` – observation operator $H$
* `initialState_mat.txt` – initial state vector

The code generates a synthetic trajectory, samples it every few steps and then applies the EnKF to produce posterior ensembles, means and confidence intervals.

# A detailed look into the code
Let's look at the code of tutorial UQ-01.

## The header files
```cpp
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
```

## The main
At the beginning, the code reads the system matrices `A`, `Aw`, and `H`, along with the initial state `x0`.
```cpp
    word outputFolder = "./ITHACAoutput/";
    int Nseeds = 1000;
    Eigen::MatrixXd A = ITHACAstream::readMatrix("A_mat.txt");
    Eigen::MatrixXd Aw = ITHACAstream::readMatrix("Awrong_mat.txt");
    Eigen::MatrixXd H = ITHACAstream::readMatrix("observation_mat.txt");
    M_Assert(ITHACAstream::readMatrix("initialState_mat.txt").cols() == 1,
             "Wrong initialState input");
    Eigen::VectorXd x0 = ITHACAstream::readMatrix("initialState_mat.txt").col(0);
    ...
```

It then defines the time discretization and simulates the true system forward in time using an explicit Euler scheme:
```cpp
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

```

During this simulation, observations are collected at regular intervals defined by `sampleDeltaStep`. The full trajectory and the sampled observations are stored for later use.

```cpp
    M_Assert(Nsamples == sampleI, "Something went wrong in the sampling");
    ITHACAstream::exportMatrix(time, "time", "eigen", outputFolder);
    ITHACAstream::exportMatrix(X, "X", "eigen", outputFolder);
```

The Bayesian filtering framework is then initialized. A Gaussian prior distribution is defined for the initial state, together with a Gaussian model-error distribution that represents uncertainty in the system dynamics.
```cpp
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
```

An ensemble of samples (`priorSamples`) is drawn from the prior distribution and used as the initial ensemble for the filter.
```cpp
    Eigen::MatrixXd posteriorSamples(stateSize, Nseeds);
    Eigen::MatrixXd priorSamples(stateSize, Nseeds);

    for (int i = 0; i < Nseeds; i++)
    {
        priorSamples.col(i) = priorDensity->Sample();
    }

    posteriorSamples = priorSamples;
```

Measurement noise is also modeled as a Gaussian distribution, and storage is allocated for the posterior mean and confidence intervals over time.
```cpp
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
```

The algorithm then enters the main EnKF loop, which iterates over all time steps.
```cpp
    for (int timeI = 0; timeI < Ntimes - 1; timeI++)
    {
        priorSamples = posteriorSamples;
        Eigen::MatrixXd forwardSamplesOld = forwardSamples;


```

At each time step, the forecast step propagates each ensemble member forward using the (possibly incorrect) system matrix and adds a realization of the model error (`modelErrorDensity`).
```cpp
        //Forecast step
        for (int i = 0; i < Nseeds; i++)
        {
            forwardSamples.col(i) = (A * deltaTime + Eigen::MatrixXd::Identity(A.rows(),
                                     A.cols())) * priorSamples.col(i) + modelErrorDensity->Sample();
        }

        sampleFlag--;
```

When a measurement is available, the update step is performed using the Ensemble Kalman Filter, which corrects the forecast ensemble based on the observed data and the observation operator.
```cpp
        if (sampleFlag == 0)
        {
            sampleFlag = sampleDeltaStep;
            Eigen::VectorXd meas = obs.col(sampleI);
            //Kalman filter
            posteriorSamples = ITHACAmuq::muq2ithaca::EnsembleKalmanFilter(forwardSamples,
                               meas, meas_cov, H * forwardSamples);
            sampleI++;
        }
```

If no measurement is available at a given step, the forecast ensemble is directly used as the posterior.
```cpp
        else
        {
            posteriorSamples = forwardSamples;
        }
```

After each update, the posterior mean is computed as the average of the ensemble, and confidence intervals are estimated using quantiles of the ensemble distribution. These quantities provide both the estimated state and a measure of uncertainty.
```cpp
        posteriorMean.col(timeI + 1) = posteriorSamples.rowwise().mean();
        minConfidence.col(timeI + 1) = ITHACAmuq::muq2ithaca::quantile(posteriorSamples,
                                       0.05);
        maxConfidence.col(timeI + 1) = ITHACAmuq::muq2ithaca::quantile(posteriorSamples,
                                       0.95);
    }
```

Finally, the code exports the time evolution of the posterior mean and the confidence bounds, allowing comparison with the true solution.
```cpp
    ITHACAstream::exportMatrix(posteriorMean, "posteriorMean", "eigen",
                               outputFolder);
    ITHACAstream::exportMatrix(minConfidence, "minConfidence", "eigen",
                               outputFolder);
    ITHACAstream::exportMatrix(maxConfidence, "maxConfidence", "eigen",
                               outputFolder);
```

## The plain code
The plain code is available [here](https://raw.githubusercontent.com/ITHACA-FV/ITHACA-FV/refs/heads/master/tutorials/UQ/01enKF/01enKF.C).