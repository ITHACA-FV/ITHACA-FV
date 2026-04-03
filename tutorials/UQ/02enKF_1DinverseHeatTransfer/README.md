# Tutorial UQ-02

## Introduction
This example demonstrates the use of the Ensemble Kalman Filter to
reconstruct the state of a one–dimensional heat conduction problem when the boundary flux is only partially known and noisy measurements of the
temperature are available inside the domain.

The mathematical model for the true problem is

$$ \rho C_p \frac{\partial T}{\partial t} + k \Delta T = 0, \quad x \in(0,1) $$

with a time‑dependent Neumann condition at $x=0$ and a fixed Dirichlet
value of 0.5 at $x=1$. The `true` flux $g_{true}(t)=5t+2$ differs from the model flux $g(t)=5t$ used in the
filter. Initial temperature $T_0=1$.

Measurements (noisy temperature values) are specified in `measurementsDict` and read by the solver. The state and error covariances are assumed Gaussian.

# A detailed look into the code
Let's look at the code of tutorial UQ-02.

## The header files

```cpp
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
```
The included headers provide OpenFOAM finite-volume tools, time and SIMPLE-control utilities, the custom Laplacian problem class, ITHACA-FV helper functions, Eigen support, the MUQ Gaussian distribution class, and the interface between MUQ and ITHACA for Ensemble Kalman filtering.
The tutorial-specific header `02enKF_1DinverseHeatTransfer.H` contains the implementation of the example class used in the driver.

## The `EnKF_1DinverseHeatTransfer` class
The `EnKF_1DinverseHeatTransfer` class is derived from `laplacianProblem` and specializes it for a 1D transient heat-transfer problem solved together with an Ensemble Kalman Filter for state reconstruction.

The constructor initializes the parent problem, creates the simpleControl object, reads the index of the left boundary patch, stores the number of ensemble members, allocates the temperature ensemble, and initializes the data structures used by the filter.
```cpp
class EnKF_1DinverseHeatTransfer: public laplacianProblem
{
    public:
        explicit EnKF_1DinverseHeatTransfer(int argc, char* argv[], int Nseeds_)
            :
            laplacianProblem(argc, argv)
        {
            fvMesh& mesh = _mesh();
            _simple = autoPtr<simpleControl>
                      (
                          new simpleControl
                          (
                              mesh
                          )
                      );
            simpleControl& simple = _simple();
            Time& runTime = _runTime();
#include "createFvOptions.H"
            left_ind = mesh.boundaryMesh().findPatchID("left");
            Nseeds = Nseeds_;
            Tensemble = PtrList<volScalarField>(Nseeds);
            volScalarField& T(_T());
            stateSize = T.size();
            priorSamples = Eigen::MatrixXd::Zero(stateSize, Nseeds);
            posteriorSamples = priorSamples;
```

It also reads the time information from the OpenFOAM runtime and prepares vectors and matrices for storing the time history, measurements, and confidence intervals at the probe location.
```cpp
            startTime = runTime.startTime().value();
            deltaTime = runTime.deltaTValue();
            endTime = runTime.endTime().value();
            Ntimes = (endTime - startTime) / deltaTime;
            Info << "startTime = " << startTime << endl;
            Info << "endTime = " << endTime << endl;
            Info << "deltaTime = " << deltaTime << endl;
            Info << "Ntimes = " << Ntimes << endl;
            timeVector.resize(Ntimes);
            probe_mean.resize(Nprobes, Ntimes);
            probe_true = probe_mean;
            probe_MaxConfidence = probe_mean;
            probe_minConfidence = probe_mean;
```

The class stores the physical parameters of the problem, namely thermal conductivity, density, and heat capacity, together with the boundary index of the left patch. 
```cpp
        autoPtr<simpleControl> _simple;
        autoPtr<fv::options> _fvOptions;
        /// Diffusivity
        double k;
        /// Density
        double rho;
        /// Heat capacity
        double Cp;
        /// Index of the left patch
        label left_ind;
        Eigen::MatrixXd measurements;
        /// Number of timesteps between each sampling
        int NtimeStepsBetweenSamps = 3;
        /// Number of seeds in the ensemble
        int Nseeds;

        /// Ensemble
        PtrList<volScalarField> Tensemble;
        /// Mean field of the ensemble
        PtrList<volScalarField> Tmean;
        int stateSize;
```

It also contains the matrices and Gaussian probability models needed by the filter: prior distribution, model-error distribution, measurement-noise distribution, prior samples, and posterior samples.
```cpp
        std::shared_ptr<muq::Modeling::Gaussian> priorDensity;
        std::shared_ptr<muq::Modeling::Gaussian> modelErrorDensity;
        std::shared_ptr<muq::Modeling::Gaussian> measNoiseDensity;
        double meas_cov;

        Eigen::MatrixXd priorSamples;
        Eigen::MatrixXd posteriorSamples;
```

Additional members are used to manage sampling in time, define the probe position, and store the true, mean, and confidence-band values of the reconstructed temperature.
```cpp
        double startTime;
        double deltaTime;
        double endTime;
        int Ntimes;

        /// Timesteps vector
        Eigen::VectorXd timeVector;
        /// Number of probes
        int Nprobes = 1;
        /// Probles position, probes are used only for visualization and postprocessing
        Foam::vector probePosition = Foam::vector(0.5, 0.1, 0.005);
        Eigen::MatrixXd probe_true;
        Eigen::MatrixXd probe_mean;
        Eigen::MatrixXd probe_MaxConfidence;
        Eigen::MatrixXd probe_minConfidence;

        int samplingTimeI = 0;
        int sampleI = 0;
        int sampleFlag = 0;
```

The method `updateBC` defines the time-dependent boundary condition applied on the left boundary. In the direct simulation it returns the true boundary heat flux:
```cpp
        double updateBC(bool reconstruction = 0)
        {
            Time& runTime = _runTime();
            scalar time = runTime.value();
            double BC = time * 5 + 2;
            ...
            return BC;
```

while in reconstruction mode it intentionally uses a wrong boundary condition:
```cpp
            if (reconstruction)
            {
                BC = time * 5 + 0;
            }

            return BC;
```

This mismatch is what makes the problem an inverse reconstruction problem rather than a direct simulation.

The method `resetSamplingCounters` resets the counters used to keep track of sampling times and measurement indices during the simulation. This is needed whenever a new direct or reconstruction phase starts.
```cpp
        void resetSamplingCounters()
        {
            samplingTimeI = 0;
            sampleI = 0;
            sampleFlag = 0;
        }
```

The method `observe` returns the values of a temperature field at the measurement locations defined in `measurementsDict`. It reads the probe positions from the dictionary, finds the corresponding mesh cells, and extracts the field values there. These sampled values represent the observations assimilated by the filter.
```cpp
        Eigen::VectorXd observe(volScalarField field)
        {
            fvMesh& mesh = _mesh();
            Time& runTime = _runTime();
            IOdictionary measurementsDict
            (
                IOobject
                (
                    "measurementsDict",
                    runTime.constant(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );
            List<vector> measurementPoints(measurementsDict.lookup("positions"));
            Eigen::VectorXd measures(measurementPoints.size());
            forAll(measurementPoints, pntI)
            {
                measures(pntI) = field[mesh.findCell(measurementPoints[pntI])];
            }
            return measures;
        }
```

The method `solve` is the core time-marching routine. It solves the transient heat equation for a given temperature field, applies the current boundary condition, advances the solution in time, and exports the field at each step. In direct mode it also stores the true probe signal and builds the matrix of observations at the sampling times.
```cpp
        /// Solve the heat transfer problem
        Eigen::MatrixXd solve(volScalarField& T, word outputFolder,
                              word outputFieldName, bool reconstruction = 0)
        {
            Eigen::MatrixXd obsMat;
            fvMesh& mesh = _mesh();
            Foam::Time& runTime = _runTime();
            simpleControl& simple = _simple();
            fv::options& fvOptions(_fvOptions());
            dimensionedScalar diffusivity("diffusivity", dimensionSet(0, 2, -1, 0, 0, 0, 0),
                                          k / (rho * Cp));
            Info << "Thermal diffusivity = " << diffusivity << " m2/s" << endl;
            int timeI = 0;

            while (runTime.loop())
            {
                scalar time = runTime.value();

                if (!reconstruction)
                {
                    sampleFlag++;
                    timeVector(timeI) = time;
                    timeI++;
                }

                // Update BC
                ITHACAutilities::assignBC(T, mesh.boundaryMesh().findPatchID("left"),
                                          updateBC(reconstruction));

                // Solve
                while (simple.correctNonOrthogonal())
                {
                    fvScalarMatrix TEqn
                    (
                        fvm::ddt(T) - fvm::laplacian(diffusivity, T)
                    );
                    fvOptions.constrain(TEqn);
                    TEqn.solve();
                    fvOptions.correct(T);
                }
```

In reconstruction mode it adds a realization of the model error to the temperature field after each solve, mimicking uncertainty in the predictive model.
```cpp
                if (reconstruction)
                {
                    //Adding model error
                    forAll(T.internalField(), cellI)
                    {
                        T.ref()[cellI] += modelErrorDensity->Sample()(cellI);
                    }
                }
                else
                {
                    probe_true.col(timeI - 1) = sampleField(T, probePosition);
                }

```
The method returns the observation matrix associated with the current run.
```cpp
if (sampleFlag == NtimeStepsBetweenSamps)
                {
                    Info << "Sampling at time = " << runTime.timeName() << nl << endl;
                    obsMat.conservativeResize(observe(T).size(), obsMat.cols() + 1);
                    obsMat.col(sampleI) = observe(T);

                    if (!reconstruction)
                    {
                        sampleFlag = 0;
                        sampleI++;
                    }
                }
               ...
               return obsMat;
```

The method `solveDirect` performs the true simulation of the heat-transfer problem. It restarts the solver, runs the direct model with the correct boundary condition, stores the resulting measurements, exports them, and resets the sampling counters. This phase generates the synthetic data that will later be assimilated during reconstruction.
```cpp
        void solveDirect()
        {
            Info << "\n****************************************\n" << endl;
            Info << "\nPerforming true solution\n" << endl;
            word outputFolder = "./ITHACAoutput/direct/";
            restart();
            volScalarField& T(_T());
            /// The measurements are obtained as observation of the true solution
            measurements = solve(T, outputFolder, "Tdirect");
            ITHACAstream::exportMatrix(measurements, "trueMeasurements", "eigen",
                                       outputFolder);
            std::cout << "Number of samples in time = " << measurements.cols() << std::endl;
            std::cout << "Number of sample in space = " << measurements.rows() << std::endl;
            resetSamplingCounters();
            Info << "\nEND true solution\n" << endl;
            Info << "\n****************************************\n" << endl;
        }
```

The method `reconstruct` performs the Ensemble Kalman Filter reconstruction. For each time step, it first calls `forecastStep`, which advances all ensemble members using the imperfect model.
```cpp
        void reconstruct()
        {
            restart();
            Time& runTime = _runTime();
            resetSamplingCounters();
            std::vector<Eigen::MatrixXd> observationTensor;
            word outputFolder = "ITHACAoutput/reconstuction";
            Tmean.resize(Ntimes);
            Info << "\n****************************************\n" << endl;
            Info << "\nStarting reconstruction\n" << endl;
            int samplesCounter = 0;
            // Read true field
            PtrList<volScalarField> Ttrue;
            ITHACAstream::read_fields(Ttrue, "Tdirect",
                                      "ITHACAoutput/direct/");

            for (int timeI = 0; timeI < Ntimes; timeI++)
            {
                double initialTime = timeI * deltaTime + startTime;
                Info << "\n\nTime " << initialTime + deltaTime << endl;
                Eigen::MatrixXd forecastOutput = forecastStep(initialTime, timeI,
                                                 initialTime + deltaTime);

```

If the current time is a sampling step, the measured data are used in the Kalman update to compute corrected posterior samples. These updated states are then assigned back to the ensemble fields.
```cpp
                if (forecastOutput.cols() > 0)
                {
                    // This is a sampling step
                    samplesCounter++;
                    observationTensor.push_back(forecastOutput);
                    sampleFlag = 0;
                    Eigen::VectorXd meas = measurements.col(samplesCounter - 1);
                    // Kalman filter
                    posteriorSamples = ITHACAmuq::muq2ithaca::EnsembleKalmanFilter(Tensemble, meas,
                                       Eigen::MatrixXd::Identity(meas.size(), meas.size()) * meas_cov, forecastOutput);

                    for (int i = 0; i < Nseeds; i++)
                    {
                        restart();
                        volScalarField& T = _T();
                        Eigen::VectorXd internalField = posteriorSamples.col(i);
                        assignIF(T, internalField);
                        Tensemble.Foam::PtrList<volScalarField>::set(i, T.clone());
                    }
                }
                ...
                Eigen::VectorXd internalField = Foam2Eigen::PtrList2Eigen(Tensemble).rowwise().mean();
                volScalarField& Tmean = _T();
                assignIF(Tmean, internalField);

```

At every time step, the method also computes the ensemble mean and the confidence bounds at the chosen probe location and exports both the reconstructed mean field and the true field for comparison.
```cpp
                volScalarField& Tmean = _T();
                assignIF(Tmean, internalField);
                ITHACAstream::exportSolution(Tmean, std::to_string(initialTime + deltaTime),
                                             outputFolder, "Tmean");
                ITHACAstream::exportSolution(Ttrue[timeI],
                                             std::to_string(initialTime + deltaTime),
                                             outputFolder, "Ttrue");
                // Save values at measurements points
                Eigen::MatrixXd ensambleProbe(Nprobes, Nseeds);
                for (int i = 0; i < Nseeds; i++)
                {
                    ensambleProbe.col(i) = sampleField(Tensemble[i], probePosition);
                }

                probe_mean.col(timeI) = ensambleProbe.rowwise().mean();
                probe_MaxConfidence.col(timeI) = ITHACAmuq::muq2ithaca::quantile(ensambleProbe,
                                                 0.95);
                probe_minConfidence.col(timeI) = ITHACAmuq::muq2ithaca::quantile(ensambleProbe,

```

The method `forecastStep` advances each ensemble member from one time interval to the next.
At the first step it initializes each field with a sample drawn from the prior distribution:
```cpp
        Eigen::MatrixXd forecastStep(double startTime_, int startIndex_,
                                     double endTime_)
        {
            word outputFolder = "ITHACAoutput/forwardSamples";
            bool reconstructionFlag = 1;
            Eigen::MatrixXd observationMat;
            sampleFlag++;

            for (int i = 0; i < Nseeds; i++)
            {
                word fieldName = "Tsample" + std::to_string(i);
                volScalarField T(_T());

                if (startIndex_ == 0)
                {
                    restart();
                    T = _T();
                    Eigen::VectorXd internalField = priorSamples.col(i);
                    assignIF(T, internalField);
                }
```

Afterward, it starts from the current ensemble state.
```cpp
               else
                {
                    T = Tensemble.set(i, nullptr)();
                }
```

Each member is then evolved by calling `solve` in reconstruction mode, and the predicted observations are collected. The returned matrix is used in the Kalman correction step.
```cpp
resetRunTime(startTime_, startIndex_, endTime_);

                if (observationMat.cols() == 0)
                {
                    observationMat = solve(T, outputFolder, fieldName, reconstructionFlag);
                }
                else
                {
                    observationMat.conservativeResize(observationMat.rows(),
                                                      observationMat.cols() + 1);
                    Eigen::MatrixXd observationForCheck = solve(T, outputFolder, fieldName,
                                                          reconstructionFlag);
                    M_Assert(observationForCheck.cols() == 1,
                             "The forecast step passed trough more than a sampling step");
                    observationMat.col(i) = observationForCheck;
                }

                Tensemble.set(i, T.clone());
            }

            return observationMat;
        }
```

The helper method `assignIF` assigns an Eigen vector to the internal field of an OpenFOAM scalar field. This is used to transfer ensemble states between the linear-algebra representation and the CFD field representation.
```cpp
        void assignIF(volScalarField& field_, Eigen::VectorXd internalField_)
        {
            for (int i = 0; i < internalField_.size(); i++)
            {
                field_.ref()[i] = internalField_(i);
            }
        }
```

The method `sampleField` extracts the value of a scalar field at a specified probe location. It is mainly used for post-processing and visualization of the reconstructed and true temperature evolution.
```cpp
        Eigen::VectorXd sampleField(volScalarField field_, Foam::vector probeLocation_)
        {
            Foam::fvMesh& mesh = _mesh();
            Eigen::VectorXd output(1);
            output(0) = field_[mesh.findCell(probeLocation_)];
            return output;
        }
```

The method `restart` resets the OpenFOAM runtime and reloads the temperature field from disk so that a new solve can start from a clean state.
```cpp
        void restart()
        {
            Time& runTime = _runTime();
            instantList Times = runTime.times();
            runTime.setTime(Times[1], 0);
            _simple.clear();
            _T.clear();
            Foam::fvMesh& mesh = _mesh();
            _simple = autoPtr<simpleControl>
                      (
                          new simpleControl
                          (
                              mesh
                          )
                      );
            _T = autoPtr<volScalarField>
                 (
                     new volScalarField
                     (
                         IOobject
                         (
                             "T",
                             runTime.timeName(),
                             mesh,
                             IOobject::MUST_READ,
                             IOobject::AUTO_WRITE
                         ),
                         mesh
                     )
                 );
        }
```

The method `resetRunTime` similarly resets the current time window used by the forecast step and rebuilds the SIMPLE controller.
```cpp
        void resetRunTime(double startTime_, int startIndex_, double endTime_)
        {
            Time& runTime = _runTime();
            instantList Times = runTime.times();
            runTime.setTime(startTime_, startIndex_);
            runTime.setEndTime(endTime_);
            _simple.clear();
            Foam::fvMesh& mesh = _mesh();
            _simple = autoPtr<simpleControl>
                      (
                          new simpleControl
                          (
                              mesh
                          )
                      );
        }
```

The methods `priorSetup`, `priorSampling`, `modelErrorSetup`, and `measNoiseSetup` configure the probabilistic ingredients of the Ensemble Kalman Filter. They define the Gaussian prior on the initial state, generate the initial ensemble, specify the model-error distribution added during forecasting, and set the measurement-noise distribution used in the update step.
```cpp
        /// Setup of the prior distribution
        void priorSetup(double mean, double cov)
        {
            volScalarField& T(_T());
            int stateSize = T.size();
            Eigen::VectorXd prior_mu = Eigen::MatrixXd::Ones(stateSize, 1) * mean;
            Eigen::MatrixXd prior_cov = Eigen::MatrixXd::Identity(stateSize,
                                        stateSize) * cov;
            priorDensity = std::make_shared<muq::Modeling::Gaussian>(prior_mu, prior_cov);
        }

        /// Samples the prior density
        void priorSampling()
        {
            for (int i = 0; i < Nseeds; i++)
            {
                priorSamples.col(i) = priorDensity->Sample();
            }

            posteriorSamples = priorSamples;
        }

        /// Setup of the model error distribution
        void modelErrorSetup(double mean, double cov)
        {
            volScalarField& T(_T());
            int stateSize = T.size();
            Eigen::VectorXd modelError_mu = Eigen::MatrixXd::Ones(stateSize, 1) * mean;
            Eigen::MatrixXd modelError_cov = Eigen::MatrixXd::Identity(stateSize,
                                             stateSize) * cov;
            modelErrorDensity = std::make_shared<muq::Modeling::Gaussian>(modelError_mu,
                                modelError_cov);
        }

        /// Setup of the measurement noise distribution
        void measNoiseSetup(double mean, double cov)
        {
            meas_cov = cov;
            int obsSize = measurements.rows();
            M_Assert(obsSize > 0, "Read measurements before setting up the noise");
            Eigen::VectorXd measNoise_mu = Eigen::MatrixXd::Ones(obsSize, 1) * mean;
            Eigen::MatrixXd measNoise_cov = Eigen::MatrixXd::Identity(obsSize,
                                            obsSize) * cov;
            measNoiseDensity = std::make_shared<muq::Modeling::Gaussian>(measNoise_mu,
                               measNoise_cov);
        }
```

## The main

The `main` routine:

1. Instantiates the example object with a specified ensemble size.
```cpp
    solverPerformance::debug = 1; //No verbose output
    // Number of elements in the ensemble
    int Nseeds = 50;
    // Create the example object of the 02enKF_1DinverseHeatTransfer type
    EnKF_1DinverseHeatTransfer example(argc, argv, Nseeds);
```

2. Reads physical parameters (`k`, `rho`, `Cp`) from the dictionary.
```cpp
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
```

3. Solves the true forward model using the correct flux.
```cpp
   example.solveDirect();
```

4. Configures prior, model error and measurement noise Gaussian densities.
```cpp
   int Ntimes = example.Ntimes;
   example.priorSetup(initialField, 0.7);
   example.modelErrorSetup(0.0, 0.7);
   example.measNoiseSetup(0.0, 0.05);
```

5. Draws prior samples, and iteratively forecasts and updates the ensemble via `ITHACAmuq::muq2ithaca::EnsembleKalmanFilter`.
```cpp
    Eigen::MatrixXd priorSamples(stateSize, Nseeds);
    example.priorSampling();
```

6. Computes posterior means and confidence intervals, which are exported to `./ITHACAoutput/`.
```cpp
    Eigen::MatrixXd posteriorSamples(stateSize, Nseeds);
    Eigen::MatrixXd posteriorMean(stateSize, Ntimes);
    Eigen::MatrixXd minConfidence = posteriorMean;
    Eigen::MatrixXd maxConfidence = minConfidence;
    posteriorMean.col(0) = posteriorSamples.rowwise().mean();
    // Performs reconstruction
    example.reconstruct();
```

Extensive documentation is also embedded in the source file.

## Post‑processing
The Python script `plots.py` in the tutorial folder visualises the reconstructed temperature at a point, along with 95% confidence bands.

## The plain code
The plain code is available [here](https://raw.githubusercontent.com/ITHACA-FV/ITHACA-FV/refs/heads/master/tutorials/UQ/02enKF_1DinverseHeatTransfer/02enKF_1DinverseHeatTransfer.C).