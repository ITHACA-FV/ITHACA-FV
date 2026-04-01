# Tutorial 16
## Introduction
This tutorial performs a sensitivity analysis using the full-order model (FOM) of the molten salt reactor cavity. The physical setup is identical to `15MSR_cavity`, but the focus here is on exploring parameter uncertainties rather than constructing a reduced order model.

### Parameters
Three random parameters are considered:

- kinematic viscosity `nu`,
- total delayed neutron fraction `betaTot` (ratios `beta_i/betaTot` fixed),
- third decay heat group constant `decLam3`.

Each parameter is assumed normally distributed and sampled at 1000 points.
The offline FOM is run for each sample, and statistics are collected.

### Figures of merit
The quantities monitored for each sample are:

- total power `P_tot`,
- mean temperature `T_m`,
both evaluated at t = 100 s. Only the final time and the initial condition
are exported since temporal evolution is not required for the analysis.

## Problem overview
The full order problem is solved by the class `msr` (see
`15MSR_cavity`), with identical geometry and boundary conditions. The
analysis uses the `LRSensitivity` utilities to compute sensitivities of the
figures of merit with respect to the parameters.

# A detailed look into the code
Let's have a look at the code of tutorial 16.

## The necessary header files
```cpp
#include "usmsrProblem.H"
#include "ReducedUnsteadyMSR.H"
#include "ITHACAstream.H"
#include "LRSensitivity.H"
#include "Tm.H"
#include "Ptot.H"
#include <Eigen/Dense>
#include <math.h>
#include <iomanip>
#include <chrono>
#include <ctime>
#include "Tm_time.H"
#include "Ptot_time.H"
```
This includes:

- the full-order MSR problem class (`usmsrProblem.H`);
- the reduced MSR class header, although it is not used in this code;
- ITHACA-FV utilities for file I/O, matrix export/import, sensitivity analysis;
- figure-of-merit classes: `Tm` for temperature-based outputs,
`Ptot` for power-based outputs;
- Eigen and standard C++ libraries for linear algebra, timing, and utilities.

## The FOM class
It defines a derived class `msr` that extends the base full-order problem `usmsrProblem`:
```cpp
class msr : public usmsrProblem
```
Then, all the flow variables, the neutron or reactor-related fields, and the thermal/material fields are initialized:
```cpp
volVectorField& U;
volScalarField& p;
volScalarField& flux;
volScalarField& prec1;
...
volScalarField& TXS;
```

Then, the number of modes used is read from the `ITHACAdict`:
```cpp
ITHACAparameters* para = ITHACAparameters::getInstance(_mesh(), _runTime());
int NUproj = para->ITHACAdict->lookupOrDefault<int>("NUproj", 2);
...
int NCproj = para->ITHACAdict->lookupOrDefault<int>("NCproj", 2);
```

Then, it defines the fraction of each delayed-neutron precursor group relative to the total delayed fraction.
```cpp
double r1 = _beta1().value() / _betaTot().value();
...
double r8 = _beta8().value() / _betaTot().value();
```
Later, when `betaTot` is changed as a sampled parameter, the code rescales each individual `beta_i` proportionally:

```cpp
_beta1().value() = r1 * mu(1, i);
...
```
So the total delayed fraction changes, but the relative distribution among the 8 groups stays fixed. That is physically consistent if you want to vary only the total amount and not the shape of the precursor spectrum.

As in all the previous tutorials, the `offlineSolve` is also defined:
```cpp
void offlineSolve(std::string dir)
```
If the database does not exist, it generates a database of full-order outputs for many sampled parameter values, updating and saving the parameters for each value:
```cpp
List<scalar> mu_now(3);
std::string folder = dir;
...
for (int i = 0; i < mu.cols(); i++)
{
    // create folders for each simulation
    folder.append(std::to_string(i));
    folder.append("/");
    // update parameters
    _nu().value() = mu(0, i);
    change_viscosity(mu(0, i));
    _betaTot().value() = mu(1, i);
    _beta1().value() = r1 * mu(1, i);
    ...
    _decLam3().value() = mu(2, i);
    // store parameters
    mu_now[0] = mu(0, i);
    mu_now[1] = mu(1, i);
    mu_now[2] = mu(2, i);
    // run full-order model
    truthSolve(mu_now, folder);
    // reset folder and restart case
    folder = dir;
    restart();
}
```
If the database exists, it just reads it:
```cpp
readMSRfields(dir);
```

Then, just as in the previous tutorial, all the figures of merit are defined in classes, and in particular the average temperature and the power density:
```cpp
class Tmlocal : public Tm {...}
class Plocal: public Ptot {...}
```

## The main function
At the beginning, the main creates the problem object:
```cpp
msr prova(argc, argv);
```
It defines which patch index is treated as the inlet for velocity-related operations:
```cpp
prova.inletIndex.resize(1, 2);
prova.inletIndex(0, 0) = 0;
prova.inletIndex(0, 1) = 0;
```
The, it sets the number of parameters and of samples:
```cpp
prova.Pnumber = 3;
prova.Tnumber = 1000;
prova.setParameters();
```
Then, it defines the distribution centers and standard deviations:
```cpp
double nu0 = 2.46E-06;
double betatot0 = 321.8E-05;
double dlam30 = 3.58e-04;
double signu = 0.1 / 3 * nu0;
double sigbeta = 0.1 / 3 * betatot0;
double sigdlam3 = 0.1 / 3 * dlam30;
```

We then read the admissible parameter ranges from file and specify that the sampling distribution is normal"
```cpp
prova.mu_range = ITHACAstream::readMatrix("./range_mat.txt");
std::string dist = {"normal"};
```

It then initializes the roes of `mu` as 1000 Monte Carlo samples for each parameter:
```cpp
prova.mu.row(0) = ITHACAsampling::samplingMC(dist, prova.mu_range(0, 0),
                  prova.mu_range(0, 1), nu0, signu, prova.Tnumber);
...
```

It saves or reload the parameters:
```cpp
double tstart = std::time(0);
Eigen::MatrixXd mu_f = prova.mu;
```

It exports the matrices if newly generated:
```cpp
if (prova.offline == false)
{
    ITHACAstream::exportMatrix(prova.mu, "mu_off", "eigen", "./ITHACAoutput");
}
```
or read existing set:
```cpp
else
{
    mu_f = ITHACAstream::readMatrix("./ITHACAoutput/mu_off_mat.txt");
}
```

Then, it runs the offline solve for all 1000 sampled parameters (this takes time):
```cpp
std::string folder = {"./ITHACAoutput/FOMoutput/"};
prova.offlineSolve(folder);
```

It measures and save the total FOM runtime:
```cpp
double tend = std::time(0);
double deltat_offline;
deltat_offline = difftime(tend, tstart);
if (prova.offline == false)
{
    std::ofstream diff_t("./ITHACAoutput/t_SA=" + name(deltat_offline));
}
```

Then, it creates the linear-regression sensitivity analysis object using the 3 parameters and 1000 samples. It computes the statistics of inputs and outputs, the regression coefficients, and a quality indicator for the linear approximation:
```cpp
int Nparameters = prova.Pnumber;
int samplingPoints = prova.Tnumber;
LRSensitivity analisi(Nparameters, samplingPoints);
```

It loads the sampled inputs into the sensitivity object `analisi` and computes mean and variance:
```cpp
analisi.MatX.col(0) = mu_f.row(0);
analisi.MatX.col(1) = mu_f.row(1);
analisi.MatX.col(2) = mu_f.row(2);
analisi.getXstats();
```
Then, the figures of merit are created as objects:
```cpp
Tmlocal TmediaSA(argc, argv, samplingPoints);
Plocal PtotSA(argc, argv, samplingPoints);
```
Then, it reads the full-order simulation results from disk and constructs the output vectors:

- one value of mean temperature per sample;
- one value of total power per sample.
```cpp
TmediaSA.buildMO(folder);
PtotSA.buildMO(folder);
```

A sensitivity analysis is performed also for the temperature:
```cpp
analisi.M = autoPtr<FofM>(& TmediaSA);
analisi.load_output(); //loads output vector into the regression machinery
```

Then, it computes output statistics (mean and variance), the regression coefficients, the linear-model quality, and saves temperature sensitivity analysis:
```cpp
analisi.getYstat();
Info << "Mean value and variance of the output:" << endl;
Info << analisi.Ey << "\t" << analisi.Vy << endl;
analisi.getBetas();
Info << "analisi regression coefficients:" << endl;
Info << analisi.betas << endl;
analisi.assessQuality();
Info << "quality: " << analisi.QI << endl;
Eigen::MatrixXd saveyout(samplingPoints, 1);
Eigen::MatrixXd saveymodel(samplingPoints, 1);
saveyout.col(0) = analisi.y;
saveymodel.col(0) = analisi.ylin;
```
It exports the real and predicted temperature:
```cpp
ITHACAstream::exportMatrix(saveyout, "T_out", "eigen", "./ITHACAoutput");
ITHACAstream::exportMatrix(saveymodel, "T_model", "eigen", "./ITHACAoutput");
```
It saves compact summary of temperature sensitivity:
```cpp
Eigen::MatrixXd coeffsT(4, 3);
coeffsT.setZero();
coeffsT(0, 0) = analisi.Ey;
coeffsT(0, 1) = analisi.Vy;
coeffsT(0, 2) = analisi.QI;
coeffsT.row(1) = analisi.EX;
coeffsT.row(2) = analisi.VX;
coeffsT.row(3) = analisi.betas;
ITHACAstream::exportMatrix(coeffsT, "coeffsT", "eigen", "./ITHACAoutput");
```
The same analysis is done for the total power:
```cpp
analisi.M = autoPtr<FofM>(& PtotSA);
analisi.load_output();
analisi.getYstat();
Info << "Mean value and variance of the output:" << endl;
Info << analisi.Ey << "\t" << analisi.Vy << endl;
analisi.getBetas();
Info << "analisi regression coefficients:" << endl;
Info << analisi.betas << endl;
analisi.assessQuality();
Info << "quality: " << analisi.QI << endl;
Eigen::MatrixXd saveyoutP(samplingPoints, 1);
Eigen::MatrixXd saveymodelP(samplingPoints, 1);
saveyoutP.col(0) = analisi.y;
saveymodelP.col(0) = analisi.ylin;
```
and the actual and predicted power outputs are also saved:
```cpp
ITHACAstream::exportMatrix(saveyoutP, "P_out", "eigen", "./ITHACAoutput");
ITHACAstream::exportMatrix(saveymodelP, "P_model", "eigen", "./ITHACAoutput");
```
Again, also the power sensitivity is saved in a compact summary form:
```cpp
Eigen::MatrixXd coeffsP(4, 3);
coeffsP.setZero();
coeffsP(0, 0) = analisi.Ey;
coeffsP(0, 1) = analisi.Vy;
coeffsP(0, 2) = analisi.QI;
coeffsP.row(1) = analisi.EX;
coeffsP.row(2) = analisi.VX;
coeffsP.row(3) = analisi.betas;
ITHACAstream::exportMatrix(coeffsP, "coeffsP", "eigen", "./ITHACAoutput");
```

## The plain code
The plain code is available [here](https://raw.githubusercontent.com/ITHACA-FV/ITHACA-FV/refs/heads/master/tutorials/CFD/16MSR_FOMSA/16MSR_FOMSA.C).