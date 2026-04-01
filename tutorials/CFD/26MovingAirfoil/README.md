# Tutorial 26
## Introduction
This tutorial covers a compressible, unsteady Reynolds-averaged Navier–Stokes (RANS) simulation of flow past a moving airfoil. The focus is on computing a reduced-order model and applying hyper-reduction to the resulting modes.

# A detailed look into the code
Let's have a look at the code of tutorial 26.

## Header files
Important headers used by this example include:
```cpp
#include "CompressibleUnSteadyRhoPimple.H"
#include "ITHACAPOD.H"
#include "ITHACAstream.H"
#include "Foam2Eigen.H"
#include "DEIM.H"
#include "HyperReducedCompressibleUnSteadyNS.H"
```
The code makes use of DEIM for hyper-reduction and includes `Eigen`, `chrono` and other utilities.

## Implementation
The `tutorial26` class derives from the compressible rhoPimple solver and holds references to the velocity `U`, pressure `p` and energy `E` fields. The `offlineSolve` method either reads precomputed snapshots or runs the full-order simulation to generate them:
```cpp
class tutorial26: public CompressibleUnSteadyRhoPimple
{
    public:
        tutorial26(int argc, char* argv[])
            : CompressibleUnSteadyRhoPimple(argc, argv),
              U(_U()),
              p(_p()),
              E(_E())
        {}

        volVectorField& U;
        volScalarField& p;
        volScalarField& E;

        void offlineSolve(word folder = "./ITHACAoutput/Offline/")
        {
            if (offline)
            {
                ITHACAstream::read_fields(Ufield, U, folder);
                ITHACAstream::read_fields(Pfield, p, folder);
                ITHACAstream::read_fields(Efield, E, folder);
            }
            else
            {
                truthSolve(folder);
            }
        }
};
```

## The main
The `main` function sets simulation parameters (time interval, time step, modes to output) from the ITHACAdict dictionary:
```cpp
    std::clock_t startOff;
    double durationOff;
    // Read some parameters from file
    ITHACAparameters* para = ITHACAparameters::getInstance(example.meshPtr(),
                             example._runTime());
    // Read the par file where the parameters are stored
    int NmodesUout  =  readInt(para->ITHACAdict->lookup("NmodesUout"));
    int NmodesPout  =  readInt(para->ITHACAdict->lookup("NmodesPout"));
    int NmodesEout  =  readInt(para->ITHACAdict->lookup("NmodesEout"));
    int NmodesUproj  = readInt(para->ITHACAdict->lookup("NmodesUproj"));
    int NmodesPproj  = readInt(para->ITHACAdict->lookup("NmodesPproj"));
    int NmodesEproj  = readInt(para->ITHACAdict->lookup("NmodesEproj"));
```
Then, it access the time parameters:
```cpp
    example.startTime  = 0;
    example.finalTime  = 0.15;
    example.timeStep   = 2e-06;
    example.writeEvery = 4e-04;
```
executes the offline stage:
```cpp
    example.offlineSolve();
```
and computes POD modes for velocity, pressure and energy (if not already saved):
```cpp
    if (example.podex == 0 )
    {
        ITHACAPOD::getModes(example.Ufield, example.Umodes, example._U().name(),
                            example.podex, 0, 0, NmodesUout);
        ITHACAPOD::getModes(example.Pfield, example.Pmodes, example._p().name(),
                            example.podex, 0, 0, NmodesPout);
        ITHACAPOD::getModes(example.Efield, example.Emodes, example.E().name(),
                            example.podex, 0, 0, NmodesEout);
    }
    else
    {
        ITHACAstream::read_fields(example.Umodes, example._U(), "./ITHACAoutput/POD/");
        ITHACAstream::read_fields(example.Pmodes, example._p(), "./ITHACAoutput/POD/");
        ITHACAstream::read_fields(example.Emodes, example._E(), "./ITHACAoutput/POD/");
    }
```
 It then constructs a `HyperReducedCompressibleUnSteadyNS` object and solves the hyper-reduced system using specified projection dimensions:
 ```cpp
    HyperReducedCompressibleUnSteadyNS hyperreduced(example);
    // Info << hyperreduced.Umodes.size() << endl;
    // Info << hyperreduced.Pmodes.size() << endl;
    // Info << hyperreduced.Emodes.size() << endl;
    hyperreduced.startTime = example.startTime;
    hyperreduced.finalTime = example.finalTime;
    hyperreduced.timeStep = example.timeStep;
    hyperreduced.writeEvery = example.writeEvery;
    /// Solving the hyper-reduced problem
    hyperreduced.SolveHyperReducedSys(NmodesUproj, NmodesPproj, NmodesEproj);
```

## Hyper-reduction
Hyper-reduction is performed on the compressible fields to accelerate online simulations. The class `HyperReducedCompressibleUnSteadyNS` handles both the
offline assembly and the online solve. Users may specify the numbers of projection modes `NmodesUproj`, `NmodesPproj` and `NmodesEproj` through the dictionary.


## The plain code
The plain code is available [here](https://raw.githubusercontent.com/ITHACA-FV/ITHACA-FV/refs/heads/master/tutorials/CFD/26MovingAirfoil/26MovingAirfoil.C).