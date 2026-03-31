# Tutorial 07

## Introduction
The problem consists of an unsteady Navier-Stokes problem coupled with heat transport in a mixing elbow geometry. The physical setup involves fluid flow with parameterized inlet velocities and temperatures, demonstrating reduced order modeling for non-isothermal flows.

The geometry is a mixing elbow where two inlets with different velocities and temperatures merge, and the outlet conditions are simulated. The tutorial showcases the application of POD-based reduction to both momentum and energy equations.

# A detailed look into the code
This section explains the main steps necessary to construct the tutorial N°7.

## The necessary header files
First of all let's have a look into the header files which have to be included, indicating what they are responsible for:
```cpp
#include "unsteadyNST.H"
#include "ITHACAPOD.H"
#include "ReducedUnsteadyNS.H"
#include "ReducedUnsteadyNST.H"
#include "ITHACAstream.H"
```
`<unsteadyNST.H>` is the base class for unsteady Navier-Stokes with temperature problems.
`<ITHACAPOD.H>` is for the computation of the POD modes.
`<ReducedUnsteadyNS.H>` and `<ReducedUnsteadyNST.H>` are for the reduced-order unsteady NS and NST problems.
`<ITHACAstream.H>` is responsible for reading and exporting the fields and other sorts of data.

Additional standard libraries:
```cpp
#include <chrono>
#include <math.h>
#include <iomanip>
```
**Chrono** is useful to compute execution times, **math.h** is used for mathematical functions, and **iomanip** for output formatting.

## Implementation of the tutorial07 class
We define the tutorial07 class as a child of the `unsteadyNST` class.
The constructor is defined with members that are the fields required to be manipulated during the resolution of the full order problem. Such fields are also initialized with the same initial conditions in the solver.
```cpp
class tutorial07: public unsteadyNST
{
    public:
        explicit tutorial07(int argc, char* argv[])
            :
            unsteadyNST(argc, argv),
            U(_U()),
            p(_p()),
            T(_T())
        {}

        // Fields To Perform
        volVectorField& U;
        volScalarField& p;
        volScalarField& T;
```
Inside the tutorial07 class we define the `offlineSolve` method. If the offline solve has been previously performed, then the method just reads the existing snapshots from the `Offline` directory. If not, it loops over all the parameter samples and performs the full-order simulations.
```cpp
        void offlineSolve()
        {
            Vector<double> inl(0, 0, 0);
            List<scalar> mu_now(1);
            Info << "here" << endl;

            if (offline)
            {
                ITHACAstream::read_fields(Ufield, U, "./ITHACAoutput/Offline/");
                ITHACAstream::read_fields(Pfield, p, "./ITHACAoutput/Offline/");
                ITHACAstream::read_fields(Tfield, T, "./ITHACAoutput/Offline/");
                mu_samples =
                    ITHACAstream::readMatrix("./ITHACAoutput/Offline/mu_samples_mat.txt");
            }
            else
            {
                for (label i = 0; i < mu.cols(); i++)
                {
                    mu_now[0] = mu(0, i);
                    truthSolve(mu_now);
                }
            }
        }
```

## Definition of the main function
The main function sets up the problem parameters, performs the offline phase, computes POD modes, and solves the online reduced problem.

First, the tutorial object is constructed and parameters are read from the ITHACA dictionary:
```cpp
tutorial07 example(argc, argv);
ITHACAparameters* para = ITHACAparameters::getInstance(example._mesh(),
                     example._runTime());
int NmodesUout = para->ITHACAdict->lookupOrDefault<int>("NmodesUout", 5);
// ... other parameters
```

Then, problem-specific settings are configured:
```cpp
example.Pnumber = 1;
example.Tnumber = 1;
example.setParameters();
example.mu_range(0, 0) = 0.1;
example.mu_range(0, 1) = 0.1;
example.genEquiPar();
example.inletIndex.resize(2, 2);
example.inletIndex << 3, 0, 2, 1;
example.inletIndexT.resize(3, 1);
example.inletIndexT << 3, 2, 0;
example.startTime = 0;
example.finalTime = 50;
example.timeStep = 0.05;
example.writeEvery = 0.1;
```

The offline solve is performed:
```cpp
example.offlineSolve();
```

Supremizer and lift functions are then computed.
Note that the lift function is assigned for both the velocity and the temperature (`liftSolveT` and `computeLiftT`).
```cpp
example.solvesupremizer();
example.liftSolve();
example.liftSolveT();
example.computeLift(example.Ufield, example.liftfield, example.Uomfield);
example.computeLiftT(example.Tfield, example.liftfieldT, example.Tomfield);
```

POD modes are extracted:
```cpp
ITHACAPOD::getModes(example.Uomfield, example.Umodes, example._U().name(),
                    example.podex, 0, 0, NmodesUout);
// ... similar for P, T, SUP
```

The reduced problem is set up and solved online:
```cpp
reducedUnsteadyNST reduced(example);
reduced.nu = 0.1;
reduced.tstart = 0;
reduced.finalTime = 50;
reduced.dt = 0.05;
reduced.DT = 1e-06;
Eigen::MatrixXd vel_now(2, 1);
vel_now << 0.6, 1.2;
Eigen::MatrixXd temp_now(3, 1);
temp_now << 60, 70, 60;
reduced.solveOnline_sup(vel_now, temp_now);
```

The reduced solution is then reconstructed and exported:
```cpp
reduced.reconstruct_sup("./ITHACAoutput/ReconstructionSUP/", 2);
reduced.reconstruct_supt("./ITHACAoutput/ReconstructionSUP/", 2);
```

This completes the tutorial, demonstrating the full workflow from full-order snapshots to reduced-order online solution.

The plain code is available [here](https://raw.githubusercontent.com/ITHACA-FV/ITHACA-FV/refs/heads/master/tutorials/CFD/07NonIsothermalMixingElbow/07NonIsothermalMixingElbow.C).