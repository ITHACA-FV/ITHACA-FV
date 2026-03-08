# Tutorial 10

## Introduction
The problem consists of an unsteady Buoyant Boussinesq (BB) 2D problem with parameterized temperature boundary conditions. The physical setup represents a differentially heated cavity with uniform temperature on the left (hot) and right (cold) sides, while other sides are adiabatic. The cavity has an aspect ratio of 1.0, laminar flow, and air as the working fluid with Pr = 0.7. Ambient temperature is 300 K, hot wall at 301.5 K, cold wall at 298.5 K.

## The necessary header files
First of all let's have a look into the header files which have to be included, indicating what they are responsible for:
```cpp
#include "UnsteadyBB.H"
#include "ITHACAPOD.H"
#include "ReducedUnsteadyBB.H"
#include "ITHACAstream.H"
```
`<UnsteadyBB.H>` is the base class for unsteady Buoyant Boussinesq problems.
`<ITHACAPOD.H>` is for the computation of the POD modes.
`<ReducedUnsteadyBB.H>` is for the reduced-order unsteady BB problem.
`<ITHACAstream.H>` is responsible for reading and exporting the fields and other sorts of data.

Additional standard libraries:
```cpp
#include <chrono>
#include <math.h>
#include <iomanip>
```
**Chrono** to compute execution times, **math.h** for mathematical functions, **iomanip** for output formatting.

## Implementation of the tutorial10 class
We define the tutorial10 class as a child of the `UnsteadyBB` class.
The constructor is defined with members that are the fields required to be manipulated during the resolution of the full order problem. Such fields are also initialized with the same initial conditions in the solver.
```cpp
class tutorial10: public UnsteadyBB
{
    public:
        explicit tutorial10(int argc, char* argv[])
            :
            UnsteadyBB(argc, argv),
            U(_U()),
            p(_p()),
            p_rgh(_p_rgh()),
            T(_T())
        {}

        // Fields To Perform
        volVectorField& U;
        volScalarField& p;
        volScalarField& p_rgh;
        volScalarField& T;
```
Inside the tutorial10 class we define the offlineSolve method with parameterized boundary conditions. If the offline solve has been previously performed then the method just reads the existing snapshots. If not, it loops over the parameter samples and boundary condition parameters to perform the full-order simulations.
```cpp
        void offlineSolve(Eigen::MatrixXd par_BC)
        {
            List<scalar> mu_now(1);

            if (offline)
            {
                ITHACAstream::read_fields(Ufield, U, "./ITHACAoutput/Offline/");
                ITHACAstream::read_fields(Tfield, T, "./ITHACAoutput/Offline/");
            }
            else
            {
                for (label k = 0; k < par_BC.rows(); k++)
                {
                    for (label j = 0; j < par_BC.cols(); j++)
                    {
                        for (label i = 0; i < mu.cols(); i++)
                        {
                            mu_now[0] = mu(0, i);
                        }

                        assignBC(T, inletIndexT(j, 0), par_BC(k, j));
                    }

                    truthSolve(mu_now);
                }
            }
        }
```

## Definition of the main function
The main function sets up the problem parameters, performs the offline phase, computes POD modes, and evaluates projection errors and online solutions.

First, the tutorial object is constructed and boundary condition parameters are read:
```cpp
tutorial10 example(argc, argv);
word par_offline_BC("./par_offline_BC");
Eigen::MatrixXd par_off_BC = ITHACAstream::readMatrix(par_offline_BC);
word par_online_BC("./par_online_BC");
Eigen::MatrixXd par_on_BC = ITHACAstream::readMatrix(par_online_BC);
ITHACAparameters* para = ITHACAparameters::getInstance(example._mesh(),
                     example._runTime());
int NmodesUproj   = para->ITHACAdict->lookupOrDefault<int>("NmodesUproj", 5);
// ... other parameters
```

Then, problem-specific settings are configured:
```cpp
example.Pnumber = 1;
example.Tnumber = 1;
example.setParameters();
example.mu_range(0, 0) = 0.00001;
example.mu_range(0, 1) = 0.00001;
example.genEquiPar();
example.inletIndexT.resize(2, 1);
example.inletIndexT << 1, 2;
example.startTime = 0.0;
example.finalTime = 10.0;
example.timeStep = 0.005;
example.writeEvery = 0.01;
```

The offline solve is performed:
```cpp
example.offlineSolve(par_off_BC);
```

Lift functions and POD modes are computed:
```cpp
example.liftSolveT();
example.computeLiftT(example.Tfield, example.liftfieldT, example.Tomfield);
ITHACAPOD::getModes(example.Ufield, example.Umodes, example._U().name(),
                    example.podex, 0, 0, NmodesOut);
ITHACAPOD::getModes(example.Tomfield, example.Tmodes, example._T().name(),
                    example.podex, 0, 0, NmodesOut);
```

Projection errors are calculated for different numbers of modes, and the reduced problem is solved online for various boundary conditions. Full-order solutions are also computed for comparison, and errors between ROM and FOM are evaluated.

This completes the tutorial, demonstrating ROM for unsteady buoyant flows with parameterized boundary conditions.