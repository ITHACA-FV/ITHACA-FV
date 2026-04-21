# Tutorial 14

## Introduction
The problem consists of an unsteady Navier-Stokes problem reproduced using Dynamic Mode Decomposition (DMD). The setup involves parameterized viscosity, and DMD is applied to extract dynamic modes from the velocity and pressure snapshots for flow reconstruction.

# A detailed look into the code
Here we document the code used for tutorial 14.

## The necessary header files
First of all, let's have a look into the header files which have to be included, indicating what they are responsible for:
`<unsteadyNS.H>` is the base class for unsteady NS problems.
`<ITHACAPOD.H>` is for the computation of the POD modes.
`<ITHACADMD.H>` is for Dynamic Mode Decomposition.
`<ReducedUnsteadyNS.H>` is for the reduced-order unsteady NS problem.
`<ITHACAstream.H>` is responsible for reading and exporting the fields and other sorts of data.

Additional libraries:
**Chrono** to compute execution times, **math.h** for mathematical functions, **iomanip** for output formatting, **redsvd** for SVD computations in DMD.

## Implementation of the tutorial14 class
We define the tutorial14 class as a child of the `unsteadyNS` class.
The constructor is defined with members that are the fields required to be manipulated during the resolution of the full order problem. Such fields are also initialized with the same initial conditions in the solver.
```cpp
class tutorial14: public unsteadyNS
{
    public:
        explicit tutorial14(int argc, char* argv[])
            :
            unsteadyNS(argc, argv),
            U(_U()),
            p(_p())
        {}

        // Fields To Perform
        volVectorField& U;
        volScalarField& p;
```
Inside the tutorial14 class we define the offlineSolve method. If the offline solve has been previously performed then the method just reads the existing snapshots. If not, it loops over the viscosity parameters, changes the viscosity, and performs the full-order simulations.
```cpp
        void offlineSolve()
        {
            List<scalar> mu_now(1);

            if (offline)
            {
                ITHACAstream::read_fields(Ufield, U, "./ITHACAoutput/Offline/");
                ITHACAstream::read_fields(Pfield, p, "./ITHACAoutput/Offline/");
            }
            else
            {
                for (label i = 0; i < mu.cols(); i++)
                {
                    mu_now[0] = mu(0, i);
                    change_viscosity( mu(0, i));
                    truthSolve(mu_now);
                }
            }
        }
```

## Definition of the main function
The main function sets up the problem parameters, performs the offline phase to generate snapshots, and applies DMD to reconstruct the flow dynamics.

First, the tutorial object is constructed and DMD parameters are read:
```cpp
tutorial14 example(argc, argv);
ITHACAparameters* para = ITHACAparameters::getInstance(example._mesh(),
                     example._runTime());
double dtDMD = readScalar(para->ITHACAdict->lookup("dtDMD"));
double finalTimeDMD = readScalar(para->ITHACAdict->lookup("finalTimeDMD"));
double startTimeDMD = readScalar(para->ITHACAdict->lookup("startTimeDMD"));
int numberOfModesDMD = readInt(para->ITHACAdict->lookup("numberOfModesDMD"));
bool exactDMD = readBool(para->ITHACAdict->lookup("exactDMD"));
bool exportDMDModes = readBool(para->ITHACAdict->lookup("exportDMDModes"));
word exportFolder = string(para->ITHACAdict->lookup("exportFolder"));
example.Pnumber = 1;
example.Tnumber = 1;
example.setParameters();
example.mu_range(0, 0) = 0.01;
example.mu_range(0, 1) = 0.01;
example.genEquiPar();
example.startTime = 60;
example.finalTime = 85;
example.timeStep = 0.01;
example.writeEvery = 0.1;
```

The offline solve is performed:
```cpp
example.offlineSolve();
```

DMD is applied to velocity and pressure fields (`DMDv` and `DMDp`, respectively):
```cpp
ITHACADMDvolVector DMDv(example.Ufield, example.writeEvery);
DMDv.getModes(numberOfModesDMD, exactDMD, exportDMDModes);
DMDv.getDynamics(startTimeDMD, finalTimeDMD, dtDMD);
DMDv.reconstruct(exportFolder, exportFieldNameU);
ITHACADMDvolScalar DMDp(example.Pfield, example.writeEvery);
DMDp.getModes(numberOfModesDMD, exactDMD, exportDMDModes);
DMDp.getDynamics(startTimeDMD, finalTimeDMD, dtDMD);
DMDp.reconstruct(exportFolder, exportFieldNameP);
```

This completes the tutorial demonstrating DMD for unsteady flow analysis.

## The plain code
The plain code is available [here](https://raw.githubusercontent.com/ITHACA-FV/ITHACA-FV/refs/heads/master/tutorials/CFD/14DMDexample/14DMDexample.C).