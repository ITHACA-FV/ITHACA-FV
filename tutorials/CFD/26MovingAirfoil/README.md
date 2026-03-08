# Tutorial 26
## Introduction
This tutorial covers a compressible, unsteady Reynolds-averaged Navier–Stokes
(RANS) simulation of flow past a moving airfoil. The focus is on computing a
reduced-order model and applying hyper-reduction to the resulting modes.

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
The code makes use of DEIM for hyper-reduction and includes Eigen, chrono and
other utilities.

## Implementation
The `tutorial26` class derives from the compressible Rho-PIMPLE solver and
holds references to the velocity `U`, pressure `p` and energy `E` fields. The
`offlineSolve` method either reads precomputed snapshots or runs the full-order
simulation to generate them:
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

## Main routine
The `main` function sets simulation parameters (time interval, time step,
modes to output) from the ITHACAdict dictionary, executes the offline stage,
and computes POD modes for velocity, pressure and energy. It then constructs a
`HyperReducedCompressibleUnSteadyNS` object and solves the hyper-reduced
system using specified projection dimensions.

## Hyper-reduction
Hyper-reduction is performed on the compressible fields to accelerate online
simulations. The class `HyperReducedCompressibleUnSteadyNS` handles both the
offline assembly and the online solve. Users may specify the numbers of
projection modes `NmodesUproj`, `NmodesPproj` and `NmodesEproj` through the
dictionary.

## Usage
Compile with the provided `Make/` target. Configure the problem via
`system/` dictionary files and run the executable to generate snapshots and
reduced models in `ITHACAoutput/`. The tutorial is designed to illustrate the
integration of compressible solvers with hyper-reduction techniques in ITHACA-FV.
