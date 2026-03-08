# Tutorial 21

## Introduction
This tutorial demonstrates reduced order modeling for an unsteady turbulent
Navier–Stokes problem with Reynolds-averaged eddy viscosity treated via
radial basis function (RBF) interpolation. The flow is driven by a parameterized
inlet velocity and turbulence is accounted for through an eddy-viscosity
field.

## The necessary header files
The code relies on the unsteady turbulent solver and the standard ITHACA
utilities:
```cpp
#include "unsteadyNS.H"
#include "UnsteadyNSTurb.H"
#include "ITHACAPOD.H"
#include "ReducedUnsteadyNS.H"
#include "ReducedUnsteadyNSTurb.H"
#include "ITHACAstream.H"
```
Additional headers provide Eigen for linear algebra, utilities for mesh/I/O,
flux options and Foam2Eigen conversions.

## Implementation of the `tutorial21` class
The tutorial class inherits from `UnsteadyNSTurb` and stores references to the
velocity, pressure and eddy viscosity (`nut`) fields. The offline solve
method reads existing snapshots if available, otherwise loops over inlet
velocity parameter samples, updates boundary conditions and calls the
full-order solver:
```cpp
class tutorial21: public UnsteadyNSTurb
{
    public:
        explicit tutorial21(int argc, char* argv[])
            :
            UnsteadyNSTurb(argc, argv),
            U(_U()),
            p(_p()),
            nut(_nut())
        {}

        volVectorField& U;
        volScalarField& p;
        volScalarField& nut;

        void offlineSolve(std::string offlinepath)
        {
            Vector<double> inl(1, 0, 0);
            List<scalar> mu_now(1);
            label BCind = 0;

            if ((offline) && (ITHACAutilities::check_folder(offlinepath) == true))
            {
                ITHACAstream::read_fields(Ufield, U, offlinepath);
                ITHACAstream::read_fields(Pfield, p, offlinepath);
                ITHACAstream::read_fields(nutFields, nut, offlinepath);
            }
            else
            {
                for (label i = 0; i < mu.cols(); i++)
                {
                    inl[0] = mu(0, i);
                    mu_now[0] = mu(0, i);
                    assignBC(U, BCind, inl);
                    truthSolve(mu_now, offlinepath);
                }
            }
        }
};
```

## Main function
The main routine configures the parameter sampling, discretization and
reduced basis extraction:
```cpp
tutorial21 example(argc, argv);
ITHACAparameters* para = ITHACAparameters::getInstance(example._mesh(),
                     example._runTime());
int NmodesU = para->ITHACAdict->lookupOrDefault<int>("NmodesU", 10);
// ... other modes for pressure, supremizer, nut, projection, etc.

example.Pnumber = 1;
example.Tnumber = 2;
example.setParameters();
example.mu_range(0, 0) = 1.0;
example.mu_range(0, 1) = 1.1;
example.genEquiPar();
example.inletIndex.resize(1, 2);
example.inletIndex(0, 0) = 0;
example.inletIndex(0, 1) = 0;

example.startTime = 0;
example.finalTime = 5;
example.timeStep = 0.001;
example.writeEvery = 0.1;

example.offlineSolve("./ITHACAoutput/Offline/");
```

After collecting offline snapshots, the script computes the RBF coefficients
from the sample matrix, solves supremizer and lifting problems, normalizes the
lift field, and calculates POD modes for velocity, pressure and eddy viscosity.
Reduced matrices are then assembled and an online object (`ReducedUnsteadyNSTurb`)
is created. The eddy viscosity is interpolated via RBF during online solves,
and the reduced model is compared to high-fidelity solutions to compute error
metrics.

## Usage
Compile using the provided `Make/` file and execute the binary. Parameter
samples are read from `./ITHACAoutput/Offline/mu_samples_mat.txt` and the
RBF interpolation vector is stored in `example.velRBF`. Online runs export
reduction coefficients and reconstructions to `ITHACAoutput/`.

## Notes
This tutorial highlights the coupling of turbulent modelling with parameter
interpolation using RBFs in a reduced context. Adjust the number of modes or
sampling range via the ITHACAdict dictionary.

## The plain code
The plain code is available [here](https://raw.githubusercontent.com/ITHACA-FV/ITHACA-FV/refs/heads/master/tutorials/CFD/21unsteadyNSTurb_RBF/21unsteadyNSTurb_RBF.C).