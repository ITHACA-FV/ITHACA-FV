# Tutorial 25
## Introduction
This tutorial shows a fluid–structure interaction (FSI) problem where a
cylinder moves in a viscous flow. The ITHACA-FV implementation leverages the
`fsiBasic` base class and constructs a reduced order model capturing both the
fluid velocity and the solid point displacement.

## Header files
Key includes for this tutorial are:
```cpp
#include "fsiBasic.H"
#include "ITHACAPOD.H"
#include "ITHACAstream.H"
#include "dynamicFvMesh.H"
#include "ReducedProblem.H"
#include "ReducedFsi.H"
```
Additional headers support mesh motion, point constraints and math constants.

## Implementation of `tutorial22` class
The `tutorial22` class inherits from `fsiBasic` and holds references to the
velocity `U`, pressure `p` and point displacement `pd` fields. The `offlineSolve`
method either loads existing snapshots or iterates over parameter samples,
updating a damping parameter and running the coupled solver. After each
simulation the code exports force and centre-of-mass data for post–processing.

```cpp
class tutorial22: public fsiBasic
{
    public:
        tutorial22(int argc, char* argv[])
            : fsiBasic(argc, argv),
              U(_U()),
              p(_p()),
              pd(_pointDisplacement())
        {}

        volVectorField& U;
        volScalarField& p;
        pointVectorField& pd;

        void offlineSolve(word folder = "./ITHACAoutput/Offline/")
        {
            if (offline )
            {
                ITHACAstream::readMiddleFields(Ufield, U, folder);
                ITHACAstream::readMiddleFields(Pfield, p, folder);
                ITHACAstream::readMiddleFields(Dfield, pd, folder);
            }
            else
            {
                for (label i = 0; i < mu.rows(); i++)
                {
                    word param_name = "damping";
                    updateStiffnessAndRebuildSolver(mu(i, 0), param_name);
                    startTime  = 0;
                    finalTime  = 30;
                    timeStep   = 0.003;
                    writeEvery = 1e-01;
                    truthSolve(i, folder);
                    word localFolder = folder + "../" + "/DataFromFoam_" + name(i + 1);
                    prepareFoamData(localFolder);
                    restart();
                    // clear time history arrays
                    fomforcex.clear();
                    fomforcey.clear();
                    centerofmassy.clear();
                }
            }
        }
};
```

## Main workflow
The `main` routine reads offline parameter files (`parsOff_mat.txt`), sets
the number of modes for velocity, pressure and displacement, and optionally
performs offline runs. POD modes are computed unless they already exist, and
additional testing data for online evaluation are generated using a second
instance of the tutorial object.

Online reduction is performed with the `ReducedFsi` class, solving the reduced
PIMPLE system for each test parameter set and exporting coefficients and
reconstructions. Error metrics comparing reduced and full-order solutions
(energy norms, centre-of-mass displacement, forces) are computed and stored
for post-processing.

## Usage
Compile via `Make/` and provide `parsOff_mat.txt`/`parsOn_mat.txt` files
containing damping parameter values. Run the executable to populate the
`ITHACAoutput` directory with offline snapshots, POD modes, and online
results. Python scripts `plotCentreMass.py`, `plotDrag.py` and `plotLift.py`
can be used to visualize the outputs.

This tutorial demonstrates building a coupled fluid-structure reduced model
with parameterized structural damping and moving mesh handling.

## The plain code
The plain code is available [here](https://raw.githubusercontent.com/ITHACA-FV/ITHACA-FV/refs/heads/master/tutorials/CFD/25Moving_Cylinder/25Moving_Cylinder.C).