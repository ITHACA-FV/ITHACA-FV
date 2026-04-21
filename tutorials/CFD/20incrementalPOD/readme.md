## Tutorial 20
### Introduction
This tutorial tests an incremental Proper Orthogonal Decomposition (POD)
algorithm, following Oxberry *et al.* [*Limited-memory adaptive snapshot
selection for POD*](https://onlinelibrary.wiley.com/doi/full/10.1002/nme.5283).
The model problem is a parameterized heat conduction equation on a square
domain subdivided into nine regions with independently varying diffusivity.

The governing equation is
$$\nabla\cdot(k\nabla T)=S$$
leading to the discretized linear system $A T = S$. The source term $S$
is defined by a hat function and the parameter vector controls the diffusivity
in each of the nine subdomains.

Snapshots of the temperature field are collected for several parameter values.
POD bases are constructed both in the classical (batch) way and incrementally,
then a new solution is projected onto each basis to compare projection
errors. The algorithm thus tests the ability of the incremental POD to adapt
while keeping memory usage limited.

### Header files
Required includes are:
```cpp
// Standard C++ I/O
#include <iostream>

// OpenFOAM utilities
#include "fvCFD.H"
#include "IOmanip.H"
#include "Time.H"

// ITHACA-FV
#include "ITHACAPOD.H"
#include "ITHACAutilities.H"

// Eigen for linear algebra
#include <Eigen/Dense>

#define _USE_MATH_DEFINES
#include <cmath>  // math constants
```

### `tutorialIPOD` class
The tutorial class inherits from `laplacianProblem` and adds fields for
temperature `T`, diffusivity `nu` and source `S`:
```cpp
class tutorialIPOD: public laplacianProblem
{
    public:
        explicit tutorialIPOD(int argc, char* argv[])
            : laplacianProblem(argc, argv),
              T(_T()),
              nu(_nu()),
              S(_S())
        {}

        volScalarField& T;
        volScalarField& nu;
        volScalarField& S;

        void offlineSolve(word folder = "./ITHACAoutput/Offline/")
        {
            if (offline)
            {
                ITHACAstream::read_fields(Tfield, "T", folder);
                mu_samples =
                    ITHACAstream::readMatrix(folder + "/mu_samples_mat.txt");
            }
            else
            {
                for (label i = 0; i < mu.rows(); i++)
                {
                    for (label j = 0; j < mu.cols() ; j++)
                    {
                        mu_now[j] = mu(i, j);
                        theta[j] = mu(i, j);
                    }
                    assignIF(T, IF);
                    truthSolve(mu_now, folder);
                }
            }
        }

        void SetSource()
        {
            volScalarField yPos = T.mesh().C().component(vector::Y).ref();
            volScalarField xPos = T.mesh().C().component(vector::X).ref();
            forAll(S, counter)
            {
                S[counter] = Foam::sin(xPos[counter] / 0.9 * M_PI) +
                             Foam::sin(yPos[counter] / 0.9 * M_PI);
            }
        }

        void compute_nu()
        {
            nu_list.resize(9);
            volScalarField nu1(nu);
            volScalarField nu2(nu);
            volScalarField nu3(nu);
            volScalarField nu4(nu);
            volScalarField nu5(nu);
            volScalarField nu6(nu);
            volScalarField nu7(nu);
            volScalarField nu8(nu);
            volScalarField nu9(nu);

            Eigen::MatrixXd Box1(2, 3);
            Box1 << 0, 0, 0, 0.3, 0.3, 0.1;
            Eigen::MatrixXd Box2(2, 3);
            Box2 << 0.3, 0, 0, 0.6, 0.3, 0.1;
            Eigen::MatrixXd Box3(2, 3);
            Box3 << 0.6, 0, 0, 0.91, 0.3, 0.1;
            Eigen::MatrixXd Box4(2, 3);
            Box4 << 0, 0.3, 0, 0.3, 0.6, 0.1;
            Eigen::MatrixXd Box5(2, 3);
            Box5 << 0.3, 0.3, 0, 0.6, 0.6, 0.1;
            Eigen::MatrixXd Box6(2, 3);
            Box6 << 0.6, 0.3, 0, 0.91, 0.6, 0.1;
            Eigen::MatrixXd Box7(2, 3);
            Box7 << 0, 0.6, 0, 0.3, 0.91, 0.1;
            Eigen::MatrixXd Box8(2, 3);
            Box8 << 0.3, 0.61, 0, 0.6, 0.91, 0.1;
            Eigen::MatrixXd Box9(2, 3);

            ITHACAutilities::setBoxToValue(nu1, Box1, 1.0);
            ITHACAutilities::setBoxToValue(nu2, Box2, 1.0);
            ITHACAutilities::setBoxToValue(nu3, Box3, 1.0);
            ITHACAutilities::setBoxToValue(nu4, Box4, 1.0);
            ITHACAutilities::setBoxToValue(nu5, Box5, 1.0);
            ITHACAutilities::setBoxToValue(nu6, Box6, 1.0);
            ITHACAutilities::setBoxToValue(nu7, Box7, 1.0);
            ITHACAutilities::setBoxToValue(nu8, Box8, 1.0);
            ITHACAutilities::setBoxToValue(nu9, Box9, 1.0);
        }
    };

```

### Main execution
The remainder of the original README contains detailed code snippets for
setting parameters, running offline solves, computing sources, building POD
bases, and performing projection tests. 
First, collect the POD modes:
```cpp
        ITHACAPOD::getModes(example.Tfield, example.Tmodes, example._T().name(),
```

Then, the incremental POD is initialized
```cpp
        scalarIncrementalPOD IPOD(example.Tfield[0], tolleranceSVD, "L2");
```
and filled
```cpp
        for (int fieldI = 1; fieldI < example.Tfield.size(); fieldI++)
        {
            IPOD.addSnapshot(example.Tfield[fieldI]);
        }
```

### Usage notes
Compile the tutorial using its `Make/` target. Provide parameter files and use
the ITHACAdict dictionary to control algorithmic options (number of modes,
incremental tolerance, etc.). The primary outputs are snapshot matrices and
error metrics comparing classical versus incremental POD projection.

## The plain code
The plain code is available [here](https://raw.githubusercontent.com/ITHACA-FV/ITHACA-FV/refs/heads/master/tutorials/CFD/20incrementalPOD/20incrementalPOD.C).