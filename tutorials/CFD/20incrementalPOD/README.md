# Tutorial 20

## Introduction
This tutorial tests an incremental Proper Orthogonal Decomposition (POD)
algorithm, following Oxberry *et al.* [*Limited-memory adaptive snapshot
selection for POD*](https://onlinelibrary.wiley.com/doi/full/10.1002/nme.5283).
The model problem is a parameterized heat conduction equation on a square
domain subdivided into nine regions with independently varying diffusivity.

The governing equation is

$$\nabla\cdot(k\nabla T)=S,$$

leading to the discretized linear system $A T = S$. The source term $S$
is defined by a hat function and the parameter vector controls the diffusivity
in each of the nine subdomains.

Snapshots of the temperature field are collected for several parameter values.
POD bases are constructed both in the classical (batch) way and incrementally,
then a new solution is projected onto each basis to compare projection
errors. The algorithm thus tests the ability of the incremental POD to adapt
while keeping memory usage limited.

# A detailed look into the code
Let's understand the code of tutorial 20.

## Header files
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

## The `tutorialIPOD` class
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
```

As usual, an offline solver is defined: it either reads the existing snapshots or runs the simulation and saves them.
```cpp
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
```
Then, it defined the source forcing term:
```cpp
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
```

As in tutorial 02, the viscosities associated with the 9 boxes of the thermal block are initialized: 
```cpp
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

Then, we construct the `operator_list` where each term of the affine decomposition is stored:
```cpp
        void assemble_operator()
        {
            for (int i = 0; i < nu_list.size(); i++)
            {
                operator_list.append(fvm::laplacian(nu_list[i], T));
            }
        }
```

It performs a full order solution for a uniform value of `mu`:
```cpp
        volScalarField solveFull(double _mu)
        {
            word folder = "./ITHACAoutput/test/";
            scalar IF = 0;
            List<scalar> mu_now(9);
            volScalarField& T = _T();

            for (label j = 0; j < mu.cols() ; j++)
            {
                mu_now[j] = _mu;
                theta[j] = _mu;
            }

            assignIF(T, IF);
            truthSolve(mu_now, folder);
            return T;
        }
```

## Main execution
First, the example object is created, the number of modes are read from file, the parameters are sampled, and the source term is computed:
```cpp
    // Create the example object of the tutorialIPOD type
    tutorialIPOD example(argc, argv);
    // Read some parameters from file
    ITHACAparameters* para = ITHACAparameters::getInstance(example._mesh(),
                             example._runTime());
    int NmodesTout = para->ITHACAdict->lookupOrDefault<int>("NmodesTout", 15);
    int NmodesTproj = para->ITHACAdict->lookupOrDefault<int>("NmodesTproj", 10);
    double tolleranceSVD =
        para->ITHACAdict->lookupOrDefault<double>("tolleranceSVD", 1);
    // Set the number of parameters
    example.Pnumber = 9;
    example.Tnumber = NmodesTout;
    // Set the parameters
    example.setParameters();
    // Set the parameter ranges, in all the subdomains the diffusivity varies between
    // 0.001 and 0.1
    example.mu_range.col(0) = Eigen::MatrixXd::Ones(9, 1) * 0.001;
    example.mu_range.col(1) = Eigen::MatrixXd::Ones(9, 1) * 0.1;
    // Generate the Parameters
    example.genRandPar(example.Tnumber);
    // Set the size of the list of values that are multiplying the affine forms
    example.theta.resize(9);
    // Set the source term
    example.SetSource();
```
Then, the diffusivity is computed in each subdomain, the operators are assembled, and the offline stage is performed:
```cpp
    example.compute_nu();
    example.assemble_operator();
    example.offlineSolve();
```
Then, collect the POD modes:
```cpp
        ITHACAPOD::getModes(example.Tfield, example.Tmodes, example._T().name(),
        example.podex, 0, 0, NmodesTout);
```

The incremental POD is initialized:
```cpp
        scalarIncrementalPOD IPOD(example.Tfield[0], tolleranceSVD, "L2");
```
and filled:
```cpp
        for (int fieldI = 1; fieldI < example.Tfield.size(); fieldI++)
        {
            IPOD.addSnapshot(example.Tfield[fieldI]);
        }
        IPOD.writeModes();
```
This final part is for computing the full order solution (for compaarison), reconstructing the reduced solution, compute the projection of the snapshots into the POD space, and compute all the relative errors (stored in a `volScalarField`):
```cpp
    word folder = "./ITHACAoutput/testReconstruction";
    // Compute new full order solution
    volScalarField Tfull(example.solveFull(0.05));
    PtrList<volScalarField> TfullList;
    TfullList.append(Tfull.clone());
    PtrList<volScalarField> Tproj;
    // Project the full order solution onto the POD space
    example.Tmodes.projectSnapshots(TfullList, Tproj, NmodesTproj);
    ITHACAstream::exportSolution(Tfull, "1", folder, "Tfull");
    ITHACAstream::exportSolution(Tproj[0], "1", folder, "Tpod");
    // Compute the relative error between POD projected field and full order snapshot
    double EPS = 1e-16;
    volScalarField relativeErrorField(Tproj[0]);

    for (label i = 0; i < relativeErrorField.internalField().size(); i++)
    {
        if (std::abs(Tfull.ref()[i]) < EPS)
        {
            relativeErrorField.ref()[i] = (std::abs(Tfull.ref()[i] - Tproj[0].ref()[i])) /
                                          EPS;
        }
        else
        {
            relativeErrorField.ref()[i] = (std::abs(Tfull.ref()[i] - Tproj[0].ref()[i])) /
                                          Tfull.ref()[i];
        }
    }

    ITHACAstream::exportSolution(relativeErrorField,
                                 "1", folder,
                                 "relativeErrorField_POD");
    Info << "Relative error L2 norm POD = " << ITHACAutilities::L2Norm(
             relativeErrorField) << endl;
    // Project the full order solution onto the incremental POD space
    IPOD.projectSnapshots(TfullList, Tproj);
    ITHACAstream::exportSolution(Tproj[0], "1", folder, "Tipod");
    volScalarField Tipod = Tproj[0];

    // Compute the relative error between incremental POD projected field and full order snapshot
    for (label i = 0; i < relativeErrorField.internalField().size(); i++)
    {
        if (std::abs(Tfull.ref()[i]) < EPS)
        {
            relativeErrorField.ref()[i] = (std::abs(Tfull.ref()[i] - Tipod.ref()[i])) / EPS;
        }
        else
        {
            relativeErrorField.ref()[i] = (std::abs(Tfull.ref()[i] - Tipod.ref()[i]))
                                          Tfull.ref()[i];
        }
    }

    ITHACAstream::exportSolution(relativeErrorField,
                                 "1", folder,
                                 "relativeErrorField_IPOD");
    Info << "\n\nRelative error L2 norm incrementalPOD = " <<
         ITHACAutilities::L2Norm(relativeErrorField) << endl;
```

## The plain code
The plain code is available [here](https://raw.githubusercontent.com/ITHACA-FV/ITHACA-FV/refs/heads/master/tutorials/CFD/20incrementalPOD/20incrementalPOD.C).