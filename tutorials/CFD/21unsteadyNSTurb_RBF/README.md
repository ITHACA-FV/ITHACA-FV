# Tutorial 21

## Introduction
This tutorial demonstrates reduced order modeling for an unsteady turbulent
Navier–Stokes problem with Reynolds-averaged eddy viscosity treated via
radial basis function (RBF) interpolation. The flow is driven by a parameterized
inlet velocity and turbulence is accounted for through an eddy-viscosity
field.

# A detailed look into the code
Let's now have a look at the code of tutorial 21.

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
The tutorial class inherits from `UnsteadyNSTurb` and stores references to the velocity, pressure and eddy viscosity (`nut`) fields. The offline solve method reads existing snapshots if available, otherwise loops over inlet velocity parameter samples, updates boundary conditions and calls the full-order solver:
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
The main routine configures the parameter sampling, discretization and performs the offline stage:
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
Then, it defines the RBF interpolation data: 
```cpp
    auto mu_mat =
        ITHACAstream::readMatrix("./ITHACAoutput/Offline/mu_samples_mat.txt");
    example.velRBF = mu_mat.col(1);
```

It solves the supremizer problem and apply the lifting function to homogenize the velocity field:
```cpp
    example.solvesupremizer();
    example.liftSolve();
    ITHACAutilities::normalizeFields(example.liftfield);
    example.computeLift(example.Ufield, example.liftfield, example.Uomfield);
    ITHACAstream::exportFields(example.Uomfield, "./ITHACAoutput/Offline", "Uofield");
```
It then performs a POD decomposition for the velocity, the pressure and the eddy viscosity:
```cpp
    ITHACAPOD::getModes(example.Uomfield, example.Umodes, example._U().name(),
                        example.podex,
                        example.supex, 0, NmodesProject);
    ITHACAPOD::getModes(example.Pfield, example.Pmodes, example.p().name(),
                        example.podex,
                        example.supex, 0, NmodesProject);
    ITHACAPOD::getModes(example.nutFields, example.nutModes, example._nut().name(),
                        example.podex,
                        example.supex, 0, NmodesProject);
```
It solves the supremizer problem based on the pressure modes (not snapshots, as before), computes the reduced order matrices with supremizer stabilization, and create the object for the eddy viscosity RBF:
```cpp
    example.solvesupremizer("modes");
    example.projectSUP("./Matrices", NmodesU, NmodesP, NmodesSUP,
                       NmodesNUT, true);
    ReducedUnsteadyNSTurb pod_rbf(
        example);
```

Then, some settings for the reduced eddy viscosity and the penalty factor are done, and the matrix storing the RBF coefficients is innitialized:
```cpp
    pod_rbf.nu = 1e-05;
    pod_rbf.tauU.resize(1, 1);
    pod_rbf.tstart = 0;
    pod_rbf.finalTime = 10.1;
    pod_rbf.dt = 0.005;
    pod_rbf.storeEvery = 0.005;
    pod_rbf.exportEvery = 0.1;
    Eigen::MatrixXd rbfCoeff;
    rbfCoeff.resize(NmodesNUT + 1, 1);
```
Finally, it performs an online solve for the new values of inlet velocities:

```cpp
    for (label k = 0; k < 1; k++)
    {
        Eigen::MatrixXd velNow(1, 1);
        velNow(0, 0) = 1.05;
        pod_rbf.tauU(0, 0) = 0;
        pod_rbf.solveOnlineSUP(velNow);
        rbfCoeff.col(k) = pod_rbf.rbfCoeffMat.col(k);
        Eigen::MatrixXd tmp_sol(pod_rbf.y.rows() + 1, 1);
        tmp_sol(0) = k + 1;
        tmp_sol.col(0).tail(pod_rbf.y.rows()) = pod_rbf.y;
        pod_rbf.online_solution.append(tmp_sol);
    }
```
The matrices of the interpolated eddy viscosity coefficients and of the online solutions are saved:

```cpp
    ITHACAstream::exportMatrix(rbfCoeff, "rbfCoeff", "python",
                               "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(rbfCoeff, "rbfCoeff", "matlab",
                               "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(pod_rbf.online_solution, "red_coeff", "python",
                               "./ITHACAoutput/red_coeff");
    ITHACAstream::exportMatrix(pod_rbf.online_solution, "red_coeff", "matlab",
                               "./ITHACAoutput/red_coeff");
    ITHACAstream::exportMatrix(pod_rbf.online_solution, "red_coeff", "eigen",
                               "./ITHACAoutput/red_coeff");
```
The remaining part is for reconstructing the solution, computing the FOM counterparts for comparison, and computing the relative L2 errors. The FOM counterpart solutions are stored in the folder `Offline_check`, and the errors are stored in the `ErrorsL2` and `ErrorsFrob` directories.

```cpp
    // Reconstruct and export the solution
    pod_rbf.reconstruct(true, "./ITHACAoutput/Reconstruction/");
    // Solve the full order problem for the online velocity values for the purpose of comparison
    tutorial21 example2(argc, argv);
    /// Set the number of parameters
    example2.Pnumber = 1;
    /// Set samples
    example2.Tnumber = 2;
    /// Set the parameters infos
    example2.setParameters();
    // Set the parameter ranges
    example2.mu_range(0, 0) = 1.05;
    example2.mu_range(0, 1) = 1.05;
    // Generate equispaced samples inside the parameter range
    example2.genEquiPar();
    // Set the inlet boundaries where we have non homogeneous boundary conditions
    example2.inletIndex.resize(1, 2);
    example2.inletIndex(0, 0) = 0;
    example2.inletIndex(0, 1) = 0;
    // Time parameters
    example2.startTime = 0;
    example2.finalTime = 5;
    example2.timeStep = 0.001;
    example2.writeEvery = 0.1;
    example2.offlineSolve("./ITHACAoutput/Offline_check/");
    Eigen::MatrixXd errFrobU = ITHACAutilities::errorFrobRel(example2.Ufield,
                               pod_rbf.uRecFields);
    Eigen::MatrixXd errFrobP =  ITHACAutilities::errorFrobRel(example2.Pfield,
                                pod_rbf.pRecFields);
    Eigen::MatrixXd errFrobNut =  ITHACAutilities::errorFrobRel(example2.nutFields,
                                  pod_rbf.nutRecFields);
    ITHACAstream::exportMatrix(errFrobU, "errFrobU", "matlab",
                               "./ITHACAoutput/ErrorsFrob/");
    ITHACAstream::exportMatrix(errFrobP, "errFrobP", "matlab",
                               "./ITHACAoutput/ErrorsFrob/");
    ITHACAstream::exportMatrix(errFrobNut, "errFrobNut", "matlab",
                               "./ITHACAoutput/ErrorsFrob/");
    Eigen::MatrixXd errL2U = ITHACAutilities::errorL2Rel(example2.Ufield,
                             pod_rbf.uRecFields);
    Eigen::MatrixXd errL2P =  ITHACAutilities::errorL2Rel(example2.Pfield,
                              pod_rbf.pRecFields);
    Eigen::MatrixXd errL2Nut =  ITHACAutilities::errorL2Rel(example2.nutFields,
                                pod_rbf.nutRecFields);
    ITHACAstream::exportMatrix(errL2U, "errL2U", "matlab",
                               "./ITHACAoutput/ErrorsL2/");
    ITHACAstream::exportMatrix(errL2P, "errL2P", "matlab",
                               "./ITHACAoutput/ErrorsL2/");
    ITHACAstream::exportMatrix(errL2Nut, "errL2Nut", "matlab",
                               "./ITHACAoutput/ErrorsL2/");
```


## Usage
Before running the `21unsteadyNSTurb_RBF` application, the mesh should be built by running `create_mesh`.

## The plain code
The plain code is available [here](https://raw.githubusercontent.com/ITHACA-FV/ITHACA-FV/refs/heads/master/tutorials/CFD/21unsteadyNSTurb_RBF/21unsteadyNSTurb_RBF.C).