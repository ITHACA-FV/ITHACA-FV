# Tutorial 10

# Introduction
The problem consists of an unsteady Buoyant Boussinesq (BB) 2D problem with parameterized temperature boundary conditions. The physical setup represents a differentially heated cavity with uniform temperature on the left (hot) and right (cold) sides, while other sides are adiabatic. The cavity has an aspect ratio of $1.0$, laminar flow, and air as the working fluid with $Pr = 0.7$. Ambient temperature is $300$ K, hot wall at $301.5$ K, cold wall at $298.5$ K.

# A detailed look into the code
Let's look at the code used for tutorial 10.

## The necessary header files
First of all let's have a look into the header files which have to be included, indicating what they are responsible for:
`<UnsteadyBB.H>` is the base class for unsteady Buoyant Boussinesq problems.
`<ITHACAPOD.H>` is for the computation of the POD modes.
`<ReducedUnsteadyBB.H>` is for the reduced-order unsteady BB problem.
`<ITHACAstream.H>` is responsible for reading and exporting the fields and other sorts of data.

Additional standard libraries:
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
Inside the tutorial10 class we define the offlineSolve method with parameterized boundary conditions. If the offline solve has been previously performed, then the method just reads the existing snapshots. If not, it loops over the parameter samples and boundary condition parameters to perform the full-order simulations.
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
Then, the same FOM solver is defined for online parameters (for comparison with the ROM):
```cpp
void onlineSolveFull(Eigen::MatrixXd par_BC, label para_set_BC, fileName folder)
{
    if (ITHACAutilities::check_folder(folder))
    {
    }
    else
    {
        mkDir(folder);
        ITHACAutilities::createSymLink(folder);
        label i = para_set_BC;

        for (label j = 0; j < par_BC.cols(); j++)
        {
            assignBC(T, inletIndexT(j, 0), par_BC(i, j));
        }

        truthSolve(folder);
    }
}
```
The `onlineSolveRead` just loads previously computed full-order solutions from disk:
```cpp
void onlineSolveRead(fileName folder)
{
    if (ITHACAutilities::check_folder(folder))
    {
        ITHACAstream::read_fields(Ufield_on, U, folder);
        ITHACAstream::read_fields(Tfield_on, T, folder);
    }
    else
    {
    }
}
```
The `liftSolveT` function is used to handle non-homogeneous boundary conditions for the temperature:

```cpp
void liftSolveT()
{
    for (label k = 0; k < inletIndexT.rows(); k++)
    {
        ...
```

It creates a copy of the temperature field:
```cpp
        volScalarField Tlift("Tlift" + name(k), T);
        ...
```
Then, it sets the boundary conditions and assigns the nonhomogeneous BCs at the inlet (T=1), while it assigns T=0 elsewhere:

```cpp
        scalar t1 = 1;
        scalar t0 = 0;
        if (j == BCind)
        {
            assignBC(Tlift, j, t1);
        }
        else if (...)
        {
            assignBC(Tlift, j, t0);
        }
```

Then, the BC is spread smoothly into the domain by solving a steady Laplace equation:
```cpp
        fvScalarMatrix TEqn
        (
            fvm::laplacian(alphaEff, Tlift)
        );
        TEqn.solve();
```

Finally, a correction for non-orthogonal cells is applied:
```cpp
        while (simple.correctNonOrthogonal())
        {
            alphat = turbulence->nut() / Prt;
            alphat.correctBoundaryConditions();
            volScalarField alphaEff("alphaEff", turbulence->nu() / Pr + alphat);
            fvScalarMatrix TEqn
            (
                fvm::laplacian(alphaEff, Tlift)
            );
            TEqn.solve();
            ...
        }
```

The lifting function is then saved:
```cpp
        Tlift.write();
        liftfieldT.append((Tlift).clone());
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

A list containing different numbers of modes is created and exported (it will be used to test projection accuracy):
```cpp
    Eigen::MatrixXd List_of_modes(NmodesOut - 5, 1);

    for (int i = 0; i < List_of_modes.rows(); i++)
    {
        List_of_modes(i, 0) = i + 1;
    }

    ITHACAstream::exportMatrix(List_of_modes, "List_of_modes", "eigen",
                               "./ITHACAoutput/l2error");
```

A temperature basis including the lifting functions is built, and the projection errors for all number of modes in the previous list are computed.
```cpp
    PtrList<volScalarField> TLmodes;

    for (label k = 0; k < example.liftfieldT.size(); k++)
    {
        TLmodes.append((example.liftfieldT[k]).clone());
    }

    for (label k = 0; k < List_of_modes.size(); k++)
    {
        TLmodes.append((example.Tmodes[k]).clone());
    }

    Eigen::MatrixXd L2errorProjMatrixU(example.Ufield.size(), List_of_modes.rows());
    Eigen::MatrixXd L2errorProjMatrixT(example.Tfield.size(), List_of_modes.rows());
```

The coefficients and $L^2$ errors and then computed and stored in a matrix for each number of modes considered.
```cpp
    for (int i = 0; i < List_of_modes.rows(); i++)
    {
        Eigen::MatrixXd coeffU = ITHACAutilities::getCoeffs(example.Ufield,
                                 example.Umodes,
                                 List_of_modes(i, 0) + example.liftfield.size() + NmodesSUPproj);
        Eigen::MatrixXd coeffT = ITHACAutilities::getCoeffs(example.Tfield, TLmodes,
                                 List_of_modes(i, 0) + example.liftfieldT.size());
        PtrList<volVectorField> rec_fieldU = ITHACAutilities::reconstructFromCoeff(
                example.Umodes, coeffU, List_of_modes(i, 0));
        PtrList<volScalarField> rec_fieldT = ITHACAutilities::reconstructFromCoeff(
                TLmodes, coeffT, List_of_modes(i, 0) + example.liftfieldT.size());
        Eigen::MatrixXd L2errorProjU = ITHACAutilities::errorL2Rel(example.Ufield,
                                       rec_fieldU);
        Eigen::MatrixXd L2errorProjT = ITHACAutilities::errorL2Rel(example.Tfield,
                                       rec_fieldT);
        L2errorProjMatrixU.col(i) = L2errorProjU;
        L2errorProjMatrixT.col(i) = L2errorProjT;
    }

    ITHACAstream::exportMatrix(L2errorProjMatrixU, "L2errorProjMatrixU", "eigen",
                               "./ITHACAoutput/l2error");
    ITHACAstream::exportMatrix(L2errorProjMatrixT, "L2errorProjMatrixT", "eigen",
                               "./ITHACAoutput/l2error");
```

Now, the ROM operators are built using the supremizer stabilization approach, and the modes are truncated:
```cpp
    example.projectSUP("./Matrices", NmodesUproj, NmodesPproj, NmodesTproj,
                       NmodesSUPproj);
    // Resize the modes for projection
    example.Tmodes.resize(NmodesTproj);
    example.Umodes.resize(NmodesUproj);
```

The ROM solver is now built and the online solve is performed, using the supremizer stabilization. Then the solution is reconstructed:
```cpp
    ReducedUnsteadyBB reduced(example);
    // Set values of the online solve
    reduced.nu = 0.00001;
    reduced.Pr = 0.71;
    reduced.tstart = 0.0;
    reduced.finalTime = 10;
    reduced.dt = 0.005;
    // No parametrization of velocity on boundary
    Eigen::MatrixXd vel_now_BC(0, 0);

    // Set the online temperature BC and solve reduced model
    for (label k = 0; k < (par_on_BC.rows()); k++)
    {
        Eigen::MatrixXd temp_now_BC(2, 1);
        temp_now_BC(0, 0) = par_on_BC(k, 0);
        temp_now_BC(1, 0) = par_on_BC(k, 1);
        reduced.solveOnline_sup(temp_now_BC, vel_now_BC, k, par_on_BC.rows());
        reduced.reconstruct_sup("./ITHACAoutput/ReconstructionSUP", 2);
    }
```

The FOM is performed for the online parameters, the saved reconstructed solution are read, and a comparison in terms of $L^2$ errors is performed:
```cpp
    // Performing full order simulation for second parameter set - temp_BC
    tutorial10 HFonline2(argc, argv);
    HFonline2.Pnumber = 1;
    HFonline2.Tnumber = 1;
    HFonline2.setParameters();
    HFonline2.mu_range(0, 0) = 0.00001;
    HFonline2.mu_range(0, 1) = 0.00001;
    HFonline2.genEquiPar();
    HFonline2.inletIndexT.resize(2, 1);
    HFonline2.inletIndexT << 1, 2;
    HFonline2.startTime = 0.0;
    HFonline2.finalTime = 10;
    HFonline2.timeStep = 0.005;
    HFonline2.writeEvery = 0.01;
    // Reconstruct the online solution
    HFonline2.onlineSolveFull(par_on_BC, 1,
                              "./ITHACAoutput/HFonline2");
    // Performing full order simulation for third parameter set - temp_BC
    tutorial10 HFonline3(argc, argv);
    HFonline3.Pnumber = 1;
    HFonline3.Tnumber = 1;
    HFonline3.setParameters();
    HFonline3.mu_range(0, 0) = 0.00001;
    HFonline3.mu_range(0, 1) = 0.00001;
    HFonline3.genEquiPar();
    HFonline3.inletIndexT.resize(2, 1);
    HFonline3.inletIndexT << 1, 2;
    HFonline3.startTime = 0.0;
    HFonline3.finalTime = 10.0;
    HFonline3.timeStep = 0.005;
    HFonline3.writeEvery = 0.01;
    // Reconstruct the online solution
    HFonline3.onlineSolveFull(par_on_BC, 2,
                              "./ITHACAoutput/HFonline3");
    // Reading in the high-fidelity solutions for the parameter set
    // for which the offline solve has been performed
    example.onlineSolveRead("./ITHACAoutput/Offline/");
    // Reading in the high-fidelity solutions for the second parameter set
    example.onlineSolveRead("./ITHACAoutput/HFonline2/");
    // Reading in the high-fidelity solutions for the second parameter set
    example.onlineSolveRead("./ITHACAoutput/HFonline3/");
    // Calculate error between online- and corresponding full order solution
    Eigen::MatrixXd L2errorMatrixU = ITHACAutilities::errorL2Rel(
                                         example.Ufield_on, reduced.UREC);
    Eigen::MatrixXd L2errorMatrixT = ITHACAutilities::errorL2Rel(
                                         example.Tfield_on, reduced.TREC);
    //Export the matrix containing the error
    ITHACAstream::exportMatrix(L2errorMatrixU, "L2errorMatrixU", "eigen",
                               "./ITHACAoutput/l2error");
    ITHACAstream::exportMatrix(L2errorMatrixT, "L2errorMatrixT", "eigen",
                               "./ITHACAoutput/l2error");
```

This completes the tutorial, demonstrating ROM for unsteady buoyant flows with parameterized boundary conditions.

## The plain code
The plain code is available [here](https://raw.githubusercontent.com/ITHACA-FV/ITHACA-FV/refs/heads/master/tutorials/CFD/10unsteadyBB_enclosed/10unsteadyBB_enclosed.C).