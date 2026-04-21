# Tutorial 11

## Introduction
The problem consists of an unsteady Buoyant Boussinesq (BB) problem for open flows where pressure cannot be neglected. The setup involves parameterized temperature and velocity boundary conditions in an open flow configuration, demonstrating ROM for buoyant flows with pressure effects.

# A detailed look into the code
This is the description of the code for tutorial 11.

## The necessary header files
First of all let's have a look into the header files which have to be included, indicating what they are responsible for:

`<UnsteadyBB.H>` is the base class for unsteady Buoyant Boussinesq problems.
`<ITHACAPOD.H>` is for the computation of the POD modes.
`<ReducedUnsteadyBB.H>` is for the reduced-order unsteady BB problem.
`<ITHACAstream.H>` is responsible for reading and exporting the fields and other sorts of data.

Additional standard libraries:
**Chrono** to compute execution times, **math.h** for mathematical functions, **iomanip** for output formatting.

# Implementation of the tutorial11 class
We define the tutorial11 class as a child of the `UnsteadyBB` class.
The constructor is defined with members that are the fields required to be manipulated during the resolution of the full order problem. Such fields are also initialized with the same initial conditions in the solver.
```cpp
class tutorial11: public UnsteadyBB
{
    public:
        explicit tutorial11(int argc, char* argv[])
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
Inside the tutorial11 class we define the offlineSolve method with parameterized boundary conditions. If the offline solve has been previously performed then the method just reads the existing snapshots including pressure fields. If not, it loops over the parameter samples and boundary condition parameters to perform the full-order simulations.
```cpp
        void offlineSolve(Eigen::MatrixXd par_BC)
        {
            List<scalar> mu_now(1);

            if (offline)
            {
                ITHACAstream::read_fields(Ufield, U, "./ITHACAoutput/Offline/");
                ITHACAstream::read_fields(Pfield, p, "./ITHACAoutput/Offline/");
                ITHACAstream::read_fields(Prghfield, p_rgh, "./ITHACAoutput/Offline/");
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
The other functions (`onlineSolveFull`, `onlineSolveRead`, `liftSolveT`) are introduced to solve the FOM problem for the novel online parameters, to read the online pre-saved solutions, and to compute the lifting function, respectively. For more details on these functions, refer to [the documentation of tutorial 10](https://github.com/ITHACA-FV/ITHACA-FV/tree/master/tutorials/CFD/10unsteadyBB_enclosed#tutorial-10). 

## Definition of the main function
The main function sets up the problem parameters, performs the offline phase including supremizer stabilization, computes POD modes for velocity, pressure, p_rgh, and temperature, evaluates projection errors, and solves the online reduced problem.

First, the tutorial object is constructed and parameters are configured. Note that the parameters are read from the files `par_offline_BC` (for the offline stage), and `par_online_BC` (for online parameters).
```cpp
tutorial11 example(argc, argv);
word par_offline_BC("./par_offline_BC");
Eigen::MatrixXd par_off_BC = ITHACAstream::readMatrix(par_offline_BC);
word par_online_BC("./par_online_BC");
Eigen::MatrixXd par_on_BC = ITHACAstream::readMatrix(par_online_BC);
ITHACAparameters* para = ITHACAparameters::getInstance(example._mesh(),
                     example._runTime());
word stabilization = para->ITHACAdict->lookupOrDefault<word>("Stabilization", "supremizer");
// ... other parameters
example.Pnumber = 1;
example.Tnumber = 1;
example.setParameters();
example.mu_range(0, 0) = 0.00001;
example.mu_range(0, 1) = 0.00001;
example.genEquiPar();
example.inletIndexT.resize(3, 1);
example.inletIndexT << 1, 2, 3;
example.inletIndex.resize(1, 2);
example.inletIndex << 3, 0;
example.startTime = 0.0;
example.finalTime = 5.0;
example.timeStep = 0.002;
example.writeEvery = 0.01;
```

The offline solve, lift functions, POD modes, and supremizer are computed:
```cpp
example.offlineSolve(par_off_BC);
example.liftSolve();
example.liftSolveT();
example.computeLift(example.Ufield, example.liftfield, example.Uomfield);
example.computeLiftT(example.Tfield, example.liftfieldT, example.Tomfield);
ITHACAPOD::getModes(example.Uomfield, example.Umodes, example._U().name(),
                    example.podex, 0, 0, NmodesOut, false);
// ... similar for P, Prgh, T
if (stabilization == "supremizer")
{
    example.solvesupremizer("modes");
}
```

Projection errors are then calculated, reduced matrices are obtained, and online solutions are performed with error evaluation against full-order solutions.
Everything repeats the same procedure of tutorial 10, for more details refer to [the documentation of tutorial 10](https://github.com/ITHACA-FV/ITHACA-FV/tree/master/tutorials/CFD/10unsteadyBB_enclosed#tutorial-10).

This completes the tutorial for open buoyant flows with pressure effects.

## The plain code
The plain code is available [here](https://raw.githubusercontent.com/ITHACA-FV/ITHACA-FV/refs/heads/master/tutorials/CFD/11UnsteadyBBOpen/11UnsteadyBBOpen.C).