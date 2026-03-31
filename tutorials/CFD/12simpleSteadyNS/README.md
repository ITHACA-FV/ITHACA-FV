# Tutorial 12

## Introduction
The problem consists of a steady Navier-Stokes problem solved using the SIMPLE algorithm. The setup involves parameterized viscosity and inlet velocities, demonstrating reduced order modeling for steady incompressible flows.

# A detailed look into the code
This is the explanation of the code used for tutorial 12.

## The necessary header files
First of all let's have a look into the header files which have to be included, indicating what they are responsible for:
`<SteadyNSSimple.H>` is the base class for steady NS problems solved with SIMPLE.
`<ITHACAstream.H>` is responsible for reading and exporting the fields and other sorts of data.
`<ITHACAPOD.H>` is for the computation of the POD modes.
`<ReducedSimpleSteadyNS.H>` is for the reduced-order steady NS problem.
`<forces.H>` and `<IOmanip.H>` are for forces computation and I/O manipulation.

## Implementation of the tutorial12 class
We define the tutorial12 class as a child of the `SteadyNSSimple` class.
The constructor is defined with members that are the fields required to be manipulated during the resolution of the full order problem. Such fields are also initialized with the same initial conditions in the solver.
```cpp
class tutorial12 : public SteadyNSSimple
{
    public:
        /// Constructor
        explicit tutorial12(int argc, char* argv[])
            :
            SteadyNSSimple(argc, argv),
            U(_U()),
            p(_p()),
            phi(_phi())
        {}

        /// Velocity field
        volVectorField& U;
        /// Pressure field
        volScalarField& p;
        ///
        surfaceScalarField& phi;
```
Inside the tutorial12 class we define the offlineSolve method. If the offline solve has been previously performed then the method just reads the existing snapshots. If not, it loops over the viscosity parameters, changes the viscosity, and performs the full-order simulations.
```cpp
        /// Perform an Offline solve
        void offlineSolve()
        {
            Vector<double> inl(0, 0, 0);
            List<scalar> mu_now(1);

            // if the offline solution is already performed read the fields
            if (offline)
            {
                ITHACAstream::read_fields(Ufield, U, "./ITHACAoutput/Offline/");
                ITHACAstream::read_fields(Pfield, p, "./ITHACAoutput/Offline/");
                mu_samples =
                    ITHACAstream::readMatrix("./ITHACAoutput/Offline/mu_samples_mat.txt");
            }
            else
            {
                Vector<double> Uinl(1, 0, 0);
                label BCind = 0;

                for (label i = 0; i < mu.cols(); i++)
                {
                    mu_now[0] = mu(0, i);
                    change_viscosity(mu(0, i));
                    assignIF(U, Uinl);
                    truthSolve2(mu_now);
                }
            }
        }
```

## Definition of the main function
The main function sets up the problem parameters, performs the offline phase, computes POD modes, and solves the online reduced problem.

First, the tutorial object is constructed and parameters are read:
```cpp
tutorial12 example(argc, argv);
ITHACAparameters* para = ITHACAparameters::getInstance(example._mesh(),
                     example._runTime());
int NmodesUout = para->ITHACAdict->lookupOrDefault<int>("NmodesUout", 15);
int NmodesPout = para->ITHACAdict->lookupOrDefault<int>("NmodesPout", 15);
int NmodesUproj = para->ITHACAdict->lookupOrDefault<int>("NmodesUproj", 10);
int NmodesPproj = para->ITHACAdict->lookupOrDefault<int>("NmodesPproj", 10);
word filename("./par");
example.mu = ITHACAstream::readMatrix(filename);
example.inletIndex.resize(1, 2);
example.inletIndex(0, 0) = 0;
example.inletIndex(0, 1) = 0;
```

The offline solve is performed:
```cpp
example.offlineSolve();
ITHACAstream::read_fields(example.liftfield, example.U, "./lift/");
example.computeLift(example.Ufield, example.liftfield, example.Uomfield);
ITHACAPOD::getModes(example.Uomfield, example.Umodes, example._U().name(),
                    example.podex, 0, 0, NmodesUout);
ITHACAPOD::getModes(example.Pfield, example.Pmodes, example._p().name(),
                    example.podex, 0, 0, NmodesPout);
```

The reduced problem is set up and solved online for different viscosities:
```cpp
reducedSimpleSteadyNS reduced(example);
word vel_file(para->ITHACAdict->lookup("online_velocities"));
Eigen::MatrixXd vel = ITHACAstream::readMatrix(vel_file);

for (label k = 0; k < (example.mu).size(); k++)
{
    scalar mu_now = example.mu(0, k);
    example.change_viscosity(mu_now);
    reduced.setOnlineVelocity(vel);
    reduced.solveOnline_Simple(mu_now, NmodesUproj, NmodesPproj);
}
```

This completes the tutorial for steady NS with SIMPLE algorithm.

## The plain code
The plain code is available [here](https://raw.githubusercontent.com/ITHACA-FV/ITHACA-FV/refs/heads/master/tutorials/CFD/12simpleSteadyNS/12simpleSteadyNS.C).