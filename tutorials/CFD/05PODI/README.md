# Tutorial 05
## Introduction to POD with Interpolation

In this tutorial the POD is applied on the backstep steady test case with parametrized viscosity.

# A detailed look into the code
In this section are explained the main steps necessary to construct the tutorial N°5.

## The necessary header files

First of all let's have a look to the header files that need to be included and what they are responsible for.

Some of the header files of ITHACA-FV necessary for this tutorial are: `<steadyNS.H>` for the full order steady NS problem, `<ITHACAPOD.H>` for the POD decomposition, `<ITHACAstream.H>` for some ITHACA input-output operations, and also `<bspline.h>` and `<rbfspline.h>` to perform the RBF interpolation in the reduced space.
```cpp
    #include "steadyNS.H"
    #include "ITHACAstream.H"
    #include "ITHACAPOD.H"
    #include "forces.H"
    #include "IOmanip.H"
    #include "bspline.h"
    #include "rbfspline.h"
```
## Implementation of the tutorial05 class

We can define the tutorial05 class as a child of the `<steadyNS>` class. The constructor is defined with members that are the fields need to be manipulated during the resolution of the full order problem using `simpleFoam`. Such fields are also initialized with the same initial conditions in the solver. 
```cpp
    class tutorial05 : public steadyNS
    {
        public:
            explicit tutorial05(int argc, char* argv[])
                : steadyNS(argc, argv), U(_U()), p(_p()) {}
```
Inside the tutorial05 class we define the `offlineSolve` method according to the specific parametrized problem that needs to be solved. If the offline solve has been previously performed, then the method just reads the existing snapshots from the `Offline` directory. Otherwise it loops over all the parameters, changes the system viscosity with the iterable parameter, and then performs the offline solve.
```cpp
        void offlineSolve()
        {
            Vector<double> inl(1, 0, 0);
            List<scalar> mu_now(1);
 
            if (offline)
            {
                ITHACAstream::read_fields(Ufield, U, "./ITHACAoutput/Offline/");
                ITHACAstream::read_fields(Pfield, p, "./ITHACAoutput/Offline/");
                mu_samples =
                ITHACAstream::readMatrix("./ITHACAoutput/Offline/mu_samples_mat.txt");
            
            }
            else
            {
                Vector<double> Uinl(0, 0, 0);
                for (label i = 0; i < mu.cols(); i++)
                {
                    mu_now[0] = mu(0, i);
                    change_viscosity(mu(0, i));
                    assignIF(U, inl);
                    truthSolve(mu_now);
                }
            }
        }
```

## Definition of the main function
In this section we show the definition of the main function. First we construct the object "example" of type tutorial05:
```cpp
    tutorial05 example(argc, argv);
```
Then we parse the `ITHACAdict` file to determine the number of modes to be written out and also the ones to be used for projection of the velocity and the pressure fields: 
```cpp
    ITHACAparameters* para = ITHACAparameters::getInstance(example._mesh(),
                             example._runTime());
    int NmodesUout = para->ITHACAdict->lookupOrDefault<int>("NmodesUout", 15);
    int NmodesPout = para->ITHACAdict->lookupOrDefault<int>("NmodesPout", 15);
```
We note that a default value can be assigned in case the parser did not find the corresponding string in the `ITHACAdict` file.

Then, the parameters for which we collect the snapshots are read from a text file (previously stored) and the offline solver is performed for these parameters.

```cpp
    // Read the par file where the parameters are stored
    word filename("./par");
    example.mu = ITHACAstream::readMatrix(filename);
    // Perform the offline solve
    example.offlineSolve();
```
Then, the velocity and pressure modes are collected:
```cpp
    // Perform the offline solve
    example.offlineSolve();
    // Perform POD on velocity and pressure, and store the first 10 modes
    ITHACAPOD::getModes(example.Ufield, example.Umodes, example._U().name(),
                        example.podex, 0, 0,
                        NmodesUout);
    ITHACAPOD::getModes(example.Pfield, example.Pmodes, example._p().name(),
                        example.podex, 0, 0,
                        NmodesPout);
```
The plain code is available [here](https://raw.githubusercontent.com/ITHACA-FV/ITHACA-FV/refs/heads/master/tutorials/CFD/05PODI/05PODI.C).