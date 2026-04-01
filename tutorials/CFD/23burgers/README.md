# Tutorial 23

## Introduction
This tutorial considers a parameterized heat-transfer–type problem based on the Burgers equation. Although the class name is `Burgers`, the implemented equations are scalar advection–diffusion; the objective is to compute a reduced basis for the velocity field using a simple offline/online split.

# A detailed look into the code
This section explains the code of tutorial 23.

## Header files
The program includes both standard and ITHACA-FV headers:
```cpp
#include <iostream>
#include "fvCFD.H"
#include "IOmanip.H"
#include "Time.H"
#include "Burgers.H"
#include "ITHACAPOD.H"
#include <Eigen/Dense>
```
`<Burgers.H>` defines the full-order problem derived from `laplacianProblem`.

## Implementation
The `tutorial23` class inherits from `Burgers` and stores a reference to the velocity field `U`:
```cpp
class tutorial23: public Burgers
{
    public:
        explicit tutorial23(int argc, char* argv[])
            : Burgers(argc, argv), U(_U()) {}

        volVectorField& U;

        void offlineSolve(word folder = "./ITHACAoutput/Offline/")
        {
            if (offline)
            {
                ITHACAstream::readMiddleFields(Ufield, U, folder);
            }
            else
            {
                truthSolve(folder);
            }
        }
};
```
During offline execution the code either reads existing snapshots or runs the full-order solver once and saves the resulting velocity fields.

## Main function
In `main` the tutorial object is created and POD mode numbers are extracted from the dictionary file:
```cpp
tutorial23 train(argc, argv);
ITHACAparameters* para = ITHACAparameters::getInstance(train._mesh(),
                     train._runTime());
int NmodesUout = para->ITHACAdict->lookupOrDefault<int>("NmodesUout", 15);
train.offlineSolve();
ITHACAPOD::getModes(train.Ufield, train.Umodes, train._U().name(),
                    train.podex, 0, 0, NmodesUout);
```
The reduced velocity basis is written to `train.Umodes` and can subsequently be used for online reduced computations.


## The plain code
The plain code is available [here](https://raw.githubusercontent.com/ITHACA-FV/ITHACA-FV/refs/heads/master/tutorials/CFD/23burgers/23burgers.C).