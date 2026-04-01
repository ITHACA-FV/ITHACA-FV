# Tutorial 25
## Introduction
This tutorial shows a fluid–structure interaction (FSI) problem where a cylinder moves in a viscous flow. The ITHACA-FV implementation leverages the `fsiBasic` base class and constructs a reduced order model capturing both the fluid velocity and the solid point displacement.

# A detailed look into the code
Let's have a look at the code of tutorial 25.

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

## Implementation of `tutorial25` class
The `tutorial25` class inherits from `fsiBasic` and holds references to the velocity `U`, pressure `p` and point displacement `pd` fields.

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
```
The `offlineSolve` method either loads existing snapshots or iterates over parameter samples, updating a damping parameter and running the coupled solver (`updateStiffnessAndRebuildSolver(mu(i, 0)`). After each simulation the code exports force and centre-of-mass data for post–processing.

```cpp
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

## The main
The `main` routine reads offline parameter files (`parsOff_mat.txt`):
```cpp
   std::ifstream exFileOff("./parsOff_mat.txt");
```
sets the number of modes for velocity, pressure and displacement:
```cpp
    int NmodesUout = para->ITHACAdict->lookupOrDefault<int>("NmodesUout", 40);
    int NmodesPout = para->ITHACAdict->lookupOrDefault<int>("NmodesPout", 40);
    int NmodesDout = para->ITHACAdict->lookupOrDefault<int>("NmodesDout", 40);
    int NmodesUproj = para->ITHACAdict->lookupOrDefault<int>("NmodesUproj", 15);
    int NmodesPproj = para->ITHACAdict->lookupOrDefault<int>("NmodesPproj", 5);
    int NmodesDproj = para->ITHACAdict->lookupOrDefault<int>("NmodesDproj", 1);
```
and optionally performs offline runs:
```cpp
    example.offlineSolve();
```
POD modes are computed unless they already exist:
```cpp
    if (example.podex == 0 )
    {
        ITHACAPOD::getModes(example.Ufield, example.Umodes, example._U().name(),
                            example.podex, 0, 0, NmodesUout);
        //ITHACAPOD::getModes(UfieldEnrich, online.Umodes, online._U().name(),
        //            example.podex, 0, 0, NmodesUout);
        ITHACAPOD::getModes(example.Pfield, example.Pmodes, example._p().name(),
                            example.podex, 0, 0, NmodesPout);
        ITHACAPOD::getModes(example.Dfield, example.Dmodes,
                            example._pointDisplacement().name(),
                            example.podex, 0, 0, NmodesDout);
    }
    else
    {
        ITHACAstream::read_fields(example.Umodes, example._U(),
                                  "./ITHACAoutput/POD/");
        ITHACAstream::read_fields(example.Pmodes, example._p(),
                                  "./ITHACAoutput/POD/");
        ITHACAstream::read_fields(example.Dmodes, example._pointDisplacement(),
                                  "./ITHACAoutput/POD/");
    }
```
and additional testing data for online evaluation (parameters sampling, FOM counterpart) are generated using a second instance of the tutorial object.

```cpp
Eigen::MatrixXd parsOn;
    std::ifstream exFileOn("./parsOn_mat.txt");

    if (exFileOn)
    {
        parsOn  = ITHACAstream::readMatrix("./parsOn_mat.txt");
    }

    /// Generate data for testing
    word test_folder = "./ITHACAoutput/TestingOff/";

    /// include Symbolinks
    if (!ITHACAutilities::check_folder(test_folder))
    {
        tutorial22 testkOff(argc, argv);
        testkOff.mu  = parsOn;
        //Perform the offline solve
        testkOff.offline = false;

        for (label k = 0; k < parsOn.rows(); k++)
        {
            word param_name = "damping";
            testkOff.updateStiffnessAndRebuildSolver(parsOn(k, 0), param_name);
            testkOff.startTime  = 0;
            testkOff.finalTime  = 30;
            testkOff.timeStep   = 0.003; //0.0025; //4e-03; // ok dt=0.001
            testkOff.writeEvery = 1e-01;
            testkOff.truthSolve(k, test_folder);
            word localFolder = test_folder + "/DataFromFoam_" + name(k + 1);
            testkOff.prepareFoamData(localFolder);
            testkOff.restart();
            testkOff.Ufield.clear();
            testkOff.Pfield.clear();
            testkOff.Dfield.clear();
            testkOff.fomforcex.clear();
            testkOff.fomforcey.clear();
            testkOff.centerofmassy.clear();
        }
        ...
    }
```

Online reduction is performed with the `ReducedFsi` class
```cpp
ReducedFsi reduced(example);
```   
It is solving the reduced PIMPLE system for each test parameter set and exporting coefficients and reconstructions.
```cpp
if (exFileOn)
    {
        parsOn  = ITHACAstream::readMatrix("./parsOn_mat.txt");
    }

    word folder = "./ITHACAoutput/Online/";
    ITHACAutilities::createSymLink("./0", folder);
    ITHACAutilities::createSymLink("./system", folder);
    ITHACAutilities::createSymLink("./constant", folder);

    for (label i = 0; i < parsOn.rows(); i++)
    {
        // Info << "Current mu = "
        //           << parsOn(i,0) << endl;
        //example.meshPtr().movePoints( example.Mesh0->points());
        word param_name = "damping";
        example.updateStiffnessAndRebuildSolver(parsOn(i, 0), param_name);
        reduced.startTime  = 0; //example.startTime;
        reduced.finalTime  = 30; //example.finalTime;;
        reduced.timeStep   = 0.003; //example.timeStep;
        reduced.writeEvery = 1e-1; //example.writeEvery;
        /// Solving the reduced problem
        reduced.solveOnline_Pimple(NmodesUproj,
                                   NmodesPproj,
                                   NmodesDproj,
                                   folder);
        word localFolder = folder +  name(i + 1);

        for (int k = 0; k < reduced.UredFields.size(); ++k)
        {
            ITHACAstream::exportSolution(reduced.UredFields[k], name(k + 1),
                                         localFolder);
            ITHACAstream::exportSolution(reduced.PredFields[k], name(k + 1),
                                         localFolder);
            ITHACAstream::exportSolution(reduced.Dfield[k],     name(k + 1),
                                         localFolder);
            ITHACAstream::writePoints(reduced.ListOfpoints[k], localFolder,
                                      name(k + 1) + "/polyMesh/");
        }

        word DataRom = folder + "../" + "/DataFromRom_" + name(i + 1);
        reduced.prepareRomData(DataRom);
        /// Restart the simulation
        example.restart();
        reduced.UredFields.clear();
        reduced.PredFields.clear();
        reduced.Dfield.clear();
        reduced.romforcey.clear();
        reduced.romforcex.clear();
        reduced.centerofmassy.clear();
    }
```
Error metrics comparing reduced and full-order solutions (energy norms, centre-of-mass displacement, forces) are computed and stored for post-processing:
```cpp
    Eigen::MatrixXd errL2U = ITHACAutilities::errorL2Rel(example.Ufield,
                             reduced.UredFields);
    Eigen::MatrixXd errL2P = ITHACAutilities::errorL2Rel(example.Pfield,
                             reduced.PredFields);
    cnpy::save(errL2U, "./ITHACAoutput/DataFromRom/errL2U_" + name(
                   NmodesUproj) + "_" + name(NmodesPproj) + ".npy");
    cnpy::save(errL2P, "./ITHACAoutput/DataFromRom/errL2P_" + name(
                   NmodesUproj) + "_" + name(NmodesPproj) + ".npy");
```

## Usage
Python scripts `plotCentreMass.py`, `plotDrag.py` and `plotLift.py` can be used to visualize the outputs.

This tutorial demonstrates building a coupled fluid-structure reduced model with parameterized structural damping and moving mesh handling.

## The plain code
The plain code is available [here](https://raw.githubusercontent.com/ITHACA-FV/ITHACA-FV/refs/heads/master/tutorials/CFD/25Moving_Cylinder/25Moving_Cylinder.C).