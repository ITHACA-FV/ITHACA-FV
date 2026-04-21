# Tutorial 24

## Introduction
This tutorial illustrates the use of hyper-reduction techniques (GappyDEIM, Empirical Cubature) to approximate a nonlinear scalar or vector function defined on a CFD mesh. The example uses analytic Gaussian fields parameterized by two variables, and demonstrates how to build a reduced representation.

# A detailed look into the code
Let's analyse the code of tutorial 24.

## Header files
The implementation depends on the hyper-reduction template and standard OpenFOAM/ITHACA utilities:
```cpp
#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include "fvMeshSubset.H"
#include "ITHACAutilities.H"
#include "Foam2Eigen.H"
#include "ITHACAstream.H"
#include "ITHACAPOD.H"
#include "hyperReduction.templates.H"
```

## Classes
Three hyper-reduction classes are defined, specializing the `HyperReduction` template for scalar (`HyperReduction_function`), vector (`HyperReduction_vectorFunction`) and combined vector+scalar (`HyperReduction_vectorScalarFunction`) fields. Each class implements `evaluate_expression` routines computing the analytic function on a field given a parameter vector `mu`, and `onlineCoeffs` which perform the RBF interpolation using precomputed pseudo-inverse matrices.

```cpp
class HyperReduction_function : public HyperReduction<PtrList<volScalarField> &>
{ ... };

class HyperReduction_vectorFunction : public HyperReduction<PtrList<volVectorField> &>
{ ... };

class HyperReduction_vectorScalarFunction : public HyperReduction<PtrList<volVectorField> &, PtrList<volScalarField> &>
{ ... };
```

The `test_scalar` function is the defined.
```cpp
void test_scalar(ITHACAparameters* para, Foam::fvMesh& mesh,
                 Foam::Time& runTime)
```
It performs the offline and online stages for a chosen field type, evaluated in a specific mesh and at a specific time.

Key steps include:

* Initialization of the field considered:
```cpp
int n_modes = para->ITHACAdict->lookupOrDefault<int>("Modes", 15);
int n_nodes = para->ITHACAdict->lookupOrDefault<int>("Nodes", 15);
simpleControl simple(mesh);
word methodName = para->ITHACAdict->lookupOrDefault<word>("HyperReduction",
                  "GappyDEIM");
#include "createFields.H"
// List of volScalarField where the snapshots are stored
PtrList<volScalarField> Sp;
// Non linear function
volScalarField S
(
    IOobject
    (
        "S",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    T.mesh(),
    dimensionedScalar("zero", dimensionSet(0, 0, -1, 1, 0, 0, 0), 0)
);
```

* Creating or loading parameter samples to train the nonlinear function(`trainingPars.npy`), for example:
```cpp
    Eigen::MatrixXd pars = ITHACAutilities::rand(100, 2, -0.5, 0.5);
```
* Generating offline snapshots by evaluating the analytic expression at each parameter and exporting them to `ITHACAoutput/Offline/`.

```cpp
    for (int i = 0; i < 100; i++)
    {
        HyperReduction_function::evaluate_expression(S, pars.row(i));
        Sp.append((S).clone());
        ITHACAstream::exportSolution(S, name(i + 1), "./ITHACAoutput/Offline/");
    }
```
* Reading the test parameters and initializing the seeds for magic points selection:
```cpp
    Eigen::MatrixXd parTest;
    cnpy::load(parTest, "testingPars.npy");
    Eigen::VectorXi initSeeds;
```
* Building a `HyperReduction` object with specified number of modes and nodes.
```cpp
HyperReduction_function c(n_modes, n_nodes, initSeeds, "Gaussian_function", Sp);
```
* Computing POD modes of the snapshots and performing the offline
  hyper-reduction, either with `GappyDEIM`:
```cpp
  if (methodName == "GappyDEIM")
    {
        HyperReduction_function c(n_modes, n_nodes, initSeeds, "Gaussian_function", Sp);
        Eigen::MatrixXd snapshotsModes;
        Eigen::VectorXd normalizingWeights;
        c.getModesSVD(c.snapshotsListTuple, snapshotsModes, normalizingWeights);
        c.offlineGappyDEIM(snapshotsModes, normalizingWeights);
        ...
        ...
    }
```
  or with `ECP`:
  ```cpp
  else if (methodName == "ECP")
    {
        HyperReduction_function c(n_modes, n_nodes, initSeeds, "Gaussian_function", Sp);
        bool saveOFModes{true};
        Eigen::MatrixXd snapshotsModes;
        Eigen::VectorXd normalizingWeights;
        c.getModesSVD(c.snapshotsListTuple, snapshotsModes, normalizingWeights,
                      saveOFModes);
        c.offlineECP(snapshotsModes, normalizingWeights);
        ...
  ```
  GappyDEIM reconstructs a nonlinear field by evaluating it only at a small set of selected points and projecting those values onto a reduced basis, significantly reducing computational cost while maintaining accuracy.
  ECP (Empirical Cubature Procedure) approximates integrals of nonlinear functions by selecting a reduced set of weighted quadrature points, enabling efficient and accurate evaluation of reduced-order models without assembling full-order integrals.

* The function generates submeshes and for each test parameter, it evaluates the nonlinear expression on the magic points and reconstructs the field. The error w.r.t. the FOM is finally computed.
  ```cpp
        // Generate the submeshes with the depth of the layer
        unsigned int layers{2};
        c.generateSubmesh(layers, mesh);
        auto sfield = c.interpolateField<volScalarField>(Sp[0]);
        Eigen::VectorXd results(parTest.rows());

        for (unsigned int idTest = 0; idTest < parTest.rows(); idTest++)
        {
            // Online evaluation of the non linear function
            Eigen::VectorXd f = c.evaluate_expression(sfield(), parTest.row(idTest),
                                c.localNodePoints);
            Eigen::VectorXd ff = Foam2Eigen::field2Eigen(c.evaluate_expression(Sp[0],
                                 parTest.row(idTest)));
            double trueIntegral = (normalizingWeights.cwiseInverse().asDiagonal() *
                                   ff).array().sum();
            double testIntegral = (c.wPU * f).array().sum();
            // Info << "Integral: " << trueIntegral << endl;
            // Info << "Reconstructed Integral: " << testIntegral << endl;
            results(idTest) = trueIntegral - testIntegral;
        }
        ...
    }
  ```

The same functionality is repeated for vector fields (`test_vector`) and for scalar+vector fields (`test_vector_scalar`).

## The main
This registers a command-line option -test with three possible values:

- scalar → run scalar test case
- vector → run vector test case
- vs → run combined vector-scalar test case

```cpp
argList::addOption("test", "scalar", "Perform scalar test");
argList::addOption("test", "vector", "Perform vector test");
argList::addOption("test", "vs", "Perform (vector, scalar) test");
```
So basically the script should be run as `24HyperReduction -test option`. 

Then, it creates a hash table listing valid values for the "test" option.
```cpp
HashTable<string> validParOptions;
validParOptions.set("test", "scalar");
validParOptions.set("test", "vector");
validParOptions.set("test", "vs");
```

Then, it initializes the Openfoam case, the simulation time object, and reads the mesh:
```cpp
#include "setRootCase.H"
Foam::Time runTime(Foam::Time::controlDictName, args);
Foam::fvMesh mesh
(
    Foam::IOobject
    (
        Foam::fvMesh::defaultRegion,
        runTime.timeName(),
        runTime,
        Foam::IOobject::MUST_READ
    )
);
```
Then, it reads the parameters from `ITHACAdict`, and runs the hyperreduction, differentiating the function depending on the `test` option given as input.

```cpp
ITHACAparameters* para = ITHACAparameters::getInstance(mesh, runTime);
if (args.get("test").match("scalar"))
{
    Info << "Init scalar testcase\n";
    test_scalar(para, mesh, runTime);
}
else if (args.get("test").match("vector"))
{
    Info << "Init vector testcase\n";
    test_vector(para, mesh, runTime);
}
else if (args.get("test").match("vs"))
{
    Info << "Init vector scalar testcase\n";
    test_vector_scalar(para, mesh, runTime);
}
else
{
    Info << "Pass '-test scalar', '-test vector', '-test vs'" << endl;
}
```

## The plain code
The plain code is available [here](https://raw.githubusercontent.com/ITHACA-FV/ITHACA-FV/refs/heads/master/tutorials/CFD/24HyperReduction/24HyperReduction.C).