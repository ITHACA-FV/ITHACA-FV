# Tutorial 24

## Introduction
This tutorial illustrates the use of hyper-reduction techniques (GappyDEIM,
Empirical Cubature) to approximate a nonlinear scalar or vector function defined
on a CFD mesh. The example uses analytic Gaussian fields parameterized by two
variables, and demonstrates how to build a reduced representation and perform
fast online evaluations.

## Header files
The implementation depends on the hyper-reduction template and standard
OpenFOAM/ITHACA utilities:
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
Three hyper-reduction classes are defined, specializing the `HyperReduction`
template for scalar, vector and combined vector+scalar fields. Each class
implements `evaluate_expression` routines computing the analytic function on a
field given a parameter vector `mu`, and `onlineCoeffs` which perform the RBF
interpolation using precomputed pseudo-inverse matrices.

```cpp
class HyperReduction_function : public HyperReduction<PtrList<volScalarField> &>
{ ... };

class HyperReduction_vectorFunction : public HyperReduction<PtrList<volVectorField> &>
{ ... };

class HyperReduction_vectorScalarFunction : public HyperReduction<PtrList<volVectorField> &, PtrList<volScalarField> &>
{ ... };
```

## Main workflow
The `test_scalar` function (called from `main`) performs the offline and online
stages for a chosen field type. Key steps include:

* Loading parameter samples (`trainingPars.npy`, `testingPars.npy`).
* Generating offline snapshots by evaluating the analytic expression at each
  parameter and exporting them to `ITHACAoutput/Offline/`.
* Building a `HyperReduction` object with specified number of modes and nodes.
* Computing POD modes of the snapshots and performing the offline
  hyper-reduction (`offlineGappyDEIM` or `offlineECP`).
* Generating submeshes and evaluating the nonlinear expression on the magic
  points.
* Online evaluation: interpolating field values at new parameter samples and
  reconstructing full fields. Relative L2 errors are computed against true
  values.

Different hyper-reduction methods can be selected via the `HyperReduction`
entry in the ITHACAdict file (`GappyDEIM` or `ECP`). The number of modes and
nodes is configurable.

## Usage
Compile the source with the provided `Make/` target. Prepare parameter files in
NumPy format (`trainingPars.npy`, `testingPars.npy`) and place them in the
current directory. Run the resulting executable; outputs (offline snapshots,
onlinearly reconstructed fields, error metrics) will be written under the
`ITHACAoutput/` directory.

Adjust hyper-reduction settings in the dictionary and regenerate snapshots if
needed. The tutorial acts as a generic template for approximating arbitrary
nonlinear functions in a reduced fashion.