# Tutorial IHTP‑01

## Introduction
This tutorial demonstrates inverse problem solving for a one–dimensional heat transfer problem. The objective is to reconstruct unknown boundary fluxes (heat flux at the left wall) and/or the internal state of the domain given noisy temperature measurements at selected points (thermocouples).

Two methods are implemented and compared:

1. **Alifanov's conjugate gradient regularization** (class `inverseLaplacianProblem_CG`)
2. **Parametrization via radial basis functions (RBFs)** (class `inverseLaplacianProblem_paramBC`)

Both are tested on a analytical benchmark problem with known solution, allowing error assessment.

## Problem formulation
The heat equation governs the domain

$$
\rho C_p \frac{\partial T}{\partial t} + k \Delta T = 0, \quad x \in (0,1)
$$

with boundary conditions:

* Left (hot side, $x=0$): heat flux $g(t)$ (unknown, to be identified);
* Right (cold side, $x=1$): fixed temperature 0.5;
* Initial: $T_0$ (spatially uniform).

Temperature measurements $y_i(t)$ are collected at locations $x_i$ inside the domain. The inverse problem is to find $g(t)$ and/or the state $T(x,t)$ that minimize the discrepancy between measurements and predictions, subject to the PDE constraint.

## Reference solution
An analytical solution is prescribed as

$$
T_{true}(x,y,z) = ax^2 + bxy + cy - az^2 + c
$$

for constants $a,b,c,d$ (set to 5, 10, 15, 20 by default). The flux $g_{true}$ is derived from this solution. During the inverse solve, a slightly different flux $g$ is assumed to intentionally introduce model
error.

## Methods overview

### Alifanov's Conjugate Gradient (CG)
This method solves

$$
\min_g \Vert H(g) - y_{meas} \Vert^2 + \lambda \Vert g \Vert^2,
$$

where $H$ is the forward operator (solving the forward heat equation) and $\lambda$ is a regularization parameter. The CG algorithm iteratively refines the flux estimate using adjoint gradients.

Implemented in `inverseLaplacianProblem_CG` class and tested via the `#include "CGtest.H"` section of the main.

### RBF Parametrization
Instead of solving for $g(t)$ directly at every boundary node, the flux is parametrized as a linear combination of basis functions centred at control points. The RBF parameters are then optimized.

Implemented in `inverseLaplacianProblem_paramBC` class; tested via `#include "parameterizedBCtest.H"`.

## Files of interest

* `IHTP01inverseLaplacian.C` – main executable; orchestrates case setup, solves reference solution, and runs inverse tests.
* `IHTP01inverseLaplacian.H` – class definitions (if needed).
* `CGtest.H` – test code for conjugate‑gradient inversion.
* `parameterizedBCtest.H` – test code for RBF‑based parametrization.
* `thermocouplesLocation_CG.H`, `thermocouplesLocation_paramBC.H` – tests
  investigating influence of thermocouple position.
* `thermocouplesNumberTest_CG.H`, `thermocouplesNumberTest_paramBC.H` – tests
  investigating influence of number of measurement points.
* `caseDir/` – case directory with `system/`, `constant/` etc.

## A detailed look into the `IHTP01inverseLaplacian.C` code 

### Header files

```cpp
#include <iostream>
#include <float.h>
#include "interpolation.H"
#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include "IOmanip.H"
#include "Time.H"
#include "inverseLaplacianProblem_CG.H"
#include "inverseLaplacianProblem_paramBC.H"
#include "ITHACAPOD.H"
#include "ITHACAutilities.H"
#include <Eigen/Dense>
#define _USE_MATH_DEFINES
#include <cmath>
#include "Foam2Eigen.H"
#include "mixedFvPatchFields.H"
#include "cellDistFuncs.H"
```

The included headers provide OpenFOAM finite-volume tools, inverse Laplacian problem classes, POD and utility routines from ITHACA-FV, interpolation and field-manipulation utilities, and the tutorial-specific header `IHTP01inverseLaplacian.H`. Together, these dependencies support both the forward heat-transfer solve and the inverse reconstruction procedures.

### Main

In main, two problem objects are created: `example_paramBC`, which uses a parameterized boundary condition, and `example_CG`, which uses a conjugate-gradient regularization method.
```cpp
    IHTP01inverseLaplacian_paramBC example_paramBC(argc, argv);
    IHTP01inverseLaplacian_CG example_CG(argc, argv);
```

The coefficients `a`, `b`, `c`, and `d` are then assigned to both objects in order to define the analytical benchmark used throughout the tutorial.
```cpp
    double a = 5;
    double b = 10;
    double c = 15;
    double d = 20;
    example_paramBC.a = a;
    example_paramBC.b = b;
    example_paramBC.c = c;
    example_paramBC.d = d;
    example_CG.a = a;
    example_CG.b = b;
    example_CG.c = c;
    example_CG.d = d;
```

The code reads from `ITHACAdict` which tests must be executed. These switches control whether the program performs the standard conjugate-gradient inversion, the parameterized boundary-condition inversion, sensitivity studies with respect to the RBF width, and tests that vary the location or number of thermocouples.
```cpp
    ITHACAparameters* para = ITHACAparameters::getInstance(example_paramBC._mesh(),
                             example_paramBC._runTime());
    label CGtest = para->ITHACAdict->lookupOrDefault<int>("CGtest", 0);
    label CGnoiseTest = para->ITHACAdict->lookupOrDefault<int>("CGnoiseTest", 0);
    label CGnoiseLevelTest =
        para->ITHACAdict->lookupOrDefault<int>("CGnoiseLevelTest", 0);
    label parameterizedBCtest =
        para->ITHACAdict->lookupOrDefault<int>("parameterizedBCtest", 0);
    label parameterizedBCtest_RBFwidth =
        para->ITHACAdict->lookupOrDefault<int>("parameterizedBCtest_RBFwidth", 0);
    label thermocouplesLocationTest_CG =
        para->ITHACAdict->lookupOrDefault<int>("thermocouplesLocationTest_CG", 0);
    label thermocouplesLocationTest_paramBC =
        para->ITHACAdict->lookupOrDefault<int>("thermocouplesLocationTest_paramBC", 0);
    label thermocouplesNumberTest_CG =
        para->ITHACAdict->lookupOrDefault<int>("thermocouplesNumberTest_CG", 0);
    label thermocouplesNumberTest_paramBC =
        para->ITHACAdict->lookupOrDefault<int>("thermocouplesNumberTest_paramBC", 0);
```

At the same time, the code reads the physical and numerical parameters needed by the inverse solvers, such as conductivity, heat-transfer coefficient, tolerances, maximum CG iterations, and the RBF shape parameter.
```cpp
    example_CG.cgIterMax = para->ITHACAdict->lookupOrDefault<int>("cgIterMax", 100);
    example_CG.Jtol =  para->ITHACAdict->lookupOrDefault<double>("Jtolerance",
                       0.000001);
    example_CG.JtolRel =
        para->ITHACAdict->lookupOrDefault<double>("JrelativeTolerance",
            0.001);
    double rbfShapePar = para->ITHACAdict->lookupOrDefault<double>("rbfShapePar",
                         0);
    M_Assert(rbfShapePar > 0, "rbfShapePar not specified");
    example_paramBC.k =
        para->ITHACAdict->lookupOrDefault<double>("thermalConductivity", 0);
    M_Assert(example_paramBC.k > 0, "thermalConductivity, k, not specified");
    example_paramBC.H =
        para->ITHACAdict->lookupOrDefault<double>("heatTranferCoeff", 0);
    M_Assert(example_paramBC.H > 0, "Heat transfer coeff, H, not specified");
    example_CG.k = example_paramBC.k;
    example_CG.H = example_paramBC.H;
    int rbfWidthTest_size =
        para->ITHACAdict->lookupOrDefault<int>("rbfWidthTest_size", 0);
```

An analytical reference temperature field `T_true` is then constructed directly from the mesh coordinates. This field acts as the exact benchmark solution. Using it, the code defines the corresponding exact boundary heat flux on the hot side through the analytical expression $g(x)=k(bx+c)$. This exact flux is assigned to both inverse-problem objects and used to generate the true forward solution.
```cpp
    volScalarField T_true(example_paramBC._T());

    for (label i = 0; i < T_true.internalField().size(); i++)
    {
        auto cx = T_true.mesh().C()[i].component(vector::X);
        auto cy = T_true.mesh().C()[i].component(vector::Y);
        auto cz = T_true.mesh().C()[i].component(vector::Z);
        T_true.ref()[i] = a * cx * cx + b * cx * cy + c * cy - a * cz * cz + c;
    }
```

The methods `solveTrue()` are then called for both objects to compute the numerical reference solution associated with the exact heat flux.
```cpp
    fvMesh& mesh = example_paramBC._mesh();
    example_paramBC.hotSide_ind = mesh.boundaryMesh().findPatchID("hotSide");
    label hotSideSize = T_true.boundaryField()[example_paramBC.hotSide_ind].size();
    example_paramBC.g.resize(hotSideSize);
    example_paramBC.gTrue.resize(hotSideSize);
    forAll(example_paramBC.g, faceI)
    {
        scalar faceX =
            mesh.boundaryMesh()[example_paramBC.hotSide_ind].faceCentres()[faceI].x();
        example_paramBC.g[faceI] = example_paramBC.k * (b * faceX + c) ;
    }
    example_paramBC.gTrue = example_paramBC.g;
    example_paramBC.solveTrue();
    example_CG.g = example_paramBC.g;
    example_CG.gTrue = example_paramBC.g;
    example_CG.solveTrue();
```

After that, the analytical field and the numerical solution are compared. The code exports the analytical solution and the numerical error field, and prints relative and absolute error norms. This verifies that the direct solver is consistent with the chosen benchmark.
```cpp
    ITHACAstream::exportSolution(T_true, "1", "./ITHACAoutput/true/",
                                 "analyticalSol");
    volScalarField error = (T_true - T).ref();
    ITHACAstream::exportSolution(error, "1", "./ITHACAoutput/true/", "error");
    ...
```

Next, the thermocouple positions are read for both problem formulations, and the synthetic measurements are generated by sampling the analytical field at those locations. These measurements become the data used in the inverse problem.
```cpp
    example_paramBC.readThermocouples();
    example_paramBC.Tmeas = example_paramBC.fieldValueAtThermocouples(T_true);
    example_CG.readThermocouples();
    example_CG.Tmeas = example_CG.fieldValueAtThermocouples(T_true);
```

The rest of the file is organized as a collection of optional tests. If `CGtest` is enabled, the inverse problem is solved with Alifanov’s regularization method. If `parameterizedBCtest` is enabled, the heat flux is reconstructed by parameterizing the boundary condition. Additional tests explore the influence of the RBF width and the sensitivity of both methods to the thermocouple placement and to the number of sensors. Each test is included through a dedicated header file, so the driver remains compact while still supporting multiple study configurations.
```cpp
    if (CGtest)
    {
#include "CGtest.H"
    }

    if (parameterizedBCtest)
    {
#include"parameterizedBCtest.H"
    }

    if (parameterizedBCtest_RBFwidth)
    {
#include"parameterizedBCtest_RBFwidth.H"
    }

    if (thermocouplesLocationTest_CG)
    {
#include"thermocouplesLocation_CG.H"
    }

    if (thermocouplesLocationTest_paramBC)
    {
#include"thermocouplesLocation_paramBC.H"
    }

    if (thermocouplesNumberTest_CG)
    {
#include"thermocouplesNumberTest_CG.H"
    }

    if (thermocouplesNumberTest_paramBC)
    {
#include"thermocouplesNumberTest_paramBC.H"
    }
```

## The plain code
The plain code is available [here](https://raw.githubusercontent.com/ITHACA-FV/ITHACA-FV/refs/heads/master/tutorials/inverseHeatTransfer/IHTP01inverseLaplacian/IHTP01inverseLaplacian.C).