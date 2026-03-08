# Tutorial IHTP‑01

## Introduction
This tutorial demonstrates inverse problem solving for a one–dimensional
heat transfer problem. The objective is to reconstruct unknown boundary
fluxes (heat flux at the left wall) and/or the internal state of the domain
given noisy temperature measurements at selected points (thermocouples).

Two methods are implemented and compared:
1. **Alifanov's conjugate gradient regularization** (class
   `inverseLaplacianProblem_CG`)
2. **Parametrization via radial basis functions (RBFs)** (class
   `inverseLaplacianProblem_paramBC`)

Both are tested on a analytical benchmark problem with known solution,
allowing error assessment.

## Problem formulation
The heat equation governs the domain

$$
\rho C_p \frac{\partial T}{\partial t} + k \Delta T = 0, \quad x\in(0,1)
$$

with boundary conditions:
* Left (hot side, $x=0$): heat flux $g(t)$ (unknown, to be identified)
* Right (cold side, $x=1$): fixed temperature 0.5
* Initial: $T_0$ (spatially uniform)

Temperature measurements $y_i(t)$ are collected at locations
$x_i$ inside the domain. The inverse problem is to find $g(t)$
and/or the state $T(x,t)$ that minimize the discrepancy between measurements
and predictions, subject to the PDE constraint.

## Reference solution
An analytical solution is prescribed as

$$
T_{true}(x,y,z) = ax^2 + bxy + cy - az^2 + c
$$

for constants $a,b,c,d$ (set to 5, 10, 15, 20 by default). The flux
$g_{true}$ is derived from this solution. During the inverse solve, a
slightly different flux $g$ is assumed to intentionally introduce model
error.

## Methods overview

### Alifanov's Conjugate Gradient (CG)
This method solves

$$
\min_g \Vert H(g) - y_{meas} \Vert^2 + \lambda \Vert g \Vert^2,
$$

where $H$ is the forward operator (solving the forward heat equation) and
$\lambda$ is a regularization parameter. The CG algorithm iteratively
refines the flux estimate using adjoint gradients.

Implemented in `inverseLaplacianProblem_CG` class and tested via the
`#include "CGtest.H"` section of the main.

### RBF Parametrization
Instead of solving for $g(t)$ directly at every boundary node, the flux is
parametrized as a linear combination of basis functions centred at control
points. The RBF parameters are then optimized.

Implemented in `inverseLaplacianProblem_paramBC` class; tested via
`#include "parameterizedBCtest.H"`.

## Files of interest
* `IHTP01inverseLaplacian.C` – main executable; orchestrates case setup,
  solves reference solution, and runs inverse tests.
* `IHTP01inverseLaplacian.H` – class definitions (if needed).
* `CGtest.H` – test code for conjugate‑gradient inversion.
* `parameterizedBCtest.H` – test code for RBF‑based parametrization.
* `thermocouplesLocation_CG.H`, `thermocouplesLocation_paramBC.H` – tests
  investigating influence of thermocouple position.
* `thermocouplesNumberTest_CG.H`, `thermocouplesNumberTest_paramBC.H` – tests
  investigating influence of number of measurement points.
* `caseDir/` – case directory with `system/`, `constant/` etc.

## Compilation and execution
Build the solver:

```bash
cd tutorials/inverseHeatTransfer/IHTP01inverseLaplacian/Make
make -j
```

Run from the tutorial root:

```bash
./IHTP01inverseLaplacian
```

The ITHACAdict controls which tests to run:

```text
CGtest              1;           // enable CG regularization
parameterizedBCtest 1;           // enable RBF parametrization
thermocouplesLocationTest_CG     1;  // vary thermocouple distance
thermocouplesLocationTest_paramBC 1;  // ... with RBF
thermocouplesNumberTest_CG       1;  // vary number of thermocouples
thermocouplesNumberTest_paramBC   1;  // ... with RBF
```

Other settings include:

```text
cgIterMax    100;             // max CG iterations
Jtolerance   1e-6;            // cost function tolerance
JrelativeTolerance 1e-3;      // relative cost tolerance
rbfShapePar  0.5;             // RBF shape parameter
thermalConductivity 1.0;      // k
heatTranferCoeff 10.0;        // H (convective heat transfer coeff)
```

## Output
Results are exported to `./ITHACAoutput/`, including:

* `true/` – reference solution (analytical) and errors.
* `CG_results/` – CG method results (reconstructed flux, state).
* `paramBC_results/` – RBF parametrization results.
* Diagnostic data (convergence, parameter values, etc.).

Compare the two methods by examining their respective error norms
(L2 relative, L∞ relative) as computed at runtime.

---
*(README summarising inverse Laplacian heat transfer tutorial.)*