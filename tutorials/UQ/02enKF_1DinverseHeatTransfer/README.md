# Tutorial UQ-02

## Introduction
This example demonstrates the use of the Ensemble Kalman Filter to
reconstruct the state of a one–dimensional heat conduction problem when the
boundary flux is only partially known and noisy measurements of the
temperature are available inside the domain.

The mathematical model for the true problem is

$$
\rho C_p \frac{\partial T}{\partial t} + k \Delta T = 0, \quad x\in(0,1)
$$

with a time‑dependent Neumann condition at $x=0$ and a fixed Dirichlet
value of 0.5 at $x=1$. The `true` flux
$g_{true}(t)=5t+2$ differs from the model flux $g(t)=5t$ used in the
filter. Initial temperature $T_0=1$.

Measurements (noisy temperature values) are specified in
`measurementsDict` and read by the solver. The state and error covariances are
assumed Gaussian.

## Code structure
The core of the tutorial is implemented in
`02enKF_1DinverseHeatTransfer.C`, which defines an `EnKF_1DinverseHeatTransfer`
class derived from `laplacianProblem`. The `main` routine:

1. Instantiates the example object with a specified ensemble size.
2. Reads physical parameters (`k`, `rho`, `Cp`) from the dictionary.
3. Solves the true forward model using the correct flux.
4. Configures prior, model error and measurement noise Gaussian densities.
5. Draws prior samples, and iteratively forecasts and updates the ensemble
   via `ITHACAmuq::muq2ithaca::EnsembleKalmanFilter`.
6. Computes posterior means and confidence intervals, which are exported to
   `./ITHACAoutput/`.

Extensive documentation is embedded in the source file, including
LaTeX‑style equations describing the algorithm and example plots.

## Compilation and execution
Build the executable using the makefile in `Make/`:

```bash
cd tutorials/UQ/02enKF_1DinverseHeatTransfer/Make
make -j
```

Run from the tutorial root:

```bash
./02enKF_1DinverseHeatTransfer
```

Dictionary entries control the mesh, time step and ensemble parameters.
Adjusting the ITHACAdict also enables or disables additional utility tests
(see the source comments for options such as `parameterizedBCtest` and
`CGtest`).

## Post‑processing
The Python script `plots.py` in the tutorial folder visualises the reconstructed
temperature at a point, along with 95% confidence bands. The `0/` subdirectory
contains example initial and boundary  condition files used by the solver.

## The plain code
The plain code is available [here](https://raw.githubusercontent.com/ITHACA-FV/ITHACA-FV/refs/heads/master/tutorials/UQ/02enKF_1DinverseHeatTransfer/02enKF_1DinverseHeatTransfer.C).