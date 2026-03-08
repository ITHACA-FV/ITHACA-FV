# Tutorial UQ-01

## Introduction
This small example tests the Ensemble Kalman Filter (EnKF) wrapper included
in ITHACA‑FV/MUQ. A simple linear dynamical system
$\dot x = A x$ is considered and only partial observations are available
through a matrix $H$. The goal is to reconstruct the full state starting
from an incorrect model matrix $A_{wrong}$.

## Problem setup
The state evolves according to

$$ \frac{dx}{dt} = A x $$

with initial condition $x(0)=x_0$. Observations at discrete times are
obtained as

$$ y = H x $$

Matrices required by the tutorial are stored in text files under the tutorial
root:

* `A_mat.txt` – true system matrix
* `Awrong_mat.txt` – incorrect matrix used for forecasting
* `observation_mat.txt` – observation operator $H$
* `initialState_mat.txt` – initial state vector

The code generates a synthetic trajectory, samples it every few steps and
then applies the EnKF to produce posterior ensembles, means and confidence
intervals.

## Implementation highlights
The entire logic resides in `01enKF.C`. After reading the matrices the
simulation is marched in time with a simple forward Euler integrator. The
ensemble Kalman update is performed using
`ITHACAmuq::muq2ithaca::EnsembleKalmanFilter`, which takes the forecasted
ensemble, measurements and noise covariance.

Results (time vector, true trajectory, posterior mean and quantiles) are
exported to `./ITHACAoutput/` in Eigen text format; a Python script
`plots.py` is included to visualise them.

## Building and running
Compile the executable with the provided makefile in `Make/`:

```bash
cd tutorials/UQ/01enKF/Make
make -j
```

Run the tutorial from its root directory:

```bash
./01enKF
```

Ensure the input matrices are present in the same folder. After execution the
`ITHACAoutput/` directory will contain the exported data, which can be
plotted using `plots.py` or any other tool.

---
*(README created to document the EnKF linear system tutorial.)*