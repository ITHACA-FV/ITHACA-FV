# Tutorial NN‑03
## Introduction
This tutorial demonstrates data‑driven turbulence modelling in a
compressible steady Navier–Stokes problem using a neural network. The problem
geometry is an open domain with a flow field parameterized by Mach number and
other parameters. The network learns to map from reduced velocity coefficients
(and parameters) to eddy‑viscosity coefficients, enabling a ROM with
inexpensive predictions.

The solver class `CompressibleSteadyNN` (in `03compSteadyNS.C`) follows the
same design as `NN‑02`: offline snapshot generation, POD mode extraction,
coefficient computation and projection, neural‑network training, and online
evaluation.

## Workflow
1. **Full‑order offline runs**: solve the compressible SIMPLE system for a
   range of parameters (e.g. Mach numbers), storing snapshots to
   `./ITHACAoutput/Offline/`.
2. **POD basis construction**: compute proper orthogonal decomposition modes
   for velocity, pressure and eddy viscosity using existing snapshots or
   incremental POD.
3. **Coefficient extraction**: project each snapshot onto the POD basis using
   `getCoeffs()` and save the resulting matrices; separate training and test
   sets.
4. **Neural network training**: `getTurbNN()` triggers network creation and
   training if the coefficient files are absent; normalisation parameters are
   automatically saved.
5. **Online reduced simulation**: load the trained TorchScript model and
   evaluate it at each time step to predict eddy‑viscosity coefficients.

Utility scripts (`compSteady.py`, `train.py`, `plots.py`, `plotsErr.py`)
perform data preprocessing, network training, and post‑processing.

## Files
* `03compSteadyNS.C` – main source with `CompressibleSteadyNN` definition.
* `compSteady.py` – Python helpers for feature engineering and data loading.
* `train.py` – standalone training script (PyTorch).
* `parsOff_mat.txt`, `parsOn_mat.txt` – sample parameter matrices.
* `Make/` – build directory; output binary: `03compSteadyNS`.

## Building and execution
Compile the executable:

```bash
cd tutorials/NN/CompressibleSteadyNS/Make
make -j
```

Run from the tutorial root:

```bash
./03compSteadyNS
```

Training is automatic if the `ITHACAoutput/NN/` directory lacks trained
coefficients. Use `train.py` to train offline or adjust POD mode counts
via the ITHACAdict.

## Output
Results are stored under `ITHACAoutput/NN/`, comprising:

* Normalisation files (`minAnglesInp_*.npy`, `scaleAnglesInp_*.npy`, etc.)
* Coefficient matrices (`coeffL2UTrain.npy`, `coeffL2UTest.npy`, etc.)
* Trained TorchScript model and optimiser checkpoints
* Error logs and diagnostic data (`loss*.dat`, `errors.dat`, etc.)

## The plain code
The plain code is available [here](https://raw.githubusercontent.com/ITHACA-FV/ITHACA-FV/refs/heads/master/tutorials/NN/CompressibleSteadyNS/03compSteadyNS.C).