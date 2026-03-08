# Tutorial NN‑02
## Introduction
This tutorial extends the reduced‑order modelling framework to a
compressible steady Navier–Stokes problem on a bump geometry. A neural network
is trained to predict the eddy viscosity (``nut``) coefficients in the reduced
system, enabling a data‑driven turbulence closure.

The C++ class `CompressibleSteadyNN` (defined in `02compBump.C`) inherits from
`CompressibleSteadyNS` and adds network architecture & training logic similar
to `NN‑01`. The network input consists of projected velocity coefficients and
optionally Mach number and Reynolds parameter values, while the output is the
reduced eddy‑viscosity coefficients.

## Workflow
1. **Offline run**: generate snapshot fields for velocity, pressure and
   turbulent viscosity using the full-order solver. Snapshots are stored under
   `./ITHACAoutput/Offline/` for later processing.
2. **Coefficient computation**: `getTurbNN()` reads middle fields from the
   offline results, projects them onto POD modes and saves the resulting
   coefficient matrices (train/test) to
   `ITHACAoutput/NN/coeffs/`.
3. **Training**: the network is constructed (linear layers with ReLU
   activations) and trained with Adam. Normalisation parameters (bias/scale)
   are saved alongside the TorchScript model.
4. **Online evaluation**: the binary loads the trained network via
   `loadNet()` and evaluates it during reduced simulations to supply eddy
   viscosity predictions.

Additional scripts (`02compBump.py`, `compBump.py`, `train.py`) support data
management and training outside the C++ code. Data files such as
`parsOff_mat.txt`/`parsOn_mat.txt` provide sample parameters, and `plots.py`
produces diagnostic figures.

## Building and execution
Build the executable from the `Make/` directory:

```bash
cd tutorials/NN/BumpCompressibleSteadyNS/Make
make -j
```

Run `./02compBump` in the tutorial root. Training occurs automatically if
coefficients are missing, or invoke the Python training script manually.

## Output
All NN‑related outputs appear under `ITHACAoutput/NN/`, including:

* `minAnglesInp_*.npy`, `scaleAnglesInp_*.npy` – input normalisation
* `minOut_*.npy`, `scaleOut_*.npy` – output normalisation
* `coeffs/coeffL2*Train.npy`, `*Test.npy` – projection coefficients
* TorchScript model file generated during training

## The plain code
The plain code is available [here](https://raw.githubusercontent.com/ITHACA-FV/ITHACA-FV/refs/heads/master/tutorials/NN/BumpCompressibleSteadyNS/02compBump.C).