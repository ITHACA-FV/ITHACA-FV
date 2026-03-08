# Tutorial NN-01

## Introduction
This tutorial couples a steady incompressible Navier–Stokes solver (SIMPLE
algorithm) with a small neural network that predicts the eddy viscosity
coefficient from reduced‑order (POD) velocity modes. The geometry is the
`simpleTurbGeomClosed` case and the objective is to learn a turbulence model
in a reduced‑order setting.

The workflow consists of two stages:
1. **Offline snapshot generation** using the parent class `SteadyNSSimple` to
   produce velocity, pressure and turbulent viscosity fields for a set of
   parameter values (e.g. angle of attack); these are written to
   `./ITHACAoutput/Offline/`.
2. **Neural‑network training**. The `SteadyNSSimpleNN` class (defined in
   `01simpleTurbGeomClosed.C`) computes L2 projection coefficients of the
   snapshots onto the POD modes and trains a fully‑connected network (two
   hidden layers of 128 and 64 neurons) using PyTorch. The network takes the
   reduced coefficients (and optionally a parameter value) as input and
   predicts the coefficients of the eddy viscosity field.

Once trained, the network is exported in TorchScript format and can be loaded
with `SteadyNSSimpleNN::loadNet`. The online solver then evaluates the
network to supply a turbulence closure while solving the reduced system.

## Files of interest
* `01simpleTurbGeomClosed.C` – contains the `SteadyNSSimpleNN` definition and
  main routine orchestrating offline/online phases.
* `train.py` – Python script for training the network outside the C++ code.
* `simpleTurbGeomClosed.py` – helper utilities for data preprocessing.
* `angOff_mat.txt`, `angOn_mat.txt`, `vel.txt` – example matrices and fields.
* `Make/` – build scripts; output binary named `01simpleTurbGeomClosed`.

## Building and running
Compile the code with:

```bash
cd tutorials/NN/01simpleTurbGeomClosed/Make
make -j
```

Use the generated executable to run offline solves and optionally train the
network (the code triggers training if output coefficients are absent). A
separate Python training script is also available:

```bash
python3 train.py
```

Prediction/inference occurs automatically when the solver is run with a
trained TorchScript model present; adjust `NmodesUproj` and `NmodesNutProj`
in the dictionary to change the number of modes.

Outputs (trained model, normalization factors, coefficient files) appear in
`ITHACAoutput/NN/`.

---
*(Tutorial README created to describe neural‑network closure example.)*