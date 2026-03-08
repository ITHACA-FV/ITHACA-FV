# Tutorial 22

This tutorial introduces data-driven correction terms into pressure-Poisson
based reduced order models (PPE-ROM) for the 2D flow around a circular
cylinder. The notebook-based workflow is available in the repository and the
Python scripts illustrate the offline data preparation and online simulation
stages.

## Offline data collection

Start from a standard ITHACA-FV offline run (use command line `offline
'poisson'` in the cylinder case). Snapshots of velocity, pressure and
supremizer fields are stored in `ITHACAoutput/Offline/` and the matrices
required for assembly of the reduced systems appear in
`ITHACAoutput/Matrices`. POD modes are saved in `ITHACAoutput/POD/`.

The offline stage is handled by the full-order OpenFOAM solver; the notebooks
explain how to read and structure the data for subsequent machine learning.

## PPE-ROM dynamical system (online stage)

The Python package `DD_PPE_ROM` defines a class `PPE_ROM` encapsulating the
reduced dynamical system for the velocity coefficient vector $a$ and the
pressure coefficients $b$. With a given number of modes $N_u, N_p$ (and
no supremizer modes for the Poisson formulation), the reduced fields are
$u_r = \sum a_i \phi_i$ and $\bar p_r = \sum b_i \chi_i$.

The standard PPE-ROM system without corrections reads:

\[
    \begin{cases}
    M \dot{a}=\nu(B+B_T)a - a^T C a - H b + \tau\Big(\sum_{k} (U_{BC,k} D^k - E^k a)\Big),\\
    D b + a^T G a - \nu N a - L =0,
    \end{cases}
\]

Matrices (M, B, C, etc.) are defined in the notebook; penalization terms
handle Dirichlet boundary conditions. Time-stepping uses a second-order
backward scheme and typical online runs last 501 steps of 0.004 s each.

The notebooks demonstrate how to compute projection error metrics and compare
reconstruction percentages of the projected fields.

## Data-driven correction terms

Exact correction terms $\tau_u^{\text{exact}},\tau_p^{\text{exact}}$ are
extracted from snapshot data; the pressure correction comprises two parts
$\tau_D, \tau_G$. An ansatz of the form
$\tilde J_A ab + (ab)^T \tilde J_B ab$ is proposed and the unknown
matrices $\tilde J_A, \tilde J_B$ are determined via least-squares
minimization against the exact terms. The notebook walks through the
derivation and shows how to select the number of retained singular values
based on pressure error metrics.

The extended system including correction terms is then solved and compared to
reference data.

## Turbulence modelling extension

The tutorial further extends to include eddy viscosity modelling. A neural
network is trained on pairs $(a,g)$ to express the reduced eddy viscosity
coefficients $g$ as a function of velocity coefficients $a$. The
reduced PDE system is modified with additional tensors accounting for
viscosity terms, and the notebook explains the required modifications.

## Running the notebooks

Open `DD-PPE-ROM_tutorial.ipynb` and `DD-SUP-ROM_tutorial.ipynb` using Jupyter
or VSCode's notebook interface. They contain both explanatory text and the
commands to execute the offline and online workflows. Python scripts
`DD_PPE_ROM.py`, `DD_SUP_ROM.py`, `TrainNet.py` provide reusable classes and
training routines.

## Output

Results are saved to the `ITHACAoutput` directory and include reduced
coefficients, prediction errors, and visualizable OpenFOAM fields for
comparison.

---
*(This README synthesizes the content of the provided notebooks into a
standalone overview.)