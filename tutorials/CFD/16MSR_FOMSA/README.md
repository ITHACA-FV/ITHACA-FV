# Tutorial 16
## Introduction
This tutorial performs a sensitivity analysis using the full-order model
(FOM) of the molten salt reactor cavity. The physical setup is identical to
`15MSR_cavity`, but the focus here is on exploring parameter uncertainties
rather than constructing a reduced order model.

### Parameters
Three random parameters are considered:
	* kinematic viscosity `nu`,
	* total delayed neutron fraction `betaTot` (ratios `beta_i/betaTot` fixed),
	* third decay heat group constant `decLam3`.

Each parameter is assumed normally distributed and sampled at 1000 points.
The offline FOM is run for each sample, and statistics are collected.

### Figures of merit
The quantities monitored for each sample are:
	* total power `P_tot`,
	* mean temperature `T_m`,
both evaluated at t = 100 s. Only the final time and the initial condition
are exported since temporal evolution is not required for the analysis.

## Problem overview
The full order problem is solved by the class `msr` (see
`15MSR_cavity`), with identical geometry and boundary conditions. The
analysis uses the `LRSensitivity` utilities to compute sensitivities of the
figures of merit with respect to the parameters.

## Usage
Compile the code with the provided `Make/` target and execute the binary.
The parameter samples are read from a normally distributed matrix generated
externally (e.g. via Python) and stored in `parsOff_mat.txt` or similar.
Outputs are written to `ITHACAoutput/`, including per-sample power and
temperature values.

Results can be post-processed using scripts or standard tools to produce
statistical distributions or sensitivity plots.

---
(*Updated README with structured description and usage information.*)

