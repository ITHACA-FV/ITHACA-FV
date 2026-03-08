# Tutorial 15
## Introduction
This tutorial considers a molten salt reactor (MSR) cavity represented by a
0.1 m × 0.1 m square domain discretized with 50 cells per side. The
physical problem is transient heat transfer coupled with neutronics and
fluid mechanics, modelled within the `usmsrProblem` class of ITHACA-FV. The
reduced order model is built for a parameterized set of quantities relevant
to reactor safety.

### Parametrization
Three parameters are varied:
	* kinematic viscosity `nu` of the fluid,
	* total delayed neutron fraction `betaTot` (with fixed ratios
		`beta_i/betaTot` for each group),
	* third decay heat group constant `decLam3`.

Sampling ranges are set to ±10 % around reference values (`nu0`,
`betaTot0`, `decLam30`) with a Reynolds number of 100 based on `nu0`.

### Sensitivity metrics
Two figures of merit are computed at 100 s of simulation time:
	* total power `P_tot`,
	* mean temperature `T_m`.

Temporal discretization during offline sampling uses Δt = 0.1 s, 31 runs,
150 s each, with snapshots every 5 s (configurable via the ITHACAdict file).

## The necessary header files
The code relies on the MSR problem class and standard ITHACA utilities:
```cpp
#include "usmsrProblem.H"
#include "ReducedUnsteadyMSR.H"
#include "ITHACAstream.H"
#include "LRSensitivity.H"
#include "Tm.H"
#include "Ptot.H"
```

Additional headers include Eigen for linear algebra and chrono for timing.

## Implementation of the `msr` class
The tutorial defines a child of `usmsrProblem` named `msr`. The constructor
initializes references to all relevant fields (velocity, pressure, flux,
precursors, temperature, decay constants, eddy viscosity, etc.) and reads
sensitivity parameters:
```cpp
class msr : public usmsrProblem
{
		public:
				explicit msr(int argc, char* argv[])
						:
						usmsrProblem(argc, argv),
						U(_U()),
						p(_p()),
						flux(_flux()),
						prec1(_prec1()),
						...
						TXS(_TXS())
				{}

				volVectorField& U;
				volScalarField& p;
				// ... other fields ...

				void offlineSolve()
				{
						if (offline == false)
						{
								List<scalar> mu_now(3);

								for (int i = 0; i < mu.cols(); i++)
								{
										// update viscosity, betas and decay constant
										mu_now[0] = mu(0, i);
										mu_now[1] = mu(1, i);
										mu_now[2] = mu(2, i);
										truthSolve(mu_now);
										restart();
								}
						}
						else
						{
								readMSRfields();
						}
				}
};
```

The offline solver cycles through parameter samples, updating the fluid
properties and performing a full-order run. If `offline` is true, precomputed
fields are simply read from disk.

## Main function
The `main` routine constructs the `msr` object, sets parameter ranges and
sampling (three parameters and 31 time samples by default), defines inlet
patch indices and simulation time settings, and then runs the offline stage.
The code also defines auxiliary classes (`Tmlocal`, `Ploct`, etc.) for
sensitivity computations of mean temperature and total power.

After offline sampling the reduced basis can be computed and the sensitivity
analysis performed using the `LRSensitivity` utilities provided by the
library.

Further details such as how to run the script, location of input files
(`parsOff_mat.txt`, etc.), and post-processing plots are available in the
original tutorial directory.

## Usage
Compile with the provided `Make/` target and run the executable. The
parameter files are read from `parsOff_mat.txt` and outputs are stored in
`ITHACAoutput/` directories. For plotting helper Python scripts
(`plotCentreMass.py`, `plotDrag.py`, `plotLift.py`) can be used.

---
(*This README replaces the previous short summary with detailed structure
and steps.*)
