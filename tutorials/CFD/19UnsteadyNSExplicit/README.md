# Tutorial 19

### Introduction
This tutorial builds a reduced order model for the classic lid-driven cavity
benchmark. The two‑dimensional square cavity has side length $L=1.0$ m and a
structured 64×64 mesh. A uniform tangential velocity $U_{lid}=1.0$ m/s is
prescribed on the top wall; no-slip boundary conditions hold on the remaining
walls. The Reynolds number based on the lid velocity is 100 (laminar flow).

An explicit Forward‑Euler time discretization is used both at the full‑order and
reduced levels. The full order model (FOM) is simulated for a constant time
step $\Delta t=0.005$ s over a total time of 1 s during the offline stage.

### Geometry sketch
![lid-driven cavity grid](https://github.com/mathLab/ITHACA-FV/blob/master/docs/images/lidDrivenCavityGrid.png)

### Problem formulation
The discretize‑then‑project approach is adopted: the fully discrete Navier–Stokes
equations are projected onto the POD basis. No supremizer modes are needed since
the time discretization already enforces pressure stability.

#### Forcing term
The forcing field used in the unsteady Navier–Stokes problem is defined as
```cpp
F = 1e-2 * (ParcelOFDriver::randomField(T.mesh()));
```

### Header files
The following headers are required:
```cpp
#include "UnsteadyNSExplicit.H"           // explicit FOM solver
#include "ITHACAPOD.H"                   // POD decomposition
#include "ReducedUnsteadyNSExplicit.H"   // reduced model
#include "ITHACAstream.H"               // I/O utilities
```

### `tutorial19` class
Derived from `UnsteadyNSExplicit`, the class manages the velocity, pressure and
flux fields used during offline/online computations.
```cpp
class tutorial19: public UnsteadyNSExplicit
{
    public:
        explicit tutorial19(int argc, char* argv[])
            :
            UnsteadyNSExplicit(argc, argv),
            U(_U()),
            p(_p()),
            phi(_phi())
        {}

        volVectorField& U;
        volScalarField& p;
        surfaceScalarField& phi;

        void offlineSolve()
        {
            Vector<double> inl(1, 0, 0);
            List<scalar> mu_now(1);

            if (offline)
            {
                ITHACAstream::read_fields(Ufield, U, "./ITHACAoutput/Offline/");
                ITHACAstream::read_fields(Pfield, p, "./ITHACAoutput/Offline/");
                if (fluxMethod == "consistent")
                {
                    ITHACAstream::read_fields(Phifield, phi, "./ITHACAoutput/Offline/");
                }
            }
            else
            {
                for (label i = 0; i < mu.cols(); i++)
                {
                    mu_now[0] = mu(0, i);
                    truthSolve(mu_now);
                }
            }
        }
};
```

### Main function overview
The main routine configures parameters, executes offline snapshot generation,
and constructs the reduced system. Key steps:

1. Parse mode counts from the dictionary, set parameter sampling (viscosity).
2. Configure inlet boundary indices and timing parameters.
3. Run `example.offlineSolve()` to generate or load snapshots.
4. Compute POD modes for velocity, pressure (and flux if consistent method).
5. Project reduced matrices via `discretizeThenProject`.
6. Build `ReducedUnsteadyNSExplicit` object and perform online simulations for
   new lid velocities. Time-stepping parameters are identical to the FOM.
7. Reconstruct reduced solution with `reduced.reconstruct(false, ...)`.

Example code snippets appear in the original README and can be copied as
needed.

### Usage
Compile with the provided `Make/` target. Configure the problem via the
`system/ITHACAdict` file (e.g. `fluxMethod` and mode counts). Run the binary to
produce `ITHACAoutput/` containing snapshots, POD modes and online results.

## The plain code
The plain code is available [here](https://raw.githubusercontent.com/ITHACA-FV/ITHACA-FV/refs/heads/master/tutorials/CFD/19UnsteadyNSExplicit/19UnsteadyNSExplicit.C).