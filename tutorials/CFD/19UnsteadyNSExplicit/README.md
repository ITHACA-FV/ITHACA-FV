# Tutorial 19

## Introduction
This tutorial builds a reduced order model for the classic lid-driven cavity
benchmark. The two‑dimensional square cavity has side length $L=1.0$ m and a
structured 64×64 mesh. A uniform tangential velocity $U_{lid}=1.0$ m/s is
prescribed on the top wall; no-slip boundary conditions hold on the remaining
walls. The Reynolds number based on the lid velocity is 100 (laminar flow).

An explicit Forward‑Euler time discretization is used both at the full‑order and
reduced levels. The full order model (FOM) is simulated for a constant time
step $\Delta t=0.005$ s over a total time of 1 s during the offline stage.

The geometry is sketched here:
![lid-driven cavity grid](https://github.com/mathLab/ITHACA-FV/blob/master/docs/images/lidDrivenCavityGrid.png)

The discretize‑then‑project approach is adopted: the fully discrete Navier–Stokes equations are projected onto the POD basis. No supremizer modes are needed since the time discretization already enforces pressure stability.

The forcing field used in the unsteady Navier–Stokes problem is defined as
```cpp
F = 1e-2 * (ParcelOFDriver::randomField(T.mesh()));
```

# A detailed look into the code

## Header files
The following headers are required:
```cpp
#include "UnsteadyNSExplicit.H"           // explicit FOM solver
#include "ITHACAPOD.H"                   // POD decomposition
#include "ReducedUnsteadyNSExplicit.H"   // reduced model
#include "ITHACAstream.H"               // I/O utilities
```

## The `tutorial19` class
Derived from `UnsteadyNSExplicit`, the class manages the velocity, pressure and
flux fields used during offline/online computations. The offline solve is defined as usual: it performs the simulations (`truthSolve`) if the data are not found; if the data already exists, it only reads the fields which will be used for the ROM.
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

1. Parse mode counts from the `ITHACAdict` dictionary, set parameter sampling (viscosity). In this case, the parameter assumes only one value (0.01).
```cpp
    ITHACAparameters* para = ITHACAparameters::getInstance(example._mesh(),
                             example._runTime());
    int NmodesUout = para->ITHACAdict->lookupOrDefault<int>("NmodesUout", 15);
    int NmodesPout = para->ITHACAdict->lookupOrDefault<int>("NmodesPout", 15);
    int NmodesSUPout = 0;
    int NmodesUproj = para->ITHACAdict->lookupOrDefault<int>("NmodesUproj", 10);
    int NmodesPproj = para->ITHACAdict->lookupOrDefault<int>("NmodesPproj", 10);
    int NmodesSUPproj = 0;
    /// Set the number of parameters
    example.Pnumber = 1;
    /// Set samples
    example.Tnumber = 1;
    /// Set the parameters infos
    example.setParameters();
    // Set the parameter ranges
    example.mu_range(0, 0) = 0.01;
    example.mu_range(0, 1) = 0.01;
    // Generate equispaced samples inside the parameter range
    example.genEquiPar();
```
2. Configure inlet boundary indices and timing parameters.
```cpp
    example.inletIndex.resize(1, 2);
    example.inletIndex(0, 0) = 0;
    example.inletIndex(0, 1) = 0;
    example.startTime = 0.0;
    example.finalTime = 1.0;
    example.timeStep = 0.005;
    example.writeEvery = 0.005;
```
3. Run `example.offlineSolve()` to generate or load snapshots.
```cpp
example.offlineSolve();
```
4. Compute POD modes for velocity, pressure (and flux if consistent method).
```cpp
    ITHACAPOD::getModes(example.Ufield, example.Umodes, example._U().name(),
                        example.podex, 0, 0,
                        NmodesUout);
    ITHACAPOD::getModes(example.Pfield, example.Pmodes, example._p().name(),
                        example.podex, 0, 0,
                        NmodesPout);

    if (example.fluxMethod == "consistent")
    {
        ITHACAPOD::getModes(example.Phifield, example.Phimodes,  example._phi().name(),
                            example.podex, 0, 0,
                            NmodesUout);
    }
```
5. Project reduced matrices via `discretizeThenProject`.
```cpp
example.discretizeThenProject("./Matrices", NmodesUproj, NmodesPproj,
                                  NmodesSUPproj);
```
6. Build `ReducedUnsteadyNSExplicit` object and perform online simulations for
   new lid velocities. Time-stepping parameters are identical to the FOM.
```cpp
    ReducedUnsteadyNSExplicit reduced(example);
    // Set values of the reduced order model
    reduced.nu = 0.01;
    reduced.tstart = 0.0;
    reduced.finalTime = 1.0;
    reduced.dt = 0.005;
    reduced.storeEvery = 0.005;
    reduced.exportEvery = 0.005;
    // Set the online velocity
    Eigen::MatrixXd vel_now(1, 1);
    vel_now(0, 0) = 1;
    reduced.solveOnline(vel_now, 1);
```
7. Reconstruct reduced solution with :
```cpp
reduced.reconstruct(false, "./ITHACAoutput/Reconstruction/");
```

## The plain code
The plain code is available [here](https://raw.githubusercontent.com/ITHACA-FV/ITHACA-FV/refs/heads/master/tutorials/CFD/19UnsteadyNSExplicit/19UnsteadyNSExplicit.C).