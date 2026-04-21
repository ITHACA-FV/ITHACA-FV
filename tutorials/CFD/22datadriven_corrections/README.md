# Tutorial 22

This tutorial introduces data-driven correction terms into supremizer-based or pressure-Poisson based reduced order models (SUP-ROM or PPE-ROM) for the 2D flow around a circular cylinder.

More details about the model can be found in [Ivagnes et al., JCP, 2023](https://doi.org/10.1016/j.jcp.2022.111904) and in [Ivagnes et al., AMC, 2023](https://doi.org/10.1016/j.amc.2023.127920).

# A detailed look into the code
This tutorial is divided into an offline stage (performed as usual in C++ using ITHACA-FV by running the `offline` script), and a second online stage performed using Python. The notebooks `DD-SUP-ROM_tutorial.ipynb` and `DD-PPE-ROM_tutorial.ipynb` illustrate the online stage using either a supremizer or a PPE stabilization approach.

## Offline data collection
The `offline` script should be run adding the specified option for the stabilization: `supremizer` or `poisson`.
In the `offline` script, the class `tutorial22` is inherited from `UnsteadyNSTurb`, and the relevant fields are innitialized:
```cpp
class tutorial22 : public UnsteadyNSTurb
{
    public:
        explicit tutorial22(int argc, char* argv[])
            :
            UnsteadyNSTurb(argc, argv),
            U(_U()),
            p(_p()),
            nut(_nut())
        {}

        // Relevant Fields
        volVectorField& U;
        volScalarField& p;
        volScalarField& nut;
```

In this class, as usual, the offline solve is defined, together with a tensor product, which will be useful in computations:
```cpp
        void offlineSolve(std::string offlinepath)
        {
            Vector<double> inl(0, 0, 0);
            List<scalar> mu_now(1);
            mu_now[0] =  1e-05;

            if ((offline) && (ITHACAutilities::check_folder(offlinepath) == true))
            {
                ITHACAstream::read_fields(Ufield, U, offlinepath);
                ITHACAstream::read_fields(Pfield, p, offlinepath);
                ITHACAstream::read_fields(nutFields, nut, offlinepath);
            }
            else
            {
                truthSolve(mu_now, offlinepath);
            }
        }

        Eigen::MatrixXd vectorTensorMult(Eigen::VectorXd g, Eigen::Tensor<double, 3> c,
                                         Eigen::VectorXd a)
        {
            int prodDim = c.dimension(0);
            Eigen::MatrixXd prod;
            prod.resize(prodDim, 1);

            for (int i = 0; i < prodDim; i++)
            {
                prod(i, 0) = g.transpose() *
                             Eigen::SliceFromTensor(c, 0, i) * a;
            }

            return prod;
        }
};
```

In the main function, two options can be given: either `supremizer` or `poisson`, depending on the chosen stabilization:

```cpp
int main(int argc, char* argv[])
{
    if (argc == 1)
    {
        Info << "Pass 'supremizer' or 'poisson' as first arguments."
                  << endl;
        return 0;
    }
    ...
    if (std::strcmp(argv[1], "supremizer") == 0)
    {
        // perform the offline stage, extracting the modes from the snapshots' dataset corresponding to parOffline
        supremizer_approach(example);
    }
    else if (std::strcmp(argv[1], "poisson") == 0)
    {
        poisson_approach(example);
    }
    else
    {
        Info << "Pass supremizer, poisson" << endl;
    }

    return 0;
}
```
Both approaches are then defined in different functions. In both cases, the physical parameters (viscosity, time stepping, initial, final time, number of modes, penalty factor for the BCs) are initialized or read from the `ITHACAdict` dictionary:

```cpp
    word filename("./par");
    Eigen::VectorXd par;
    example.inletIndex.resize(1, 2);
    example.inletIndex << 0, 0;
    example.inletIndexT.resize(1, 1);
    example.inletIndexT << 1;
    ITHACAparameters* para = ITHACAparameters::getInstance(example._mesh(),
                             example._runTime());
    int NmodesU = para->ITHACAdict->lookupOrDefault<int>("NmodesU", 5);
    int NmodesP = para->ITHACAdict->lookupOrDefault<int>("NmodesP", 5);
    int NmodesSUP = para->ITHACAdict->lookupOrDefault<int>("NmodesSUP", 5);
    int NmodesNUT = para->ITHACAdict->lookupOrDefault<int>("NmodesNUT", 5);
    int NmodesProject = para->ITHACAdict->lookupOrDefault<int>("NmodesProject", 5);
    int NmodesMatrixRec = para->ITHACAdict->lookupOrDefault<int>("NmodesMatrixRec",
                          5);
    double penaltyFactor =
        para->ITHACAdict->lookupOrDefault<double>("penaltyFactor", 5);
    double U_BC = para->ITHACAdict->lookupOrDefault<double>("U_BC", 0.001);
    double romStartTime = para->ITHACAdict->lookupOrDefault<double>("romStartTime",
                          0);
    double romEndTime = para->ITHACAdict->lookupOrDefault<double>("romEndTime", 3);
    double romTimeStep = para->ITHACAdict->lookupOrDefault<double>("romTimeStep",
                         0.001);
    double e = para->ITHACAdict->lookupOrDefault<double>("RBFradius", 1);
    example.startTime = 20;
    example.finalTime = 40;
    example.timeStep = 0.0002;
    example.writeEvery = 0.004;
```

The FOM simulation is then performed (together with the supremizer solve, if the supremizer approach is used):

```cpp
example.offlineSolve("./ITHACAoutput/Offline");
example.solvesupremizer();
```
The POD is performed and the modes are computed:

```cpp
    ITHACAPOD::getModes(example.Ufield, example.Umodes, example._U().name(),
                        example.podex, 0, 0, NmodesProject);
    ITHACAPOD::getModes(example.Pfield, example.Pmodes, example._p().name(),
                        example.podex, 0, 0, NmodesProject);
    ITHACAPOD::getModes(example.supfield, example.supmodes, example._U().name(),
                        example.podex,
                        example.supex, 1, NmodesProject);
```
Then, depending on the stabilization, either `projectSUP` or `projectPPE` are used to compute the reduced order operators to assemble the online system.
Finally, the projected coefficients are stored and will be used in the online stage in the Python scripts:

```cpp
    Eigen::MatrixXd coeefs = ITHACAutilities::getCoeffs(example.Ufield,
                             example.L_U_SUPmodes);
    Eigen::MatrixXd coeefsNut = ITHACAutilities::getCoeffs(example.nutFields,
                                example.nutModes);
    Eigen::MatrixXd coeefsP = ITHACAutilities::getCoeffs(example.Pfield,
                              example.Pmodes);
    cnpy::save(coeefs, "./ITHACAoutput/Matrices/coeefs.npy");
    cnpy::save(coeefsNut, "./ITHACAoutput/Matrices/coeefsNut.npy");
    cnpy::save(coeefsP, "./ITHACAoutput/Matrices/coeefsP.npy");
```
Note that in the case of supremizer approach, the supremizer coefficients are appended in the `coeefs` object.


## Online stage: SUP-ROM or PPE-ROM dynamical system

The Python files `DD_SUP_ROM` and `DD_PPE_ROM` define two classes (`SUP_ROM` and `PPE_ROM`, respectively) encapsulating the
reduced dynamical system for the velocity coefficient vector $a$ and the pressure coefficients $b$.
With a given number of modes $N_u, N_p$ (and the supremizer modes for the SUP formulation), the reduced SUP or PPE systems can be assembled.
The tutorials `DD-SUP-ROM_tutorial.ipynb` and `DD-PPE-ROM_tutorial.ipynb` show how to use the `SUP_ROM` and `PPE_ROM` classes to enhance the dynamical systems with data-driven corrections.

In order to do so, the data-driven extra terms are built using data by solving a least-squares problem, and then added into the ROM system.
The extended system including correction terms is then solved and its solution is compared to the reference data and to the standard ROM approach.

### Turbulence modelling extension

The tutorial further extends to include eddy viscosity modelling. A neural network is trained on pairs $(a,g)$ to express the reduced eddy viscosity coefficients $g$ as a function of velocity coefficients $a$. The reduced PDE system is modified with additional tensors accounting for viscosity terms, and the notebook explains the required modifications.

