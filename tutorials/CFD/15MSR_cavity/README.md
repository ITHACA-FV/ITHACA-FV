# Tutorial 15
## Introduction
This tutorial considers a molten salt reactor (MSR) cavity represented by a 0.1 m × 0.1 m square domain discretized with 50 cells per side. The physical problem is transient heat transfer coupled with neutronics and fluid mechanics, modelled within the `usmsrProblem` class of ITHACA-FV. The reduced order model is built for a parameterized set of quantities relevant to reactor safety.

### Parametrization
Three parameters are varied:

- kinematic viscosity `nu` of the fluid,
- total delayed neutron fraction `betaTot` (with fixed ratios `beta_i/betaTot` for each group),
- third decay heat group constant `decLam3`.

Sampling ranges are set to ±10 % around reference values (`nu0`, `betaTot0`, `decLam30`) with a Reynolds number of 100 based on `nu0`.

### Sensitivity metrics
Two figures of merit are computed at 100 s of simulation time:

- total power `P_tot`,
- mean temperature `T_m`.

Temporal discretization during offline sampling uses $\Delta t = 0.1$ s, 31 runs,
150 s each, with snapshots every 5 s (which are configurable via the `system/ITHACAdict` file).

# A detailed look into the code
Let's now have a look into the code of tutorial 15.

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

Additional headers include `Eigen` for linear algebra and `chrono` for timing.

## Implementation of the `msr` class
The tutorial defines a child of `usmsrProblem` named `msr`. The constructor initializes references to all relevant fields (velocity, pressure, flux,
precursors, temperature, decay constants, eddy viscosity, etc.) and reads sensitivity parameters:
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

The offline solver cycles through parameter samples, updating the fluid properties and performing a full-order run. If `offline` is true, precomputed
fields are simply read from disk.

The code also defines auxiliary classes (`Tmlocal`, `Ploct`, etc.) for sensitivity computations of mean temperature and total power.
In particular, `Tmlocal` is for the domain-averaged temperature:
```cpp
class Tmlocal : public Tm
{
    public:
        explicit Tmlocal(int argc, char* argv[], int Nsampled)
            :
            Tm(argc, argv, Nsampled),
            T(_T())
        {}

        volScalarField& T;
};
```
`Ploct` computes the time-dependent counterpart of the total power, inherited from `Ptot_time`. It stores a reference to the power density field:
```cpp
class Ploct: public Ptot_time
{
    public:
        explicit Ploct(int argc, char* argv[], int Nsampled)
            :
            Ptot_time(argc, argv, Nsampled),
            powerDens(_powerDens())
        {}

        volScalarField& powerDens;
};
```
`Plocal` is similar to `Ploct`, but likely for the non-time-indexed version of the total power figure of merit:
```cpp
class Plocal: public Ptot
{
    public:
        explicit Plocal(int argc, char* argv[], int Nsampled)
            :
            Ptot(argc, argv, Nsampled),
            powerDens(_powerDens())
        {}

        volScalarField& powerDens;
};
```
Finally, `Tmloct` is the time-dependent version of the mean-temperature output object.

```cpp
class Tmloct: public Tm_time
{
    public:
        explicit Tmloct(int argc, char* argv[], int Nsampled)
            :
            Tm_time(argc, argv, Nsampled),
            T(_T())
        {}

        volScalarField& T;
};
```


## Main function
The `main` routine constructs the `msr` object, sets parameter ranges and sampling (three parameters and 31 time samples by default), defines inlet patch indices and the simulation time settings:
```cpp
	//object initialization
    msr prova(argc, argv);
    //inletIndex set
    prova.inletIndex.resize(1, 2);
    prova.inletIndex(0, 0) = 0;
    prova.inletIndex(0, 1) = 0;
    prova.inletIndexT.resize(4, 1);
    prova.inletIndexT(0, 0) = 0;
    prova.inletIndexT(1, 0) = 1;
    prova.inletIndexT(2, 0) = 2;
    prova.inletIndexT(3, 0) = 3;
	//set the number of parameters and values of each
    prova.Pnumber = 3;
    prova.Tnumber = 31;
    prova.setParameters();
    int central = static_cast<int>(prova.Tnumber / 2);
    double nu0 = 2.46E-06;
    double betatot0 = 3.218E-05;
    double dlam30 = 3.58e-04;
    //+-10% from default values
    prova.mu_range(0, 0) = 1.1 * nu0; //nu range
    prova.mu_range(0, 1) = 0.9 * nu0;
    prova.mu_range(1, 0) = 1.1 * betatot0; //betaTot range
    prova.mu_range(1, 1) = 0.9 * betatot0;
    prova.mu_range(2, 0) = 1.1 * dlam30; //decLam3 range
    prova.mu_range(2, 1) = 0.9 * dlam30;
```
The parameters are sampled in a random way, and a central reference point is forced:
```cpp
    prova.genRandPar();
	//override central column of mu with reference values of the parameters
    prova.mu(0, central) = nu0;
    prova.mu(1, central) = betatot0;
    prova.mu(2, central) = dlam30;
    double tstart = std::time(0);
    Eigen::MatrixXd arange(prova.Pnumber, 2);
    Eigen::MatrixXd mu_f = prova.mu;
```
It then runs the offline stage and collects information on the elapsed time:
```cpp
if (prova.offline == false)
    {
        for (int i = 0; i < prova.Pnumber; i++)
        {
            arange(i, 0) = prova.mu.row(i).minCoeff();
            arange(i, 1) = prova.mu.row(i).maxCoeff();
        }

        ITHACAstream::exportMatrix(arange, "range", "eigen", "./ITHACAoutput");
        ITHACAstream::exportMatrix(prova.mu, "mu_off", "eigen", "./ITHACAoutput");
    }
    else
    {
        mu_f = ITHACAstream::readMatrix("./ITHACAoutput/mu_off_mat.txt");
        arange = ITHACAstream::readMatrix("./ITHACAoutput/range_mat.txt");
    }

    //perform the offline step
    prova.offlineSolve();
    prova._nu().value() = nu0;
    prova.change_viscosity(nu0);
    double tend = std::time(0);
    double deltat_offline;
    deltat_offline = difftime(tend, tstart);

    //save the time needed to perform offline stage
    if (prova.offline == false)
    {
        std::ofstream diff_t("./ITHACAoutput/t_offline=" + name(deltat_offline));
    }

    Info << "Offline stage time= " << deltat_offline << endl;
```


After offline sampling, the lifting function for the velocity and the pressure is computed, then the reduced basis can be computed and the reduced PPE operators are extrated:
```cpp
	//compute the liftfunction for the velocity
    prova.liftSolve();
    //compute the liftfuncion for the temperature
    prova.liftSolveT();
    //homogenize the velocity field
    prova.homogenizeU();
    //homogenize the temperature field
    prova.homogenizeT();
    //compute the spatial modes
    prova.msrgetModesEVD();
    // perform the projection
    Eigen::VectorXi NPrec(8);
    NPrec(0) = prova.NPrecproj1;
    NPrec(1) = prova.NPrecproj2;
    NPrec(2) = prova.NPrecproj3;
    NPrec(3) = prova.NPrecproj4;
    NPrec(4) = prova.NPrecproj5;
    NPrec(5) = prova.NPrecproj6;
    NPrec(6) = prova.NPrecproj7;
    NPrec(7) = prova.NPrecproj8;
    Eigen::VectorXi NDec(3);
    NDec(0) = prova.NDecproj1;
    NDec(1) = prova.NDecproj2;
    NDec(2) = prova.NDecproj3;
    prova.projectPPE("./Matrices", prova.NUproj, prova.NPproj, prova.NFluxproj,
                     NPrec, prova.NTproj, NDec, prova.NCproj);
```
Then, the BC conditions for the velocity and temperature are set, the ROM object is initialized: 
```cpp
    //define the moving wall bc condition for U
    double umedio = 2.46E-03;
    //define the temperature at the moving wall
    Eigen::MatrixXd Twall(4, 2);
    Twall.setZero();
    Twall(0, 0) = 850;
    Twall(1, 0) = 950;
    Twall(2, 0) = 950;
    Twall(3, 0) = 950;
    //define the ROM object and its time settings
    reducedusMSR ridotto(prova);
    ridotto.tstart = 0;
    ridotto.finalTime = 10;
    ridotto.dt = 0.1;
```
Now the velocity and temperature inputs are initialized and passed online to the reduce order solver:
```cpp
    //vel_now for the ROM object
    Eigen::MatrixXd vel_now(1, 1);
    vel_now(0, 0) = umedio;
    //temp_now for the ROM object
    Eigen::MatrixXd temp_now = Twall;
```

The online ROM parameters using the central reference are now assigned to the ROM:
```cpp
    // set the online value of mu equal to mu(:,central)
    Eigen::VectorXd mu_on(3);
    ridotto.nu = mu_f(0, central);
    ridotto.btot = mu_f(1, central);
    ridotto.b1 = prova.r1 * ridotto.btot;
    ridotto.b2 = prova.r2 * ridotto.btot;
    ridotto.b3 = prova.r3 * ridotto.btot;
    ridotto.b4 = prova.r4 * ridotto.btot;
    ridotto.b5 = prova.r5 * ridotto.btot;
    ridotto.b6 = prova.r6 * ridotto.btot;
    ridotto.b7 = prova.r7 * ridotto.btot;
    ridotto.b8 = prova.r8 * ridotto.btot;
    ridotto.dl3 = mu_f(2, central);
    // index corresponding at central solution
    int ref_start = 465;
    mu_on(0) = ridotto.nu;
    mu_on(1) = ridotto.btot;
    mu_on(2) = ridotto.dl3;
```

The ROM is then solved, and the solution is reconstructed and exported.
```cpp
    ridotto.solveOnline(vel_now, temp_now, mu_on, ref_start);
    ridotto.reconstructAP("./ITHACAoutput/Reconstruction/", 10);
    //ridotto.reconstructAP("./ITHACAoutput/Reconstruction/",50);
    //final snapshots for ROM and FOM respectively
    int tf = 21;
    int tf_fom = 31;
    //Some post process calculations
    int c = 0;
```

Some figures of merit are computed for the FOM database (average temperature and total power):
```cpp
    Eigen::MatrixXd err(10, tf);
    Eigen::MatrixXd TmFOM(prova.Tnumber, tf);
    Eigen::MatrixXd PtotFOM(prova.Tnumber, tf);
    Eigen::MatrixXd TmROM(prova.Tnumber, tf);
    Eigen::MatrixXd PtotROM(prova.Tnumber, tf);
```
The errors are then computed:
```cpp
    Info << "Computing errors on reference case..." << endl;
	for (int i = ref_start; i < ref_start + tf; i++)
    {
        err(0, c) = ITHACAutilities::errorL2Rel(prova.Ufield[i], ridotto.UREC[c]);
        err(1, c) = ITHACAutilities::errorL2Rel(prova.Fluxfield[i],
                                                ridotto.FLUXREC[c]);
        err(2, c) = ITHACAutilities::errorL2Rel(prova.Prec1field[i],
                                                ridotto.PREC1REC[c]);
        err(3, c) = ITHACAutilities::errorL2Rel(prova.Prec4field[i],
                                                ridotto.PREC4REC[c]);
        err(4, c) = ITHACAutilities::errorL2Rel(prova.Prec8field[i],
                                                ridotto.PREC8REC[c]);
        err(5, c) = ITHACAutilities::errorL2Rel(prova.Tfield[i], ridotto.TREC[c]);
        err(6, c) = ITHACAutilities::errorL2Rel(prova.Dec1field[i],
                                                ridotto.DEC1REC[c]);
        err(7, c) = ITHACAutilities::errorL2Rel(prova.Dec3field[i],
                                                ridotto.DEC3REC[c]);
        err(8, c) = ITHACAutilities::errorL2Rel(prova.PowerDensfield[i],
                                                ridotto.POWERDENSREC[c]);
        err(9, c) = ITHACAutilities::errorL2Rel(prova.TXSFields[i],
                                                ridotto.TXSREC[c]);
        c++;
    }
```

The same figures of merit are also computed for the ROM are then collected using the above-defined functions:
```cpp
	//compute the figures of merit for every "snapshot" time
    Tmloct Tmedia(argc, argv, prova.Tnumber);
    Ploct Ptotale(argc, argv, prova.Tnumber);
    c = 0;

    for (int i = 0; i < tf; i++)
    {
        Tmedia.buildMO(folder, i);
        Ptotale.buildMO(folder, i);
        TmROM.col(c) = Tmedia.modelOutput;
        PtotROM.col(c) = Ptotale.modelOutput;
        c++;
    }

    ITHACAstream::exportMatrix(TmROM, "rom_resT", "eigen",
                               "./ITHACAoutput");
    ITHACAstream::exportMatrix(PtotROM, "rom_resP", "eigen",
                               "./ITHACAoutput");
    
```

Then, the sensitivity analysis performed using the `LRSensitivity` utilities provided by the library. Here, the sampling points are defined and the ROM is re-initialized:
```cpp
    //define the folder where save the output
    folder = {"./ITHACAoutput/Sensitivity/ModelOutput/"};
    int Nparameters = prova.Pnumber;
    std::vector<std::string> pdflist = {"normal", "normal", "normal"};
    int samplingPoints = 1000;
    LRSensitivity analisi(Nparameters, samplingPoints);
    //define the training range, values sampled outside are rejected
    analisi.trainingRange = arange;
    //set the parameters of each distributions
    Eigen::MatrixXd Mp(Nparameters, 2);
    Mp(0, 0) = nu0;
    Mp(0, 1) = nu0 * 0.1 / 3;
    Mp(1, 0) = betatot0;
    Mp(1, 1) = betatot0 * 0.1 / 3;
    Mp(2, 0) = dlam30;
    Mp(2, 1) = dlam30 * 0.1 / 3;
    //build and export the sampling sets
    analisi.buildSamplingSet(pdflist, Mp);
    ITHACAstream::exportMatrix(analisi.MatX, "sampled", "eigen", "./ITHACAoutput");
    //compute parameters statistics
    analisi.getXstats();
    //perform the simulations with the reduced object
    int counter = 0;
    tstart = std::time(0);
    ridotto.clearFields();
```
Then, the ROM is solved in a loop over all sampled parameters for sensitivity purpose, and the solution is reconstructed:
```cpp
	for (int j = 0; j < samplingPoints; j++)
    {
        ridotto.nu = analisi.MatX(j, 0);
        ridotto.btot = analisi.MatX(j, 1);
        ridotto.b1 = prova.r1 * ridotto.btot;
        ridotto.b2 = prova.r2 * ridotto.btot;
        ridotto.b3 = prova.r3 * ridotto.btot;
        ridotto.b4 = prova.r4 * ridotto.btot;
        ridotto.b5 = prova.r5 * ridotto.btot;
        ridotto.b6 = prova.r6 * ridotto.btot;
        ridotto.b7 = prova.r7 * ridotto.btot;
        ridotto.b8 = prova.r8 * ridotto.btot;
        ridotto.dl3 = analisi.MatX(j, 2);
        mu_on(0) = ridotto.nu;
        mu_on(1) = ridotto.btot;
        mu_on(2) = ridotto.dl3;
        folder.append(std::to_string(counter));
        ridotto.solveOnline(vel_now, temp_now, mu_on, ref_start);
        ridotto.reconstructAP(folder, 1000);
        folder = {"./ITHACAoutput/Sensitivity/ModelOutput/"};
        ridotto.clearFields();
        counter++;
    }
```
Finally, the figures of merit, timings info, and LRSensitivity analysis objects are collected:
```cpp
	tend = std::time(0);
    int Ntot = samplingPoints;
    double deltat_online = difftime(tend, tstart);
    // save SA online time
    std::ofstream diff_ton("./ITHACAoutput/t_onlineSA_" + name(Ntot) + "=" + name(
                               deltat_online));
    //initialize the figure of merit objects
    Tmlocal TmediaSA(argc, argv, Ntot);
    Plocal PtotSA(argc, argv, Ntot);
    //build the Model Output
    TmediaSA.buildMO(folder);
    PtotSA.buildMO(folder);
    //assign FofM to LRSensitivity object, first figure of merit considered is the average temperature
    analisi.M = autoPtr<FofM>(& TmediaSA);
    //load the model output in the LRSensitivity object
    analisi.load_output();
    //compute output statistics
    analisi.getYstat();
    Info << "Mean value and variance of the output:" << endl;
    Info << analisi.Ey << "\t" << analisi.Vy << endl;
    Info << "-------------" << endl;
    // compute linear regression coefficients
    analisi.getBetas();
    Info << "analisi regression coefficients:" << endl;
    Info << analisi.betas << endl;
    Info << "-------------" << endl;
    // assess the quality of the linear regression, i.e. R^2
    analisi.assessQuality();
    Info << "quality: " << analisi.QI << endl;
    // save the ROM output, linear output, statistics
    Eigen::MatrixXd saveyout(samplingPoints, 1);
    Eigen::MatrixXd saveymodel(samplingPoints, 1);
    saveyout.col(0) = analisi.y;
    saveymodel.col(0) = analisi.ylin;
    ITHACAstream::exportMatrix(saveyout, "T_out", "eigen", "./ITHACAoutput");
    ITHACAstream::exportMatrix(saveymodel, "T_model", "eigen", "./ITHACAoutput");
    Eigen::MatrixXd coeffsT(4, 3);
    coeffsT.setZero();
    coeffsT(0, 0) = analisi.Ey;
    coeffsT(0, 1) = analisi.Vy;
    coeffsT(0, 2) = analisi.QI;
    coeffsT.row(1) = analisi.EX;
    coeffsT.row(2) = analisi.VX;
    coeffsT.row(3) = analisi.betas;
    ITHACAstream::exportMatrix(coeffsT, "coeffsT", "eigen", "./ITHACAoutput");
    //assign FofM object to LRSensitivity object
    analisi.M = autoPtr<FofM>(& PtotSA);
    //load the model output in the LRSensitivity object, repeat all the previous step
    // for the total power
    analisi.load_output();
    analisi.getYstat();
    Info << "Mean value and variance of the output:" << endl;
    Info << analisi.Ey << "\t" << analisi.Vy << endl;
    Info << "-------------" << endl;
    analisi.getBetas();
    Info << "analisi regression coefficients:" << endl;
    Info << analisi.betas << endl;
    Info << "-------------" << endl;
    analisi.assessQuality();
    Info << "quality: " << analisi.QI << endl;
    Eigen::MatrixXd saveyoutP(samplingPoints, 1);
    Eigen::MatrixXd saveymodelP(samplingPoints, 1);
    saveyoutP.col(0) = analisi.y;
    saveymodelP.col(0) = analisi.ylin;
    ITHACAstream::exportMatrix(saveyoutP, "P_out", "eigen", "./ITHACAoutput");
    ITHACAstream::exportMatrix(saveymodelP, "P_model", "eigen", "./ITHACAoutput");
    Eigen::MatrixXd coeffsP(4, 3);
    coeffsP.setZero();
    coeffsP(0, 0) = analisi.Ey;
    coeffsP(0, 1) = analisi.Vy;
    coeffsP(0, 2) = analisi.QI;
    coeffsP.row(1) = analisi.EX;
    coeffsP.row(2) = analisi.VX;
    coeffsP.row(3) = analisi.betas;
    ITHACAstream::exportMatrix(coeffsP, "coeffsP", "eigen", "./ITHACAoutput");
    return 0;
```

## The plain code
The plain code is available [here](https://raw.githubusercontent.com/ITHACA-FV/ITHACA-FV/refs/heads/master/tutorials/CFD/15MSR_cavity/15MSR_cavity.C).
