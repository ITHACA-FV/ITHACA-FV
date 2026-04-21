# Tutorial NN‑03
## Introduction
This tutorial demonstrates data‑driven turbulence modelling in a compressible steady Navier–Stokes problem using a neural network. The problem geometry is an open domain with a flow field parameterized by Mach number and other parameters. The network learns to map from reduced velocity coefficients (and parameters) to eddy‑viscosity coefficients, enabling a ROM with
inexpensive predictions.

## The script 03compSteadyNS.C

The solver class `CompressibleSteadyNN` (in `03compSteadyNS.C`) follows the
same design as `NN‑02`: offline snapshot generation, POD mode extraction, coefficient computation and projection, neural‑network training, and online evaluation.

The class `CompressibleSteadyNN` extends `CompressibleSteadyNS` by adding the neural-network interface.

```cpp
class CompressibleSteadyNN : public CompressibleSteadyNS
{
    public:
        /// Constructor
        CompressibleSteadyNN(int argc, char* argv[])
```

It stores the numbers of projected velocity and nut modes, the normalization matrices used for scaling network inputs and outputs, and the loaded `TorchScript` model.

```cpp
        {
            ITHACAparameters* para = ITHACAparameters::getInstance();
            NUmodes = para->ITHACAdict->lookupOrDefault<label>("NmodesUproj", 10);
            NNutModes = para->ITHACAdict->lookupOrDefault<label>("NmodesNutProj", 10);
        };

        label NUmodes;
        label NNutModes;

        Eigen::MatrixXd bias_inp;
        Eigen::MatrixXd scale_inp;
        Eigen::MatrixXd bias_out;
        Eigen::MatrixXd scale_out;

        Eigen::MatrixXd coeffL2Nut;
        Eigen::MatrixXd coeffL2U;

        torch::Tensor coeffL2U_tensor;
        torch::Tensor coeffL2Nut_tensor;

        torch::nn::Sequential Net;
        torch::optim::Optimizer* optimizer;
        torch::jit::script::Module netTorchscript;
```

The method `loadNet` reads the trained network and its normalization data from disk:
```cpp
        void loadNet(word filename)
        {
            std::string Msg = filename +
                              " is not existing, please run the training stage of the net with the correct number of modes for U and Nut";
            M_Assert(ITHACAutilities::check_file(filename), Msg.c_str());
            netTorchscript = torch::jit::load(filename);
            cnpy::load(bias_inp, "ITHACAoutput/NN/minAnglesInp_" + name(
                           NUmodes) + "_" + name(NNutModes) + ".npy");
            cnpy::load(scale_inp, "ITHACAoutput/NN/scaleAnglesInp_" + name(
                           NUmodes) + "_" + name(NNutModes) + ".npy");
            cnpy::load(bias_out, "ITHACAoutput/NN/minOut_" + name(NUmodes) + "_" + name(
                           NNutModes) + ".npy");
            cnpy::load(scale_out, "ITHACAoutput/NN/scaleOut_" + name(NUmodes) + "_" + name(
                           NNutModes) + ".npy");
            netTorchscript.eval();
        }
```

`getTurbNN` builds the reduced training and test datasets by projecting offline and validation snapshots onto the POD bases and exporting the resulting coefficients as NumPy arrays.
```cpp
        // This function computes the coefficients which are later used for training
        void getTurbNN()
        {
            if (!ITHACAutilities::check_folder("ITHACAoutput/NN/coeffs"))
            {
                mkDir("ITHACAoutput/NN/coeffs");
                // Read Fields for Train
                PtrList<volVectorField> UfieldTrain;
                PtrList<volScalarField> PfieldTrain;
                PtrList<volScalarField> nutFieldsTrain;
                ITHACAstream::readMiddleFields(UfieldTrain, _U(),
                                               "./ITHACAoutput/Offline/");
                ITHACAstream::readMiddleFields(PfieldTrain, _p(),
                                               "./ITHACAoutput/Offline/");
                auto nut_train = _mesh().lookupObject<volScalarField>("nut");
                ITHACAstream::readMiddleFields(nutFieldsTrain, nut_train,
                                               "./ITHACAoutput/Offline/");
                /// Compute the coefficients for train
                Info << "Computing the coefficients for U train" << endl;
                Eigen::MatrixXd coeffL2U_train = ITHACAutilities::getCoeffs(UfieldTrain,
                                                 Umodes,
                                                 0, true);
                Info << "Computing the coefficients for p train" << endl;
                Eigen::MatrixXd coeffL2P_train = ITHACAutilities::getCoeffs(PfieldTrain,
                                                 Pmodes,
                                                 0, true);
                Info << "Computing the coefficients for nuT train" << endl;
                Eigen::MatrixXd coeffL2Nut_train = ITHACAutilities::getCoeffs(nutFieldsTrain,
                                                   nutModes,
                                                   0, true);
                coeffL2U_train.transposeInPlace();
                coeffL2P_train.transposeInPlace();
                coeffL2Nut_train.transposeInPlace();
                cnpy::save(coeffL2U_train, "ITHACAoutput/NN/coeffs/coeffL2UTrain.npy");
                cnpy::save(coeffL2P_train, "ITHACAoutput/NN/coeffs/coeffL2PTrain.npy");
                cnpy::save(coeffL2Nut_train, "ITHACAoutput/NN/coeffs/coeffL2NutTrain.npy");
                // Read Fields for Test
                PtrList<volVectorField> UfieldTest;
                PtrList<volScalarField> PfieldTest;
                PtrList<volScalarField> nutFieldsTest;
                /// Compute the coefficients for test
                ITHACAstream::readMiddleFields(UfieldTest, _U(),
                                               "./ITHACAoutput/checkOff/");
                ITHACAstream::readMiddleFields(PfieldTest, _p(),
                                               "./ITHACAoutput/checkOff/");
                auto nut_test = _mesh().lookupObject<volScalarField>("nut");
                ITHACAstream::readMiddleFields(nutFieldsTest, nut_test,
                                               "./ITHACAoutput/checkOff/");
                // Compute the coefficients for test
                Eigen::MatrixXd coeffL2U_test = ITHACAutilities::getCoeffs(UfieldTest,
                                                Umodes,
                                                0, true);
                Eigen::MatrixXd coeffL2P_test = ITHACAutilities::getCoeffs(PfieldTest,
                                                Pmodes,
                                                0, true);
                Eigen::MatrixXd coeffL2Nut_test = ITHACAutilities::getCoeffs(nutFieldsTest,
                                                  nutModes,
                                                  0, true);
                coeffL2U_test.transposeInPlace();
                coeffL2P_test.transposeInPlace();
                coeffL2Nut_test.transposeInPlace();
                cnpy::save(coeffL2U_test, "ITHACAoutput/NN/coeffs/coeffL2UTest.npy");
                cnpy::save(coeffL2P_test, "ITHACAoutput/NN/coeffs/coeffL2PTest.npy");
                cnpy::save(coeffL2Nut_test, "ITHACAoutput/NN/coeffs/coeffL2NutTest.npy");
            }
        }
```

The method `evalNet` assembles the parameter vector and reduced velocity coefficients into a single input, normalizes it, performs inference, and returns the predicted reduced nut coefficients.
```cpp
        Eigen::MatrixXd evalNet(Eigen::MatrixXd a, Eigen::MatrixXd mu_now)
        {
            Eigen::MatrixXd xpred(a.rows() + mu_now.rows(), 1);

            if (xpred.cols() != 1)
            {
                throw std::runtime_error("Predictions for more than one sample not supported yet.");
            }

            xpred.bottomRows(a.rows()) = a;
            xpred.topRows(mu_now.rows()) = mu_now;
            xpred = xpred.array() * scale_inp.array() + bias_inp.array() ;
            xpred.transposeInPlace();
            torch::Tensor xTensor = eigenMatrix2torchTensor(xpred);
            torch::Tensor out;
            std::vector<torch::jit::IValue> XTensorInp;
            XTensorInp.push_back(xTensor);
            out = netTorchscript.forward(XTensorInp).toTensor();
            Eigen::MatrixXd g = torchTensor2eigenMatrix<double>(out);
            g.transposeInPlace();
            g = (g.array() - bias_out.array()) / scale_out.array();
            return g;
        }
```

The class `ReducedCompressibleSteadyNN` extends `ReducedCompressibleSteadyNS` and implements the online reduced solver.

```cpp
class ReducedCompressibleSteadyNN : public ReducedCompressibleSteadyNS
{
    public:
        /// Constructor
        explicit ReducedCompressibleSteadyNN(CompressibleSteadyNN& FOMproblem)
            :
            problem(&FOMproblem),
            ReducedCompressibleSteadyNS(FOMproblem)
        {}

        /// Full problem.
        CompressibleSteadyNN* problem;
```

The method `projectReducedOperators` precomputes the reduced pressure-gradient contribution:
```cpp
        void projectReducedOperators(int NmodesUproj, int NmodesPproj, int NmodesEproj)
        {
            PtrList<volVectorField> gradModP;

            for (label i = 0; i < NmodesPproj; i++)
            {
                gradModP.append(fvc::grad(problem->Pmodes[i]));
            }

            projGradModP = problem->Umodes.project(gradModP,
                                                   NmodesUproj); // Modes without lifting
        }
```

The main routine, `solveOnlineCompressible`, initializes the reduced coefficients for velocity, pressure, and energy from the current fields:
```cpp
        void solveOnlineCompressible(int NmodesUproj, int NmodesPproj, int NmodesEproj,
                                     int NmodesNutProj, Eigen::MatrixXd mu_now,
                                     word Folder = "./ITHACAoutput/Online/")
        {
            counter++;
            // Residuals initialization
            scalar residualNorm(1);
            scalar residualJump(1);
            Eigen::MatrixXd uResidualOld = Eigen::MatrixXd::Zero(1, NmodesUproj);
            Eigen::MatrixXd eResidualOld = Eigen::MatrixXd::Zero(1, NmodesEproj);
            Eigen::MatrixXd pResidualOld = Eigen::MatrixXd::Zero(1, NmodesPproj);
            Eigen::VectorXd uResidual(Eigen::Map<Eigen::VectorXd>(uResidualOld.data(),
                                      NmodesUproj));
            Eigen::VectorXd eResidual(Eigen::Map<Eigen::VectorXd>(eResidualOld.data(),
                                      NmodesEproj));
            Eigen::VectorXd pResidual(Eigen::Map<Eigen::VectorXd>(pResidualOld.data(),
                                      NmodesPproj));
            ...
            Eigen::MatrixXd e = ITHACAutilities::getCoeffs(E, problem->Emodes, NmodesEproj,
                                true);
            Eigen::MatrixXd u = ITHACAutilities::getCoeffs(U, problem->Umodes, NmodesUproj,
                                true);
            Eigen::MatrixXd p = ITHACAutilities::getCoeffs(P, problem->Pmodes, NmodesPproj,
                                true);
```

It predicts the turbulent-viscosity coefficients with the neural network, reconstructs the `nut` field:
```cpp
            Eigen::MatrixXd  nutCoeff = problem->evalNet(u, mu_now);
            problem->nutModes.reconstruct(nut, nutCoeff, "nut");
```

It then enters a SIMPLE-like reduced iteration loop. In each iteration, the momentum, energy, and pressure equations are assembled in full-order form, projected onto the reduced bases, solved in reduced coordinates, and reconstructed back into the full fields:

```cpp
            while ((residualJump > residualJumpLim
                    || residualNorm > normalizedResidualLim) && csolve < maxIter)
            {
                csolve++;
                Info << "csolve:" << csolve << endl;
#if OFVER == 6
                problem->_simple().loop(runTime);
#else
                problem->_simple().loop();
#endif
                uResidualOld = uResidual;
                eResidualOld = eResidual;
                pResidualOld = pResidual;
                //Momentum equation phase
                List<Eigen::MatrixXd> RedLinSysU;
                ITHACAutilities::assignBC(U, idInl, Uinlet);
                /// nuEff calculation.
                nueff = nut + problem->turbulence->nu();
                fvVectorMatrix UEqnR
                (
                    fvm::div(phi, U)
                    - fvc::div((rho * nueff) * dev2(T(fvc::grad(U))))
                    - fvm::laplacian(rho * nueff, U)
                    ==
                    fvOptions(rho, U)
                );
                UEqnR.relax();
                fvOptions.constrain(UEqnR);
                RedLinSysU = problem->Umodes.project(UEqnR,
                                                     NmodesUproj); // Modes without lifting
                Eigen::MatrixXd projGradP = projGradModP * p;
                RedLinSysU[1] = RedLinSysU[1] - projGradP; // projGradP;
                u = reducedProblem::solveLinearSys(RedLinSysU, u,
                                                   uResidual);//, "fullPivLu");//"bdcSvd");
                problem->Umodes.reconstruct(U, u, "U");
                ITHACAutilities::assignBC(U, idInl, Uinlet);
                fvOptions.correct(U);
                //Energy equation phase
                fvScalarMatrix EEqnR
                (
                    fvm::div(phi, E)
                    + fvc::div(phi, volScalarField("Ekp", 0.5 * magSqr(U) + P / rho))
                    - fvm::laplacian(problem->turbulence->alphaEff(), E)
                    ==
                    fvOptions(rho, E)
                );
                EEqnR.relax();
                fvOptions.constrain(EEqnR);
                List<Eigen::MatrixXd> RedLinSysE = problem->Emodes.project(EEqnR, NmodesEproj);
                e = reducedProblem::solveLinearSys(RedLinSysE, e, eResidual);
                problem->Emodes.reconstruct(E, e, "e");
                //EEqnR.solve(); //For debug purposes only
                fvOptions.correct(E);
                thermo.correct(); // Here are calculated both temperature and density based on P,U and he.
                // Pressure equation phase
                constrainPressure(P, rho, U, problem->getPhiHbyA(UEqnR, U, P),
                                  problem->getRhorAUf(
                                      UEqnR));// Update the pressure BCs to ensure flux consistency
                surfaceScalarField phiHbyACalculated = problem->getPhiHbyA(UEqnR, U, P);
                closedVolume = adjustPhi(phiHbyACalculated, U, P);
                List<Eigen::MatrixXd> RedLinSysP;

                while (problem->_simple().correctNonOrthogonal())
                {
                    volScalarField rAU(1.0 /
                                       UEqnR.A()); // Inverse of the diagonal part of the U equation matrix
                    volVectorField HbyA(constrainHbyA(rAU * UEqnR.H(), U,
                                                      P)); // H is the extra diagonal part summed to the r.h.s. of the U equation
                    surfaceScalarField phiHbyA("phiHbyA", fvc::interpolate(rho)*fvc::flux(HbyA));
                    surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho * rAU));
                    fvScalarMatrix PEqnR
                    (
                        fvc::div(phiHbyA)
                        - fvm::laplacian(rhorAUf, P)
                        ==
                        fvOptions(psi, P, rho.name())
                    );
                    PEqnR.setReference
                    (
                        problem->_pressureControl().refCell(),
                        problem->_pressureControl().refValue()
                    );
                    RedLinSysP = problem->Pmodes.project(PEqnR, NmodesPproj);
                    p = reducedProblem::solveLinearSys(RedLinSysP, p, pResidual);
                    problem->Pmodes.reconstruct(P, p, "p");

                    if (problem->_simple().finalNonOrthogonalIter())
                    {
                        phi = problem->getPhiHbyA(UEqnR, U, P) + PEqnR.flux();
                    }
                }

                //#include "continuityErrs.H"
#include "incompressible/continuityErrs.H"
                P.relax();// Explicitly relax pressure for momentum corrector
                U = problem->HbyA() - (1.0 / UEqnR.A()) * problem->getGradP(P);
                U.correctBoundaryConditions();
                fvOptions.correct(U);
                bool pLimited = problem->_pressureControl().limit(P);

                // For closed-volume cases adjust the pressure and density levels to obey overall mass continuity
                if (closedVolume)
                {
                    P += (problem->_initialMass() - fvc::domainIntegrate(psi * P))
                         / fvc::domainIntegrate(psi);
                }

                if (pLimited || closedVolume)
                {
                    P.correctBoundaryConditions();
                }

                rho = thermo.rho(); // Here rho is calculated as p*psi = p/(R*T)
                rho.relax();
```
Residuals are monitored through both normalized magnitude and residual jump, and the neural-network closure for nut is updated during the iterations. At convergence, the final fields are exported.
```cpp
                residualNorm = max(max((uResidual.cwiseAbs()).sum() /
                                       (RedLinSysU[1].cwiseAbs()).sum(),
                                       (pResidual.cwiseAbs()).sum() / (RedLinSysP[1].cwiseAbs()).sum()),
                                   (eResidual.cwiseAbs()).sum() / (RedLinSysE[1].cwiseAbs()).sum());
                residualJump = max(max(((uResidual - uResidualOld).cwiseAbs()).sum() /
                                       (RedLinSysU[1].cwiseAbs()).sum(),
                                       ((pResidual - pResidualOld).cwiseAbs()).sum() /
                                       (RedLinSysP[1].cwiseAbs()).sum()),
                                   ((eResidual - eResidualOld).cwiseAbs()).sum() /
                                   (RedLinSysE[1].cwiseAbs()).sum());
                //problem->turbulence->correct(); // Resolution of the full turbulence (debug purposes only)
                nutCoeff = problem->evalNet(u, mu_now);
                problem->nutModes.reconstruct(nut, nutCoeff, "nut");
            }

            ITHACAstream::exportSolution(U, name(counter), Folder);
            ITHACAstream::exportSolution(P, name(counter), Folder);
            ITHACAstream::exportSolution(E, name(counter), Folder);
            ITHACAstream::exportSolution(nut, name(counter), Folder);
```

The class `tutorial03` provides the specific tutorial setup.
```cpp
class tutorial03 : public CompressibleSteadyNN
{
    public:
        /// Constructor
        explicit tutorial03(int argc, char* argv[])
            :
            CompressibleSteadyNN(argc, argv)
        {
            /// Export intermediate steps
            middleExport = para->ITHACAdict->lookupOrDefault<bool>("middleExport", true);
        }

        ITHACAparameters* para = ITHACAparameters::getInstance();
```

Its `offlineSolve` method either reads previously computed snapshots or performs the full-order offline simulations over the parameter samples, where the varying parameter is the viscosity.
```cpp
        /// Perform an Offline solve
        void offlineSolve(word folder = "./ITHACAoutput/Offline/")
        {
            //std::ofstream cpuTimes;
            //double durationOff;
            //cpuTimes.open(folder + "/cpuTimes", std::ios_base::app);
            /// Velocity field
            volVectorField& U = _U();
            /// Pressure field
            volScalarField& p = _p();
            /// Energy field
            volScalarField& E = _E();

            // if the offline solution is already performed read the fields
            if (offline && !ITHACAutilities::check_folder("./ITHACAoutput/POD/1"))
            {
                ITHACAstream::readMiddleFields(Ufield, U, folder);
                ITHACAstream::readMiddleFields(Efield, E, folder);
                ITHACAstream::readMiddleFields(Pfield, p, folder);
                /// Eddy viscosity field
                auto nut = _mesh().lookupObject<volScalarField>("nut");
                ITHACAstream::readMiddleFields(nutFields, nut, folder);
                mu_samples = ITHACAstream::readMatrix("./parsOff_mat.txt");
            }
            // else perform offline stage
            else if (!offline)
            {
                double UIFinit = para->ITHACAdict->lookupOrDefault<double>("UIFinit", 250);
                Vector<double> Uinl(UIFinit, 0, 0);
                //Vector<double> Uinl(250, 0, 0);
            ...
```

For each sample, the code updates the viscosity, assigns the inlet condition, and runs the truth solve, storing velocity, pressure, energy, and turbulent-viscosity snapshots.
```cpp
                for (label i = 0; i < mu.rows(); i++)
                {
                    Info << "Current mu = " << mu(i, 0) << endl;
                    changeViscosity(mu(i, 0));
                    assignIF(_U(), Uinl);
                    truthSolve(folder);
                }
```

### The main
In main, the tutorial object is created and the offline and online parameter sets are either loaded from file or generated automatically.
```cpp
int main(int argc, char* argv[])
{
    // Construct the tutorial object
    tutorial03 example(argc, argv);
    //Info << example.pThermo().p() << endl;
    //Info << example._p() << endl;
    // exit(0);
    //std::clock_t startOff, startOn;
    //double durationOn, durationOff;
    std::cerr << "debug point 1" << endl;
    ITHACAparameters* para = ITHACAparameters::getInstance();
    //Eigen::MatrixXd parOff;
    std::ifstream exFileOff("./parsOff_mat.txt");

    if (exFileOff)
    {
        example.mu  = ITHACAstream::readMatrix("./parsOff_mat.txt");
    }
    else
    {
        //example.mu  = ITHACAutilities::rand(20, 1, 1.00e-05, 1.00e-2);
        label OffNum = para->ITHACAdict->lookupOrDefault<label>("OffNum", 25);
        example.mu  = Eigen::VectorXd::LinSpaced(OffNum, 1.00e-05, 1.00e-02);
        ITHACAstream::exportMatrix(example.mu, "parsOff", "eigen", "./");
    }

    Eigen::MatrixXd parsOn;
    std::ifstream exFileOn("./parsOn_mat.txt");

    if (exFileOn)
    {
        parsOn = ITHACAstream::readMatrix("./parsOn_mat.txt");
    }
    else
    {
        label OnNum = para->ITHACAdict->lookupOrDefault<label>("OnNum", 20);
        parsOn = ITHACAutilities::rand(OnNum, 1, 1.00e-05, 1.00e-02);
        ITHACAstream::exportMatrix(parsOn, "parsOn", "eigen", "./");
    }
```

The code then reads the numbers of POD and projection modes from the dictionary, performs the offline stage:

```cpp
    int NmodesUout = para->ITHACAdict->lookupOrDefault<int>("NmodesUout", 15);
    int NmodesPout = para->ITHACAdict->lookupOrDefault<int>("NmodesPout", 15);
    int NmodesEout = para->ITHACAdict->lookupOrDefault<int>("NmodesEout", 15);
    int NmodesNutOut = para->ITHACAdict->lookupOrDefault<int>("NmodesNutOut", 15);
    int NmodesUproj = para->ITHACAdict->lookupOrDefault<int>("NmodesUproj", 10);
    int NmodesPproj = para->ITHACAdict->lookupOrDefault<int>("NmodesPproj", 10);
    int NmodesEproj = para->ITHACAdict->lookupOrDefault<int>("NmodesEproj", 10);
    int NmodesNutProj = para->ITHACAdict->lookupOrDefault<int>("NmodesNutProj", 10);
    example.offlineSolve();
```

It then computes POD modes for velocity, pressure, energy, and `nut`:
```cpp
    ITHACAPOD::getModes(example.Ufield, example.Umodes, example._U().name(),
                        example.podex, 0, 0,
                        NmodesUout);
    ITHACAPOD::getModes(example.Pfield, example.Pmodes, example._p().name(),
                        example.podex, 0, 0,
                        NmodesPout);
    ITHACAPOD::getModes(example.Efield, example.Emodes, example._E().name(),
                        example.podex, 0, 0,
                        NmodesEout);
    ITHACAPOD::getModes(example.nutFields, example.nutModes, "nut", example.podex,
                        0, 0,
                        NmodesNutOut);
```

If the validation folder `checkOff` is missing, a second full-order run is performed on the online parameter set to generate reference solutions.
```cpp
if (!ITHACAutilities::check_folder("./ITHACAoutput/checkOff"))
    {
        tutorial03 checkOff(argc, argv);
        checkOff.mu  = ITHACAstream::readMatrix("./parsOn_mat.txt");
        //Perform the offline solve
        //checkOff.middleExport = para->ITHACAdict->lookupOrDefault<bool>("middleExport", true);
        checkOff.offline = false;
        checkOff.restart();
        checkOff.offlineSolve("./ITHACAoutput/checkOff/");
        //durationOff = (std::clock() - startOff);
        checkOff.offline = true;
    }
```

Otherwise, the stored validation snapshots are reorganized into a one-sample-per-case structure for direct comparison.
```cpp
    else 
    {
        PtrList<volVectorField> UfieldCheck;
        PtrList<volScalarField> PfieldCheck;
        PtrList<volScalarField> EfieldCheck;
        PtrList<volScalarField> nutFieldsCheck;
        ITHACAstream::readMiddleFields(UfieldCheck, example._U(),
                                       "./ITHACAoutput/checkOff/");
        ITHACAstream::readMiddleFields(PfieldCheck, example._p(),
                                       "./ITHACAoutput/checkOff/");
        ITHACAstream::readMiddleFields(EfieldCheck, example._E(),
                                       "./ITHACAoutput/checkOff/");
        auto nutCheck = example._mesh().lookupObject<volScalarField>("nut");
        ITHACAstream::readMiddleFields(nutFieldsCheck, nutCheck,
                                       "./ITHACAoutput/checkOff/");
        Eigen::MatrixXd snapsCheck =
            ITHACAstream::readMatrix("./ITHACAoutput/checkOff/snaps");
        label fieldNum = 0;

        for (label k = 0; k < snapsCheck.rows(); k++)
        {
            fieldNum = fieldNum + snapsCheck(k, 0);
            //Info << "fieldNum" << fieldNum << endl;
            ITHACAstream::exportSolution(UfieldCheck[fieldNum - 1], name(k + 1),
                                         "./ITHACAoutput/checkOffSingle/");
            ITHACAstream::exportSolution(PfieldCheck[fieldNum - 1], name(k + 1),
                                         "./ITHACAoutput/checkOffSingle/");
            ITHACAstream::exportSolution(EfieldCheck[fieldNum - 1], name(k + 1),
                                         "./ITHACAoutput/checkOffSingle/");
            ITHACAstream::exportSolution(nutFieldsCheck[fieldNum - 1], name(k + 1),
                                         "./ITHACAoutput/checkOffSingle/");
        }
    }
```

The reduced datasets for neural-network training are then generated, the trained `TorchScript` model is loaded, and the reduced solver is created.
```cpp
    example.getTurbNN();
    example.loadNet("ITHACAoutput/NN/Net_" + name(NmodesUproj) + "_" + name(
                        NmodesNutProj) + ".pt");
    ReducedCompressibleSteadyNN reduced(example);
```

If the online results are not yet available, the code loops over all online parameters, updates the viscosity, projects the reduced operators, resets the solver state, validates the turbulence model, and performs the online reduced solve while storing CPU times.
```cpp
if (!ITHACAutilities::check_folder("./ITHACAoutput/Online_" + name(
                                           NmodesUproj) + "_" + name(NmodesNutProj) + "/" ) )
    {
        std::ofstream cpuTimes;
        word OnlineFolder = "./ITHACAoutput/Online_" + name(NmodesUproj) + "_" + name(
                                NmodesNutProj) + "/";
        cpuTimes.open(OnlineFolder + "/cpuTimes", std::ios_base::app);

        for (label k = 0; k < parsOn.rows(); k++)
        {
            std::clock_t startOn;
            startOn = std::clock();
            double durationOn;
            Eigen::MatrixXd mu_now = parsOn.row(k);
            mu_now.transposeInPlace();
            example.changeViscosity(mu_now(0, 0));
            reduced.projectReducedOperators(NmodesUproj, NmodesPproj, NmodesEproj);
            example.restart();
            example.turbulence->validate();
            reduced.solveOnlineCompressible(NmodesUproj, NmodesPproj, NmodesEproj,
                                            NmodesNutProj, mu_now, "./ITHACAoutput/Online_" + name(NmodesUproj) + "_" +
                                            name(NmodesNutProj) + "/");
            durationOn = std::clock() - startOn;
            cpuTimes << durationOn << endl;
        }
    }
```

If the online folder already exists, the code reads both reduced and full-order fields, computes normalized error fields for velocity, pressure, energy, and turbulent viscosity, exports them, and also computes relative L2 error matrices for all variables.
```cpp
   else if (ITHACAutilities::check_folder("./ITHACAoutput/Online_" + name(
            NmodesUproj) + "_" + name(NmodesNutProj) + "/" ) )
    {
        PtrList<volVectorField> offlineU, onlineU;
        PtrList<volScalarField> offlineP, onlineP;
        PtrList<volScalarField> offlineE, onlineE;
        PtrList<volScalarField> offlineNut, onlineNut;
        ITHACAstream::read_fields(onlineU, example._U(),
                                  "./ITHACAoutput/Online_" + name(NmodesUproj) + "_" + name(NmodesNutProj) + "/");
        ITHACAstream::read_fields(onlineP, example._p(),
                                  "./ITHACAoutput/Online_" + name(NmodesUproj) + "_" + name(NmodesNutProj) + "/");
        ITHACAstream::read_fields(onlineE, example._E(),
                                  "./ITHACAoutput/Online_" + name(NmodesUproj) + "_" + name(NmodesNutProj) + "/");
        auto nut = example._mesh().lookupObject<volScalarField>("nut");
        /// Save the online fields
        ITHACAstream::read_fields(onlineNut, nut,
                                  "./ITHACAoutput/Online_" + name(NmodesUproj) + "_" + name(NmodesNutProj) + "/");
        /// Save the errors
        ITHACAstream::read_fields(offlineU, example._U(),
                                  "./ITHACAoutput/checkOffSingle/");
        ITHACAstream::read_fields(offlineP, example._p(),
                                  "./ITHACAoutput/checkOffSingle/");
        ITHACAstream::read_fields(offlineE, example._E(),
                                  "./ITHACAoutput/checkOffSingle/");
        ITHACAstream::read_fields(offlineNut, nut, "./ITHACAoutput/checkOffSingle/");
        for (label j = 0; j < parsOn.rows(); j++)
        {
            volVectorField Ue = offlineU[j] - onlineU[j];
            auto u = ITHACAutilities::L2Norm(offlineU[j]);
            Ue /= u;
            volScalarField Pe = offlineP[j] - onlineP[j];
            auto p = ITHACAutilities::L2Norm(offlineP[j]);
            Pe /= p;
            volScalarField Ee = offlineE[j] - onlineE[j];
            auto e = ITHACAutilities::L2Norm(offlineE[j]);
            Ee /= e;
            volScalarField Nute = offlineNut[j] - onlineNut[j];
            auto n = ITHACAutilities::L2Norm(offlineNut[j]);
            Nute /= n;
            Ue.rename("Ue");
            Pe.rename("Pe");
            Ee.rename("Ee");
            Nute.rename("Nute");
            ITHACAstream::exportSolution(Ue,   name(j + 1),
                                         "./ITHACAoutput/ErrorFields_" + name(NmodesUproj) + "_" + name(
                                             NmodesNutProj) + "/");
            ITHACAstream::exportSolution(Pe,   name(j + 1),
                                         "./ITHACAoutput/ErrorFields_" + name(NmodesUproj) + "_" + name(
                                             NmodesNutProj) + "/");
            ITHACAstream::exportSolution(Ee,   name(j + 1),
                                         "./ITHACAoutput/ErrorFields_" + name(NmodesUproj) + "_" + name(
                                             NmodesNutProj) + "/");
            ITHACAstream::exportSolution(Nute, name(j + 1),
                                         "./ITHACAoutput/ErrorFields_" + name(NmodesUproj) + "_" + name(
                                             NmodesNutProj) + "/");
        }

        Eigen::MatrixXd errorU = ITHACAutilities::errorL2Rel(offlineU, onlineU);
        Eigen::MatrixXd errorP = ITHACAutilities::errorL2Rel(offlineP, onlineP);
        Eigen::MatrixXd errorE = ITHACAutilities::errorL2Rel(offlineE, onlineE);
        Eigen::MatrixXd errorNut = ITHACAutilities::errorL2Rel(offlineNut, onlineNut);
        ///
        ITHACAstream::exportMatrix(errorU,
                                   "errorU" + name(NmodesUproj) + "_" + name(NmodesNutProj),     "python",
                                   "./ITHACAoutput/ErrorFields_" + name(NmodesUproj) + "_" + name(
                                       NmodesNutProj) + "/");
        ITHACAstream::exportMatrix(errorP,
                                   "errorP" + name(NmodesUproj) + "_" + name(NmodesNutProj),     "python",
                                   "./ITHACAoutput/ErrorFields_" + name(NmodesUproj) + "_" + name(
                                       NmodesNutProj) + "/");
        ITHACAstream::exportMatrix(errorE,
                                   "errorE" + name(NmodesUproj) + "_" + name(NmodesNutProj),     "python",
                                   "./ITHACAoutput/ErrorFields_" + name(NmodesUproj) + "_" + name(
                                       NmodesNutProj) + "/");
        ITHACAstream::exportMatrix(errorNut,
                                   "errorNut" + name(NmodesUproj) + "_" + name(NmodesNutProj), "python",
                                   "./ITHACAoutput/ErrorFields_" + name(NmodesUproj) + "_" + name(
                                       NmodesNutProj) + "/");
    }
```
## Other files
* `03compSteadyNS.C` – main source with `CompressibleSteadyNN` definition.
* `compSteady.py` – Python helpers for feature engineering and data loading.
* `train.py` – standalone training script (PyTorch).
* `parsOff_mat.txt`, `parsOn_mat.txt` – sample parameter matrices.
* `Make/` – build directory; output binary: `03compSteadyNS`.

## The plain code
The plain code is available [here](https://raw.githubusercontent.com/ITHACA-FV/ITHACA-FV/refs/heads/master/tutorials/NN/CompressibleSteadyNS/03compSteadyNS.C).