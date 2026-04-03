# Tutorial NN‑02
## Introduction
This tutorial extends the reduced‑order modelling framework to a compressible steady Navier–Stokes problem on a bump geometry. A neural network is trained to predict the eddy viscosity (`nut`) coefficients in the reduced system, enabling a data‑driven turbulence closure.

## The script `02compBump.C`

The included headers provide the compressible full-order solver, the reduced-order solver, POD and I/O utilities, mesh-motion tools, Eigen/Torch conversion functions, and OpenFOAM support for forces and field manipulation. The Torch interface is used to load and evaluate the trained neural network during the online stage.

```cpp
#include <torch/script.h>
#include <torch/torch.h>
#include "torch2Eigen.H"
#include "CompressibleSteadyNS.H"
#include "ReducedCompressibleSteadyNS.H"
#include "RBFMotionSolver.H"
#include "ITHACAstream.H"
#include "ITHACAPOD.H"
#include "forces.H"
#include "IOmanip.H"
#include "Foam2Eigen.H"
```

The class `CompressibleSteadyNN` inherits from `CompressibleSteadyNS` and adds network architecture & training logic similar to `NN‑01` tutorial. In the constructor, the number of projected velocity and turbulent-viscosity modes is read from the dictionary, and a fully connected neural network with two hidden layers is defined. The class also stores the normalization matrices required for preprocessing the network input and postprocessing its output, together with the `TorchScript` model used for inference.

```cpp
class CompressibleSteadyNN : public CompressibleSteadyNS
{
    public:
        /// Constructor
        CompressibleSteadyNN(int argc, char* argv[])
            :
            CompressibleSteadyNS(argc, argv)
        {
            ITHACAparameters* para = ITHACAparameters::getInstance();
            NUmodes = para->ITHACAdict->lookupOrDefault<label>("NmodesUproj", 10);
            NNutModes = para->ITHACAdict->lookupOrDefault<label>("NmodesNutProj", 10);
            Net->push_back(torch::nn::Linear(NUmodes, 128));
            Net->push_back(torch::nn::ReLU());
            Net->push_back(torch::nn::Linear(128, 64));
            Net->push_back(torch::nn::ReLU());
            Net->push_back(torch::nn::Linear(64, NNutModes));
            optimizer = new torch::optim::Adam(Net->parameters(),
                                               torch::optim::AdamOptions(2e-2));
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

The method `loadNet` loads the trained network from disk and reads the corresponding normalization arrays. These arrays are required because the network is not evaluated directly on raw reduced coefficients, but on scaled inputs, and its predictions are mapped back to the reduced physical space after inference.

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

The method `getTurbNN` builds the reduced datasets for neural-network training and testing. It reads velocity, pressure, and turbulent-viscosity snapshots from the offline and validation folders, projects them onto the POD bases, transposes the resulting coefficient matrices into sample-wise form, and stores them as NumPy files. This prepares the data used by the external Python script that trains the closure model.

```cpp
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
                ...
            }
        }
```

The two overloads of `evalNet` perform neural-network evaluation. One version takes reduced velocity coefficients together with the current parameter vector:

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

while the other uses a preassembled input:

```cpp
        // Function to eval the NN once the input is provided
        Eigen::MatrixXd evalNet(Eigen::MatrixXd a)
        {
            netTorchscript.eval();
            a.transposeInPlace();
            a = (a.array() - bias_inp.array()) / scale_inp.array();
            torch::Tensor aTensor = eigenMatrix2torchTensor(a);
            torch::Tensor out;// = Net->forward(aTensor);
            std::vector<torch::jit::IValue> XTensorInp;
            XTensorInp.push_back(aTensor);
            out = netTorchscript.forward(XTensorInp).toTensor();
            Eigen::MatrixXd g = torchTensor2eigenMatrix<double>(out);
            g = g.array() * scale_out.array() + bias_out.array();
            g.transposeInPlace();
            return g;
        }
```

In both cases, the input is normalized, converted into a Torch tensor, passed through the TorchScript model, and converted back into denormalized reduced coefficients for `nut`.

The class `ReducedCompressibleSteadyNN` extends `ReducedCompressibleSteadyNS` and contains the online reduced solver.

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

The method `projectReducedOperators` precomputes the projection of the pressure-gradient contribution onto the velocity reduced space. This avoids rebuilding the same operator structure repeatedly during the online iterations.

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

The main online routine is `solveOnlineCompressible`. It performs one reduced simulation for a given parameter vector and initializes the reduced coefficients for velocity, pressure, energy, and turbulent viscosity, together with the convergence monitors.

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
            // Parameters definition
            ITHACAparameters* para = ITHACAparameters::getInstance();
            float residualJumpLim =
                para->ITHACAdict->lookupOrDefault<float>("residualJumpLim", 1e-5);
            float normalizedResidualLim =
                para->ITHACAdict->lookupOrDefault<float>("normalizedResidualLim", 1e-5);
            int maxIter =
                para->ITHACAdict->lookupOrDefault<float>("maxIter", 2000);
            bool closedVolume = false;
            label csolve = 0;
            // Full variables initialization
            fluidThermo& thermo = problem->pThermo();
            volVectorField& U = problem->_U();
            volScalarField& P = problem->_p();
            volScalarField& E = problem->pThermo->he();
            volScalarField& nut = const_cast<volScalarField&>
                                  (problem->_mesh().lookupObject<volScalarField>("nut"));
            volScalarField& rho = problem->_rho();
            volScalarField& psi = problem->_psi();
            surfaceScalarField& phi = problem->_phi();
            Time& runTime = problem->_runTime();
            fvMesh& mesh = problem->_mesh();
            fv::options& fvOptions = problem->_fvOptions();
            scalar cumulativeContErr = problem->cumulativeContErr;
            // Reduced variables initialization
            Eigen::MatrixXd u = Eigen::MatrixXd::Zero(NmodesUproj, 1);
            Eigen::MatrixXd e = Eigen::MatrixXd::Zero(NmodesEproj, 1);
            Eigen::MatrixXd p = Eigen::MatrixXd::Zero(NmodesPproj, 1);
            Eigen::MatrixXd nutCoeff = ITHACAutilities::getCoeffs(nut, problem->nutModes,
                                       NmodesNutProj, true);
            //vector Uinlet(170,0,0); // Vector for the inlet boundary condition
            label idInl =
                problem->_mesh().boundaryMesh().findPatchID("inlet"); // ID of the inlet patch
            vector Uinlet(problem->_U().boundaryFieldRef()[idInl][0][0], 0, 0);
            P.storePrevIter();
```

It then enters an iterative SIMPLE-like loop in which the momentum, energy, and pressure equations are assembled in full-order form, projected onto the corresponding reduced bases, and solved in reduced coordinates.

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
```

In the momentum step, the compressible turbulent momentum equation is projected onto the velocity basis, and the pressure-gradient contribution is added through the precomputed reduced operator. The resulting reduced system is solved and the full velocity field is reconstructed. 

```cpp
            //Momentum equation phase
                List<Eigen::MatrixXd> RedLinSysU;
                ITHACAutilities::assignBC(U, idInl, Uinlet);
                fvVectorMatrix UEqnR
                (
                    fvm::div(phi, U)
                    - fvc::div((rho * problem->turbulence->nuEff()) * dev2(T(fvc::grad(U))))
                    - fvm::laplacian(rho * problem->turbulence->nuEff(), U)
                    ==
                    fvOptions(rho, U)
                );
                UEqnR.relax();
                fvOptions.constrain(UEqnR);
                RedLinSysU = problem->Umodes.project(UEqnR,
                                                     NmodesUproj); // Modes without lifting
                Eigen::MatrixXd projGradP = projGradModP * p;
                RedLinSysU[1] = RedLinSysU[1] - projGradP;
                u = reducedProblem::solveLinearSys(RedLinSysU, u,
                                                   uResidual);//, "fullPivLu");//"bdcSvd");
                problem->Umodes.reconstruct(U, u, "U");
                ITHACAutilities::assignBC(U, idInl, Uinlet);
                //solve(UEqnR == -problem->getGradP(P)); //For debug purposes only, second part only useful when using uEqn_global==-getGradP
                fvOptions.correct(U);
```
In the energy step, the total-energy equation is assembled, projected, solved in reduced form, and reconstructed into the full field. After that, the thermodynamic state is updated so that temperature and density remain consistent with pressure and energy.

```cpp
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
                thermo.correct(); 
```
The pressure step follows the standard compressible pressure-correction structure. The pressure equation is assembled inside the non-orthogonal correction loop, projected onto the pressure basis, and solved in reduced form.

```cpp
                constrainPressure(P, rho, U, problem->getPhiHbyA(UEqnR, U, P),
                                  problem->getRhorAUf(
                                      UEqnR));// Update the pressure BCs to ensure flux consistency
                surfaceScalarField phiHbyACalculated = problem->getPhiHbyA(UEqnR, U, P);
                closedVolume = adjustPhi(phiHbyACalculated, U, P);
                List<Eigen::MatrixXd> RedLinSysP;
```

The pressure field is reconstructed, the flux is corrected, and continuity consistency is enforced.

```cpp
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
```

At the end of each global iteration, pressure, velocity, and density are updated, and convergence is checked through both normalized residuals and residual jumps.


```cpp
#include "continuityErrs.H"
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
                Info << "Ures = " << (uResidual.cwiseAbs()).sum() /
                          (RedLinSysU[1].cwiseAbs()).sum() << endl;
                Info << "Eres = " << (eResidual.cwiseAbs()).sum() /
                          (RedLinSysE[1].cwiseAbs()).sum() << endl;
                Info << "Pres = " << (pResidual.cwiseAbs()).sum() /
                          (RedLinSysP[1].cwiseAbs()).sum() << endl;
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
            }
```

Once convergence is reached, the neural network is evaluated using the final reduced velocity coefficients and the current parameter vector. The predicted reduced eddy-viscosity coefficients are reconstructed into the full nut field, and the solution fields `U`, `P`, `E`, and `nut` are exported. In this implementation, the turbulence closure is therefore applied after the reduced compressible system has converged, using the neural model as a low-dimensional reconstruction of the turbulent viscosity field.


```cpp

            nutCoeff = problem->evalNet(u, mu_now);
            problem->nutModes.reconstruct(nut, nutCoeff, "nut");
            label k = 1;
            // U.rename("Ur");
            // P.rename("Pr");
            // E.rename("Er");
            // nut.rename("nutR");
            ITHACAstream::exportSolution(U, name(counter), Folder);
            ITHACAstream::exportSolution(P, name(counter), Folder);
            ITHACAstream::exportSolution(E, name(counter), Folder);
            ITHACAstream::exportSolution(nut, name(counter), Folder);
```

The class `tutorial02` specializes the problem by introducing geometric parametrization through mesh motion. In the constructor, it loads the RBF mesh-motion dictionary, extracts the points belonging to the upper and lower deformable patches, initializes the RBF motion solver:

```cpp
class tutorial02 : public CompressibleSteadyNN
{
    public:
        /// Constructor
        explicit tutorial02(int argc, char* argv[])
            :
            CompressibleSteadyNN(argc, argv)
        {
            dyndict = new IOdictionary
            (
                IOobject
                (
                    "dynamicMeshDictRBF",
                    "./constant",
                    _mesh(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );
            ITHACAutilities::getPointsFromPatch(_mesh(), 0, top0, top0_ind);
            ITHACAutilities::getPointsFromPatch(_mesh(), 1, bot0, bot0_ind);
            // Info << _mesh().points().size() << endl;
            ms = new RBFMotionSolver(_mesh(), *dyndict);
```

Then, it stores the original mesh coordinates, and reads the option controlling intermediate export:

```cpp
            movingIDs = ms->movingIDs();
            x0 = ms->movingPoints();
            curX = x0;
            point0 = ms->curPoints();
            /// Export intermediate steps
            middleExport = para->ITHACAdict->lookupOrDefault<bool>("middleExport", true);
        }
        ...
```

The helper function `f1` defines the scalar shape used for the deformation. It produces a localized bump-like variation along the streamwise direction. 

```cpp
        double f1(double chord, double x)
        {
            double res = chord * (std::pow((x) / chord,
                                           0.5) * (1 - (x) / chord)) / (std::exp(15 * (x) / chord));
            return res;
        }
```

The method `moveBasis` applies this shape to a list of points by modifying their vertical coordinate according to the parameter amplitude.

```cpp
        List<vector> moveBasis(const List<vector>& originalPoints, double par)
        {
            List<vector> movedPoints(originalPoints);

            for (int i = 0; i < originalPoints.size(); i++)
            {
                movedPoints[i][2] += par * f1(1, movedPoints[i][0]);
            }

            return movedPoints;
        }
```

The method `updateMesh` resets the mesh to the reference configuration and, if nonzero parameters are provided, deforms the top and bottom boundaries using the RBF motion solver. This gives a smooth geometry update for each parameter sample.

```cpp
        void updateMesh(double parTop = 0, double parBot = 0)
        {
            _mesh().movePoints(point0);

            if (parTop != 0 || parBot != 0)
            {
                // Info << parTop << endl;
                List<vector> top0_cur = moveBasis(top0, parTop);
                List<vector> bot0_cur = moveBasis(bot0, parBot);
                ITHACAutilities::setIndices2Value(top0_ind, top0_cur, movingIDs, curX);
                ITHACAutilities::setIndices2Value(bot0_ind, bot0_cur, movingIDs, curX);
                ms->setMotion(curX - x0);
                point = ms->curPoints();
                _mesh().movePoints(point);
            }
        }
```

The method `offlineSolve` manages the offline stage. If offline snapshots already exist, it simply reads the velocity, energy, pressure, and turbulent-viscosity fields from disk, together with the stored parameter samples. 

```cpp
        /// Perform an Offline solve
        void offlineSolve(word folder = "./ITHACAoutput/Offline/")
        {
            /// Velocity field
            volVectorField& U = _U();
            /// Pressure field
            volScalarField& p = _p();
            /// Energy field
            volScalarField& E = _E();

            // If the offline solution is already performed but POD modes are not present, then read the fields
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
```

Otherwise, it loops over all offline parameter values, updates the mesh, writes the deformed geometry, assigns the inlet initial condition, and runs the full-order truth solve. For each sample, the generated mesh is linked into the corresponding time directories so that the stored snapshots remain consistent with the deformed geometry:

```cpp
            else if (!offline)
            {
                double UIFinit = para->ITHACAdict->lookupOrDefault<double>("UIFinit", 170);
                Vector<double> Uinl(UIFinit, 0, 0);

                for (label i = 0; i < mu.rows(); i++)
                {
                    updateMesh(mu(i, 0), mu(i, 1));
                    ITHACAstream::writePoints(_mesh().points(), folder, name(i + 1) + "/polyMesh/");
                    assignIF(_U(), Uinl);
                    truthSolve(folder);
                    label j = 1;
                    word polyMesh2beLinked = folder + name(i + 1) + "/" + "polyMesh/";

                    while (ITHACAutilities::check_folder(folder + name(i + 1) + "/" + name(j)))
                    {
                        word folderContLink = folder + name(i + 1) + "/" + name(j) + "/";
                        system("ln -s  $(readlink -f " + polyMesh2beLinked + ") " + folderContLink +
                               " >/dev/null 2>&1");
                        j++;
                    }
                }
            }
        }
};
```


### The main
1. **Offline run**: read or generate the parameters sampling:
```cpp
// Construct the tutorial object
    tutorial02 example(argc, argv);
    ITHACAparameters* para = ITHACAparameters::getInstance();
    std::ifstream exFileOff("./parsOff_mat.txt");

    if (exFileOff)
    {
        example.mu  = ITHACAstream::readMatrix("./parsOff_mat.txt");
    }
    else
    {
        int OffNum = para->ITHACAdict->lookupOrDefault<int>("OffNum", 100);
        double BumpAmp = para->ITHACAdict->lookupOrDefault<double>("BumpAmp", 0.1);
        example.mu.resize(OffNum, 2);
        Eigen::MatrixXd parTop = ITHACAutilities::rand(example.mu.rows(), 1, 0,
                                 BumpAmp);
        Eigen::MatrixXd parBot = ITHACAutilities::rand(example.mu.rows(), 1, -BumpAmp,
                                 0);
        example.mu.leftCols(1) = parTop;
        example.mu.rightCols(1) = parBot;
        ITHACAstream::exportMatrix(example.mu, "parsOff", "eigen", "./");
    }
    ...
```
Read or generate the snapshot fields for velocity, pressure and turbulent viscosity using the full-order solver. Snapshots are stored under `./ITHACAoutput/Offline/` for later processing.

```cpp
    example.updateMesh();
    //Perform the offline solve
    example.offlineSolve();
    // Move the mesh to the original geometry to get the modes into a mid mesh
    example.updateMesh();
```

2. **Coefficient computation**: reads middle fields from the offline results, projects them onto POD modes, eventually train the network, and saves the resulting coefficient matrices (train/test):
```cpp
    // Perform POD on velocity and pressure
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

3. **Training**: the network is constructed with `getTurbNN` (linear layers with ReLU activations) and trained with Adam. Normalisation parameters (bias/scale) are saved alongside the TorchScript model.

4. **ROM and Online evaluation**: the binary loads the trained network via `loadNet()` and evaluates it during reduced simulations to supply eddy viscosity predictions, just as in the previous tutorial NN-01.

```cpp
    example.loadNet("ITHACAoutput/NN/Net_" + name(example.NUmodes) + "_" + name(
                        example.NNutModes) + ".pt");
```

The ROM with efficient `nut` evaluation through the network is performed via the online compressible solver:

```cpp
    ReducedCompressibleSteadyNN reduced(example);

    //Perform the online solutions
    for (label k = 0; k < parsOn.rows(); k++)
    {
        example.updateMesh();
        example.updateMesh(parsOn(k, 0), parsOn(k, 1));
        ITHACAstream::writePoints(example._mesh().points(),
                                  "./ITHACAoutput/Online_" + name(NmodesUproj) + "_" + name(NmodesNutProj) + "/",
                                  name(k + 1) + "/polyMesh/");
        //Info << example.inletIndex.rows() << endl;
        //reduced.setOnlineVelocity(vel);
        reduced.projectReducedOperators(NmodesUproj, NmodesPproj, NmodesEproj);
        example.restart();
        example.turbulence->validate(); ///////////////////////////////////////////////// Is it needed to validate the nut field?
        Eigen::MatrixXd mu_now = parsOn.row(k);
        mu_now.transposeInPlace();
        reduced.solveOnlineCompressible(NmodesUproj, NmodesPproj, NmodesEproj,
                                        NmodesNutProj, mu_now, "./ITHACAoutput/Online_" + name(NmodesUproj) + "_" +
                                        name(NmodesNutProj) + "/");
    }
```
The same procedure is repeated on the validation set (saved in the `ITHACAoutput/checkOff` folder).

Additional scripts (`02compBump.py`, `compBump.py`, `train.py`) support data management and training outside the C++ code. Data files such as `parsOff_mat.txt`/`parsOn_mat.txt` provide sample parameters, and `plots.py` produces diagnostic figures.


## The plain code
The plain code is available [here](https://raw.githubusercontent.com/ITHACA-FV/ITHACA-FV/refs/heads/master/tutorials/NN/BumpCompressibleSteadyNS/02compBump.C).