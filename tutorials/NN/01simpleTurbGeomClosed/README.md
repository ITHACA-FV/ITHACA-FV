# Tutorial NN-01

## Introduction
This tutorial couples a steady incompressible Navier–Stokes solver (SIMPLE algorithm) with a small neural network that predicts the eddy viscosity coefficient from reduced‑order (POD) velocity modes. The geometry is the `simpleTurbGeomClosed` case and the objective is to learn a turbulence model in a reduced‑order setting.

The workflow consists of two stages:
1. **Offline snapshot generation** using the parent class `SteadyNSSimple` to produce velocity, pressure and turbulent viscosity fields for a set of parameter values (e.g. angle of attack); these are written to `./ITHACAoutput/Offline/`.
2. **Neural‑network training**. The `SteadyNSSimpleNN` class (defined in `01simpleTurbGeomClosed.C`) computes L2 projection coefficients of the snapshots onto the POD modes and trains a fully‑connected network (two hidden layers of 128 and 64 neurons) using PyTorch. The network takes the reduced coefficients (and optionally a parameter value) as input and predicts the coefficients of the eddy viscosity field.

Once trained, the network is exported in TorchScript format and can be loaded with `SteadyNSSimpleNN::loadNet`. The online solver then evaluates the network to supply a turbulence closure while solving the reduced system.

# A detailed look into the code
Let's have a look at all the useful scripts for tutorial NN-01.

## Files of interest

* `01simpleTurbGeomClosed.C` – contains the `SteadyNSSimpleNN` definition and main routine orchestrating offline and online phases: we will comment this script, which is the core of the tutorial, in the next section.
* `train.py` – Python script for training the network outside the C++ code.
* `simpleTurbGeomClosed.py` – helper utilities for data preprocessing.
* `angOff_mat.txt`, `angOn_mat.txt`, `vel.txt` – example matrices and fields.

## Script 01simpleTurbGeomClosed.C
Useful libraries are first included:
```cpp
#include <torch/script.h>
#include <torch/torch.h>
#include "torch2Eigen.H"
#include "SteadyNSSimple.H"
#include "ITHACAstream.H"
#include "ITHACAPOD.H"
#include "ReducedSimpleSteadyNS.H"
#include "forces.H"
#include "IOmanip.H"
#include "RBFMotionSolver.H"
```
among which we have torch, which we will use to load the pre-trained neural network and to infer new outputs (the eddy viscosity coefficients) from new inputs (the velocity coefficients).

The statement `using namespace ITHACAtorch::torch2Eigen;` makes the tensor/matrix conversion utilities directly accessible. This is important because the code frequently moves data from the reduced basis representation stored as Eigen objects into Torch tensors for neural network inference, and then converts the predicted output back to Eigen form.

A class `SteadyNSSimpleNN` is here defined. This class extends `SteadyNSSimple`. Its purpose is to augment the standard steady SIMPLE full-order problem with a neural network interface. In other words, this class is responsible for everything related to the coupling between the CFD reduced-order workflow and the trained Torch model.

```cpp
class SteadyNSSimpleNN : public SteadyNSSimple
{
    public:
        /// Constructor
        SteadyNSSimpleNN(int argc, char* argv[])
            :
            SteadyNSSimple(argc, argv)
        {
            Net->push_back(torch::nn::Linear(NUmodes, 128));
            Net->push_back(torch::nn::ReLU());
            Net->push_back(torch::nn::Linear(128, 64));
            Net->push_back(torch::nn::ReLU());
            Net->push_back(torch::nn::Linear(64, NNutModes));
            optimizer = new torch::optim::Adam(Net->parameters(),
                                               torch::optim::AdamOptions(2e-2));
        };

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
The constructor first initializes the parent class and then defines a fully connected feedforward neural network. The network maps the reduced velocity coefficients to the reduced turbulent viscosity coefficients. The architecture consists of an input linear layer from `NUmodes` to 128 neurons, followed by a ReLU activation, then a second linear layer from 128 to 64 neurons, another ReLU, and finally an output layer from 64 to `NNutModes`. This is a straightforward multilayer perceptron used as a nonlinear regression model, using Adam as optimizer in the training process.
The matrices `bias_inp`, `scale_inp`, `bias_out`, and `scale_out` store the normalization parameters used for input and output preprocessing. These variables are essential because the neural network is not evaluated on raw reduced coefficients. Instead, the input is normalized before inference, and the output is denormalized afterward. This helps stabilize training and makes the model more robust numerically.
The matrices `coeffL2Nut` and `coeffL2U`, together with their tensor counterparts, are intended to store reduced coefficients for the turbulent viscosity and velocity fields. These quantities are part of the dataset preparation process used for neural network training. The members `Net`, `optimizer`, and `netTorchscript` represent the neural network in two different forms. `Net` is the trainable C++ Torch module, optimizer manages parameter updates during training, and `netTorchscript` is the serialized, pre-trained model loaded from disk for inference. In this file, the `TorchScript` model is the one actually used during the online phase.

The method `loadNet(word filename)` loads the trained neural network and all associated normalization data. First, it checks whether the file exists and aborts with a clear error message if it does not. This is a protective step because the reduced turbulent closure cannot work unless the network has already been trained for the chosen numbers of velocity and turbulent-viscosity modes. Then the `TorchScript` model is loaded from disk, and the normalization arrays are read from .npy files. Finally, the model is switched to evaluation mode, which disables any training-specific behavior and ensures deterministic inference.

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

The method `getTurbNN()` prepares the training and testing datasets for the neural network. It only runs if the folder `ITHACAoutput/NN/coeffs` does not already exist. The function first reads the offline snapshots for velocity, pressure, and turbulent viscosity from `./ITHACAoutput/Offline/`. Note that the snapshots are read with `readMiddleFields`: since the case is steady, we also include inside the snapshots matrix the middle states. This will improve the accuracy of the ROM.

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
```

It then projects these full fields onto their respective reduced bases using `ITHACAutilities::getCoeffs`. The result is a set of reduced coordinates for the training data. These coefficient matrices are transposed so that each row corresponds to a sample, which is the natural layout for machine learning datasets, and then saved as NumPy files.

```cpp
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
```

The second half of `getTurbNN()` performs the same procedure on the test snapshots stored in `./ITHACAoutput/checkOff/`. This generates reduced coefficients for velocity, pressure, and turbulent viscosity on the validation set.

```cpp
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

The first overload of `evalNet`, namely `Eigen::MatrixXd evalNet(Eigen::MatrixXd a, double mu_now)`, evaluates the neural network when both the reduced velocity coefficients and the current parameter value are part of the input. The function creates an input vector `xpred` with one additional entry beyond the reduced velocity coefficients.
The first element is the current parameter `mu_now`, while the remaining entries are the reduced velocity coefficients. After assembling the input, the code applies an affine transformation using `scale_inp` and `bias_inp`, converting the physical input into the normalized space expected by the network. The normalized input is then transposed, converted into a Torch tensor, and passed through the `TorchScript` network. The output tensor is converted back to an Eigen matrix, transposed into column format, and finally denormalized using `bias_out` and `scale_out`. The returned matrix therefore contains the predicted reduced coefficients for the turbulent viscosity field.

```cpp
        Eigen::MatrixXd evalNet(Eigen::MatrixXd a, double mu_now)
        {
            Eigen::MatrixXd xpred(a.rows() + 1, 1);

            if (xpred.cols() != 1)
            {
                throw std::runtime_error("Predictions for more than one sample not supported yet.");
            }

            xpred.bottomRows(a.rows()) = a;
            xpred(0, 0) = mu_now;
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

The second overload of `evalNet`, namely `Eigen::MatrixXd evalNet(Eigen::MatrixXd a)`, is a simpler variant in which the network input consists only of the reduced coefficients already packed into `a`. It normalizes the input, performs inference, denormalizes the output, and returns the predicted coefficients. Compared with the previous overload, this version assumes that the parameter dependence has already been embedded in the provided input matrix.

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

The class `reducedSimpleSteadyNN` extends `reducedSimpleSteadyNS`. Its role is to define the online reduced-order solver that uses the neural network closure. While the previous class handles the full-order problem plus the machine learning interface, this class actually performs the reduced SIMPLE iterations in the online stage.

```cpp
class reducedSimpleSteadyNN : public reducedSimpleSteadyNS
{
    public:
        /// Constructor
        explicit reducedSimpleSteadyNN(SteadyNSSimpleNN& FOMproblem)
            :
            problem(& FOMproblem)
        {}

        /// Full problem.
        SteadyNSSimpleNN* problem;
```

The core routine is `solveOnline_Simple`. This method performs one online reduced simulation for a given parameter value `mu_now`. It initializes reduced coefficients for velocity, pressure, and turbulent viscosity:

```cpp
        // Function to perform the online phase
        void solveOnline_Simple(scalar mu_now,
                                int NmodesUproj, int NmodesPproj, int NmodesNut = 0,
                                word Folder = "./ITHACAoutput/Reconstruct/")
        {
            counter++;

            ...

            // Initializations
            Eigen::VectorXd uresidualOld = Eigen::VectorXd::Zero(NmodesUproj);
            Eigen::VectorXd presidualOld = Eigen::VectorXd::Zero(NmodesPproj);
            Eigen::VectorXd uresidual = Eigen::VectorXd::Zero(NmodesUproj);
            Eigen::VectorXd presidual = Eigen::VectorXd::Zero(NmodesPproj);
            scalar U_norm_res(1);
            scalar P_norm_res(1);
            Eigen::MatrixXd a = Eigen::VectorXd::Zero(NmodesUproj);
            Eigen::MatrixXd a0 = Eigen::VectorXd::Zero(NmodesUproj);
            Eigen::MatrixXd b = Eigen::VectorXd::Zero(NmodesPproj);
            Eigen::MatrixXd b0 = Eigen::VectorXd::Zero(NmodesPproj);
            Eigen::MatrixXd bOld = b;
            Eigen::MatrixXd nutCoeff = Eigen::VectorXd::Zero(NmodesNut);
            Eigen::MatrixXd nutCoeffOld = Eigen::VectorXd::Zero(NmodesNut);
            float residualJumpLim =
                problem->para->ITHACAdict->lookupOrDefault<float>("residualJumpLim", 1e-5);
            float normalizedResidualLim =
                problem->para->ITHACAdict->lookupOrDefault<float>("normalizedResidualLim",
                    1e-5);
            scalar residual_jump(1 + residualJumpLim);
            volScalarField& P = problem->_p();
            volVectorField& U = problem->_U();
            volScalarField& nut = const_cast<volScalarField&>
                                  (problem->_mesh().lookupObject<volScalarField>("nut"));
            volVectorField u2 = U;
            a0 = ITHACAutilities::getCoeffs(U, problem->Umodes, NmodesUproj, true);
            b = ITHACAutilities::getCoeffs(P, problem->Pmodes, NmodesPproj, true);
            nutCoeff = ITHACAutilities::getCoeffs(nut, problem->nutModes, NmodesNut, true);
            a(0) = a0(0); // Do not remove: it is not working without this condition
            b = b0;
            fvMesh& mesh = problem->_mesh();
            Time& runTime = problem->_runTime();
            P.rename("p");
            surfaceScalarField& phi(problem->_phi());
```

It then reconstructs the corresponding full fields, and then enters a SIMPLE iteration loop.

```cpp
            problem->Umodes.reconstruct(U, a, "U");
            problem->Pmodes.reconstruct(P, b, "p");
            vector v(1, 0, 0);
            ITHACAutilities::assignBC(U, 0, v);
            //problem->nutModes.reconstruct(nut, nutCoeff, "nut");
            phi = fvc::flux(U);
            int iter = 0;
            simpleControl& simple = problem->_simple();
            std::ofstream res_os_U;
            std::ofstream res_os_P;
            res_os_U.open(Folder + name(counter) + "/residualsU", std::ios_base::app);
            res_os_P.open(Folder + name(counter) + "/residualsP", std::ios_base::app);

            // SIMPLE algorithm starts here
            while ((residual_jump > residualJumpLim
                    || std::max(U_norm_res, P_norm_res) > normalizedResidualLim) &&
                    iter < maxIterOn)
            {
                iter++;
                Info << iter << endl;
#if defined(OFVER) && (OFVER == 6)
                simple.loop(runTime);
#else
                problem->_simple().loop();
#endif
```
Inside the SIMPLE loop, if the case is turbulent, the neural network is evaluated to predict the reduced turbulent-viscosity coefficients from the current reduced velocity state and parameter. These coefficients are then reconstructed into the full nut field.

```cpp
                if (ITHACAutilities::isTurbulent())
                {
                    nutCoeff = problem->evalNet(a, mu_now);
                    //nutCoeff = nutCoeffOld + 0.7*(nutCoeff - nutCoeffOld);
                    volScalarField& nut = const_cast<volScalarField&>
                                          (problem->_mesh().lookupObject<volScalarField>("nut"));
                    problem->nutModes.reconstruct(nut, nutCoeff, "nut");
                    //ITHACAstream::exportSolution(nut, name(counter), Folder);
                }
```

The momentum equation is assembled, projected onto the velocity POD space, and solved in reduced form.

```cpp
                volScalarField nueff = problem->turbulence->nuEff();
                vector v(1, 0, 0);
                ITHACAutilities::assignBC(U, 0, v);
                // Momentum equation
                fvVectorMatrix UEqn
                (
                    fvm::div(phi, U)
                    - fvm::laplacian(nueff, U)
                    - fvc::div(nueff * dev2(T(fvc::grad(U))))
                );
                UEqn.relax();
                //solve(UEqn == - fvc::grad(P));
                List<Eigen::MatrixXd> RedLinSysU = problem->Umodes.project(UEqn, NmodesUproj);
                volVectorField gradPfull = -fvc::grad(P);
                Eigen::MatrixXd projGrad = problem->Umodes.project(gradPfull, NmodesUproj);
                RedLinSysU[1] = RedLinSysU[1] + projGrad;
                a = reducedProblem::solveLinearSys(RedLinSysU, a, uresidual);
                problem->Umodes.reconstruct(U, a, "U");
                ITHACAutilities::assignBC(U, 0, v);
                volScalarField rAU(1.0 / UEqn.A());
                volVectorField HbyA(constrainHbyA(1.0 / UEqn.A() * UEqn.H(), U, P));
                surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
                adjustPhi(phiHbyA, U, P);
                tmp<volScalarField> rAtU(rAU);

                if (simple.consistent())
                {
                    rAtU = 1.0 / (1.0 / rAU - UEqn.H1());
                    phiHbyA +=
                        fvc::interpolate(rAtU() - rAU) * fvc::snGrad(P) * mesh.magSf();
                    HbyA -= (rAU - rAtU()) * fvc::grad(P);
                }

                List<Eigen::MatrixXd> RedLinSysP;
                bOld = b;
```

Then the pressure-correction equation is assembled, projected onto the pressure POD space, and solved in reduced form through the standard SIMPLE logic.

```cpp
                while (simple.correctNonOrthogonal())
                {
                    // Continuity equation
                    fvScalarMatrix pEqn
                    (
                        fvm::laplacian(rAtU(), P) == fvc::div(phiHbyA)
                    );
                    RedLinSysP = problem->Pmodes.project(pEqn, NmodesPproj);
                    //pEqn.solve();
                    b = reducedProblem::solveLinearSys(RedLinSysP, b, presidual);
                    problem->Pmodes.reconstruct(P, b, "p");

                    if (simple.finalNonOrthogonalIter())
                    {
                        phi = phiHbyA - pEqn.flux();
                    }
                }

                b = bOld + mesh.fieldRelaxationFactor("p") * (b - bOld);
                problem->Pmodes.reconstruct(P, b, "p");
                nutCoeffOld = nutCoeff;
                // P.relax();
                U = HbyA - rAtU() * fvc::grad(P);
                U.correctBoundaryConditions();
```

Residuals are monitored using both their normalized magnitude and their jump between consecutive iterations, and the loop stops when the chosen thresholds are satisfied or the maximum number of iterations is reached. At convergence, the reconstructed fields are exported:

```cpp
                uresidualOld = uresidualOld - uresidual;
                presidualOld = presidualOld - presidual;
                uresidualOld = uresidualOld.cwiseAbs();
                presidualOld = presidualOld.cwiseAbs();
                residual_jump = std::max(uresidualOld.sum(), presidualOld.sum());
                uresidualOld = uresidual;
                presidualOld = presidual;
                uresidual = uresidual.cwiseAbs();
                presidual = presidual.cwiseAbs();
                U_norm_res = uresidual.sum() / (RedLinSysU[1].cwiseAbs()).sum();
                P_norm_res = presidual.sum() / (RedLinSysP[1].cwiseAbs()).sum();
                ...
                res_os_U << U_norm_res << endl;
                res_os_P << P_norm_res << endl;
            }

            res_os_U.close();
            res_os_P.close();
            ...
            problem->Umodes.reconstruct(U, a, "Uaux");
            problem->Pmodes.reconstruct(P, b, "Paux");

            if (ITHACAutilities::isTurbulent())
            {
                volScalarField& nut = const_cast<volScalarField&>
                                      (problem->_mesh().lookupObject<volScalarField>("nut"));
                nut.rename("nutAux");
                ITHACAstream::exportSolution(nut, name(counter), Folder);
            }

            ITHACAstream::exportSolution(U, name(counter), Folder);
            ITHACAstream::exportSolution(P, name(counter), Folder);
            runTime.setTime(runTime.startTime(), 0);
        }
```

The class `tutorial01cl` specializes the problem by adding geometric parameterization. It stores the original mesh points and the current moved coordinates, and defines how the mesh is deformed as a function of the input parameter, which in this example is an angle.

```cpp
class tutorial01cl : public SteadyNSSimpleNN
{
    public:
        /// Constructor
        explicit tutorial01cl(int argc, char* argv[])
            :
            SteadyNSSimpleNN(argc, argv)
        {
            curX = _mesh().points();
            point0 = _mesh().points();
        }

        /// Initial coordinates of the grid points
        vectorField point0;
        /// List to store the moved coordinates
        List<vector> curX;
```

The method `offlineSolve` either reads previously generated snapshots or performs the full-order simulations over the training parameters. 

```cpp
        void offlineSolve(Eigen::MatrixXd Box, List<label> patches,
                          word folder = "./ITHACAoutput/Offline/")
        {
            Vector<double> inl(0, 0, 0);
            List<scalar> mu_now(1);
            volVectorField& U = _U();
            volScalarField& p = _p();

            // if the offline solution is already performed read the fields
            if (offline && ITHACAutilities::isTurbulent() &&
                    !ITHACAutilities::check_folder("./ITHACAoutput/POD/1"))
            {
                ITHACAstream::readMiddleFields(Ufield, U, folder);
                ITHACAstream::readMiddleFields(Pfield, p, folder);
                auto nut = _mesh().lookupObject<volScalarField>("nut");
                ITHACAstream::readMiddleFields(nutFields, nut, folder);
            }
            else if (offline && !ITHACAutilities::check_folder("./ITHACAoutput/POD/1"))
            {
                ITHACAstream::readMiddleFields(Ufield, U, folder);
                ITHACAstream::readMiddleFields(Pfield, p, folder);
            }
            else if (!offline)
            {
                ...
```

For each parameter value, the mesh is reset, the selected points are moved:
```cpp
            
                for (label i = 0; i < mu.rows(); i++)
                {
                    _mesh().movePoints(point0);
                    List<vector> points2Move;
                    labelList boxIndices = ITHACAutilities::getIndicesFromBox(_mesh(), patches, Box,
                                           points2Move);
                    mu_now[0] = mu(i, 0);
                    linearMovePts(mu_now[0], points2Move);

                    for (int j = 0; j < boxIndices.size(); j++)
                    {
                        curX[boxIndices[j]] = points2Move[j];
                    }

                    Field<vector> curXV(curX);
                    _mesh().movePoints(curXV);
                    ITHACAstream::writePoints(_mesh().points(), folder,
                                              name(i + 1) + "/polyMesh/");
                    //assignIF(U, Uinl);
```
The truth problem is then solved, and the resulting snapshots are stored.
```cpp
                    truthSolve2(mu_now, folder);
                    int middleF = 1;

                    while (ITHACAutilities::check_folder(folder + name(
                            i + 1) + "/" + name(middleF)))
                    {
                        word command("ln -s  $(readlink -f " + folder + name(
                                         i + 1) + "/polyMesh/) " + folder + name(i + 1) + "/" + name(
                                         middleF) + "/" + " >/dev/null 2>&1");
                        std::cout.setstate(std::ios_base::failbit);
                        system(command);
                        std::cout.clear();
                        middleF++;
                    }

                    restart();
                }
            }
        }
```

The method `linearMovePts` defines the actual deformation law, applying a smooth position-dependent horizontal displacement controlled by the angle.

```cpp
        void linearMovePts(double angle, List<vector>& points2Move)
        {
            double sMax;
            scalarList x;
            scalarList y;

            for (label i = 0; i < points2Move.size(); i++)
            {
                x.append(points2Move[i].component(0));
                y.append(points2Move[i].component(1));
            }

            double xMin = min(x);
            double xMax = max(x);
            double yMin = min(y);
            double yMax = max(y);
            sMax = (yMax - yMin) * std::tan(M_PI * angle / 180);

            for (label i = 0; i < points2Move.size(); i++)
            {
                points2Move[i].component(0) = points2Move[i].component(0) +
                                              (yMax - points2Move[i].component(1)) / (yMax - yMin) * (xMax -
                                                  points2Move[i].component(0)) / (xMax - xMin) * sMax;
            }
        }
};
```

In `main`, the tutorial object is created and the parameter values for offline and online stages are loaded from file or generated if missing.

```cpp
int main(int argc, char* argv[])
{
    // Construct the tutorial object
    tutorial01cl example(argc, argv);
    // Read some parameters from file
    ITHACAparameters* para = ITHACAparameters::getInstance(example._mesh(),
                             example._runTime());
    // Read the files where the parameters are stored
    std::ifstream exFileOff("./angOff_mat.txt");

    if (exFileOff)
    {
        example.mu  = ITHACAstream::readMatrix("./angOff_mat.txt");
    }
    else
    {
        example.mu  = Eigen::VectorXd::LinSpaced(50, 0, 75);
        ITHACAstream::exportMatrix(example.mu, "angOff", "eigen", "./");
    }

    Eigen::MatrixXd angOn;
    std::ifstream exFileOn("./angOn_mat.txt");

    if (exFileOn)
    {
        angOn = ITHACAstream::readMatrix("./angOn_mat.txt");
    }
    else
    {
        //angOn = Eigen::VectorXd::LinSpaced(20, 5, 70);
        angOn = ITHACAutilities::rand(10, 1, 0, 70);
        ITHACAstream::exportMatrix(angOn, "angOn", "eigen", "./");
    }
```

A geometric box and a set of moving patches are defined to identify the portion of the mesh that must be deformed.

```cpp
    //Set the box including the step
    Eigen::MatrixXd Box(2, 3);
    Box << 1.98, 0.01, 0.11,
        7.02, -0.666669, -0.01;
    //Select the patches to be moved
    List<label> movPat;
    movPat.append(3);
    movPat.append(5);
    // Set the maximum iterations number for the offline phase
    example.maxIter = para->ITHACAdict->lookupOrDefault<int>("maxIter", 2000);

```

The offline phase is then executed to generate or load the training snapshots:

```cpp
    example.offlineSolve(Box, movPat);
    List<vector> points2Move;
    labelList boxIndices = ITHACAutilities::getIndicesFromBox(example._mesh(),
                           movPat, Box,
                           points2Move);
    example.linearMovePts((example.mu.maxCoeff() + example.mu.minCoeff()) / 2,
                          points2Move);

    for (int j = 0; j < boxIndices.size(); j++)
    {
        example.curX[boxIndices[j]] = points2Move[j];
    }

    Field<vector> curXV(example.curX);
    example._mesh().movePoints(curXV);
```

After the offline stage, POD modes are computed for velocity and pressure, and stored.

```cpp
    ITHACAPOD::getModes(example.Ufield, example.Umodes, example._U().name(),
                        example.podex, 0, 0,
                        example.NUmodesOut);
    ITHACAPOD::getModes(example.Pfield, example.Pmodes, example._p().name(),
                        example.podex, 0, 0,
                        example.NPmodesOut);
```
A second full-order problem is then used to generate reference solutions on the online parameter set, which will later be used for error analysis.

```cpp
    tutorial01cl checkOff(argc, argv);
    std::clock_t startOff;
    double durationOff;
    startOff = std::clock();

    if (!ITHACAutilities::check_folder("./ITHACAoutput/checkOff"))
    {
        checkOff.restart();
        ITHACAparameters* para = ITHACAparameters::getInstance(checkOff._mesh(),
                                 checkOff._runTime());
        checkOff.offline = false;
        checkOff.mu = angOn;
        checkOff.offlineSolve(Box, movPat, "./ITHACAoutput/checkOff/");
        checkOff.offline = true;
    }
```

If the case is turbulent, POD modes for turbulent viscosity are also computed, and `getTurbNN` is called to export the reduced datasets needed for neural-network training:

```cpp
    if (ITHACAutilities::isTurbulent())
    {
        ITHACAPOD::getModes(example.nutFields, example.nutModes, "nut",
                            example.podex, 0, 0, example.NNutModesOut);
        // Create the RBF for turbulence
        example.getTurbNN();
    }
```

The code then checks that the trained network file exists, loads it, and creates the reduced-order problem.

```cpp
    std::string wrng =
        "The network for this problem has not been calculated yet. Perform the Python training (see train.py).";
    M_Assert(ITHACAutilities::check_file("ITHACAoutput/NN/Net_" + name(
            example.NUmodes) + "_" + name(
                    example.NNutModes) + ".pt"), wrng.c_str());
    example.loadNet("ITHACAoutput/NN/Net_" + name(example.NUmodes) + "_" + name(
                        example.NNutModes) + ".pt");
    example._mesh().movePoints(example.point0);
    // Create the reduced object
    reducedSimpleSteadyNN reduced(example);
    ...
```

For each online parameter, the mesh is deformed consistently with the training stage, the reduced SIMPLE solver is executed, and the reduced solution is written to disk:

```cpp
   for (label k = 0; k < angOn.rows(); k++)
      {
        Info << "Solving online for parameter number " + name(k + 1) << endl;
        scalar mu_now = angOn(k, 0);
        example.restart();
        reduced.vel_now = vel;

        if (ITHACAutilities::isTurbulent())
        {
            example._mesh().movePoints(example.point0);
            List<vector> points2Move;
            labelList boxIndices = ITHACAutilities::getIndicesFromBox(example._mesh(),
                                   movPat, Box, points2Move);
            example.linearMovePts(mu_now, points2Move);

            for (int j = 0; j < boxIndices.size(); j++)
            {
                example.curX[boxIndices[j]] = points2Move[j];
            }

            Field<vector> curXV(example.curX);
            example._mesh().movePoints(curXV);
            ITHACAstream::writePoints(example._mesh().points(),
                                      "./ITHACAoutput/Reconstruct_" + name(example.NUmodes) + "_" + name(
                                          example.NPmodes) + "/", name(k + 1) + "/polyMesh/");
            reduced.solveOnline_Simple(mu_now, example.NUmodes, example.NPmodes,
                                       example.NNutModes, "./ITHACAoutput/Reconstruct_" + name(
                                           example.NUmodes) + "_" + name(example.NPmodes) + "/");
        }
        else
        {
            example._mesh().movePoints(example.point0);
            List<vector> points2Move;
            labelList boxIndices = ITHACAutilities::getIndicesFromBox(example._mesh(),
                                   movPat, Box, points2Move);
            example.linearMovePts(mu_now, points2Move);

            for (int j = 0; j < boxIndices.size(); j++)
            {
                example.curX[boxIndices[j]] = points2Move[j];
            }

            Field<vector> curXV(example.curX);
            example._mesh().movePoints(curXV);
            ITHACAstream::writePoints(example._mesh().points(),
                                      "./ITHACAoutput/Reconstruct_" + name(example.NUmodes) + "_" + name(
                                          example.NPmodes) + "/", name(k + 1) + "/polyMesh/");
            reduced.solveOnline_Simple(mu_now, example.NUmodes, example.NPmodes,
                                       example.NNutModes, "./ITHACAoutput/Reconstruct_" + name(
                                           example.NUmodes) + "_" + name(example.NPmodes) + "/");
        }
    }
```

In the final section, the code reads the full-order reference fields and the reconstructed reduced fields, computes relative errors for velocity and pressure, exports these error matrices, and prints the offline and online computational times.

```cpp
    PtrList<volVectorField> Ufull;
    PtrList<volScalarField> Pfull;
    PtrList<volVectorField> Ured;
    PtrList<volScalarField> Pred;
    volVectorField Uaux("Uaux", checkOff._U());
    volScalarField Paux("Paux", checkOff._p());
    ITHACAstream::readConvergedFields(Ufull, checkOff._U(),
                                      "./ITHACAoutput/checkOff/");
    ITHACAstream::readConvergedFields(Pfull, checkOff._p(),
                                      "./ITHACAoutput/checkOff/");
    ITHACAstream::read_fields(Ured, Uaux,
                              "./ITHACAoutput/Reconstruct_" + name(example.NUmodes) + "_" + name(
                                  example.NPmodes) + "/");
    ITHACAstream::read_fields(Pred, Paux,
                              "./ITHACAoutput/Reconstruct_" + name(example.NUmodes) + "_" + name(
                                  example.NPmodes) + "/");
    Eigen::MatrixXd relErrorU(Ufull.size(), 1);
    Eigen::MatrixXd relErrorP(Pfull.size(), 1);
    dimensionedVector U_fs("U_fs", dimVelocity, vector(1, 0, 0));

    for (label k = 0; k < Ufull.size(); k++)
    {
        volVectorField errorU = (Ufull[k] - Ured[k]).ref();
        volVectorField devU = (Ufull[k] - U_fs).ref();
        volScalarField errorP = (Pfull[k] - Pred[k]).ref();
        relErrorU(k, 0) = ITHACAutilities::frobNorm(errorU) /
                          ITHACAutilities::frobNorm(devU);
        relErrorP(k, 0) = ITHACAutilities::frobNorm(errorP) /
                          ITHACAutilities::frobNorm(Pfull[k]);
    }

    ITHACAstream::exportMatrix(relErrorU,
                               "errorU_" + name(example.NUmodes) + "_" + name(example.NPmodes) + "_" + name(
                                   example.NNutModes), "python", ".");
    ITHACAstream::exportMatrix(relErrorP,
                               "errorP_" + name(example.NUmodes) + "_" + name(example.NPmodes) + "_" + name(
                                   example.NNutModes), "python", ".");
    ...
```

## The plain code
The plain code is available [here](https://raw.githubusercontent.com/ITHACA-FV/ITHACA-FV/refs/heads/master/tutorials/NN/01simpleTurbGeomClosed/01simpleTurbGeomClosed.C).