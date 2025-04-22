/*---------------------------------------------------------------------------*\
     ██╗████████╗██╗  ██╗ █████╗  ██████╗ █████╗       ███████╗██╗   ██╗
     ██║╚══██╔══╝██║  ██║██╔══██╗██╔════╝██╔══██╗      ██╔════╝██║   ██║
     ██║   ██║   ███████║███████║██║     ███████║█████╗█████╗  ██║   ██║
     ██║   ██║   ██╔══██║██╔══██║██║     ██╔══██║╚════╝██╔══╝  ╚██╗ ██╔╝
     ██║   ██║   ██║  ██║██║  ██║╚██████╗██║  ██║      ██║      ╚████╔╝
     ╚═╝   ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝      ╚═╝       ╚═══╝

 * In real Time Highly Advanced Computational Applications for Finite Volumes
 * Copyright (C) 2017 by the ITHACA-FV authors
-------------------------------------------------------------------------------
License
    This file is part of ITHACA-FV
    ITHACA-FV is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    ITHACA-FV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Lesser General Public License for more details.
    You should have received a copy of the GNU Lesser General Public License
    along with ITHACA-FV. If not, see <http://www.gnu.org/licenses/>.
Description
    Example of steady NS Reduction Problem
SourceFiles
    03steadyNS.C
\*---------------------------------------------------------------------------*/

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

using namespace ITHACAtorch::torch2Eigen;

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
            // Net->push_back(torch::nn::Linear(NUmodes, 128));
            // Net->push_back(torch::nn::ReLU());
            // Net->push_back(torch::nn::Linear(128, 64));
            // Net->push_back(torch::nn::ReLU());
            // Net->push_back(torch::nn::Linear(64, NNutModes));
            // optimizer = new torch::optim::Adam(Net->parameters(),
            //                                    torch::optim::AdamOptions(2e-2));
        };

        label NUmodes;
        label NNutModes;

        ////////////////////////////////////
        // Eddy viscosity initialization //
        //////////////////////////////////

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
                std::cout << "Computing the coefficients for U train" << std::endl;
                Eigen::MatrixXd coeffL2U_train = ITHACAutilities::getCoeffs(UfieldTrain,
                                                 Umodes,
                                                 0, true);
                std::cout << "Computing the coefficients for p train" << std::endl;
                Eigen::MatrixXd coeffL2P_train = ITHACAutilities::getCoeffs(PfieldTrain,
                                                 Pmodes,
                                                 0, true);
                std::cout << "Computing the coefficients for nuT train" << std::endl;
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

        // // This function computes the coefficients which are later used for training
        // void getTurbNN()
        // {
        //     if (!ITHACAutilities::check_folder("ITHACAoutput/NN/coeffs"))
        //     {
        //         mkDir("ITHACAoutput/NN/coeffs");
        //         // Read Fields for Train
        //         PtrList<volVectorField> UfieldTrain;
        //         PtrList<volScalarField> PfieldTrain;
        //         PtrList<volScalarField> nutFieldsTrain;
        //         ITHACAstream::readMiddleFields(UfieldTrain, _U(),
        //                                        "./ITHACAoutput/Offline/");
        //         ITHACAstream::readMiddleFields(PfieldTrain, _p(),
        //                                        "./ITHACAoutput/Offline/");
        //         auto nut_train = _mesh().lookupObject<volScalarField>("nut");
        //         ITHACAstream::readMiddleFields(nutFieldsTrain, nut_train,
        //                                        "./ITHACAoutput/Offline/");
        //         /// Compute the coefficients for train
        //         std::cout << "Computing the coefficients for U train" << std::endl;
        //         Eigen::MatrixXd coeffL2U_train = ITHACAutilities::getCoeffs(UfieldTrain,
        //                                          Umodes,
        //                                          0, true);
        //         std::cout << "Computing the coefficients for p train" << std::endl;
        //         Eigen::MatrixXd coeffL2P_train = ITHACAutilities::getCoeffs(PfieldTrain,
        //                                          Pmodes,
        //                                          0, true);
        //         std::cout << "Computing the coefficients for nuT train" << std::endl;
        //         Eigen::MatrixXd coeffL2Nut_train = ITHACAutilities::getCoeffs(nutFieldsTrain,
        //                                            nutModes,
        //                                            0, true);
        //         coeffL2U_train.transposeInPlace();
        //         coeffL2P_train.transposeInPlace();
        //         coeffL2Nut_train.transposeInPlace();
        //         cnpy::save(coeffL2U_train, "ITHACAoutput/NN/coeffs/coeffL2UTrain.npy");
        //         cnpy::save(coeffL2P_train, "ITHACAoutput/NN/coeffs/coeffL2PTrain.npy");
        //         cnpy::save(coeffL2Nut_train, "ITHACAoutput/NN/coeffs/coeffL2NutTrain.npy");
        //         // Read Fields for Test
        //         PtrList<volVectorField> UfieldTest;
        //         PtrList<volScalarField> PfieldTest;
        //         PtrList<volScalarField> nutFieldsTest;
        //         /// Compute the coefficients for test
        //         ITHACAstream::readMiddleFields(UfieldTest, _U(),
        //                                        "./ITHACAoutput/checkOff/");
        //         ITHACAstream::readMiddleFields(PfieldTest, _p(),
        //                                        "./ITHACAoutput/checkOff/");
        //         auto nut_test = _mesh().lookupObject<volScalarField>("nut");
        //         ITHACAstream::readMiddleFields(nutFieldsTest, nut_test,
        //                                        "./ITHACAoutput/checkOff/");

        //         Info << "UfieldTest.size" << UfieldTest.size() << endl;
        //         //Info << "PfieldCheck.size" << PfieldCheck.size() << endl;
        //     //Info << "EfieldCheck.size" << EfieldCheck.size() << endl;
        //     //Info << "nutFieldsCheck.size" << nutFieldsCheck.size() << endl;
        //     //exit(0);

        //         // Compute the coefficients for test
        //         Eigen::MatrixXd coeffL2U_test = ITHACAutilities::getCoeffs(UfieldTest,
        //                                         Umodes,
        //                                         0, true);
        //         Eigen::MatrixXd coeffL2P_test = ITHACAutilities::getCoeffs(PfieldTest,
        //                                         Pmodes,
        //                                         0, true);
        //         Eigen::MatrixXd coeffL2Nut_test = ITHACAutilities::getCoeffs(nutFieldsTest,
        //                                           nutModes,
        //                                           0, true);
        //         coeffL2U_test.transposeInPlace();
        //         coeffL2P_test.transposeInPlace();
        //         coeffL2Nut_test.transposeInPlace();
        //         cnpy::save(coeffL2U_test, "ITHACAoutput/NN/coeffs/coeffL2UTest.npy");
        //         cnpy::save(coeffL2P_test, "ITHACAoutput/NN/coeffs/coeffL2PTest.npy");
        //         cnpy::save(coeffL2Nut_test, "ITHACAoutput/NN/coeffs/coeffL2NutTest.npy");
        //     }
        // }

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
};

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
            //volScalarField& P = problem->pThermo().p();
            volScalarField& P = problem->_p();
            volScalarField& E = problem->pThermo->he();
            /// nueff
            volScalarField nueff = problem->turbulence->nuEff();
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
            //Eigen::MatrixXd u = Eigen::MatrixXd::Zero(NmodesUproj, 1);
            //Eigen::MatrixXd e = Eigen::MatrixXd::Zero(NmodesEproj, 1);
            Eigen::MatrixXd e = ITHACAutilities::getCoeffs(E, problem->Emodes, NmodesEproj,
                                true);
            Eigen::MatrixXd u = ITHACAutilities::getCoeffs(U, problem->Umodes, NmodesUproj,
                                true);
            //Eigen::MatrixXd p = Eigen::MatrixXd::Zero(NmodesPproj, 1);
            Eigen::MatrixXd p = ITHACAutilities::getCoeffs(P, problem->Pmodes, NmodesPproj,
                                true);
            //Eigen::MatrixXd nutCoeff = ITHACAutilities::getCoeffs(nut, problem->nutModes, NmodesNutProj, true);
            //problem->nutModes.reconstruct(nut, nutCoeff, "nut");
            Eigen::MatrixXd  nutCoeff = problem->evalNet(u, mu_now);
            problem->nutModes.reconstruct(nut, nutCoeff, "nut");
            //vector Uinlet(170,0,0); // Vector for the inlet boundary condition
            label idInl =
                problem->_mesh().boundaryMesh().findPatchID("inlet"); // ID of the inlet patch
            vector Uinlet(problem->_U().boundaryFieldRef()[idInl][0][0], 0, 0);
            P.storePrevIter();

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
                //volVectorField gradpfull = -fvc::grad(P);
                //Eigen::MatrixXd projGrad = problem->Umodes.project(gradpfull, NmodesUproj);
                Eigen::MatrixXd projGradP = projGradModP * p;
                RedLinSysU[1] = RedLinSysU[1] - projGradP; // projGradP;
                u = reducedProblem::solveLinearSys(RedLinSysU, u,
                                                   uResidual);//, "fullPivLu");//"bdcSvd");
                problem->Umodes.reconstruct(U, u, "U");
                ITHACAutilities::assignBC(U, idInl, Uinlet);
                //solve(UEqnR == -problem->getGradP(P)); //For debug purposes only, second part only useful when using uEqn_global==-getGradP
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
                std::cout << "Ures = " << (uResidual.cwiseAbs()).sum() /
                          (RedLinSysU[1].cwiseAbs()).sum() << std::endl;
                std::cout << "Eres = " << (eResidual.cwiseAbs()).sum() /
                          (RedLinSysE[1].cwiseAbs()).sum() << std::endl;
                std::cout << "Pres = " << (pResidual.cwiseAbs()).sum() /
                          (RedLinSysP[1].cwiseAbs()).sum() << std::endl;
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

            //label k = 1;
            // U.rename("Ur");
            // P.rename("Pr");
            // E.rename("Er");
            // nut.rename("nutR");
            ITHACAstream::exportSolution(U, name(counter), Folder);
            ITHACAstream::exportSolution(P, name(counter), Folder);
            ITHACAstream::exportSolution(E, name(counter), Folder);
            ITHACAstream::exportSolution(nut, name(counter), Folder);
        }

};


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

                for (label i = 0; i < mu.rows(); i++)
                {
                    //std::clock_t startOff;
                    //startOff = std::clock();
                    std::cout << "Current mu = " << mu(i, 0) << std::endl;
                    changeViscosity(mu(i, 0));
                    assignIF(_U(), Uinl);
                    truthSolve(folder);
                    //durationOff = (std::clock() - startOff);
                    //cpuTimes << durationOff << std::endl;
                }
            }

            //cpuTimes.close();
        }

};

int main(int argc, char* argv[])
{
    // Construct the tutorial object
    tutorial03 example(argc, argv);
    //Info << example.pThermo().p() << endl;
    //Info << example._p() << endl;
    // exit(0);
    //std::clock_t startOff, startOn;
    //double durationOn, durationOff;
    std::cerr << "debug point 1" << std::endl;
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

    std::cerr << "debug point 2" << std::endl;
    // Read some parameters from file
    int NmodesUout = para->ITHACAdict->lookupOrDefault<int>("NmodesUout", 15);
    int NmodesPout = para->ITHACAdict->lookupOrDefault<int>("NmodesPout", 15);
    int NmodesEout = para->ITHACAdict->lookupOrDefault<int>("NmodesEout", 15);
    int NmodesNutOut = para->ITHACAdict->lookupOrDefault<int>("NmodesNutOut", 15);
    int NmodesUproj = para->ITHACAdict->lookupOrDefault<int>("NmodesUproj", 10);
    int NmodesPproj = para->ITHACAdict->lookupOrDefault<int>("NmodesPproj", 10);
    int NmodesEproj = para->ITHACAdict->lookupOrDefault<int>("NmodesEproj", 10);
    int NmodesNutProj = para->ITHACAdict->lookupOrDefault<int>("NmodesNutProj", 10);
    //Set the inlet boundaries patch 0 directions x and y
    // example.inletIndex.resize(1, 2);
    // example.inletIndex(0, 0) = 1;
    // example.inletIndex(0, 1) = 0;
    //Perform the offline solve
    std::cerr << "debug point 3" << std::endl;
    example.offlineSolve();
    std::cerr << "debug point 4" << std::endl;
    //Read the lift field
    // ITHACAstream::read_fields(example.liftfield, example._U(), "./lift/");
    // ITHACAutilities::normalizeFields(example.liftfield);
    // Homogenize the snapshots
    // example.computeLift(example.Ufield, example.liftfield, example.Uomfield);
    // Perform POD on velocity and pressure and store the first 10 modes
    // ITHACAPOD::getModes(example.Uomfield, example.Umodes, example._U().name(), example.podex, 0, 0,
    //                     NmodesUout);
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
    // Info << Foam::mag(example.Umodes[10] ) << endl;
    // Info << Foam::max(Foam::mag(example.Umodes[10] )) << endl;
    // exit(0);

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
    else //(ITHACAutilities::check_folder("./ITHACAoutput/checkOff"))
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
        // Info << "UfieldCheck.size" << UfieldCheck.size() << endl;
        // Info << "PfieldCheck.size" << PfieldCheck.size() << endl;
        // Info << "EfieldCheck.size" << EfieldCheck.size() << endl;
        // Info << "nutFieldsCheck.size" << nutFieldsCheck.size() << endl;
        // exit(0);
        //////////
        Eigen::MatrixXd snapsCheck =
            ITHACAstream::readMatrix("./ITHACAoutput/checkOff/snaps");
        label fieldNum = 0;

        for (label k = 0; k < snapsCheck.rows(); k++)
        {
            //Info << "snapsCheck(" <<  k <<  ",0)=" << snapsCheck(k,0)<< endl;
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
            //ITHACAutilities::createSymLink("./ITHACAoutput/checkOff/"+name(k+1)+"/polyMesh", "./ITHACAoutput/checkOffSingle/"+name(k+1)+"/");
        }
    }

    // Create the coefficients to train the net
    example.getTurbNN();
    //Before loading the net, it has to be created through the python script
    example.loadNet("ITHACAoutput/NN/Net_" + name(NmodesUproj) + "_" + name(
                        NmodesNutProj) + ".pt");
    // Create the reduced object
    ReducedCompressibleSteadyNN reduced(example);

    if (!ITHACAutilities::check_folder("./ITHACAoutput/Online_" + name(
                                           NmodesUproj) + "_" + name(NmodesNutProj) + "/" ) )
    {
        //Perform the online solutions
        std::ofstream cpuTimes;
        word OnlineFolder = "./ITHACAoutput/Online_" + name(NmodesUproj) + "_" + name(
                                NmodesNutProj) + "/";
        cpuTimes.open(OnlineFolder + "/cpuTimes", std::ios_base::app);

        for (label k = 0; k < parsOn.rows(); k++)
        {
            //scalar mu_now = parOn(k, 0);
            std::clock_t startOn;
            startOn = std::clock();
            double durationOn;
            Eigen::MatrixXd mu_now = parsOn.row(k);
            mu_now.transposeInPlace();
            example.changeViscosity(mu_now(0, 0));
            // reduced.setOnlineVelocity(vel);
            reduced.projectReducedOperators(NmodesUproj, NmodesPproj, NmodesEproj);
            //std::cout << "############################" << std::endl;
            example.restart();
            example.turbulence->validate();
            //std::cout << "##############################################################" << std::endl;
            // reduced.solveOnlineCompressible(mu_now, NmodesUproj, NmodesPproj, NmodesEproj);
            reduced.solveOnlineCompressible(NmodesUproj, NmodesPproj, NmodesEproj,
                                            NmodesNutProj, mu_now, "./ITHACAoutput/Online_" + name(NmodesUproj) + "_" +
                                            name(NmodesNutProj) + "/");
            durationOn = std::clock() - startOn;
            cpuTimes << durationOn << std::endl;
        }
    }
    //// Read the files
    else if (ITHACAutilities::check_folder("./ITHACAoutput/Online_" + name(
            NmodesUproj) + "_" + name(NmodesNutProj) + "/" ) )
    {
        PtrList<volVectorField> offlineU, onlineU;
        PtrList<volScalarField> offlineP, onlineP;
        PtrList<volScalarField> offlineE, onlineE;
        PtrList<volScalarField> offlineNut, onlineNut;
        //////
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
        //PtrList<volVectorField> Uerrflds;
        //PtrList<volScalarField> Perrflds,Eerrflds,NuTerrflds;

        for (label j = 0; j < parsOn.rows(); j++)
        {
            volVectorField Ue = offlineU[j] - onlineU[j];
            //auto err = foam2Eigen::foam2Eigen
            auto u = ITHACAutilities::L2Norm(offlineU[j]);
            Ue /= u;
            //std::cout << "Ue = " << ITHACAutilities::L2Norm(Ue)/u << std::endl;
            //////////
            volScalarField Pe = offlineP[j] - onlineP[j];
            auto p = ITHACAutilities::L2Norm(offlineP[j]);
            Pe /= p;
            //std::cout << "Pe = " <<  ITHACAutilities::L2Norm(Pe)/p << std::endl;
            ////////
            volScalarField Ee = offlineE[j] - onlineE[j];
            auto e = ITHACAutilities::L2Norm(offlineE[j]);
            Ee /= e;
            //std::cout << "Ee = " <<  ITHACAutilities::L2Norm(Ee)/e << std::endl;
            volScalarField Nute = offlineNut[j] - onlineNut[j];
            auto n = ITHACAutilities::L2Norm(offlineNut[j]);
            Nute /= n;
            //std::cout << "Nute = " <<  ITHACAutilities::L2Norm(Nute)/n << std::endl;
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

    // if(ITHACAutilities::check_folder("./ITHACAoutput/checkOff"))
    // {
    //         PtrList<volVectorField> UfieldCheck;
    //         PtrList<volScalarField> PfieldCheck;
    //         PtrList<volScalarField> EfieldCheck;
    //         PtrList<volScalarField> nutFieldsCheck;
    //         ITHACAstream::readMiddleFields(UfieldCheck, example._U(),"./ITHACAoutput/checkOff/");
    //         ITHACAstream::readMiddleFields(PfieldCheck, example._p(),"./ITHACAoutput/checkOff/");
    //         ITHACAstream::readMiddleFields(EfieldCheck, example._E(),"./ITHACAoutput/checkOff/");
    //         auto nutCheck = example._mesh().lookupObject<volScalarField>("nut");
    //         ITHACAstream::readMiddleFields(nutFieldsCheck, nutCheck, "./ITHACAoutput/checkOff/");
    //         // Info << "UfieldCheck.size" << UfieldCheck.size() << endl;
    //         // Info << "PfieldCheck.size" << PfieldCheck.size() << endl;
    //         // Info << "EfieldCheck.size" << EfieldCheck.size() << endl;
    //         // Info << "nutFieldsCheck.size" << nutFieldsCheck.size() << endl;
    //         // exit(0);
    //         //////////
    //         Eigen::MatrixXd snapsCheck = ITHACAstream::readMatrix("./ITHACAoutput/checkOff/snaps");
    //         label fieldNum = 0;
    //         for(label k=0; k<snapsCheck.rows(); k++)
    //         {
    //             Info << "snapsCheck(" <<  k <<  ",0)=" << snapsCheck(k,0)<< endl;
    //             fieldNum = fieldNum + snapsCheck(k,0);
    //             Info << "fieldNum" << fieldNum << endl;
    //             ITHACAstream::exportSolution(UfieldCheck[fieldNum-1], name(k+1), "./ITHACAoutput/checkOffSingle/");
    //             ITHACAstream::exportSolution(PfieldCheck[fieldNum-1], name(k+1), "./ITHACAoutput/checkOffSingle/");
    //             ITHACAstream::exportSolution(EfieldCheck[fieldNum-1], name(k+1), "./ITHACAoutput/checkOffSingle/");
    //             ITHACAstream::exportSolution(nutFieldsCheck[fieldNum-1], name(k+1), "./ITHACAoutput/checkOffSingle/");
    //             //ITHACAutilities::createSymLink("./ITHACAoutput/checkOff/"+name(k+1)+"/polyMesh", "./ITHACAoutput/checkOffSingle/"+name(k+1)+"/");
    //         }
    //         exit(0);
    //         ITHACAutilities::createSymLink("./0", "./ITHACAoutput/checkOffSingle/");
    //         ITHACAutilities::createSymLink("./system", "./ITHACAoutput/checkOffSingle/");
    //         ITHACAutilities::createSymLink("./constant", "./ITHACAoutput/checkOffSingle/");
    //         if(!ITHACAutilities::check_folder("./ITHACAoutput/Online_"+name(NmodesUproj)+"_"+name(NmodesNutProj)+"/" ) )
    //         {
    //             //Perform the online solutions
    //             std::ofstream cpuTimes;
    //             word OnlineFolder = "./ITHACAoutput/Online_" + name(NmodesUproj) + "_" + name(NmodesNutProj) + "/";
    //             cpuTimes.open(OnlineFolder + "/cpuTimes", std::ios_base::app);
    //             for (label k = 0; k < parsOn.rows(); k++)
    //             {
    //                 //scalar mu_now = parOn(k, 0);
    //                 std::clock_t startOn;
    //                 startOn = std::clock();
    //                 double durationOn;
    //                 Eigen::MatrixXd mu_now = parsOn.row(k);
    //                 mu_now.transposeInPlace();
    //                 example.changeViscosity(mu_now(0,0));
    //                 // reduced.setOnlineVelocity(vel);
    //                 reduced.projectReducedOperators(NmodesUproj, NmodesPproj, NmodesEproj);
    //                 //std::cout << "############################" << std::endl;
    //                 example.restart();
    //                 example.turbulence->validate();
    //                 std::cout << "##############################################################" << std::endl;
    //                 // reduced.solveOnlineCompressible(mu_now, NmodesUproj, NmodesPproj, NmodesEproj);
    //                 reduced.solveOnlineCompressible(mu_now, NmodesUproj, NmodesPproj, NmodesEproj, NmodesNutProj, "./ITHACAoutput/Online_" + name(NmodesUproj) + "_"+name(NmodesNutProj) + "/");
    //                 durationOn = std::clock() - startOn;
    //                 cpuTimes << durationOn << std::endl;
    //             }
    //         }
    //         //// Read the files
    //         else if (ITHACAutilities::check_folder("./ITHACAoutput/Online_" + name(NmodesUproj) + "_" + name(NmodesNutProj) + "/" ) )
    //         {
    //             PtrList<volVectorField> offlineU, onlineU;
    //             PtrList<volScalarField> offlineP, onlineP;
    //             PtrList<volScalarField> offlineE, onlineE;
    //             PtrList<volScalarField>offlineNut, onlineNut;
    //             //////
    //             ITHACAstream::read_fields(onlineU,example._U(),"./ITHACAoutput/Online_" + name(NmodesUproj) + "_" + name(NmodesNutProj) + "/");
    //             ITHACAstream::read_fields(onlineP,example._p(),"./ITHACAoutput/Online_" + name(NmodesUproj) + "_" + name(NmodesNutProj) + "/");
    //             ITHACAstream::read_fields(onlineE,example._E(),"./ITHACAoutput/Online_" + name(NmodesUproj) + "_" + name(NmodesNutProj) + "/");
    //             auto nut = example._mesh().lookupObject<volScalarField>("nut");
    //             /// Save the online fields
    //             ITHACAstream::read_fields(onlineNut,nut,"./ITHACAoutput/Online_" + name(NmodesUproj) + "_" + name(NmodesNutProj) + "/");
    //             /// Save the errors
    //             ITHACAstream::read_fields(offlineU,example._U(),"./ITHACAoutput/checkOffSingle/");
    //             ITHACAstream::read_fields(offlineP,example._p(),"./ITHACAoutput/checkOffSingle/");
    //             ITHACAstream::read_fields(offlineE,example._E(),"./ITHACAoutput/checkOffSingle/");
    //             ITHACAstream::read_fields(offlineNut,nut,       "./ITHACAoutput/checkOffSingle/");
    //             //PtrList<volVectorField> Uerrflds;
    //             //PtrList<volScalarField> Perrflds,Eerrflds,NuTerrflds;
    //             for(label j=0; j<parsOn.rows(); j++)
    //             {
    //                 volVectorField Ue = offlineU[j] - onlineU[j];
    //                 //auto err = foam2Eigen::foam2Eigen
    //                 volScalarField Pe = offlineP[j] - onlineP[j];
    //                 volScalarField Ee = offlineE[j] - onlineE[j];
    //                 volScalarField Nute = offlineNut[j] - onlineNut[j];
    //                 Ue.rename("Ue");
    //                 Pe.rename("Pe");
    //                 Ee.rename("Ee");
    //                 Nute.rename("Nute");
    //                 ITHACAstream::exportSolution(Ue,   name(j+1), "./ITHACAoutput/ErrorFields_" + name(NmodesUproj) + "_" + name(NmodesNutProj) + "/");
    //                 ITHACAstream::exportSolution(Pe,   name(j+1), "./ITHACAoutput/ErrorFields_" + name(NmodesUproj) + "_" + name(NmodesNutProj) + "/");
    //                 ITHACAstream::exportSolution(Ee,   name(j+1), "./ITHACAoutput/ErrorFields_" + name(NmodesUproj) + "_" + name(NmodesNutProj) + "/");
    //                 ITHACAstream::exportSolution(Nute, name(j+1), "./ITHACAoutput/ErrorFields_" + name(NmodesUproj) + "_" + name(NmodesNutProj) + "/");
    //             }
    //             Eigen::MatrixXd errorU = ITHACAutilities::errorL2Rel(offlineU, onlineU);
    //             Eigen::MatrixXd errorP = ITHACAutilities::errorL2Rel(offlineP, onlineP);
    //             Eigen::MatrixXd errorE = ITHACAutilities::errorL2Rel(offlineE, onlineE);
    //             Eigen::MatrixXd errorNut = ITHACAutilities::errorL2Rel(offlineNut, onlineNut);
    //             ///
    //             ITHACAstream::exportMatrix(errorU,"errorU" + name(NmodesUproj) + "_" + name(NmodesNutProj),     "python", "./ITHACAoutput/ErrorFields_" + name(NmodesUproj) + "_" + name(NmodesNutProj) + "/");
    //             ITHACAstream::exportMatrix(errorP,"errorP" + name(NmodesUproj) + "_" + name(NmodesNutProj),     "python", "./ITHACAoutput/ErrorFields_" + name(NmodesUproj) + "_" + name(NmodesNutProj) + "/");
    //             ITHACAstream::exportMatrix(errorE,"errorE" + name(NmodesUproj) + "_" + name(NmodesNutProj),     "python", "./ITHACAoutput/ErrorFields_" + name(NmodesUproj) + "_" + name(NmodesNutProj) + "/");
    //             ITHACAstream::exportMatrix(errorNut,"errorNut" + name(NmodesUproj) + "_" + name(NmodesNutProj), "python", "./ITHACAoutput/ErrorFields_" + name(NmodesUproj) + "_" + name(NmodesNutProj) + "/");
    //         }
    // }
    // else
    // {
    //     std::cerr << "CheckOff folder is missing, error analysis cannot be performed." << std::endl;
    // }
    exit(0);
}
