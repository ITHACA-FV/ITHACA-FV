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
    Example of turbulent steady NS Reduction Problem solved by the use of the SIMPLE algorithm
SourceFiles
    019simpleTurbGeom.C
\*---------------------------------------------------------------------------*/

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


using namespace ITHACAtorch::torch2Eigen;

class SteadyNSSimpleNN : public SteadyNSSimple
{
    public:
        /// Constructor
        SteadyNSSimpleNN(int argc, char* argv[])
            :
            SteadyNSSimple(argc, argv)
        {};

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

        // This function creates the NN, trains it and saves the weights
        void getTurbNN()
        {
            Eigen::MatrixXd NNdata(2, 1);

            if (ITHACAutilities::check_file("./ITHACAoutput/Offline/NN/NNdata.txt"))
            {
                NNdata = ITHACAstream::readMatrix("./ITHACAoutput/Offline/NN/NNdata.txt");
            }
            // If the NN does not exist or input/output dimensions are changed, a new NN is constructed
            if (!ITHACAutilities::check_file("./ITHACAoutput/Offline/NN/Net.pt") ||
                    NNdata(0, 0) != NUmodes || NNdata(1, 0) != NNutModes)
            {
                // Since offline velocity and eddy viscosity are used to train the NN, they are read jyst in case
                // it has not been done during the offline phase
                if (nutFields.size() == 0 || Ufield.size() == 0)
                {
                    ITHACAstream::readMiddleFields(Ufield, _U(), "./ITHACAoutput/Offline/");
                    auto nut = _mesh().lookupObject<volScalarField>("nut");
                    ITHACAstream::readMiddleFields(nutFields, nut, "./ITHACAoutput/Offline/");
                }

                coeffL2Nut = ITHACAutilities::getCoeffs(nutFields, nutModes,
                                                        NNutModes, true);
                coeffL2Nut.transposeInPlace();
                coeffL2U = ITHACAutilities::getCoeffs(Ufield, Umodes,
                                                      NUmodes, true);
                coeffL2U.transposeInPlace();
                bias_inp  = coeffL2U.colwise().minCoeff();
                scale_inp = coeffL2U.colwise().maxCoeff() - coeffL2U.colwise().minCoeff();
                bias_out  = coeffL2Nut.colwise().minCoeff();
                scale_out = coeffL2Nut.colwise().maxCoeff() - coeffL2Nut.colwise().minCoeff();

                for (unsigned int i = 0; i < coeffL2U.rows(); i++)
                {
                    coeffL2U.row(i) = (coeffL2U.row(i).array() - bias_inp.array()).array()
                                      / scale_inp.array();
                }

                for (unsigned int i = 0; i < coeffL2Nut.rows(); i++)
                {
                    coeffL2Nut.row(i) = (coeffL2Nut.row(i).array() - bias_out.array()).array()
                                        / scale_out.array();
                }

                coeffL2U_tensor = eigenMatrix2torchTensor(coeffL2U);
                coeffL2Nut_tensor = eigenMatrix2torchTensor(coeffL2Nut);
                Net->push_back(torch::nn::Linear(NUmodes, 128));
                Net->push_back(torch::nn::ReLU());
                Net->push_back(torch::nn::Linear(128, 64));
                Net->push_back(torch::nn::ReLU());
                Net->push_back(torch::nn::Linear(64, NNutModes));
                optimizer = new torch::optim::Adam(Net->parameters(),
                                                   torch::optim::AdamOptions(2e-2));
                ITHACAparameters* para = ITHACAparameters::getInstance();
                int epochs = para->ITHACAdict->lookupOrDefault<int>("epochs", 20000);

                for (int64_t epoch = 1; epoch <= epochs; ++epoch)
                {
                    std::cerr << epoch << std::endl;
                    Net->zero_grad();
                    torch::Tensor y_r = Net->forward(coeffL2U_tensor);
                    torch::Tensor d_loss_real = torch::nn::functional::mse_loss(y_r,
                                                coeffL2Nut_tensor);
                    std::cout << d_loss_real.item<float>() << std::endl;
                    d_loss_real.backward();
                    optimizer->step();
                }

                // The NN weights are saved
                mkDir("./ITHACAoutput/Offline/NN");
                torch::save(Net, "./ITHACAoutput/Offline/NN/Net.pt");
                ITHACAstream::SaveDenseMatrix(bias_inp, "./ITHACAoutput/Offline/NN/",
                                              "bias_inp");
                ITHACAstream::SaveDenseMatrix(scale_inp, "./ITHACAoutput/Offline/NN/",
                                              "scale_inp");
                ITHACAstream::SaveDenseMatrix(bias_out, "./ITHACAoutput/Offline/NN/",
                                              "bias_out");
                ITHACAstream::SaveDenseMatrix(scale_out, "./ITHACAoutput/Offline/NN/",
                                              "scale_out");
                std::ofstream NNFile;
                NNFile.open ("./ITHACAoutput/Offline/NN/NNdata.txt");
                NNFile << std::to_string(NUmodes) + "\n";
                NNFile << std::to_string(NNutModes) + "\n";
                NNFile.close();
            }

            // If the NN already exists, all the weights are just read
            else
            {
                Net->push_back(torch::nn::Linear(NUmodes, 128));
                Net->push_back(torch::nn::ReLU());
                Net->push_back(torch::nn::Linear(128, 64));
                Net->push_back(torch::nn::ReLU());
                Net->push_back(torch::nn::Linear(64, NNutModes));
                torch::load(Net, "./ITHACAoutput/Offline/NN/Net.pt");
 
                ITHACAstream::ReadDenseMatrix(bias_inp, "./ITHACAoutput/Offline/NN/",
                                              "bias_inp");
                ITHACAstream::ReadDenseMatrix(scale_inp, "./ITHACAoutput/Offline/NN/",
                                              "scale_inp");
                ITHACAstream::ReadDenseMatrix(bias_out, "./ITHACAoutput/Offline/NN/",
                                              "bias_out");
                ITHACAstream::ReadDenseMatrix(scale_out, "./ITHACAoutput/Offline/NN/",
                                              "scale_out");
            }
        }

        // Function to eval the NN once the input is provided
        Eigen::MatrixXd evalNet(Eigen::MatrixXd a)
        {
            a.transposeInPlace();
            a = (a.array() - bias_inp.array()) / scale_inp.array();
            torch::Tensor aTensor = eigenMatrix2torchTensor(a);
            torch::Tensor out = Net->forward(aTensor);
            Eigen::MatrixXd g = torchTensor2eigenMatrix<double>(out);
            g = g.array() * scale_out.array() + bias_out.array();
            g.transposeInPlace();
            return g;
        }
};

class reducedSimpleSteadyNN : public reducedSimpleSteadyNS
{
    public:
        /// Constructor
        explicit reducedSimpleSteadyNN(SteadyNSSimpleNN& FOMproblem)
            :
            problem(&FOMproblem)
        {}

        /// Full problem.
        SteadyNSSimpleNN* problem;

        // Function to perform the online phase
        void solveOnline_Simple(scalar mu_now,
                                int NmodesUproj, int NmodesPproj, int NmodesNut = 0,
                                word Folder = "./ITHACAoutput/Reconstruct/")
        {
            counter++;

            // For all the variables, in case one wants to use all the available modes it is just necessary to set
            // the requested number into the ITHACAdict to zero
            if (NmodesUproj == 0)
            {
                NmodesUproj = problem->Umodes.size();
            }

            if (NmodesPproj == 0)
            {
                NmodesPproj = problem->Pmodes.size();
            }

            if (NmodesNut == 0)
            {
                NmodesNut = problem->nutModes.size();
            }

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
            Eigen::MatrixXd nutCoeff = Eigen::VectorXd::Zero(NmodesNut);
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
            problem->Umodes.reconstruct(U, a, "U");
            Info << ITHACAutilities::errorL2Rel(U, u2) << endl;
            problem->Pmodes.reconstruct(P, b, "p");
            problem->nutModes.reconstruct(nut, nutCoeff, "nut");
            phi = fvc::flux(U);
            int iter = 0;
            simpleControl& simple = problem->_simple();

            if (ITHACAutilities::isTurbulent())
            {
                Eigen::MatrixXd nutCoeff = problem->evalNet(a);
                problem->nutModes.reconstruct(nut, nutCoeff, "nut");
                ITHACAstream::exportSolution(nut, name(counter), Folder);
            }

            PtrList<volVectorField> gradModP;

            for (label i = 0; i < NmodesPproj; i++)
            {
                gradModP.append(fvc::grad(problem->Pmodes[i]));
            }

            projGradModP = (problem->Umodes).project(gradModP, NmodesUproj);

            // SIMPLE algorithm starts here
            while ((residual_jump > residualJumpLim
                    || std::max(U_norm_res, P_norm_res) > normalizedResidualLim) &&
                    iter < maxIterOn)
            {
                iter++;
                std::cout << iter << std::endl;
#if OFVER == 6
                simple.loop(runTime);
#else
                simple.loop();
#endif
                volScalarField nueff = problem->turbulence->nuEff();
                vector v(1, 0, 0);
                ITHACAutilities::assignBC(U, 0, v);

                // If the case is turbulent, then the network is evaluated
                if (ITHACAutilities::isTurbulent())
                {
                    Eigen::MatrixXd nutCoeff = problem->evalNet(a);
                    volScalarField& nut = const_cast<volScalarField&>
                                          (problem->_mesh().lookupObject<volScalarField>("nut"));
                    problem->nutModes.reconstruct(nut, nutCoeff, "nut");
                    ITHACAstream::exportSolution(nut, name(counter), Folder);
                }

                // Momentum equation
                fvVectorMatrix UEqn
                (
                    fvm::div(phi, U)
                    - fvm::laplacian(nueff, U)
                    - fvc::div(nueff * dev2(T(fvc::grad(U))))
                );
                UEqn.relax();
                List<Eigen::MatrixXd> RedLinSysU = problem->Umodes.project(UEqn, NmodesUproj);
                RedLinSysU[1] = RedLinSysU[1] - projGradModP * b;
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

                while (simple.correctNonOrthogonal())
                {
                    // Continuity equation
                    fvScalarMatrix pEqn
                    (
                        fvm::laplacian(rAtU(), P) == fvc::div(phiHbyA)
                    );
                    RedLinSysP = problem->Pmodes.project(pEqn, NmodesPproj);
                    b = reducedProblem::solveLinearSys(RedLinSysP, b, presidual);
                    problem->Pmodes.reconstruct(P, b, "p");

                    if (simple.finalNonOrthogonalIter())
                    {
                        phi = phiHbyA - pEqn.flux();
                    }
                }

                P.relax();
                U = HbyA - rAtU() * fvc::grad(P);
                U.correctBoundaryConditions();
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

                if (problem->para->debug)
                {
                    std::cout << "Residual jump = " << residual_jump << std::endl;
                    std::cout << "Normalized residual = " << std::max(U_norm_res,
                              P_norm_res) << std::endl;
                    std::cout << "Final normalized residual for velocity: " << U_norm_res <<
                              std::endl;
                    std::cout << "Final normalized residual for pressure: " << P_norm_res <<
                              std::endl;
                }
            }

            std::cout << "Solution " << counter << " converged in " << iter <<
                      " iterations." << std::endl;
            std::cout << "Final normalized residual for velocity: " << U_norm_res <<
                      std::endl;
            std::cout << "Final normalized residual for pressure: " << P_norm_res <<
                      std::endl;
            problem->Umodes.reconstruct(U, a, "Uaux");
            problem->Pmodes.reconstruct(P, b, "Paux");
            ITHACAstream::exportSolution(U, name(counter), Folder);
            ITHACAstream::exportSolution(P, name(counter), Folder);
            runTime.setTime(runTime.startTime(), 0);
        }
};

class tutorial19cl : public SteadyNSSimpleNN
{
    public:
        /// Constructor
        explicit tutorial19cl(int argc, char* argv[])
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

        /// Perform an Offline solve
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
                std::cout << "giusto" << std::endl;
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
                Vector<double> Uinl(1, 0, 0);

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
                    assignIF(U, Uinl);
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
                }
            }
        }

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

int main(int argc, char* argv[])
{
    // Construct the tutorial object
    tutorial19cl example(argc, argv);
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
        example.mu  = Eigen::VectorXd::LinSpaced(50, 0, 50);
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
        angOn = Eigen::VectorXd::LinSpaced(20, 1, 48);
        ITHACAstream::exportMatrix(angOn, "angOn", "eigen", "./");
    }

    // Set the inlet boundaries patch 0 directions x and y
    example.inletIndex.resize(1, 2);
    example.inletIndex(0, 0) = 0;
    example.inletIndex(0, 1) = 0;
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
    // Perform the offline solve
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

    // Perform POD on velocity and pressure and store the first 10 modes
    ITHACAPOD::getModes(example.Ufield, example.Umodes, example._U().name(),
                        example.podex, 0, 0,
                        example.NUmodesOut);
    ITHACAPOD::getModes(example.Pfield, example.Pmodes, example._p().name(),
                        example.podex, 0, 0,
                        example.NPmodesOut);

    if (ITHACAutilities::isTurbulent())
    {
        ITHACAPOD::getModes(example.nutFields, example.nutModes, "nut",
                            example.podex, 0, 0, example.NNutModesOut);
        // Create the RBF for turbulence
        example.getTurbNN();
    }

    PtrList<volVectorField> UfullTest;
    PtrList<volScalarField> PfullTest;
    PtrList<volVectorField> projUfullTest;
    PtrList<volScalarField> projPfullTest;
    volVectorField& U = example._U();
    volScalarField& p = example._p();
    ITHACAstream::read_fields(UfullTest, U, "./ITHACAoutput/Offline/1/");
    ITHACAstream::read_fields(PfullTest, p, "./ITHACAoutput/Offline/1/");
    example.Umodes.projectSnapshots(UfullTest, projUfullTest, example.NUmodes,
                                    "L2");
    example.Pmodes.projectSnapshots(PfullTest, projPfullTest, example.NPmodes,
                                    "L2");
    Eigen::MatrixXd projUErr = ITHACAutilities::errorL2Rel(UfullTest,
                               projUfullTest);
    Eigen::MatrixXd projPErr = ITHACAutilities::errorL2Rel(PfullTest,
                               projPfullTest);
    std::cout << "U errors: " << projUErr << std::endl;
    std::cout << "P errors: " << projPErr << std::endl;

    for (int i = 0; i < example.inletIndex.rows(); i++)
    {
        volVectorField ul(example._U());
        vector v(0, 0, 0);
        vector v0(0, 0, 0);
        v[example.inletIndex(i, 1)] = 1;
        ITHACAutilities::assignIF(ul, v0);
        ITHACAutilities::assignBC(ul, example.inletIndex(i, 0), v);

        for (int j = 0; j < example.Umodes.size(); j++)
        {
            ITHACAutilities::assignBC(example.Umodes[j], example.inletIndex(i, 0), v0);
        }

        for (int k = 0; k < ul.boundaryField().size(); k++)
        {
            if (k != example.inletIndex(i, 0))
            {
                ITHACAutilities::changeBCtype(ul, "fixedValue", k);
                ITHACAutilities::assignBC(ul, k, v0);
            }
        }

        example.liftfield.append(ul);
    }

    example._mesh().movePoints(example.point0);
    // Create the reduced object
    reducedSimpleSteadyNN reduced(example);
    PtrList<volVectorField> U_rec_list;
    PtrList<volScalarField> P_rec_list;
    // Reads inlet volocities boundary conditions.
    word vel_file(para->ITHACAdict->lookup("online_velocities"));
    Eigen::MatrixXd vel = ITHACAstream::readMatrix(vel_file);
    // Set the maximum iterations number for the online phase
    reduced.maxIterOn = para->ITHACAdict->lookupOrDefault<int>("maxIterOn", 2000);

    //Perform the online solutions
    for (label k = 0; k < angOn.rows(); k++)
    {
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
                                      "./ITHACAoutput/Reconstruct/", name(k + 1) + "/polyMesh/");
            reduced.solveOnline_Simple(mu_now, example.NUmodes, example.NPmodes,
                                       example.NNutModes);
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
                                      "./ITHACAoutput/Reconstruct/", name(k + 1) + "/polyMesh/");
            reduced.solveOnline_Simple(mu_now, example.NUmodes, example.NPmodes,
                                       example.NNutModes);
        }
    }

    // Error analysis
    tutorial19cl checkOff(argc, argv);

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
    ITHACAstream::read_fields(Ured, Uaux, "./ITHACAoutput/Reconstruct/");
    ITHACAstream::read_fields(Pred, Paux, "./ITHACAoutput/Reconstruct/");
    Eigen::MatrixXd relErrorU(Ufull.size(), 1);
    Eigen::MatrixXd relErrorP(Pfull.size(), 1);
    dimensionedVector U_fs("U_fs", dimVelocity, vector(1, 0, 0));

    for (label k = 0; k < Ufull.size(); k++)
    {
        volVectorField errorU = Ufull[k] - Ured[k];
        volVectorField devU = Ufull[k] - U_fs;
        volScalarField errorP = Pfull[k] - Pred[k];
        relErrorU(k, 0) = ITHACAutilities::frobNorm(errorU) /
                          ITHACAutilities::frobNorm(devU);
        relErrorP(k, 0) = ITHACAutilities::frobNorm(errorP) /
                          ITHACAutilities::frobNorm(Pfull[k]);
    }

    ITHACAstream::exportMatrix(relErrorU,
                               "errorU_" + name(example.NUmodes) + "_" + name(example.NPmodes), "python", ".");
    ITHACAstream::exportMatrix(relErrorP,
                               "errorP_" + name(example.NUmodes) + "_" + name(example.NPmodes), "python", ".");
    exit(0);
}
