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
    Example of model reduction problem using the DEIM for a Heat Transfer Problem
SourceFiles
    09DEIM_ROM.C
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IOmanip.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include "ITHACAutilities.H"
#include "Foam2Eigen.H"
#include <chrono>
#include "fvMeshSubset.H"
#include "ITHACAstream.H"
#include "ITHACAPOD.H"
#include "EigenFunctions.H"
#include "DEIM.H"
#include <chrono>
#include <Eigen/SVD>
#include <Eigen/SparseLU>
#include "laplacianProblem.H"

class DEIM_function : public DEIM<fvScalarMatrix>
{
        using DEIM::DEIM;
    public:
        static fvScalarMatrix evaluate_expression(volScalarField& T, Eigen::MatrixXd mu)
        {
            volScalarField yPos = T.mesh().C().component(vector::Y);
            volScalarField xPos = T.mesh().C().component(vector::X);
            volScalarField nu(T);

            for (auto i = 0; i < nu.size(); i++)
            {
                nu[i] = std::exp( - 2 * std::pow(xPos[i] - mu(0) - 1,
                                                 2) - 2 * std::pow(yPos[i] - mu(1) - 0.5, 2)) + 1;
            }

            nu.correctBoundaryConditions();
            fvScalarMatrix TiEqn22
            (
                fvm::laplacian(nu, T, "Gauss linear")
            );
            return TiEqn22;
        }

        Eigen::MatrixXd onlineCoeffsA(Eigen::MatrixXd mu)
        {
            Eigen::MatrixXd theta(magicPointsAcol().size(), 1);
            fvScalarMatrix Aof = evaluate_expression(fieldA(), mu);
            Eigen::SparseMatrix<double> Mr;
            Eigen::VectorXd br;
            Foam2Eigen::fvMatrix2Eigen(Aof, Mr, br);

            for (int i = 0; i < magicPointsAcol().size(); i++)
            {
                int ind_row = localMagicPointsArow[i] + xyz_Arow()[i] *
                              fieldA().size();
                int ind_col = localMagicPointsAcol[i] + xyz_Acol()[i] *
                              fieldA().size();
                theta(i) = Mr.coeffRef(ind_row, ind_col);
            }

            return theta;
        }

        Eigen::MatrixXd onlineCoeffsB(Eigen::MatrixXd mu)
        {
            Eigen::MatrixXd theta(magicPointsB().size(), 1);
            fvScalarMatrix Aof = evaluate_expression(fieldB(), mu);
            Eigen::SparseMatrix<double> Mr;
            Eigen::VectorXd br;
            Foam2Eigen::fvMatrix2Eigen(Aof, Mr, br);

            for (int i = 0; i < magicPointsB().size(); i++)
            {
                int ind_row = localMagicPointsB[i] + xyz_B()[i] * fieldB().size();
                theta(i) = br(ind_row);
            }

            return theta;
        }

        PtrList<volScalarField> fieldsA;
        autoPtr<volScalarField> fieldA;
        autoPtr<volScalarField> fieldB;
        PtrList<volScalarField> fieldsB;
};

class DEIMLaplacian: public laplacianProblem
{
    public:
        explicit DEIMLaplacian(int argc, char* argv[])
            :
            laplacianProblem(argc, argv),
            nu(_nu()),
            S(_S()),
            T(_T())
        {
            fvMesh& mesh = _mesh();
            ITHACAdict = new IOdictionary
            (
                IOobject
                (
                    "ITHACAdict",
                    "./system",
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );
            NTmodes = readInt(ITHACAdict->lookup("N_modes_T"));
            NmodesDEIMA = readInt(ITHACAdict->lookup("N_modes_DEIM_A"));
            NmodesDEIMB = readInt(ITHACAdict->lookup("N_modes_DEIM_B"));
        }

        volScalarField& nu;
        volScalarField& S;
        volScalarField& T;


        DEIM_function* DEIMmatrice;
        PtrList<fvScalarMatrix> Mlist;
        Eigen::MatrixXd ModesTEig;
        std::vector<Eigen::MatrixXd> ReducedMatricesA;
        std::vector<Eigen::MatrixXd> ReducedVectorsB;

        int NTmodes;
        int NmodesDEIMA;
        int NmodesDEIMB;

        double time_full;
        double time_rom;

        void OfflineSolve(Eigen::MatrixXd par, word Folder)
        {
            if (offline)
            {
                ITHACAstream::read_fields(Tfield, T, "./ITHACAoutput/Offline/");
            }
            else
            {
                for (int i = 0; i < par.rows(); i++)
                {
                    fvScalarMatrix Teqn = DEIMmatrice->evaluate_expression(T, par.row(i));
                    Teqn.solve();
                    Mlist.append((Teqn).clone());
                    ITHACAstream::exportSolution(T, name(i + 1), "./ITHACAoutput/" + Folder);
                    Tfield.append((T).clone());
                }
            }
        };

        void OnlineSolveFull(Eigen::MatrixXd par, word Folder)
        {
            auto t1 = std::chrono::high_resolution_clock::now();
            auto t2 = std::chrono::high_resolution_clock::now();
            auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>
                             (t2 - t1);
            time_full = 0;

            for (int i = 0; i < par.rows(); i++)
            {
                fvScalarMatrix Teqn = DEIMmatrice->evaluate_expression(T, par.row(i));
                t1 = std::chrono::high_resolution_clock::now();
                Teqn.solve();
                t2 = std::chrono::high_resolution_clock::now();
                time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
                time_full += time_span.count();
                ITHACAstream::exportSolution(T, name(i + 1), "./ITHACAoutput/" + Folder);
                Tfield.append((T).clone());
            }
        };

        void PODDEIM()
        {
            PODDEIM(NTmodes, NmodesDEIMA, NmodesDEIMB);
        }

        void PODDEIM(int NmodesT, int NmodesDEIMA, int NmodesDEIMB)
        {
            DEIMmatrice = new DEIM_function(Mlist, NmodesDEIMA, NmodesDEIMB, "T_matrix");
            fvMesh& mesh  =  const_cast<fvMesh&>(T.mesh());
            // Differential Operator
            DEIMmatrice->fieldA = autoPtr<volScalarField>(new volScalarField(
                                      DEIMmatrice->generateSubmeshMatrix(2, mesh, T)));
            DEIMmatrice->fieldB = autoPtr<volScalarField>(new volScalarField(
                                      DEIMmatrice->generateSubmeshVector(2, mesh, T)));
            // Source Terms
            ModesTEig = Foam2Eigen::PtrList2Eigen(Tmodes);
            ModesTEig.conservativeResize(ModesTEig.rows(), NmodesT);
            ReducedMatricesA.resize(NmodesDEIMA);
            ReducedVectorsB.resize(NmodesDEIMB);

            for (int i = 0; i < NmodesDEIMA; i++)
            {
                ReducedMatricesA[i] = ModesTEig.transpose() * DEIMmatrice->MatrixOnlineA[i] *
                                      ModesTEig;
            }

            for (int i = 0; i < NmodesDEIMB; i++)
            {
                ReducedVectorsB[i] = ModesTEig.transpose() * DEIMmatrice->MatrixOnlineB;
            }
        };

        void OnlineSolve(Eigen::MatrixXd par_new, word Folder)
        {
            auto t1 = std::chrono::high_resolution_clock::now();
            auto t2 = std::chrono::high_resolution_clock::now();
            auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>
                             (t2 - t1);
            time_rom = 0;

            for (int i = 0; i < par_new.rows(); i++)
            {
                // solve
                t1 = std::chrono::high_resolution_clock::now();
                Eigen::MatrixXd thetaonA = DEIMmatrice->onlineCoeffsA(par_new.row(i));
                Eigen::MatrixXd thetaonB = DEIMmatrice->onlineCoeffsB(par_new.row(i));
                Eigen::MatrixXd A = EigenFunctions::MVproduct(ReducedMatricesA, thetaonA);
                Eigen::VectorXd B = EigenFunctions::MVproduct(ReducedVectorsB, thetaonB);
                Eigen::VectorXd x = A.fullPivLu().solve(B);
                Eigen::VectorXd full = ModesTEig * x;
                t2 = std::chrono::high_resolution_clock::now();
                time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
                time_rom += time_span.count();
                // Export
                volScalarField Tred("Tred", T);
                Tred = Foam2Eigen::Eigen2field(Tred, full);
                ITHACAstream::exportSolution(Tred, name(i + 1), "./ITHACAoutput/" + Folder);
                Tonline.append((Tred).clone());
            }
        }
};

int main(int argc, char* argv[])
{
    // Construct the case
    DEIMLaplacian example(argc, argv);
    // Read some parameters from file
    ITHACAparameters* para = ITHACAparameters::getInstance(example._mesh(),
                             example._runTime());
    // Create the offline parameters for the solve
    example.mu = ITHACAutilities::rand(100, 2, -0.5, 0.5);
    // Solve the offline problem to compute the snapshots for the projections
    example.OfflineSolve(example.mu, "Offline");
    // Compute the POD modes
    ITHACAPOD::getModes(example.Tfield, example.Tmodes, example._T().name(),
                        example.podex, 0, 0, 20);
    // Compute the offline part of the DEIM procedure
    example.PODDEIM();
    // Construct a new set of parameters
    Eigen::MatrixXd par_new1 = ITHACAutilities::rand(100, 2, -0.5, 0.5);
    // Solve the online problem with the new parameters
    example.OnlineSolve(par_new1, "Online_red");
    // Solve a new full problem with the new parameters (necessary to compute speed up and error)
    DEIMLaplacian example_new(argc, argv);
    example_new.OnlineSolveFull(par_new1, "Online_full");
    // Output some infos
    std::cout << std::endl << "The FOM Solve took: " << example_new.time_full  <<
              " seconds." << std::endl;
    std::cout << std::endl << "The ROM Solve took: " << example.time_rom  <<
              " seconds." << std::endl;
    std::cout << std::endl << "The Speed-up is: " << example_new.time_full /
              example.time_rom  << std::endl << std::endl;
    Eigen::MatrixXd error = ITHACAutilities::errorL2Abs(example_new.Tfield,
                            example.Tonline);
    std::cout << "The mean L2 error is: " << error.mean() << std::endl;
    exit(0);
}

/// \dir 09DEIM_ROM Folder of the turorial 9
/// \file
/// \brief Implementation of tutorial 9 which presents DEIM for a Heat Conduction Problem

/// \example 09DEIM_ROM.C
/// \section intro_09DEIM_ROM Introduction to tutorial 9
/// In this tutorial we implement a test where we use the Discrete Empirical Interpolation
/// Method for a case where we have a non-linear dependency with respect to the
/// input parameters.
///
/// The following image illustrates the computational domain which is the same as the previous example
/// \image html domain_deim.png
///
/// The physical problem is given by a heat transfer problem which is described by the Poisson equation:
///
/// \f[ \nabla \cdot (\nu \nabla T) = S \f]
///
/// The parametric diffusivity is described by a parametric Gaussian function:
/// \f[
/// \nu(\mathbf{x},\mathbf{\mu}) = e^{-2(x-\mu_x-1)^2 - 2(y-\mu_y-0.5)^2},
///  \f]
///
/// The problem is then discretized as:
///
/// \f[ A(\mu)T = b \f]
///
/// In this case, even if the problem is linear, due to non-linearity with respect to the
/// input parameter of the conductivity constant it is not possible to have an affine decomposition
/// of the discretized differential operator.
///
/// We seek therefore an approximate affine expansion of the differential operator of this type:
///
/// \f[ A(\mu) = \sum_{i = 1}^{N_D} \theta_i(\mu) A_i  \f]
///
/// using the Discrete Empirical Interpolation Method
///
/// \section plaincode The plain program
/// Here there's the plain code
///


