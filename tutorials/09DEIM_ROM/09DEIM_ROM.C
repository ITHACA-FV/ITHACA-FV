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
#include <GenEigsSolver.h>
#include <Eigen/SparseLU>
#include "laplacianProblem.H"

class DEIM_function : public DEIM<PtrList<fvScalarMatrix>, volScalarField >
{
    public:
        using DEIM::DEIM;
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
            dimensionedScalar correct
            (
                "correct",
                dimensionSet(0, 1, 0, 0, 0, 0, 0),
                scalar(3.0)
            );
            fvScalarMatrix TiEqn22
            (
                fvm::laplacian(nu, T, "Gauss linear")
            );
            return TiEqn22;
        }

        Eigen::MatrixXd onlineCoeffsA(Eigen::MatrixXd mu)
        {
            Eigen::MatrixXd theta(fieldsA.size(), 1);

            for (int i = 0; i < fieldsA.size(); i++)
            {
                Eigen::SparseMatrix<double> Mr;
                Eigen::VectorXd br;
                fvScalarMatrix Aof = evaluate_expression(fieldsA[i], mu);
                Foam2Eigen::fvMatrix2Eigen(Aof, Mr, br);
                int ind_row = localMagicPointsA[i].first() + xyz_A[i].first() *
                              fieldsA[i].size();
                int ind_col = localMagicPointsA[i].second() + xyz_A[i].second() *
                              fieldsA[i].size();
                theta(i) = Mr.coeffRef(ind_row, ind_col);
            }

            return theta;
        }

        Eigen::MatrixXd onlineCoeffsB(Eigen::MatrixXd mu)
        {
            Eigen::MatrixXd theta(fieldsB.size(), 1);

            for (int i = 0; i < fieldsB.size(); i++)
            {
                Eigen::SparseMatrix<double> Mr;
                Eigen::VectorXd br;
                fvScalarMatrix Aof = evaluate_expression(fieldsB[i], mu);
                Foam2Eigen::fvMatrix2Eigen(Aof, Mr, br);
                int ind_row = localMagicPointsB[i] + xyz_B[i] * fieldsB[i].size();
                theta(i) = br(ind_row);
            }

            return theta;
        }
};

class DEIMlaplacian: public laplacianProblem
{
    public:
        explicit DEIMlaplacian(int argc, char* argv[])
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
                    Mlist.append(Teqn);
                    ITHACAstream::exportSolution(T, "./ITHACAoutput/", Folder, name(i + 1));
                    Tfield.append(T);
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
                ITHACAstream::exportSolution(T, "./ITHACAoutput/", Folder, name(i + 1));
                Tfield.append(T);
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
            DEIMmatrice->generateSubmeshesMatrix(2, mesh, T);
            DEIMmatrice->generateSubmeshesVector(2, mesh, T);
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
                Eigen::MatrixXd A = EigenFunctions::MVproduct(ReducedMatricesA, thetaonA);
                Eigen::VectorXd B = ReducedVectorsB[0];
                Eigen::VectorXd x = A.ldlt().solve(B);
                Eigen::VectorXd full = ModesTEig * x;
                t2 = std::chrono::high_resolution_clock::now();
                time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
                time_rom += time_span.count();
                // Export
                volScalarField Tred("Tred", T);
                Tred = Foam2Eigen::Eigen2field(Tred, full);
                ITHACAstream::exportSolution(Tred, "./ITHACAoutput/", Folder, name(i + 1));
                Tonline.append(Tred);
            }
        }
};

int main(int argc, char* argv[])
{
    // Construct the case
    DEIMlaplacian example(argc, argv);
    // Create the offline parameters for the solve
    example.mu = ITHACAutilities::rand(100, 2, -0.5, 0.5);
    // Solve the offline problem to compute the snapshots for the projections
    example.OfflineSolve(example.mu, "Offline");
    // Compute the POD modes
    ITHACAPOD::getModes(example.Tfield, example.Tmodes, example.podex, 0, 0, 20);
    // Compute the offline part of the DEIM procedure
    example.PODDEIM();
    // Construct a new set of parameters
    Eigen::MatrixXd par_new1 = ITHACAutilities::rand(100, 2, -0.5, 0.5);
    // Solve the online problem with the new parameters
    example.OnlineSolve(par_new1, "Online_red");
    // Solve a new full problem with the new parameters (necessary to compute speed up and error)
    DEIMlaplacian example_new(argc, argv);
    example_new.OnlineSolveFull(par_new1, "Online_full");
    // Output some infos
    std::cout << std::endl << "The FOM Solve took: " << example_new.time_full  <<
              " seconds." << std::endl;
    std::cout << std::endl << "The ROM Solve took: " << example.time_rom  <<
              " seconds." << std::endl;
    std::cout << std::endl << "The Speed-up is: " << example_new.time_full /
              example.time_rom  << std::endl << std::endl;
    Eigen::MatrixXd error = ITHACAutilities::error_listfields(example_new.Tfield,
                            example.Tonline);
    std::cout << "The mean L2 error is: " << error.mean() << std::endl;
    exit(0);
}

/// \dir 09DEIM_ROM Folder of the turorial 9
/// \file
/// \brief Implementation of tutorial 9 for an unsteady Navier-Stokes problem

/// \example 09DEIM_ROM.C
/// \section intro_09DEIM_ROM Introduction to tutorial 9
/// In this tutorial we implement test
///
/// The following image illustrates blabla
/// \image html cylinder.png
///
/// \section code09 A detailed look into the code
///
/// In this section we explain the main steps necessary to construct the tutorial N°9
///
/// \subsection header ITHACA-FV header files
///
/// First of all let's have a look at the header files that need to be included and what they are responsible for.
///
/// The header files of ITHACA-FV necessary for this tutorial are: <unsteadyNS.H> for the full order unsteady NS problem,
/// <ITHACAPOD.H> for the POD decomposition, <reducedUnsteadyNS.H> for the construction of the reduced order problem,
/// and finally <ITHACAstream.H> for some ITHACA input-output operations.
///
/// \section plaincode The plain program
/// Here there's the plain code
///


