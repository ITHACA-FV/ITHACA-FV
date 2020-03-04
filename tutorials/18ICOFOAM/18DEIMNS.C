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

#include "SteadyNSSimple.H"
#include "ITHACAstream.H"
#include "ITHACAPOD.H"
#include "ReducedSimpleSteadyNS.H"
#include "forces.H"
#include "IOmanip.H"
#include "DEIM.H"


class DEIM_functionU : public DEIM<fvVectorMatrix>
{
    public:
        using DEIM::DEIM;



        PtrList<volVectorField> fieldsA;
        PtrList<volVectorField> fieldsB;
        PtrList<surfaceScalarField> phi;
};

class DEIM_functionP : public DEIM<fvScalarMatrix>
{
    public:
        using DEIM::DEIM;

};

class tutorial18 : public SteadyNSSimple
{
    public:
        /// Constructor
        explicit tutorial18(int argc, char* argv[])
            :
            SteadyNSSimple(argc, argv),
            U(_U()),
            p(_p()),
            phi(_phi())
        {}

        /// Velocity field
        volVectorField& U;
        /// Pressure field
        volScalarField& p;
        ///
        surfaceScalarField& phi;



        /// Perform an Offline solve
        void offlineSolve()
        {
            Vector<double> inl(0, 0, 0);
            List<scalar> mu_now(1);

            // if the offline solution is already performed read the fields
            if (offline)
            {
                ITHACAstream::read_fields(Ufield, U, "./ITHACAoutput/Offline/");
                ITHACAstream::read_fields(Pfield, p, "./ITHACAoutput/Offline/");
                mu_samples =
                    ITHACAstream::readMatrix("./ITHACAoutput/Offline/mu_samples_mat.txt");
            }

            else
            {
                label BCind = 0;

                for (label i = 0; i < mu.rows(); i++)
                {
                    mu_now[0] = mu(i, 0);
                    change_viscosity(mu(i, 0));
                    truthSolve2(mu_now);
                }
            }
        }



};

class tutorial18red : public reducedSimpleSteadyNS
{
    public:

        tutorial18red(SteadyNSSimple& FOMproblem)
            :
            reducedSimpleSteadyNS(FOMproblem)
        {}

        DEIM_functionU* DEIMmatriceU;
        DEIM_functionP* DEIMmatriceP;
        std::vector<Eigen::MatrixXd> ReducedMatricesA_U;
        Eigen::MatrixXd ReducedVectorsB_U;
        std::vector<Eigen::MatrixXd> ReducedMatricesA_p;
        Eigen::MatrixXd ReducedVectorsB_p;
        PtrList<volScalarField> pFieldsA_U;
        PtrList<volVectorField> UFieldsA_U;
        PtrList<surfaceScalarField> phiFieldsA_U;
        PtrList<volScalarField> pFieldsB_U;
        PtrList<volVectorField> UFieldsB_U;
        PtrList<surfaceScalarField> phiFieldsB_U;
        PtrList<volScalarField> pFieldsA_p;
        PtrList<volVectorField> UFieldsA_p;
        PtrList<surfaceScalarField> phiFieldsA_p;
        PtrList<volScalarField> pFieldsB_p;
        PtrList<volVectorField> UFieldsB_p;
        PtrList<surfaceScalarField> phiFieldsB_p;

        void PODDEIM(int NmodesU, int NmodesP, int NmodesDEIMAU, int NmodesDEIMBU,
                     int NmodesDEIMAP,
                     int NmodesDEIMBP)
        {
            problem->Pmodes.toEigen();
            fvMesh& mesh  =  const_cast<fvMesh&>(problem->_U().mesh());
            DEIMmatriceU = new DEIM_functionU(problem->Ulist, NmodesDEIMAU, NmodesDEIMBU,
                                              "U_matrix");
            DEIMmatriceP = new DEIM_functionP(problem->Plist, NmodesDEIMAP, NmodesDEIMBP,
                                              "P_matrix");
            DEIMmatriceU->generateSubmeshesMatrix(1, mesh,
                                                  problem->_U());
            DEIMmatriceU->generateSubmeshesVector(1, mesh,
                                                  problem->_U());
            DEIMmatriceP->generateSubmeshesMatrix(1, mesh,
                                                  problem->_p());
            DEIMmatriceP->generateSubmeshesVector(1, mesh,
                                                  problem->_p());
            // U Eqn
            pFieldsA_U = DEIMmatriceU->generateSubFieldsMatrix(problem->_p());
            UFieldsA_U = DEIMmatriceU->generateSubFieldsMatrix(problem->_U());
            phiFieldsA_U = DEIMmatriceU->generateSubFieldsMatrix(problem->_phi());
            pFieldsB_U = DEIMmatriceU->generateSubFieldsVector(problem->_p());
            UFieldsB_U = DEIMmatriceU->generateSubFieldsVector(problem->_U());
            phiFieldsB_U = DEIMmatriceU->generateSubFieldsVector(problem->_phi());
            // P Eqn
            pFieldsA_p = DEIMmatriceP->generateSubFieldsMatrix(problem->_p());
            UFieldsA_p = DEIMmatriceP->generateSubFieldsMatrix(problem->_U());
            phiFieldsA_p = DEIMmatriceP->generateSubFieldsMatrix(problem->_phi());
            pFieldsB_p = DEIMmatriceP->generateSubFieldsVector(problem->_p());
            UFieldsB_p = DEIMmatriceP->generateSubFieldsVector(problem->_U());
            phiFieldsB_p = DEIMmatriceP->generateSubFieldsVector(problem->_phi());
            // reduced operators
            ReducedMatricesA_U.resize(NmodesDEIMAU);
            ReducedMatricesA_p.resize(NmodesDEIMAP);
            

            // U
            for (int i = 0; i < NmodesDEIMAU; i++)
            {
                ReducedMatricesA_U[i] = ULmodes.EigenModes[0].leftCols(NmodesU).transpose() *
                                        DEIMmatriceU->MatrixOnlineA[i] *
                                        ULmodes.EigenModes[0].leftCols(NmodesU);
            }

            ReducedVectorsB_U = ULmodes.EigenModes[0].leftCols(NmodesU).transpose() *
                                DEIMmatriceU->MatrixOnlineB;

            // P
            for (int i = 0; i < NmodesDEIMAP; i++)
            {
                ReducedMatricesA_p[i] = problem->Pmodes.EigenModes[0].leftCols(
                                            NmodesP).transpose() *
                                        DEIMmatriceP->MatrixOnlineA[i] *
                                        problem->Pmodes.EigenModes[0].leftCols(NmodesP);
            }

            ReducedVectorsB_p = problem->Pmodes.EigenModes[0].leftCols(
                                    NmodesP).transpose() *
                                DEIMmatriceP->MatrixOnlineB;
        };

        Eigen::MatrixXd onlineCoeffsAU(volVectorField& U, surfaceScalarField& phi,
                                       volScalarField& A, volVectorField& H)
        {
            Eigen::MatrixXd theta(UFieldsA_U.size(), 1);
            volScalarField nueff = problem->_laminarTransport().nu();
            fvVectorMatrix ueqn
            (
                fvm::div(phi, U) - fvm::laplacian(nueff,
                                                  U) - fvc::div(nueff * dev2(T(fvc::grad(U))))
            );
            ueqn.relax();
            A = ueqn.A();
            H = ueqn.H();
            Eigen::SparseMatrix<double> Mr;
            Eigen::VectorXd br;
            Foam2Eigen::fvMatrix2Eigen(ueqn, Mr, br);

            for (int i = 0; i < UFieldsA_U.size(); i++)
            {
                int ind_row = DEIMmatriceU->magicPointsA[i].first() +
                              DEIMmatriceU->xyz_A[i].first() *
                              problem->_U().size();
                int ind_col = DEIMmatriceU->magicPointsA[i].second() +
                              DEIMmatriceU->xyz_A[i].second() *
                              problem->_U().size();
                theta(i) = Mr.coeffRef(ind_row, ind_col);
            }

            return theta;
        }

        Eigen::MatrixXd onlineCoeffsBU(volVectorField& U, surfaceScalarField& phi)
        {
            Eigen::MatrixXd theta(UFieldsB_U.size(), 1);
            volScalarField nueff = problem->_laminarTransport().nu();
            fvVectorMatrix ueqn
            (
                fvm::div(phi, U) - fvm::laplacian(nueff,
                                                  U) - fvc::div(nueff * dev2(T(fvc::grad(U))))
            );
            ueqn.relax();
            Eigen::SparseMatrix<double> Mr;
            Eigen::VectorXd br;
            Foam2Eigen::fvMatrix2Eigen(ueqn, Mr, br);

            for (int i = 0; i < UFieldsB_U.size(); i++)
            {
                int ind_row = DEIMmatriceU->magicPointsB[i] +
                              DEIMmatriceU->xyz_B[i] *
                              problem->_U().size();
                theta(i) = br(ind_row);
            }

            return theta;
        }



        Eigen::MatrixXd onlineCoeffsAP(volVectorField& U, volScalarField& P,
                                       volScalarField& A, volVectorField& H)
        {
            volVectorField HbyA(constrainHbyA(1.0 / A * H, U, P));
            surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
            Eigen::MatrixXd theta(UFieldsA_p.size(), 1);
            fvScalarMatrix pEqn
            (
                fvm::laplacian(1 / A, P) == fvc::div(phiHbyA)
            );
            Eigen::SparseMatrix<double> Mr;
            Eigen::VectorXd br;
            Foam2Eigen::fvMatrix2Eigen(pEqn, Mr, br);

            for (int i = 0; i < pFieldsA_p.size(); i++)
            {
                int ind_row = DEIMmatriceP->magicPointsA[i].first() +
                              DEIMmatriceP->xyz_A[i].first() *
                              problem->_p().size();
                int ind_col = DEIMmatriceP->magicPointsA[i].second() +
                              DEIMmatriceP->xyz_A[i].second() *
                              problem->_p().size();
                theta(i) = Mr.coeffRef(ind_row, ind_col);
            }

            return theta;
        }

        Eigen::MatrixXd onlineCoeffsBP(volVectorField& U, volScalarField& P,
                                       volScalarField& A, volVectorField& H)
        {
            volVectorField HbyA(constrainHbyA(1.0 / A * H, U, P));
            surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
            Eigen::MatrixXd theta(pFieldsB_p.size(), 1);
            fvScalarMatrix pEqn
            (
                fvm::laplacian(1 / A, P) == fvc::div(phiHbyA)
            );
            Eigen::SparseMatrix<double> Mr;
            Eigen::VectorXd br;
            Foam2Eigen::fvMatrix2Eigen(pEqn, Mr, br);

            for (int i = 0; i < pFieldsB_p.size(); i++)
            {
                int ind_row = DEIMmatriceP->magicPointsB[i] +
                              DEIMmatriceP->xyz_B[i] *
                              problem->_p().size();
                theta(i) = br(ind_row);
            }

            return theta;
        }



        fvVectorMatrix Ueqn(volVectorField& U, surfaceScalarField& phi)
        {
            volScalarField nueff = problem->_laminarTransport().nu();
            fvVectorMatrix ueqn
            (
                fvm::div(phi, U)
            );
            return ueqn;
        }

        void solveOnline_Simple(scalar mu_now,
                                int NmodesUproj, int NmodesPproj, int NmodesSup, word Folder)
        {
            ULmodes.resize(0);

            for (int i = 0; i < problem->inletIndex.rows(); i++)
            {
                ULmodes.append(problem->liftfield[i]);
            }

            for (int i = 0; i < NmodesUproj; i++)
            {
                ULmodes.append(problem->Umodes.toPtrList()[i]);
            }

            for (int i = 0; i < NmodesSup; i++)
            {
                ULmodes.append(problem->supmodes.toPtrList()[i]);
            }

            counter++;
            scalar UprojN;
            scalar PprojN;

            if (NmodesUproj == 0)
            {
                UprojN = ULmodes.size();
            }

            else
            {
                UprojN = NmodesUproj + NmodesSup;
            }

            if (NmodesPproj == 0)
            {
                PprojN = problem->Pmodes.size();
            }

            else
            {
                PprojN = NmodesPproj;
            }

            Eigen::VectorXd uresidualOld = Eigen::VectorXd::Zero(UprojN);
            Eigen::VectorXd presidualOld = Eigen::VectorXd::Zero(PprojN);
            Eigen::VectorXd uresidual = Eigen::VectorXd::Zero(UprojN);
            Eigen::VectorXd presidual = Eigen::VectorXd::Zero(PprojN);
            scalar U_norm_res(1);
            scalar P_norm_res(1);
            Eigen::MatrixXd a = Eigen::VectorXd::Zero(UprojN);
            Eigen::MatrixXd b = Eigen::VectorXd::Zero(PprojN);
            a(0) = vel_now(0, 0);
            ITHACAparameters para;
            float residualJumpLim =
                para.ITHACAdict->lookupOrDefault<float>("residualJumpLim", 1e-5);
            float normalizedResidualLim =
                para.ITHACAdict->lookupOrDefault<float>("normalizedResidualLim", 1e-5);
            scalar residual_jump(1 + residualJumpLim);
            problem->restart();
            volScalarField& P = problem->_p();
            volVectorField& U = problem->_U();
            P.rename("p");
            surfaceScalarField& phi = problem->_phi();
            phi = fvc::interpolate(U) & U.mesh().Sf();
            int iter = 0;
            simpleControl& simple = problem->_simple();

            while (residual_jump > residualJumpLim
                    || std::max(U_norm_res, P_norm_res) > normalizedResidualLim)
            {
                iter++;
                P.storePrevIter();
                volScalarField nu = problem->_laminarTransport().nu();
                fvVectorMatrix UEqn1
                (
                    fvm::div(phi, U) - fvm::laplacian(nu, U) - fvc::div(nu * dev2(T(fvc::grad(U))))
                );
                UEqn1.relax();
                List<Eigen::MatrixXd> RedLinSysU = ULmodes.project(UEqn1, UprojN);
                volVectorField H = UEqn1.H();
                volScalarField A = UEqn1.A();
                Eigen::MatrixXd thetaonAU = onlineCoeffsAU(U, phi, A, H);
                Eigen::MatrixXd thetaonBU = onlineCoeffsBU(U, phi);
                RedLinSysU[0] = EigenFunctions::MVproduct(ReducedMatricesA_U, thetaonAU);
                RedLinSysU[1] = ReducedVectorsB_U*thetaonBU;
                //RedLinSysU[1] = ReducedVectorsB_U*thetaonBU;
                // for (int i = 0; i < ReducedOperators.size(); i++)
                // {
                //     RedLinSysU[0] += onlineViscosity * ReducedOperators[i][0];
                //     RedLinSysU[1] += onlineViscosity * ReducedOperators[i][1];
                // }
                //RedLinSysU[0] += onlineViscosity * redDivDev2;
                RedLinSysU[1] += redGradP * b;
                a = reducedProblem::solveLinearSys(RedLinSysU, a, uresidual, vel_now);
                ULmodes.reconstruct(U, a, "U");
                volVectorField HbyA(constrainHbyA(1.0 / A * H, U, P));
                surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
                List<Eigen::MatrixXd> RedLinSysP;

                while (simple.correctNonOrthogonal())
                {
                    Eigen::MatrixXd thetaonAP = onlineCoeffsAP(U, P, A, H);
                    Eigen::MatrixXd thetaonBP = onlineCoeffsBP(U, P, A, H);
                    fvScalarMatrix pEqn
                    (
                        fvm::laplacian(1 / A, P) == fvc::div(phiHbyA)
                    );
                    RedLinSysP = problem->Pmodes.project(pEqn, PprojN);
                    // RedLinSysP[0] = EigenFunctions::MVproduct(ReducedMatricesA_p, thetaonAP);
                    // RedLinSysP[1] = ReducedVectorsB_p*thetaonBP;
                    b = reducedProblem::solveLinearSys(RedLinSysP, b, presidual);
                    problem->Pmodes.reconstruct(P, b, "p");

                    if (simple.finalNonOrthogonalIter())
                    {
                        phi = phiHbyA - pEqn.flux();
                    }
                }

                P.relax();
                U = HbyA - 1.0 / A * fvc::grad(P);
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

                if (para.debug)
                {
                    std::cout << "Residual jump = " << residual_jump << std::endl;
                    std::cout << "Normalized residual = " << std::max(U_norm_res,
                              P_norm_res) << std::endl;
                }
            }

            std::cout << "Solution " << counter << " converged in " << iter <<
                      " iterations." << std::endl;
            std::cout << "Final normalized residual for velocity: " << U_norm_res <<
                      std::endl;
            std::cout << "Final normalized residual for pressure: " << P_norm_res <<
                      std::endl;
            ULmodes.reconstruct(U, a, "Uaux");
            P.rename("Paux");
            problem->Pmodes.reconstruct(P, b, "Paux");
            ITHACAstream::exportSolution(U, name(counter), Folder);
            ITHACAstream::exportSolution(P, name(counter), Folder);
        }

};

int main(int argc, char* argv[])
{
    // Construct the tutorial object
    tutorial18 example(argc, argv);
    // Read some parameters from file
    ITHACAparameters para;
    int NmodesUout = para.ITHACAdict->lookupOrDefault<int>("NmodesUout", 15);
    int NmodesPout = para.ITHACAdict->lookupOrDefault<int>("NmodesPout", 15);
    int NmodesUproj = para.ITHACAdict->lookupOrDefault<int>("NmodesUproj", 10);
    int NmodesPproj = para.ITHACAdict->lookupOrDefault<int>("NmodesPproj", 10);
    // Read the par file where the parameters are stored
    example.mu = Eigen::VectorXd::LinSpaced(100, 0.01, 0.001);
    // Set the inlet boundaries patch 0 directions x and y
    example.inletIndex.resize(1, 2);
    example.inletIndex(0, 0) = 0;
    example.inletIndex(0, 1) = 0;
    // Perform the offline solve
    example.offlineSolve();
    //example.liftSolve();
    example.liftfield.append(ITHACAutilities::computeAverage(example.Ufield));
    ITHACAutilities::normalizeFields(example.liftfield);
    example.liftfield[0].rename("Ulift0");
    //ITHACAutilities::changeBCtype(example.liftfield[0],"fixedValue",0);
    ITHACAstream::exportSolution(example.liftfield[0], "0", "./");
    // Homogenize the snapshots
    example.computeLift(example.Ufield, example.liftfield, example.Uomfield);
    // Perform POD on velocity and pressure and store the first 10 modes
    ITHACAPOD::getModes(example.Uomfield, example.Umodes, example.podex, 0, 0,
                        NmodesUout);
    ITHACAPOD::getModes(example.Pfield, example.Pmodes, example.podex, 0, 0,
                        NmodesPout);
    // Create the reduced object
    tutorial18red reduced(example);
    // Compute the offline part of the DEIM procedure
    reduced.PODDEIM(10, 10, 20, 1, 20, 20);
    reduced.project(10, 10);
    //reduced.project(NmodesUproj, NmodesPproj);
    PtrList<volVectorField> U_rec_list;
    PtrList<volScalarField> P_rec_list;
    Eigen::MatrixXd vel(1, 1);
    vel(0, 0) = 1.0;

    //Perform the online solutions
    for (label k = 0; k < (example.mu).size(); k++)
    {
        scalar mu_now = example.mu(k);
        example.restart();
        example.change_viscosity(mu_now);
        reduced.onlineViscosity = mu_now;
        reduced.setOnlineVelocity(vel);
        reduced.solveOnline_Simple(mu_now, 10, 10, 0,
                                   "./ITHACAoutput/Offline/");
    }

    exit(0);
}