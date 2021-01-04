/*---------------------------------------------------------------------------*\
     ██╗████████╗██╗  ██╗ █████╗  ██████╗ █████╗       ███████╗██╗   ██╗
     ██║╚══██╔══╝██║  ██║██╔══██╗██╔════╝██╔══██╗      ██╔════╝██║   ██║
     ██║   ██║   ███████║███████║██║     ███████║█████╗█████╗  ██║   ██║
     ██║   ██║   ██╔══██║██╔══██║██║     ██╔══██║╚════╝██╔══╝  ╚██╗ ██╔╝
     ██║   ██║   ██║  ██║██║  ██║╚██████╗██║  ██║      ██║      ╚████╔╝
     ╚═╝   ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝      ╚═╝       ╚═══╝

 * In real Time Highly Advanced Computational Applications for Finite Volumes
 * Copyright (C) 2020 by the ITHACA-FV authors
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

\*---------------------------------------------------------------------------*/

/// \file
/// Source file of the ReducedUnsteadyNSExplicit class


#include "ReducedUnsteadyNSExplicit.H"


// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructor initialization
ReducedUnsteadyNSExplicit::ReducedUnsteadyNSExplicit()
{
}

ReducedUnsteadyNSExplicit::ReducedUnsteadyNSExplicit(UnsteadyNSExplicit&
        FOMproblem)
    :
    problem(&FOMproblem)
{
    N_BC = problem->inletIndex.rows();
    Nphi_u = problem->B_matrix.rows();
    Nphi_p = problem->K_matrix.cols();
}

// * * * * * * * * * * * * * Solve Functions  * * * * * * * * * * * //

void ReducedUnsteadyNSExplicit::solveOnline(Eigen::MatrixXd vel,
        label startSnap)
{
    if (problem->fluxMethod == "inconsistent")
    {
        // Create and resize the solution vectors
        Eigen::VectorXd a_o = Eigen::VectorXd::Zero(Nphi_u);
        Eigen::VectorXd a_n = a_o;
        Eigen::MatrixXd b = Eigen::VectorXd::Zero(Nphi_p);
        Eigen::MatrixXd x = Eigen::VectorXd::Zero(Nphi_p);
        Eigen::VectorXd presidual = Eigen::VectorXd::Zero(Nphi_p);
        Eigen::VectorXd RHS  = Eigen::VectorXd::Zero(Nphi_p);
        // Counting variable
        int counter = 0;
        // Set the initial time
        time = tstart;

        // Determine number of time steps
        while (time < finalTime - 0.5 * dt)
        {
            time = time + dt;
            counter ++;
        }

        // Set the initial time
        time = tstart;
        // Initial conditions / guesses
        a_o = ITHACAutilities::getCoeffs(problem->Ufield[0],
                                         problem->Umodes);
        b   = ITHACAutilities::getCoeffs(problem->Pfield[0],
                                         problem->Pmodes);
        // Set size of online solution
        online_solution.resize(counter + 1);
        // Create vector to store temporal solution and save initial condition as first solution
        Eigen::MatrixXd tmp_sol(Nphi_u + Nphi_p + 1, 1);
        tmp_sol(0) = time;
        tmp_sol.col(0).segment(1, Nphi_u) = a_o;
        tmp_sol.col(0).tail(b.rows()) = b;
        online_solution[0] = tmp_sol;

        for (label i = 1; i < online_solution.size(); i++)
        {
            time = time + dt;
            std::cout << " ################## time =   " << time <<
                      " ##################" << std::endl;
            // Pressure Poisson Equation
            // Diffusion Term
            Eigen::VectorXd M1 = problem->BP_matrix * a_o * nu ;
            // Convection Term
            Eigen::MatrixXd cf(1, 1);
            // Divergence term
            Eigen::MatrixXd M2 = problem->P_matrix * a_o;

            for (label l = 0; l < Nphi_p; l++)
            {
                cf = a_o.transpose() * Eigen::SliceFromTensor(problem->Cf_tensor, 0,
                        l) * a_o;
                RHS(l) = (1 / dt) * M2(l, 0) - cf(0, 0) + M1(l, 0);
            }

            List<Eigen::MatrixXd> RedLinSysP = problem->LinSysDiv;
            RedLinSysP[1] = (1 / dt) * problem->LinSysDiv[1] + nu * problem->LinSysDiff[1] +
                            problem->LinSysConv[1] + RHS;
            b = reducedProblem::solveLinearSys(RedLinSysP, x, presidual);
            // Momentum Equation
            // Convective term
            Eigen::MatrixXd cc(1, 1);
            // Diffusion Term
            Eigen::VectorXd M5 = problem->B_matrix * a_o * nu ;
            // Pressure Gradient Term
            Eigen::VectorXd M3 = problem->K_matrix * b;
            // Boundary Vector - Diffusion
            Eigen::VectorXd M4 = problem->RD_matrix * nu;
            // Boundary Vector - Convection
            Eigen::VectorXd M6 = problem->RC_matrix;

            for (label l = 0; l < Nphi_u; l++)
            {
                cc = a_o.transpose() * Eigen::SliceFromTensor(problem->C_tensor, 0,
                        l) * a_o;
                a_n(l) = a_o(l) + (M5(l) - cc(0, 0) - M3(l) + M4(l, 0) + M6(l, 0)) * dt;
            }

            tmp_sol(0) = time;
            tmp_sol.col(0).segment(1, Nphi_u) = a_n;
            tmp_sol.col(0).tail(b.rows()) = b;
            online_solution[i] = tmp_sol;
            a_o = a_n;
        }
    }
    else
    {
        std::cout << " Only the inconsistent flux method is implemented. " << std::endl;
        exit(0);
    }
}


void ReducedUnsteadyNSExplicit::reconstruct(bool exportFields, fileName folder)
{
    if (exportFields)
    {
        mkDir(folder);
        ITHACAutilities::createSymLink(folder);
    }

    int counter = 0;
    int nextwrite = 0;
    List < Eigen::MatrixXd> CoeffU;
    List < Eigen::MatrixXd> CoeffP;
    List <double> tValues;
    CoeffU.resize(0);
    CoeffP.resize(0);
    tValues.resize(0);
    int exportEveryIndex = round(exportEvery / storeEvery);

    for (label i = 0; i < online_solution.size(); i++)
    {
        if (counter == nextwrite)
        {
            Eigen::MatrixXd currentUCoeff;
            Eigen::MatrixXd currentPCoeff;
            currentUCoeff = online_solution[i].block(1, 0, Nphi_u, 1);
            currentPCoeff = online_solution[i].bottomRows(Nphi_p);
            CoeffU.append(currentUCoeff);
            CoeffP.append(currentPCoeff);
            nextwrite += exportEveryIndex;
            double timeNow = online_solution[i](0, 0);
            tValues.append(timeNow);
        }

        counter++;
    }

    volVectorField uRec("uRec", problem->Umodes[0]);
    volScalarField pRec("pRec", problem->Pmodes[0]);
    uRecFields = problem->Umodes.reconstruct(uRec, CoeffU, "uRec");
    pRecFields = problem->Pmodes.reconstruct(pRec, CoeffP, "pRec");

    if (exportFields)
    {
        ITHACAstream::exportFields(uRecFields, folder,
                                   "uRec");
        ITHACAstream::exportFields(pRecFields, folder,
                                   "pRec");
    }
}


//************************************************************************* //
