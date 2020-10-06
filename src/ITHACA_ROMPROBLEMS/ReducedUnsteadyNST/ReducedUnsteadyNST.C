/*---------------------------------------------------------------------------*\
v/*---------------------------------------------------------------------------*\
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

\*---------------------------------------------------------------------------*/

/// \file
/// Source file of the reducedUnsteadyNST class


#include "ReducedUnsteadyNST.H"


// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructor
reducedUnsteadyNST::reducedUnsteadyNST()
{
}

reducedUnsteadyNST::reducedUnsteadyNST(unsteadyNST& FOMproblem)
//problem(&FOMproblem)
{
    problem   = &FOMproblem;
    N_BC      = problem->inletIndex.rows();
    N_BC_t    = problem->inletIndexT.rows();
    Nphi_u    = problem->B_matrix.rows();
    Nphi_p    = problem->K_matrix.cols();
    Nphi_t    = problem->Y_matrix.rows();

    // Create locally the velocity modes
    for (label k = 0; k < problem->liftfield.size(); k++)
    {
        Umodes.append(problem->liftfield[k]);
    }

    for (label k = 0; k < problem->NUmodes; k++)
    {
        Umodes.append(problem->Umodes[k]);
    }

    for (label k = 0; k < problem->NSUPmodes; k++)
    {
        Umodes.append(problem->supmodes[k]);
    }

    // Create locally the pressure modes
    for (label k = 0; k < problem->NPmodes; k++)
    {
        Pmodes.append(problem->Pmodes[k]);
    }

    for (label k = 0; k < problem->liftfieldT.size(); k++)
    {
        Tmodes.append(problem->liftfieldT[k]);
    }

    // Create locally the temperature modes
    for (label k = 0; k < problem->NTmodes; k++)
    {
        Tmodes.append(problem->Tmodes[k]);
    }

    // Store locally the snapshots for projections
    for (label k = 0; k < problem->Ufield.size(); k++)
    {
        Usnapshots.append(problem->Ufield[k]);
        Psnapshots.append(problem->Pfield[k]);
    }

    for (label k = 0; k < problem->Tfield.size(); k++)
    {
        Tsnapshots.append(problem->Tfield[k]);
    }

    newton_object_sup = newton_unsteadyNST_sup(Nphi_u + Nphi_p, Nphi_u + Nphi_p,
                        FOMproblem);
    newton_object_sup_t = newton_unsteadyNST_sup_t(Nphi_t, Nphi_t, FOMproblem);
}


// * * * * * * * * * * * * * * * Operators supremizer  * * * * * * * * * * * * * //
//Operator to evaluate the residual for the supremizer approach
label newton_unsteadyNST_sup::operator()(const Eigen::VectorXd& x,
        Eigen::VectorXd& fvec) const
{
    Eigen::VectorXd a_dot(Nphi_u);
    Eigen::VectorXd a_tmp(Nphi_u);
    Eigen::VectorXd b_tmp(Nphi_p);
    a_tmp = x.head(Nphi_u);
    b_tmp = x.tail(Nphi_p);
    a_dot = (x.head(Nphi_u) - y_old.head(Nphi_u)) / dt;
    // Convective term
    Eigen::MatrixXd cc(1, 1);
    // Momentum Term
    Eigen::VectorXd M1 = problem->B_matrix * a_tmp * nu;
    // Gradient of pressure
    Eigen::VectorXd M2 = problem->K_matrix * b_tmp;
    // Mass Term Velocity
    Eigen::VectorXd M5 = problem->M_matrix * a_dot;
    // Pressure Term
    Eigen::VectorXd M3 = problem->P_matrix * a_tmp;

    for (label i = 0; i < Nphi_u; i++)
    {
        cc = a_tmp.transpose() * problem->C_matrix[i] * a_tmp;
        fvec(i) = - M5(i) + M1(i) - cc(0, 0) - M2(i);
    }

    for (label j = 0; j < Nphi_p; j++)
    {
        label k = j + Nphi_u;
        fvec(k) = M3(j);
    }

    for (label j = 0; j < N_BC; j++)
    {
        fvec(j) = x(j) - BC(j);
    }

    return 0;
}

// Operator to evaluate the Jacobian for the supremizer approach
label newton_unsteadyNST_sup::df(const Eigen::VectorXd& x,
                                 Eigen::MatrixXd& fjac) const
{
    Eigen::NumericalDiff<newton_unsteadyNST_sup> numDiff(*this);
    numDiff.df(x, fjac);
    return 0;
}

label newton_unsteadyNST_sup_t::operator()(const Eigen::VectorXd& t,
        Eigen::VectorXd& fvect) const
{
    Eigen::VectorXd c_dot(Nphi_t);
    Eigen::VectorXd c_tmp(Nphi_t);
    c_tmp = t.head(Nphi_t);
    c_dot = (t.head(Nphi_t) - z_old.head(Nphi_t)) / dt;
    // Convective term temperature
    Eigen::MatrixXd qq(1, 1);
    // diffusive term temperature
    Eigen::VectorXd M6 = problem->Y_matrix * c_tmp * DT;
    // Mass Term Temperature
    Eigen::VectorXd M8 = problem->MT_matrix * c_dot;

    for (label i = 0; i < Nphi_t; i++)
    {
        qq = a_tmp.transpose() * problem-> Q_matrix[i] * c_tmp;
        fvect(i) = -M8(i) + M6(i) - qq(0, 0);
    }

    for (label j = 0; j < N_BC_t; j++)
    {
        fvect(j) = t(j) - BC_t(j);
    }

    return 0;
}
label newton_unsteadyNST_sup_t::df(const Eigen::VectorXd& t,
                                   Eigen::MatrixXd& fjact) const
{
    Eigen::NumericalDiff<newton_unsteadyNST_sup_t> numDiff(*this);
    numDiff.df(t, fjact);
    return 0;
}

// * * * * * * * * * * * * * * * Solve Functions  * * * * * * * * * * * * * //
void reducedUnsteadyNST::solveOnline_sup(Eigen::MatrixXd& vel_now,
        Eigen::MatrixXd& temp_now, label startSnap)
{
    // Create and resize the solution vector
    y.resize(Nphi_u + Nphi_p, 1);
    y.setZero();
    z.resize(Nphi_t, 1);
    z.setZero();
    // Set Initial Conditions
    y.head(Nphi_u) = ITHACAutilities::getCoeffs(Usnapshots[startSnap], Umodes);
    y.tail(Nphi_p) = ITHACAutilities::getCoeffs(Psnapshots[startSnap], Pmodes);
    z.head(Nphi_t) = ITHACAutilities::getCoeffs(Tsnapshots[startSnap], Tmodes);

    // Change initial condition for the lifting function
    for (label j = 0; j < N_BC; j++)
    {
        y(j) = vel_now(j, 0);
    }

    for (label j = 0; j < N_BC_t; j++)
    {
        z(j) = temp_now(j, 0);
    }

    // Set some properties of the newton object
    newton_object_sup.nu = nu;
    newton_object_sup.y_old = y;
    newton_object_sup.dt = dt;
    newton_object_sup_t.DT = DT;
    newton_object_sup_t.z_old = z;
    newton_object_sup_t.dt = dt;
    newton_object_sup.BC.resize(N_BC);
    newton_object_sup_t.BC_t.resize(N_BC_t);

    for (label j = 0; j < N_BC_t; j++)
    {
        newton_object_sup_t.BC_t(j) = temp_now(j, 0);
    }

    for (label j = 0; j < N_BC; j++)
    {
        newton_object_sup.BC(j) = vel_now(j, 0);
    }

    // Set number of online solutions
    label Ntsteps = static_cast<int>((finalTime - tstart) / dt);
    online_solution.resize(Ntsteps);
    online_solutiont.resize(Ntsteps);
    // Set the initial time
    time = tstart;
    // Counting variable
    label counter = 0;
    // Create vector to store temporal solution and save initial condition as first solution
    Eigen::MatrixXd tmp_sol(Nphi_u + Nphi_p + 1, 1);
    tmp_sol(0) = time;
    tmp_sol.col(0).tail(y.rows()) = y;
    online_solution[counter] = tmp_sol;
    Eigen::MatrixXd tmp_solt(Nphi_t + 1, 1);
    tmp_solt(0) = time;
    tmp_solt.col(0).tail(z.rows()) = z;
    online_solutiont[counter] = tmp_solt;
    counter ++;
    // Create nonlinear solver object
    Eigen::HybridNonLinearSolver<newton_unsteadyNST_sup> hnls(newton_object_sup);
    Eigen::HybridNonLinearSolver<newton_unsteadyNST_sup_t> hnlst(
        newton_object_sup_t);
    // Set output colors for fancy output
    Color::Modifier red(Color::FG_RED);
    Color::Modifier green(Color::FG_GREEN);
    Color::Modifier def(Color::FG_DEFAULT);

    // Start the time loop
    while (time < finalTime)
    {
        time = time + dt;
        Eigen::VectorXd res(y);
        Eigen::VectorXd rest(z);
        res.setZero();
        rest.setZero();
        hnls.solve(y);

        for (label j = 0; j < N_BC; j++)
        {
            y(j) = vel_now(j, 0);
        }

        // set the a_temp
        // solve for temperature
        newton_object_sup_t.a_tmp = y.head(Nphi_u);
        hnlst.solve(z);

        for (label j = 0; j < N_BC_t; j++)
        {
            z(j) = temp_now(j, 0);
        }

        newton_object_sup.operator()(y, res);
        newton_object_sup.y_old = y;
        newton_object_sup_t.operator()(z, rest);
        newton_object_sup_t.z_old = z;
        std::cout << "################## Online solve N° " << count_online_solve <<
                  " ##################" << std::endl;
        Info << "Time = " << time << endl;
        std::cout << "Solving for the parameter: " << vel_now << std::endl;
        std::cout << "Solving for the parameter: " << temp_now << std::endl;

        if (res.norm() < 1e-5)
        {
            std::cout << green << "|F(x)| = " << res.norm() << " - Minimun reached in " <<
                      hnls.iter << " iterations " << def << std::endl << std::endl;
        }
        else
        {
            std::cout << red << "|F(x)| = " << res.norm() << " - Minimun reached in " <<
                      hnls.iter << " iterations " << def << std::endl << std::endl;
        }

        if (rest.norm() < 1e-5)
        {
            std::cout << green << "|F(x)| = " << rest.norm() << " - Minimun reached in " <<
                      hnlst.iter << " iterations " << def << std::endl << std::endl;
        }
        else
        {
            std::cout << red << "|F(x)| = " << rest.norm() << " - Minimun reached in " <<
                      hnlst.iter << " iterations " << def << std::endl << std::endl;
        }

        count_online_solve += 1;
        tmp_sol(0) = time;
        tmp_sol.col(0).tail(y.rows()) = y;

        if (counter >= online_solution.size())
        {
            online_solution.append(tmp_sol);
        }
        else
        {
            online_solution[counter] = tmp_sol;
        }

        tmp_solt(0) = time;
        tmp_solt.col(0).tail(z.rows()) = z;

        if (counter >= online_solutiont.size())
        {
            online_solutiont.append(tmp_solt);
        }
        else
        {
            online_solutiont[counter] = tmp_solt;
        }

        counter ++;
    }

    // Save the solution
    ITHACAstream::exportMatrix(online_solution, "red_coeff", "python",
                               "./ITHACAoutput/red_coeff");
    ITHACAstream::exportMatrix(online_solution, "red_coeff", "matlab",
                               "./ITHACAoutput/red_coeff");
    ITHACAstream::exportMatrix(online_solutiont, "red_coeff", "python",
                               "./ITHACAoutput/red_coeff_t");
    ITHACAstream::exportMatrix(online_solutiont, "red_coeff", "matlab",
                               "./ITHACAoutput/red_coeff_t");
    count_online_solve += 1;
}


void reducedUnsteadyNST::reconstruct_sup(fileName folder, label printevery)
{
    mkDir(folder);
    ITHACAutilities::createSymLink(folder);
    label counter = 0;
    label nextwrite = 0;
    label counter2 = 1;

    for (label i = 0; i < online_solution.size(); i++)
    {
        if (counter == nextwrite)
        {
            volVectorField U_rec("U_rec", Umodes[0] * 0);

            for (label j = 0; j < Nphi_u; j++)
            {
                U_rec += Umodes[j] * online_solution[i](j + 1, 0);
            }

            //problem.exportSolution(U_rec, name(online_solution[i](0, 0)), folder);
            ITHACAstream::exportSolution(U_rec,  name(counter2), folder);
            volScalarField P_rec("P_rec", Pmodes[0] * 0);

            for (label j = 0; j < Nphi_p; j++)
            {
                P_rec += Pmodes[j] * online_solution[i](j + Nphi_u + 1, 0);
            }

            //problem.exportSolution(P_rec, name(online_solution[i](0, 0)), folder);
            ITHACAstream::exportSolution(P_rec, name(counter2), folder);
            nextwrite += printevery;
            counter2 ++;
            UREC.append(U_rec);
            PREC.append(P_rec);
        }

        counter++;
    }
}



void reducedUnsteadyNST::reconstruct_supt(fileName folder, label printevery)
{
    mkDir(folder);
    ITHACAutilities::createSymLink(folder);
    label counter = 0;
    label nextwrite = 0;
    label counter2 = 1;

    for (label i = 0; i < online_solutiont.size(); i++)
    {
        if (counter == nextwrite)
        {
            volScalarField T_rec("T_rec", Tmodes[0] * 0);

            for (label j = 0; j < Nphi_t; j++)
            {
                T_rec += Tmodes[j] * online_solutiont[i](j + 1, 0);
            }

            ITHACAstream::exportSolution(T_rec,  name(counter2), folder);
            nextwrite += printevery;
            counter2 ++;
            TREC.append(T_rec);
        }

        counter++;
    }
}

// ************************************************************************* //
