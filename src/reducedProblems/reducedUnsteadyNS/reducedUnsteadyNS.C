/*---------------------------------------------------------------------------*\
Copyright (C) 2017 by the ITHACA-FV authors

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
/// Source file of the reducedUnsteadyNS class


#include "reducedUnsteadyNS.H"


// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructor
reducedUnsteadyNS::reducedUnsteadyNS()
{

}

reducedUnsteadyNS::reducedUnsteadyNS(unsteadyNS& problem, word tipo)
{
    B_matrix = problem.B_matrix;
    C_matrix = problem.C_matrix;
    M_matrix = problem.M_matrix;
    K_matrix = problem.K_matrix;
    N_BC = problem.inletIndex.rows();

    Nphi_u = B_matrix.rows();
    Nphi_p = K_matrix.cols();

    // Pressure Poisson equation approach
    if (tipo == "PPE")
    {
        D_matrix = problem.D_matrix;
        G_matrix = problem.G_matrix;
        BC1_matrix = problem.BC1_matrix;
        BC2_matrix = problem.BC2_matrix;
        BC3_matrix = problem.BC3_matrix;
        newton_object_PPE = newton_unsteadyNS_PPE(Nphi_u + Nphi_p , Nphi_u + Nphi_p, problem);
        newton_object_sup.nu = nu;
    }

    // Supremizer approach
    if (tipo == "SUP")
    {
        P_matrix = problem.P_matrix;
        newton_object_sup = newton_unsteadyNS_sup(Nphi_u + Nphi_p , Nphi_u + Nphi_p, problem);
        newton_object_sup.nu = nu;
    }

    // Create locally the velocity modes
    for (label k = 0; k < problem.liftfield.size(); k++)
    {
        Umodes.append(problem.liftfield[k]);
    }
    for (label k = 0; k < problem.NUmodes; k++)
    {
        Umodes.append(problem.Umodes[k]);
    }
    for (label k = 0; k < problem.NSUPmodes; k++)
    {
        Umodes.append(problem.supmodes[k]);
    }

    // Create locally the pressure modes
    for (label k = 0; k < problem.NPmodes; k++)
    {
        Pmodes.append(problem.Pmodes[k]);
    }

    // Store locally the snapshots for prokections
    for (label k = 0; k < problem.Ufield.size(); k++)
    {
        Usnapshots.append(problem.Ufield[k]);
        Psnapshots.append(problem.Pfield[k]);
    }
}


// * * * * * * * * * * * * * * * Operators supremizer  * * * * * * * * * * * * * //

// Operator to evaluate the residual for the supremizer approach
int newton_unsteadyNS_sup::operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const
{
    Eigen::VectorXd a_dot(Nphi_u);
    Eigen::VectorXd a_tmp(Nphi_u);
    Eigen::VectorXd b_tmp(Nphi_p);
    a_tmp = x.head(Nphi_u);
    b_tmp = x.tail(Nphi_p);
    a_dot = (x.head(Nphi_u) - y_old.head(Nphi_u)) / dt;

    // Convective term
    Eigen::MatrixXd cc(1, 1);
    // Mom Term
    Eigen::VectorXd M1 = B_matrix * a_tmp * nu;
    // Gradient of pressure
    Eigen::VectorXd M2 = K_matrix * b_tmp;
    // Mass Term
    Eigen::VectorXd M5 = M_matrix * a_dot;
    // Pressure Term
    Eigen::VectorXd M3 = P_matrix * a_tmp;

    for (label i = 0; i < Nphi_u; i++)
    {
        cc = a_tmp.transpose() * C_matrix[i] * a_tmp;
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
int newton_unsteadyNS_sup::df(const Eigen::VectorXd &x,  Eigen::MatrixXd &fjac) const
{
    Eigen::NumericalDiff<newton_unsteadyNS_sup> numDiff(*this);
    numDiff.df(x, fjac);

    return 0;
}

// * * * * * * * * * * * * * * * Operators PPE * * * * * * * * * * * * * * * //

// Operator to evaluate the residual for the supremizer approach
int newton_unsteadyNS_PPE::operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const
{
    Eigen::VectorXd a_dot(Nphi_u);
    Eigen::VectorXd a_tmp(Nphi_u);
    Eigen::VectorXd b_tmp(Nphi_p);
    a_tmp = x.head(Nphi_u);
    b_tmp = x.tail(Nphi_p);
    a_dot = (x.head(Nphi_u) - y_old.head(Nphi_u)) / dt;

    // Convective terms
    Eigen::MatrixXd cc(1, 1);
    Eigen::MatrixXd gg(1, 1);
    Eigen::MatrixXd bb(1, 1);
    // Mom Term
    Eigen::VectorXd M1 = B_matrix * a_tmp * nu;
    // Gradient of pressure
    Eigen::VectorXd M2 = K_matrix * b_tmp;
    // Mass Term
    Eigen::VectorXd M5 = M_matrix * a_dot;
    // Pressure Term
    Eigen::VectorXd M3 = D_matrix * b_tmp;

    // BC PPE
    Eigen::VectorXd M6 = BC1_matrix * a_tmp * nu;

    // BC PPE
    Eigen::VectorXd M7 = BC3_matrix * a_tmp * nu;



    for (label i = 0; i < Nphi_u; i++)
    {
        cc = a_tmp.transpose() * C_matrix[i] * a_tmp;
        fvec(i) = - M5(i) + M1(i) - cc(0, 0) - M2(i);
    }
    for (label j = 0; j < Nphi_p; j++)
    {
        label k = j + Nphi_u;
        gg = a_tmp.transpose() * G_matrix[j] * a_tmp;
        bb = a_tmp.transpose() * BC2_matrix[j] * a_tmp;
        //fvec(k) = M3(j, 0) - gg(0, 0) - M6(j, 0) + bb(0, 0);
        fvec(k) = M3(j, 0) + gg(0, 0) - M7(j, 0);
    }
    for (label j = 0; j < N_BC; j++)
    {
        fvec(j) = x(j) - BC(j);
    }
    return 0;
}

// Operator to evaluate the Jacobian for the supremizer approach
int newton_unsteadyNS_PPE::df(const Eigen::VectorXd &x,  Eigen::MatrixXd &fjac) const
{
    Eigen::NumericalDiff<newton_unsteadyNS_PPE> numDiff(*this);
    numDiff.df(x, fjac);
    return 0;
}



// * * * * * * * * * * * * * * * Solve Functions  * * * * * * * * * * * * * //
void reducedUnsteadyNS::solveOnline_sup(Eigen::MatrixXd vel_now, label startSnap)
{
    // Create and resize the solution vector
    y.resize(Nphi_u + Nphi_p, 1);
    y.setZero();

    // Set Initial Conditions
    // if(this->tstart!=0)
    // {
    y.head(Nphi_u) = ITHACAutilities::get_coeffs(Usnapshots[startSnap], Umodes);
    y.tail(Nphi_p) = ITHACAutilities::get_coeffs(Psnapshots[startSnap], Pmodes);
    // }



    // Change initial condition for the lifting function
    for (label j = 0; j < N_BC; j++)
    {
        y(j) = vel_now(j, 0);
    }


    // Set some properties of the newton object
    newton_object_sup.nu = nu;
    newton_object_sup.y_old = y;
    newton_object_sup.dt = dt;
    newton_object_sup.B_matrix = B_matrix;
    newton_object_sup.BC.resize(N_BC);
    for (label j = 0; j < N_BC; j++)
    {
        newton_object_sup.BC(j) = vel_now(j, 0);
    }


    // Set number of online solutions
    int Ntsteps = (int) ((finalTime - tstart) / dt);
    online_solution.resize(Ntsteps);

    // Set the initial time
    time = tstart;

    // Counting variable
    int counter = 0;

    // Create vectpr to store temporal solution and save initial condition as first solution
    Eigen::MatrixXd tmp_sol(Nphi_u + Nphi_p + 1, 1);
    tmp_sol(0) = time;
    tmp_sol.col(0).tail(y.rows()) = y;
    online_solution[counter] = tmp_sol;
    counter ++;

    // Create nonlinear solver object
    Eigen::HybridNonLinearSolver<newton_unsteadyNS_sup> hnls(newton_object_sup);

    // Set output colors for fancy output
    Color::Modifier red(Color::FG_RED);
    Color::Modifier green(Color::FG_GREEN);
    Color::Modifier def(Color::FG_DEFAULT);

    // Start the time loop
    while (time < finalTime)
    {
        time = time + dt;
        Eigen::VectorXd res(y);
        res.setZero();
        hnls.solve(y);
        for (label j = 0; j < N_BC; j++)
        {
            y(j) = vel_now(j, 0);
        }
        newton_object_sup.operator()(y, res);
        newton_object_sup.y_old = y;


        std::cout << "################## Online solve N° " << count_online_solve << " ##################" << std::endl;
        Info << "Time = " << time << endl;
        std::cout << "Solving for the parameter: " << vel_now << std::endl;
        if (res.norm() < 1e-5)
        {
            std::cout << green << "|F(x)| = " << res.norm() << " - Minimun reached in " << hnls.iter << " iterations " << def << std::endl << std::endl;
        }
        else
        {
            std::cout << red << "|F(x)| = " << res.norm() << " - Minimun reached in " << hnls.iter << " iterations " << def << std::endl << std::endl;
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
        counter ++;
    }
    // Save the solution
    ITHACAstream::exportMatrix(online_solution, "red_coeff", "python", "./ITHACAoutput/red_coeff");
    ITHACAstream::exportMatrix(online_solution, "red_coeff", "matlab", "./ITHACAoutput/red_coeff");
    count_online_solve += 1;
}

// * * * * * * * * * * * * * * * Solve Functions  * * * * * * * * * * * * * //
void reducedUnsteadyNS::solveOnline_PPE(Eigen::MatrixXd vel_now, label startSnap)
{
    // Create and resize the solution vector
    y.resize(Nphi_u + Nphi_p, 1);
    y.setZero();

    // Set Initial Conditions
    y.head(Nphi_u) = ITHACAutilities::get_coeffs(Usnapshots[startSnap], Umodes);;
    y.tail(Nphi_p) = ITHACAutilities::get_coeffs(Psnapshots[startSnap], Pmodes);;

    // Change initial condition for the lifting function
    for (label j = 0; j < N_BC; j++)
    {
        y(j) = vel_now(j, 0);
    }

    // Set some properties of the newton object
    newton_object_PPE.nu = nu;
    newton_object_PPE.y_old = y;
    newton_object_PPE.dt = dt;
    newton_object_PPE.B_matrix = B_matrix;
    newton_object_PPE.BC.resize(N_BC);
    for (label j = 0; j < N_BC; j++)
    {
        newton_object_PPE.BC(j) = vel_now(j, 0);
    }

    // Set number of online solutions
    int Ntsteps = (int) ((finalTime - tstart) / dt);
    online_solution.resize(Ntsteps);

    // Set the initial time
    time = tstart;

    // Counting variable
    int counter = 0;

    // Create vectpr to store temporal solution and save initial condition as first solution
    Eigen::MatrixXd tmp_sol(Nphi_u + Nphi_p + 1, 1);
    tmp_sol(0) = time;
    tmp_sol.col(0).tail(y.rows()) = y;
    online_solution[counter] = tmp_sol;
    counter ++;

    // Create nonlinear solver object
    Eigen::HybridNonLinearSolver<newton_unsteadyNS_PPE> hnls(newton_object_PPE);

    // Set output colors for fancy output
    Color::Modifier red(Color::FG_RED);
    Color::Modifier green(Color::FG_GREEN);
    Color::Modifier def(Color::FG_DEFAULT);

    // Start the time loop
    while (time < finalTime + dt)
    {
        time = time + dt;
        Eigen::VectorXd res(y);
        res.setZero();
        hnls.solve(y);
        for (label j = 0; j < N_BC; j++)
        {
            y(j) = vel_now(j, 0);
        }
        newton_object_PPE.operator()(y, res);
        newton_object_PPE.y_old = y;

        std::cout << "################## Online solve N° " << count_online_solve << " ##################" << std::endl;
        Info << "Time = " << time << endl;
        std::cout << "Solving for the parameter: " << vel_now << std::endl;
        if (res.norm() < 1e-5)
        {
            std::cout << green << "|F(x)| = " << res.norm() << " - Minimun reached in " << hnls.iter << " iterations " << def << std::endl << std::endl;
        }
        else
        {
            std::cout << red << "|F(x)| = " << res.norm() << " - Minimun reached in " << hnls.iter << " iterations " << def << std::endl << std::endl;
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
        counter ++;
    }

    // Save the solution
    ITHACAstream::exportMatrix(online_solution, "red_coeff", "python", "./ITHACAoutput/red_coeff");
    ITHACAstream::exportMatrix(online_solution, "red_coeff", "matlab", "./ITHACAoutput/red_coeff");
    count_online_solve += 1;
}

void reducedUnsteadyNS::reconstruct_PPE(unsteadyNS& problem, fileName folder, int printevery)
{
    mkDir(folder);
    system("ln -s ../../constant " + folder + "/constant");
    system("ln -s ../../0 " + folder + "/0");
    system("ln -s ../../system " + folder + "/system");

    int counter = 0;
    int nextwrite = 0;
    int counter2 = 1;

    for (label i = 0; i < online_solution.size(); i++)
    {
        if (counter == nextwrite)
        {
            volVectorField U_rec("U_rec", Umodes[0] * 0);
            for (label j = 0; j < Nphi_u; j++)
            {
                U_rec += Umodes[j] * online_solution[i](j + 1, 0);
            }
            problem.exportSolution(U_rec,  name(counter2), folder);

            volScalarField P_rec("P_rec", Pmodes[0] * 0);
            for (label j = 0; j < Nphi_p; j++)
            {
                P_rec += Pmodes[j] * online_solution[i](j + Nphi_u + 1, 0);
            }
            problem.exportSolution(P_rec, name(counter2), folder);
            nextwrite += printevery;
            counter2 ++;
        }
        counter++;
    }
}

void reducedUnsteadyNS::reconstruct_sup(unsteadyNS& problem, fileName folder, int printevery)
{
    mkDir(folder);
    system("ln -s ../../constant " + folder + "/constant");
    system("ln -s ../../0 " + folder + "/0");
    system("ln -s ../../system " + folder + "/system");

    int counter = 0;
    int nextwrite = 0;
    int counter2 = 1;

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
            problem.exportSolution(U_rec,  name(counter2), folder);

            volScalarField P_rec("P_rec", Pmodes[0] * 0);
            for (label j = 0; j < Nphi_p; j++)
            {
                P_rec += Pmodes[j] * online_solution[i](j + Nphi_u + 1, 0);
            }
            //problem.exportSolution(P_rec, name(online_solution[i](0, 0)), folder);
            problem.exportSolution(P_rec, name(counter2), folder);
            nextwrite += printevery;
            counter2 ++;
        }
        counter++;
    }
}
// ************************************************************************* //

