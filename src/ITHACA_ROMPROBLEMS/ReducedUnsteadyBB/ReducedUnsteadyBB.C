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
/// Source file of the ReducedUnsteadyBB class


#include "ReducedUnsteadyBB.H"


// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructor
ReducedUnsteadyBB::ReducedUnsteadyBB()
{
}

ReducedUnsteadyBB::ReducedUnsteadyBB(UnsteadyBB& FOMproblem)
    :
    problem(&FOMproblem)
{
    N_BC_t    = problem->inletIndexT.rows();
    N_BC      = problem->inletIndex.rows();
    Nphi_u    = problem->B_matrix.rows();
    Nphi_prgh = problem->K_matrix.cols();
    Nphi_t    = problem->Y_matrix.rows();

    // Create locally the velocity modes
    for (int k = 0; k < problem->liftfield.size(); k++)
    {
        LUmodes.append((problem->liftfield[k]).clone());
    }

    for (int k = 0; k < problem->NUmodes; k++)
    {
        LUmodes.append((problem->Umodes[k]).clone());
    }

    for (int k = 0; k < problem->NSUPmodes; k++)
    {
        LUmodes.append((problem->supmodes[k]).clone());
    }

    // Create locally the temperature modes including BC with liftfield
    for (int k = 0; k < problem->liftfieldT.size(); k++)
    {
        LTmodes.append((problem->liftfieldT[k]).clone());
    }

    for (int k = 0; k < problem->NTmodes; k++)
    {
        LTmodes.append((problem->Tmodes[k]).clone());
    }

    newton_object_sup = newton_unsteadyBB_sup(Nphi_u + Nphi_prgh + Nphi_t,
                        Nphi_u + Nphi_prgh + Nphi_t,
                        FOMproblem);
}


// * * * * * * * * * * * * * * * Operators supremizer  * * * * * * * * * * * * * //
//Operator to evaluate the residual for the supremizer approach
int newton_unsteadyBB_sup::operator()(const Eigen::VectorXd& x,
                                      Eigen::VectorXd& fvec) const
{
    Eigen::VectorXd a_dot(Nphi_u);
    Eigen::VectorXd a_tmp(Nphi_u);
    Eigen::VectorXd b_tmp(Nphi_prgh);
    Eigen::VectorXd c_dot(Nphi_t);
    Eigen::VectorXd c_tmp(Nphi_t);
    a_tmp = x.head(Nphi_u);
    a_dot = (x.head(Nphi_u) - y_old.head(Nphi_u)) / dt;
    b_tmp = x.segment(Nphi_u, Nphi_prgh);
    c_tmp = x.tail(Nphi_t);
    c_dot = (x.tail(Nphi_t) - y_old.tail(Nphi_t)) / dt;
    // Convective term
    Eigen::MatrixXd cc(1, 1);
    // Diffusive Term
    Eigen::VectorXd M1 = problem->B_matrix * a_tmp * nu;
    // Mass Term Velocity
    Eigen::VectorXd M5 = problem->M_matrix * a_dot;
    // Gradient of pressure
    Eigen::VectorXd M2 = problem->K_matrix * b_tmp;
    // Continuity
    Eigen::VectorXd M3 = problem->P_matrix * a_tmp;
    // Buoyancy Term
    Eigen::VectorXd M10 = problem->H_matrix * c_tmp;
    // Convective term temperature
    Eigen::MatrixXd qq(1, 1);
    // diffusive term temperature
    Eigen::VectorXd M6 = problem->Y_matrix * c_tmp * (nu / Pr);
    // Mass Term Temperature
    Eigen::VectorXd M8 = problem->W_matrix * c_dot;

    for (int i = 0; i < Nphi_u; i++)
    {
        cc = a_tmp.transpose() * Eigen::SliceFromTensor(problem->C_tensor, 0,
                i) * a_tmp;
        fvec(i) = - M5(i) + M1(i) - cc(0, 0) - M10(i) - M2(i);
    }

    for (int j = 0; j < Nphi_prgh; j++)
    {
        int k = j + Nphi_u;
        fvec(k) = M3(j);
    }

    for (int j = 0; j < Nphi_t; j++)
    {
        int k = j + Nphi_u + Nphi_prgh;
        qq = a_tmp.transpose() * problem->Q_matrix[j]  * c_tmp;
        fvec(k) =   -M8(j) + M6(j) - qq(0, 0);
    }

    for (int j = 0; j < N_BC; j++)
    {
        fvec(j) = x(j) - BC(j);
    }

    for (int j = 0; j < N_BC_t; j++)
    {
        int k = j + Nphi_u + Nphi_prgh;
        fvec(k) = x(k) - BC_t(j);
    }

    return 0;
}

// Operator to evaluate the Jacobian for the supremizer approach
int newton_unsteadyBB_sup::df(const Eigen::VectorXd& x,
                              Eigen::MatrixXd& fjac) const
{
    Eigen::NumericalDiff<newton_unsteadyBB_sup> numDiff(*this);
    numDiff.df(x, fjac);
    return 0;
}

// * * * * * * * * * * * * * * * Operators PPE * * * * * * * * * * * * * * * //

// Operator to evaluate the residual for the supremizer approach
int newton_unsteadyBB_PPE::operator()(const Eigen::VectorXd& x,
                                      Eigen::VectorXd& fvec) const
{
    Eigen::VectorXd a_dot(Nphi_u);
    Eigen::VectorXd a_tmp(Nphi_u);
    Eigen::VectorXd b_tmp(Nphi_prgh);
    Eigen::VectorXd c_dot(Nphi_t);
    Eigen::VectorXd c_tmp(Nphi_t);
    a_tmp = x.head(Nphi_u);
    b_tmp = x.segment(Nphi_u, Nphi_prgh);
    a_dot = (x.head(Nphi_u) - y_old.head(Nphi_u)) / dt;
    c_tmp = x.tail(Nphi_t);
    c_dot = (x.tail(Nphi_t) - y_old.tail(Nphi_t)) / dt;
    // Convective terms
    Eigen::MatrixXd cc(1, 1);
    Eigen::MatrixXd gg(1, 1);
    Eigen::MatrixXd bb(1, 1);
    // Convective term temperature
    Eigen::MatrixXd qq(1, 1);
    // Eigen::MatrixXd st(1, 1);
    // Mom Term
    Eigen::VectorXd M1 = problem->B_matrix * a_tmp * nu;
    // Gradient of pressure
    Eigen::VectorXd M2 = problem->K_matrix * b_tmp;
    // Mass Term
    Eigen::VectorXd M5 = problem->M_matrix * a_dot;
    // Pressure Term
    Eigen::VectorXd M3 = problem->D_matrix * b_tmp;
    // BC PPE
    Eigen::VectorXd M6 = problem->BC1_matrix * a_tmp * nu;
    // BC PPE
    // Buoyancy Term
    Eigen::VectorXd M10 = problem->H_matrix * c_tmp;
    Eigen::VectorXd M11 = problem->HP_matrix * c_tmp;
    Eigen::VectorXd M7 = problem->BC3_matrix * a_tmp * nu;
    // diffusive term temperature
    Eigen::VectorXd M9 = problem->Y_matrix * c_tmp * (nu / Pr);
    // Mass Term Temperature
    Eigen::VectorXd M8 = problem->W_matrix * c_dot;

    for (int i = 0; i < Nphi_u; i++)
    {
        cc = a_tmp.transpose() * Eigen::SliceFromTensor(problem->C_tensor, 0,
                i) * a_tmp;
        fvec(i) = - M5(i) + M1(i) - cc(0, 0) - M10(i) - M2(i);
    }

    for (int j = 0; j < Nphi_prgh; j++)
    {
        int k = j + Nphi_u;
        gg = a_tmp.transpose() * problem->G_matrix[j] * a_tmp;
        bb = a_tmp.transpose() * problem->BC2_matrix[j] * a_tmp;
        //fvec(k) = M3(j, 0) - gg(0, 0) - M6(j, 0) + bb(0, 0);
        fvec(k) = M3(j, 0) + gg(0, 0) + M11(j, 0) - M7(j, 0);
    }

    for (int j = 0; j < Nphi_t; j++)
    {
        int k = j + Nphi_u + Nphi_prgh;
        qq = a_tmp.transpose() * problem->Q_matrix[j] * c_tmp;
        fvec(k) = -M8(j) + M9(j) - qq(0, 0);
    }

    // for (int j = 0; j < N_BC; j++)
    //{
    //     fvec(j) = x(j) - BC(j);
    // }
    return 0;
}

// Operator to evaluate the Jacobian for the supremizer approach
int newton_unsteadyBB_PPE::df(const Eigen::VectorXd& x,
                              Eigen::MatrixXd& fjac) const
{
    Eigen::NumericalDiff<newton_unsteadyBB_PPE> numDiff(*this);
    numDiff.df(x, fjac);
    return 0;
}

// * * * * * * * * * * * * * * * Solve Functions  * * * * * * * * * * * * * //
Eigen::MatrixXd ReducedUnsteadyBB::solveOnline_sup(Eigen::MatrixXd& temp_now_BC,
        Eigen::MatrixXd& vel_now_BC, int NParaSet, int startSnap)
{
    std::cout << "################## Online solve N° " << NParaSet <<
              " ##################" << std::endl;
    std::cout << "Solving for the parameter: " << temp_now_BC << std::endl;
    // Count number of time steps
    int counter = 0;
    time = tstart;

    while (time < finalTime - 0.5 * dt)
    {
        time = time + dt;
        counter ++;
    }

    // Set size of online solution
    online_solutiont.resize(Nphi_u + Nphi_prgh + Nphi_t + 1, counter + 1);
    // Set initial condition for online solve
    volScalarField T_IC("T_IC", problem->Tfield[0]);

    for (int j = 0; j < T_IC.boundaryField().size(); j++)
    {
        for (int i = 0; i < N_BC_t; i++)
        {
            if (j == problem->inletIndexT(i, 0))
            {
                T_IC.boundaryFieldRef()[problem->inletIndexT(i, 0)][j] = temp_now_BC(i, 0);
            }
            else
            {
            }
        }
    }

    // Create and resize the solution vector
    y.resize(Nphi_u + Nphi_prgh + Nphi_t, 1);
    y.setZero();
    // Calculate the time-dependent coefficients
    y.head(Nphi_u) = ITHACAutilities::getCoeffs(problem->Ufield[startSnap],
                     LUmodes);

    if  (Nphi_prgh != 0)
    {
        y.segment(Nphi_u, Nphi_prgh) =  ITHACAutilities::getCoeffs(
                                            problem->Prghfield[2],
                                            problem->Prghmodes);
    }

    y.tail(Nphi_t) = ITHACAutilities::getCoeffs(T_IC, LTmodes);
    // Set some properties of the newton object
    newton_object_sup.nu = nu;
    newton_object_sup.y_old = y;
    newton_object_sup.dt = dt;
    newton_object_sup.Pr = Pr;
    newton_object_sup.BC_t.resize(N_BC_t);
    newton_object_sup.BC.resize(N_BC);

    // Change initial condition for the lifting function
    for (int j = 0; j < N_BC_t; j++)
    {
        newton_object_sup.BC_t(j) = temp_now_BC(j, 0);
    }

    for (int j = 0; j < N_BC; j++)
    {
        newton_object_sup.BC(j) = vel_now_BC(j, 0);
    }

    // Set the initial time
    time = tstart;
    // Create vector to store temporal solution and save initial condition as first solution
    Eigen::MatrixXd tmp_sol(Nphi_u + Nphi_prgh + Nphi_t + 1, 1);
    tmp_sol(0) = time;
    tmp_sol.col(0).tail(y.rows()) = y;
    online_solutiont.col(0) = tmp_sol;
    // Create nonlinear solver object
    Eigen::HybridNonLinearSolver<newton_unsteadyBB_sup> hnls(newton_object_sup);
    // Set output colors for fancy output
    Color::Modifier red(Color::FG_RED);
    Color::Modifier green(Color::FG_GREEN);
    Color::Modifier def(Color::FG_DEFAULT);

    // Start the time loop
    for (int i = 1; i < online_solutiont.cols(); i++)
    {
        time = time + dt;
        Eigen::VectorXd res(y);
        res.setZero();
        hnls.solve(y);

        for (int j = 0; j < N_BC; j++)
        {
            y(j) = vel_now_BC(j, 0);
        }

        for (int j = 0; j < N_BC_t; j++)
        {
            int k = j + Nphi_prgh + Nphi_u;
            y(k) = temp_now_BC(j, 0);
        }

        newton_object_sup.operator()(y, res);
        newton_object_sup.y_old = y;
        Info << "Time = " << time << endl;

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

        tmp_sol(0) = time;
        tmp_sol.col(0).tail(y.rows()) = y;
        online_solutiont.col(i) = tmp_sol;
    }

    // Save the current solution
    ITHACAstream::exportMatrix(online_solutiont, "red_coeff", "python",
                               "./ITHACAoutput/red_coeff/" + name(NParaSet) + "/");
    ITHACAstream::exportMatrix(online_solutiont, "red_coeff", "matlab",
                               "./ITHACAoutput/red_coeff/" + name(NParaSet)  + "/");
    ITHACAstream::exportMatrix(online_solutiont, "red_coeff", "eigen",
                               "./ITHACAoutput/red_coeff/" + name(NParaSet)  + "/");
    return online_solutiont;
}

// * * * * * * * * * * * * * * * Solve Functions  * * * * * * * * * * * * * //
Eigen::MatrixXd ReducedUnsteadyBB::solveOnline_PPE(Eigen::MatrixXd&
        temp_now_BC,
        Eigen::MatrixXd& vel_now_BC, int NParaSet, int startSnap)
{
    std::cout << "################## Online solve N° " << NParaSet <<
              " ##################" << std::endl;
    std::cout << "Solving for the parameter: " << temp_now_BC << std::endl;
    // Count number of time steps
    int counter = 0;
    time = tstart;

    while (time < finalTime - 0.5 * dt)
    {
        time = time + dt;
        counter ++;
    }

    // Set size of online solution
    online_solutiont.resize(Nphi_u + Nphi_prgh + Nphi_t + 1, counter + 1);
    //Average Method
    //if (problem->AveMethod == "mean")
    //{
    // for (int j = 0; j < temp_now_BC.cols(); j++)
    //{
    //    for (int i = 0; i < N_BC_t; i++)
    //        {
    //int patche = problem->inletIndexT(i, 0);
    //           temp_now_BC(i, j)= temp_now_BC(i, j)-problem->Tsub[0].boundaryField()[problem->inletIndexT(i, 0)][0];
    //    }
    // }
    //volScalarField T_IC("T_IC", problem->Tsub[0]);
    //for (int j = 0; j < T_IC.boundaryField().size(); j++)
    //   {
    //     for (int i = 0; i < N_BC_t; i++)
    //      {
    //       T_IC.boundaryFieldRef()[problem->inletIndexT(i, 0)][j] =
    //        temp_now_BC(i, 0)-problem->Tsub[0].boundaryField()[problem->inletIndexT(i, 0)][0];
    //
    //   }
    // }
    volScalarField T_IC("T_IC", problem->Tfield[0]);

    for (int j = 0; j < T_IC.boundaryField().size(); j++)
    {
        for (int i = 0; i < N_BC_t; i++)
        {
            if (j == problem->inletIndexT(i, 0))
            {
                T_IC.boundaryFieldRef()[problem->inletIndexT(i, 0)][j] = temp_now_BC(i, 0);
            }
            else
            {
            }
        }
    }

    // Create and resize the solution vector
    y.resize(Nphi_u + Nphi_prgh + Nphi_t, 1);
    y.setZero();
    // Calculate the time-dependent coefficients
    y.head(Nphi_u) = ITHACAutilities::getCoeffs(problem->Ufield[startSnap],
                     LUmodes);

    if  (Nphi_prgh != 0)
    {
        y.segment(Nphi_u, Nphi_prgh) =  ITHACAutilities::getCoeffs(
                                            problem->Prghfield[2], problem->Prghmodes);
    }

    y.tail(Nphi_t) = ITHACAutilities::getCoeffs(T_IC, LTmodes);
    // Set some properties of the newton object
    newton_object_PPE.nu = nu;
    newton_object_PPE.y_old = y;
    newton_object_PPE.dt = dt;
    newton_object_PPE.Pr = Pr;
    newton_object_PPE.BC_t.resize(N_BC_t);
    newton_object_PPE.BC.resize(N_BC);

    //Eigen::MatrixXd Ncoeff = ITHACAutilities::getCoeffs(problem->nutFields
    //        , problem->nuTmodes);
    // std::vector<double> tv;
    //    tv.resize(1);
    //    tv[0] = time;

    //for (int l = 0; l < Nphi_nut; l++)
    //{
    //        newton_object_PPE.nu_c(l) = problem->rbfsplines[l]->eval(tv);
    //}

    //Change initial condition for the lifting function
    for (int j = 0; j < N_BC_t; j++)
    {
        newton_object_PPE.BC_t(j) = temp_now_BC(j, 0);
    }

    for (int j = 0; j < N_BC; j++)
    {
        newton_object_PPE.BC(j) = vel_now_BC(j, 0);
    }

    // Set the initial time
    time = tstart;
    // Create vector to store temporal solution and save initial condition as first solution
    Eigen::MatrixXd tmp_sol(Nphi_u + Nphi_prgh + Nphi_t + 1, 1);
    tmp_sol(0) = time;
    tmp_sol.col(0).tail(y.rows()) = y;
    online_solutiont.col(0) = tmp_sol;
    // Create nonlinear solver object
    Eigen::HybridNonLinearSolver<newton_unsteadyBB_PPE> hnls(newton_object_PPE);
    // Set output colors for fancy output
    Color::Modifier red(Color::FG_RED);
    Color::Modifier green(Color::FG_GREEN);
    Color::Modifier def(Color::FG_DEFAULT);

    // Start the time loop
    for (int i = 1; i < online_solutiont.cols(); i++)
    {
        time = time + dt;
        Eigen::VectorXd res(y);
        res.setZero();
        hnls.solve(y);

        //if (problem->BCmethod == "lift")
        // {

        //}

        for (int j = 0; j < N_BC; j++)
        {
            y(j) = vel_now_BC(j, 0);
        }

        for (int j = 0; j < N_BC_t; j++)
        {
            int k = j + Nphi_prgh + Nphi_u;
            y(k) = temp_now_BC(j, 0);
        }

        newton_object_PPE.operator()(y, res);
        newton_object_PPE.y_old = y;
        Info << "Time = " << time << endl;

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

        tmp_sol(0) = time;
        tmp_sol.col(0).tail(y.rows()) = y;
        online_solutiont.col(i) = tmp_sol;
    }

    // Save the current solution
    ITHACAstream::exportMatrix(online_solutiont, "red_coeff", "python",
                               "./ITHACAoutput/red_coeff/" + name(NParaSet) + "/");
    ITHACAstream::exportMatrix(online_solutiont, "red_coeff", "python",
                               ".");
    ITHACAstream::exportMatrix(online_solutiont, "red_coeff", "matlab",
                               "./ITHACAoutput/red_coeff/" + name(NParaSet)  + "/");
    ITHACAstream::exportMatrix(online_solutiont, "red_coeff", "eigen",
                               "./ITHACAoutput/red_coeff/" + name(NParaSet)  + "/");
    return online_solutiont;
}

void ReducedUnsteadyBB::reconstruct_sup(fileName folder, int printevery)
{
    if (ITHACAutilities::check_folder(folder))
    {
    }
    else
    {
        mkDir(folder);
        ITHACAutilities::createSymLink(folder);
    }

    int counter = 0;
    int nextwrite = 0;
    int counter2 = 1 + TREC.size();

    for (int i = 0; i < online_solutiont.cols(); i++)
    {
        if (counter == nextwrite)
        {
            volVectorField U_rec("U_rec", LUmodes[0] * 0);

            for (int j = 0; j < Nphi_u; j++)
            {
                U_rec += LUmodes[j] * online_solutiont(j + 1, i);
            }

            ITHACAstream::exportSolution(U_rec,  name(counter2), folder);

            if  (Nphi_prgh != 0)
            {
                volScalarField P_rec("P_rec", problem->Prghmodes[0] * 0);

                for (int j = 0; j < Nphi_prgh; j++)
                {
                    P_rec += problem->Prghmodes[j] * online_solutiont(j + Nphi_u + 1, i);
                }

                ITHACAstream::exportSolution(P_rec,  name(counter2), folder);
                PREC.append((P_rec).clone());
            }

            volScalarField T_rec("T_rec", LTmodes[0] * 0);

            for (int j = 0; j < Nphi_t; j++)
            {
                T_rec += LTmodes[j] * online_solutiont(j + Nphi_prgh + Nphi_u + 1, i);
            }

            ITHACAstream::exportSolution(T_rec, name(counter2), folder);
            nextwrite += printevery;
            double timenow = online_solutiont(0, i);
            std::ofstream of(folder + name(counter2) + "/" + name(timenow));
            counter2 ++;
            UREC.append((U_rec).clone());
            TREC.append((T_rec).clone());
        }

        counter++;
    }
}


// ************************************************************************* //
