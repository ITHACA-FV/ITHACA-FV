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

\*---------------------------------------------------------------------------*/

/// \file
/// Source file of the ReducedUnsteadyNSTurbIntrusive class


#include "ReducedUnsteadyNSTurbIntrusive.H"


// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructor
ReducedUnsteadyNSTurbIntrusive::ReducedUnsteadyNSTurbIntrusive()
{
}

ReducedUnsteadyNSTurbIntrusive::ReducedUnsteadyNSTurbIntrusive(
    UnsteadyNSTurbIntrusive& fomProblem)
{
    problem = &fomProblem;
    N_BC = problem->inletIndex.rows();
    Nphi_u = problem->B_matrix.rows();
    Nphi_p = problem->K_matrix.cols();

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

    // Store locally the snapshots for projections
    for (label k = 0; k < problem->Ufield.size(); k++)
    {
        Usnapshots.append(problem->Ufield[k]);
        Psnapshots.append(problem->Pfield[k]);
    }

    newtonObjectSUP = newtonUnsteadyNSTurbSUPIntrusive(Nphi_u + Nphi_p,
                      Nphi_u + Nphi_p,
                      fomProblem);
    newtonObjectPPE = newtonUnsteadyNSTurbPPEIntrusive(Nphi_u + Nphi_p,
                      Nphi_u + Nphi_p,
                      fomProblem);
}


// * * * * * * * * * * * * * * * Operators supremizer  * * * * * * * * * * * * * //

// Operator to evaluate the residual for the supremizer approach
int newtonUnsteadyNSTurbSUPIntrusive::operator()(const Eigen::VectorXd& x,
        Eigen::VectorXd& fvec) const
{
    Eigen::VectorXd a_dot(Nphi_u);
    Eigen::VectorXd aTmp(Nphi_u);
    Eigen::VectorXd bTmp(Nphi_p);
    aTmp = x.head(Nphi_u);
    bTmp = x.tail(Nphi_p);

    // Choose the order of the numerical difference scheme for approximating the time derivative
    if (problem->timeDerivativeSchemeOrder == "first")
    {
        a_dot = (x.head(Nphi_u) - y_old.head(Nphi_u)) / dt;
    }
    else
    {
        a_dot = (1.5 * x.head(Nphi_u) - 2 * y_old.head(Nphi_u) + 0.5 * yOldOld.head(
                     Nphi_u)) / dt;
    }

    // Convective term
    Eigen::MatrixXd cc(1, 1);
    // Mom Term
    Eigen::VectorXd m1 = problem->bTotalMatrix * aTmp * nu;
    // Gradient of pressure
    Eigen::VectorXd m2 = problem->K_matrix * bTmp;
    // Mass Term
    Eigen::VectorXd m5 = problem->M_matrix * a_dot;
    // Pressure Term
    Eigen::VectorXd m3 = problem->P_matrix * aTmp;

    for (label i = 0; i < Nphi_u; i++)
    {
        cc = aTmp.transpose() * Eigen::SliceFromTensor(problem->cTotalTensor, 0,
                i) * aTmp;
        fvec(i) = - m5(i) + m1(i) - cc(0, 0) - m2(i);
    }

    for (label j = 0; j < Nphi_p; j++)
    {
        label k = j + Nphi_u;
        fvec(k) = m3(j);
    }

    for (label j = 0; j < N_BC; j++)
    {
        fvec(j) = x(j) - bc(j);
    }

    return 0;
}

// Operator to evaluate the Jacobian for the supremizer approach
int newtonUnsteadyNSTurbSUPIntrusive::df(const Eigen::VectorXd& x,
        Eigen::MatrixXd& fjac) const
{
    Eigen::NumericalDiff<newtonUnsteadyNSTurbSUPIntrusive> numDiff(*this);
    numDiff.df(x, fjac);
    return 0;
}

// * * * * * * * * * * * * * * * Operators PPE * * * * * * * * * * * * * * * //

// Operator to evaluate the residual for the supremizer approach
int newtonUnsteadyNSTurbPPEIntrusive::operator()(const Eigen::VectorXd& x,
        Eigen::VectorXd& fvec) const
{
    Eigen::VectorXd a_dot(Nphi_u);
    Eigen::VectorXd aTmp(Nphi_u);
    Eigen::VectorXd bTmp(Nphi_p);
    aTmp = x.head(Nphi_u);
    bTmp = x.tail(Nphi_p);

    // Choose the order of the numerical difference scheme for approximating the time derivative
    if (problem->timeDerivativeSchemeOrder == "first")
    {
        a_dot = (x.head(Nphi_u) - y_old.head(Nphi_u)) / dt;
    }
    else
    {
        a_dot = (1.5 * x.head(Nphi_u) - 2 * y_old.head(Nphi_u) + 0.5 * yOldOld.head(
                     Nphi_u)) / dt;
    }

    // Convective terms
    Eigen::MatrixXd cc(1, 1);
    Eigen::MatrixXd gg(1, 1);
    Eigen::MatrixXd bb(1, 1);
    // Mom Term
    Eigen::VectorXd m1 = problem->B_matrix * aTmp * nu;
    // Gradient of pressure
    Eigen::VectorXd m2 = problem->K_matrix * bTmp;
    // Mass Term
    Eigen::VectorXd m5 = problem->M_matrix * a_dot;
    // Pressure Term
    Eigen::VectorXd m3 = problem->D_matrix * bTmp;
    // BC PPE
    Eigen::VectorXd m6 = problem->BC1_matrix * aTmp * nu;
    // BC PPE
    Eigen::VectorXd m7 = problem->BC3_matrix * aTmp * nu;

    for (label i = 0; i < Nphi_u; i++)
    {
        cc = aTmp.transpose() * Eigen::SliceFromTensor(problem->cTotalTensor, 0,
                i) * aTmp;
        fvec(i) = - m5(i) + m1(i) - cc(0, 0) - m2(i);
    }

    for (label j = 0; j < Nphi_p; j++)
    {
        label k = j + Nphi_u;
        gg = aTmp.transpose() * Eigen::SliceFromTensor(problem->gTensor, 0,
                j) * aTmp;
        bb = aTmp.transpose() * Eigen::SliceFromTensor(problem->bc2Tensor, 0,
                j) * aTmp;
        //fvec(k) = m3(j, 0) - gg(0, 0) - m6(j, 0) + bb(0, 0);
        fvec(k) = m3(j, 0) + gg(0, 0) - m7(j, 0);
    }

    for (label j = 0; j < N_BC; j++)
    {
        fvec(j) = x(j) - bc(j);
    }

    return 0;
}

// Operator to evaluate the Jacobian for the supremizer approach
int newtonUnsteadyNSTurbPPEIntrusive::df(const Eigen::VectorXd& x,
        Eigen::MatrixXd& fjac) const
{
    Eigen::NumericalDiff<newtonUnsteadyNSTurbPPEIntrusive> numDiff(*this);
    numDiff.df(x, fjac);
    return 0;
}



// * * * * * * * * * * * * * * * Solve Functions  * * * * * * * * * * * * * //
void ReducedUnsteadyNSTurbIntrusive::solveOnlineSUP(Eigen::MatrixXd vel,
        label startSnap)
{
    M_Assert(exportEvery >= dt,
             "The time step dt must be smaller than exportEvery.");
    M_Assert(storeEvery >= dt,
             "The time step dt must be smaller than storeEvery.");
    M_Assert(ITHACAutilities::isInteger(storeEvery / dt) == true,
             "The variable storeEvery must be an integer multiple of the time step dt.");
    M_Assert(ITHACAutilities::isInteger(exportEvery / dt) == true,
             "The variable exportEvery must be an integer multiple of the time step dt.");
    M_Assert(ITHACAutilities::isInteger(exportEvery / storeEvery) == true,
             "The variable exportEvery must be an integer multiple of the variable storeEvery.");
    int numberOfStores = round(storeEvery / dt);
    vel_now = setOnlineVelocity(vel);
    // Create and resize the solution vector
    y.resize(Nphi_u + Nphi_p, 1);
    y.setZero();
    y.head(Nphi_u) = ITHACAutilities::get_coeffs(Usnapshots[startSnap], Umodes);
    y.tail(Nphi_p) = ITHACAutilities::get_coeffs_ortho(Psnapshots[startSnap],
                     Pmodes);
    int nextStore = 0;
    int counter2 = 0;

    // Change initial condition for the lifting function
    for (label j = 0; j < N_BC; j++)
    {
        y(j) = vel_now(j, 0);
    }

    // Set some properties of the newton object
    newtonObjectSUP.nu = nu;
    newtonObjectSUP.y_old = y;
    newtonObjectSUP.yOldOld = newtonObjectSUP.y_old;
    newtonObjectSUP.dt = dt;
    newtonObjectSUP.bc.resize(N_BC);

    for (label j = 0; j < N_BC; j++)
    {
        newtonObjectSUP.bc(j) = vel_now(j, 0);
    }

    // Set number of online solutions
    int Ntsteps = static_cast<int>((finalTime - tstart) / dt);
    int onlineSize = static_cast<int>(Ntsteps / numberOfStores);
    online_solution.resize(onlineSize);
    // Set the initial time
    time = tstart;
    // Counting variable
    int counter = 0;
    // Create vector to store temporal solution and save initial condition as first solution
    Eigen::MatrixXd tmp_sol(Nphi_u + Nphi_p + 1, 1);
    tmp_sol(0) = time;
    tmp_sol.col(0).tail(y.rows()) = y;

    if (time != 0)
    {
        online_solution[counter] = tmp_sol;
        counter ++;
        counter2++;
        nextStore += numberOfStores;
    }

    // Create nonlinear solver object
    Eigen::HybridNonLinearSolver<newtonUnsteadyNSTurbSUPIntrusive> hnls(
        newtonObjectSUP);
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

        newtonObjectSUP.operator()(y, res);
        newtonObjectSUP.yOldOld = newtonObjectSUP.y_old;
        newtonObjectSUP.y_old = y;
        std::cout << "################## Online solve N° " << count_online_solve <<
                  " ##################" << std::endl;
        Info << "Time = " << time << endl;
        std::cout << "Solving for the parameter: " << vel_now << std::endl;

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

        count_online_solve += 1;
        tmp_sol(0) = time;
        tmp_sol.col(0).tail(y.rows()) = y;

        if (counter == nextStore)
        {
            if (counter2 >= online_solution.size())
            {
                online_solution.append(tmp_sol);
            }
            else
            {
                online_solution[counter2] = tmp_sol;
            }

            nextStore += numberOfStores;
            counter2 ++;
        }

        counter ++;
    }

    // Save the solution
    ITHACAstream::exportMatrix(online_solution, "red_coeff", "python",
                               "./ITHACAoutput/red_coeff");
    ITHACAstream::exportMatrix(online_solution, "red_coeff", "matlab",
                               "./ITHACAoutput/red_coeff");
    count_online_solve += 1;
}

// * * * * * * * * * * * * * * * Solve Functions  * * * * * * * * * * * * * //
void ReducedUnsteadyNSTurbIntrusive::solveOnlinePPE(Eigen::MatrixXd vel,
        label startSnap)
{
    M_Assert(exportEvery >= dt,
             "The time step dt must be smaller than exportEvery.");
    M_Assert(storeEvery >= dt,
             "The time step dt must be smaller than storeEvery.");
    M_Assert(ITHACAutilities::isInteger(storeEvery / dt) == true,
             "The variable storeEvery must be an integer multiple of the time step dt.");
    M_Assert(ITHACAutilities::isInteger(exportEvery / dt) == true,
             "The variable exportEvery must be an integer multiple of the time step dt.");
    M_Assert(ITHACAutilities::isInteger(exportEvery / storeEvery) == true,
             "The variable exportEvery must be an integer multiple of the variable storeEvery.");
    int numberOfStores = round(storeEvery / dt);
    vel_now = setOnlineVelocity(vel);
    // Create and resize the solution vector
    y.resize(Nphi_u + Nphi_p, 1);
    y.setZero();
    // Set Initial Conditions
    y.head(Nphi_u) = ITHACAutilities::get_coeffs(Usnapshots[startSnap], Umodes);
    y.tail(Nphi_p) = ITHACAutilities::get_coeffs(Psnapshots[startSnap], Pmodes);
    int nextStore = 0;
    int counter2 = 0;

    // Change initial condition for the lifting function
    for (label j = 0; j < N_BC; j++)
    {
        y(j) = vel_now(j, 0);
    }

    // Set some properties of the newton object
    newtonObjectPPE.nu = nu;
    newtonObjectPPE.y_old = y;
    newtonObjectPPE.yOldOld = newtonObjectPPE.y_old;
    newtonObjectPPE.dt = dt;
    newtonObjectPPE.bc.resize(N_BC);

    for (label j = 0; j < N_BC; j++)
    {
        newtonObjectPPE.bc(j) = vel_now(j, 0);
    }

    // Set number of online solutions
    int Ntsteps = static_cast<int>((finalTime - tstart) / dt);
    int onlineSize = static_cast<int>(Ntsteps / numberOfStores);
    online_solution.resize(onlineSize);
    // Set the initial time
    time = tstart;
    // Counting variable
    int counter = 0;
    // Create vectpr to store temporal solution and save initial condition as first solution
    Eigen::MatrixXd tmp_sol(Nphi_u + Nphi_p + 1, 1);
    tmp_sol(0) = time;
    tmp_sol.col(0).tail(y.rows()) = y;

    if (time != 0)
    {
        online_solution[counter] = tmp_sol;
        counter ++;
    }

    // Create nonlinear solver object
    Eigen::HybridNonLinearSolver<newtonUnsteadyNSTurbPPEIntrusive> hnls(
        newtonObjectPPE);
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

        newtonObjectPPE.operator()(y, res);
        newtonObjectPPE.yOldOld = newtonObjectPPE.y_old;
        newtonObjectPPE.y_old = y;
        std::cout << "################## Online solve N° " << count_online_solve <<
                  " ##################" << std::endl;
        Info << "Time = " << time << endl;
        std::cout << "Solving for the parameter: " << vel_now << std::endl;

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

        count_online_solve += 1;
        tmp_sol(0) = time;
        tmp_sol.col(0).tail(y.rows()) = y;

        if (counter == nextStore)
        {
            if (counter2 >= online_solution.size())
            {
                online_solution.append(tmp_sol);
            }
            else
            {
                online_solution[counter2] = tmp_sol;
            }

            nextStore += numberOfStores;
            counter2 ++;
        }

        counter ++;
    }

    // Save the solution
    ITHACAstream::exportMatrix(online_solution, "red_coeff", "python",
                               "./ITHACAoutput/red_coeff");
    ITHACAstream::exportMatrix(online_solution, "red_coeff", "matlab",
                               "./ITHACAoutput/red_coeff");
    count_online_solve += 1;
}

void ReducedUnsteadyNSTurbIntrusive::reconstructPPE(fileName folder)
{
    mkDir(folder);
    ITHACAutilities::createSymLink(folder);
    int counter = 0;
    int nextWrite = 0;
    int counter2 = 1;
    int exportEveryIndex = round(exportEvery / storeEvery);

    for (label i = 0; i < online_solution.size(); i++)
    {
        if (counter == nextWrite)
        {
            volVectorField uRec("uRec", Umodes[0] * 0);

            for (label j = 0; j < Nphi_u; j++)
            {
                uRec += Umodes[j] * online_solution[i](j + 1, 0);
            }

            ITHACAstream::exportSolution(uRec,  name(counter2), folder);
            volScalarField pRec("pRec", Pmodes[0] * 0);

            for (label j = 0; j < Nphi_p; j++)
            {
                pRec += Pmodes[j] * online_solution[i](j + Nphi_u + 1, 0);
            }

            ITHACAstream::exportSolution(pRec, name(counter2), folder);
            nextWrite += exportEveryIndex;
            double timeNow = online_solution[i](0, 0);
            std::ofstream of(folder + name(counter2) + "/" + name(timeNow));
            counter2 ++;
        }

        counter++;
    }
}

void ReducedUnsteadyNSTurbIntrusive::reconstructSUP(fileName folder)
{
    mkDir(folder);
    ITHACAutilities::createSymLink(folder);
    int counter = 0;
    int nextWrite = 0;
    int counter2 = 1;
    int exportEveryIndex = round(exportEvery / storeEvery);

    for (label i = 0; i < online_solution.size(); i++)
    {
        if (counter == nextWrite)
        {
            volVectorField uRec("uRec", Umodes[0] * 0);

            for (label j = 0; j < Nphi_u; j++)
            {
                uRec += Umodes[j] * online_solution[i](j + 1, 0);
            }

            ITHACAstream::exportSolution(uRec,  name(counter2), folder);
            volScalarField pRec("pRec", Pmodes[0] * 0);

            for (label j = 0; j < Nphi_p; j++)
            {
                pRec += Pmodes[j] * online_solution[i](j + Nphi_u + 1, 0);
            }

            ITHACAstream::exportSolution(pRec, name(counter2), folder);
            nextWrite += exportEveryIndex;
            double timeNow = online_solution[i](0, 0);
            std::ofstream of(folder + name(counter2) + "/" + name(timeNow));
            counter2 ++;
        }

        counter++;
    }
}

Eigen::MatrixXd ReducedUnsteadyNSTurbIntrusive::setOnlineVelocity(
    Eigen::MatrixXd vel)
{
    assert(problem->inletIndex.rows() == vel.rows()
           && "Imposed boundary conditions dimensions do not match given values matrix dimensions");
    Eigen::MatrixXd vel_scal;
    vel_scal.resize(vel.rows(), vel.cols());

    for (int k = 0; k < problem->inletIndex.rows(); k++)
    {
        label p = problem->inletIndex(k, 0);
        label l = problem->inletIndex(k, 1);
        scalar area = gSum(problem->liftfield[0].mesh().magSf().boundaryField()[p]);
        scalar u_lf = gSum(problem->liftfield[k].mesh().magSf().boundaryField()[p] *
                           problem->liftfield[k].boundaryField()[p]).component(l) / area;
        vel_scal(k, 0) = vel(k, 0) / u_lf;
    }

    return vel_scal;
}
// ************************************************************************* //

