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
/// Source file of the reducedUnsteadyNS class


#include "ReducedUnsteadyNSTurb.H"


// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructor
ReducedUnsteadyNSTurb::ReducedUnsteadyNSTurb()
{
}

ReducedUnsteadyNSTurb::ReducedUnsteadyNSTurb(UnsteadyNSTurb& fomProblem)
{
    problem = &fomProblem;
    N_BC = problem->inletIndex.rows();
    Nphi_u = problem->B_matrix.rows();
    Nphi_p = problem->K_matrix.cols();
    nphiNut = problem->cTotalTensor.dimension(1);

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

    newtonObjectSUP = newtonUnsteadyNSTurbSUP(Nphi_u + Nphi_p, Nphi_u + Nphi_p,
                      fomProblem);
    newtonObjectPPE = newton_UnsteadyNSTurb_PPE(Nphi_u + Nphi_p, Nphi_u + Nphi_p,
                      fomProblem);
}


// * * * * * * * * * * * * * * * Operators supremizer  * * * * * * * * * * * * * //

// Operator to evaluate the residual for the supremizer approach
int newtonUnsteadyNSTurbSUP::operator()(const Eigen::VectorXd& x,
                                        Eigen::VectorXd& fvec) const
{
    Eigen::VectorXd a_dot(Nphi_u);
    Eigen::VectorXd aTmp(Nphi_u);
    Eigen::VectorXd bTmp(Nphi_p);
    aTmp = x.head(Nphi_u);
    bTmp = x.tail(Nphi_p);
    a_dot = (x.head(Nphi_u) - y_old.head(Nphi_u)) / dt;
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
        cc = aTmp.transpose() * Eigen::SliceFromTensor(problem->C_tensor, 0,
                i) * aTmp - gNut.transpose() *
             Eigen::SliceFromTensor(problem->cTotalTensor, 0, i) * aTmp;
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
int newtonUnsteadyNSTurbSUP::df(const Eigen::VectorXd& x,
                                Eigen::MatrixXd& fjac) const
{
    Eigen::NumericalDiff<newtonUnsteadyNSTurbSUP> numDiff(*this);
    numDiff.df(x, fjac);
    return 0;
}

// * * * * * * * * * * * * * * * Operators PPE * * * * * * * * * * * * * * * //

// Operator to evaluate the residual for the supremizer approach
int newton_UnsteadyNSTurb_PPE::operator()(const Eigen::VectorXd& x,
        Eigen::VectorXd& fvec) const
{
    Eigen::VectorXd a_dot(Nphi_u);
    Eigen::VectorXd aTmp(Nphi_u);
    Eigen::VectorXd bTmp(Nphi_p);
    aTmp = x.head(Nphi_u);
    bTmp = x.tail(Nphi_p);
    a_dot = (x.head(Nphi_u) - y_old.head(Nphi_u)) / dt;
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
        cc = aTmp.transpose() * problem->C_matrix[i] * aTmp - gNut.transpose() *
             problem->cTotalMatrix[i] * aTmp;
        fvec(i) = - m5(i) + m1(i) - cc(0, 0) - m2(i);
    }

    for (label j = 0; j < Nphi_p; j++)
    {
        label k = j + Nphi_u;
        gg = aTmp.transpose() * problem->G_matrix[j] * aTmp;
        bb = aTmp.transpose() * problem->BC2_matrix[j] * aTmp;
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
int newton_UnsteadyNSTurb_PPE::df(const Eigen::VectorXd& x,
                                  Eigen::MatrixXd& fjac) const
{
    Eigen::NumericalDiff<newton_UnsteadyNSTurb_PPE> numDiff(*this);
    numDiff.df(x, fjac);
    return 0;
}



// * * * * * * * * * * * * * * * Solve Functions  * * * * * * * * * * * * * //
void ReducedUnsteadyNSTurb::solveOnlineSUP(Eigen::MatrixXd& velNow,
        label startSnap)
{
    // Create and resize the solution vector
    y.resize(Nphi_u + Nphi_p, 1);
    y.setZero();
    y.head(Nphi_u) = ITHACAutilities::get_coeffs(Usnapshots[startSnap], Umodes);
    y.tail(Nphi_p) = ITHACAutilities::get_coeffs_ortho(Psnapshots[startSnap],
                     Pmodes);

    // Change initial condition for the lifting function
    for (label j = 0; j < N_BC; j++)
    {
        y(j) = velNow(j, 0);
    }

    // Set some properties of the newton object
    newtonObjectSUP.nu = nu;
    newtonObjectSUP.y_old = y;
    newtonObjectSUP.dt = dt;
    newtonObjectSUP.bc.resize(N_BC);

    for (label j = 0; j < N_BC; j++)
    {
        newtonObjectSUP.bc(j) = velNow(j, 0);
    }

    // Set number of online solutions
    int Ntsteps = static_cast<int>((finalTime - tstart) / dt);
    online_solution.resize(Ntsteps);
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
    }

    // Create nonlinear solver object
    Eigen::HybridNonLinearSolver<newtonUnsteadyNSTurbSUP> hnls(newtonObjectSUP);
    // Set output colors for fancy output
    Color::Modifier red(Color::FG_RED);
    Color::Modifier green(Color::FG_GREEN);
    Color::Modifier def(Color::FG_DEFAULT);

    // This part up to the while loop is just to compute the eddy viscosity field at time = startTime if time is not equal to zero
    if (time != 0)
    {
        volScalarField nutTmp("nutRec", problem->nutModes[0] * 0);
        Eigen::VectorXd gNut0;
        gNut0.resize(nphiNut);
        std::vector<double> tv0;
        tv0.resize(1);
        tv0[0] = time;

        for (label i = 0; i < nphiNut; i++)
        {
            gNut0(i) = problem->rbfSplines[i]->eval(tv0);
        }

        for (label j = 0; j < nphiNut; j++)
        {
            nutTmp += problem->nutModes[j] * gNut0(j);
        }

        nutRec.append(nutTmp);
    }

    // Start the time loop
    while (time < finalTime)
    {
        time = time + dt;
        std::vector<double> tv;
        tv.resize(1);
        tv[0] = time;

        for (label i = 0; i < nphiNut; i++)
        {
            newtonObjectSUP.gNut(i) = problem->rbfSplines[i]->eval(tv);
        }

        volScalarField nutTmp("nutRec", problem->nutModes[0] * 0);

        for (label j = 0; j < nphiNut; j++)
        {
            nutTmp += problem->nutModes[j] * newtonObjectSUP.gNut(j);
        }

        nutRec.append(nutTmp);
        Eigen::VectorXd res(y);
        res.setZero();
        hnls.solve(y);

        for (label j = 0; j < N_BC; j++)
        {
            y(j) = velNow(j, 0);
        }

        newtonObjectSUP.operator()(y, res);
        newtonObjectSUP.y_old = y;
        std::cout << "################## Online solve N° " << count_online_solve <<
                  " ##################" << std::endl;
        Info << "Time = " << time << endl;
        std::cout << "Solving for the parameter: " << velNow << std::endl;

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
    ITHACAstream::exportMatrix(online_solution, "red_coeff", "python",
                               "./ITHACAoutput/red_coeff");
    ITHACAstream::exportMatrix(online_solution, "red_coeff", "matlab",
                               "./ITHACAoutput/red_coeff");
    count_online_solve += 1;
}

// * * * * * * * * * * * * * * * Solve Functions  * * * * * * * * * * * * * //
void ReducedUnsteadyNSTurb::solveOnlinePPE(Eigen::MatrixXd& velNow,
        label startSnap)
{
    // Create and resize the solution vector
    y.resize(Nphi_u + Nphi_p, 1);
    y.setZero();
    // Set Initial Conditions
    y.head(Nphi_u) = ITHACAutilities::get_coeffs(Usnapshots[startSnap], Umodes);
    y.tail(Nphi_p) = ITHACAutilities::get_coeffs(Psnapshots[startSnap], Pmodes);

    // Change initial condition for the lifting function
    for (label j = 0; j < N_BC; j++)
    {
        y(j) = velNow(j, 0);
    }

    // Set some properties of the newton object
    newtonObjectPPE.nu = nu;
    newtonObjectPPE.y_old = y;
    newtonObjectPPE.dt = dt;
    newtonObjectPPE.bc.resize(N_BC);

    for (label j = 0; j < N_BC; j++)
    {
        newtonObjectPPE.bc(j) = velNow(j, 0);
    }

    // Set number of online solutions
    int Ntsteps = static_cast<int>((finalTime - tstart) / dt);
    online_solution.resize(Ntsteps);
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
    Eigen::HybridNonLinearSolver<newton_UnsteadyNSTurb_PPE> hnls(newtonObjectPPE);
    // Set output colors for fancy output
    Color::Modifier red(Color::FG_RED);
    Color::Modifier green(Color::FG_GREEN);
    Color::Modifier def(Color::FG_DEFAULT);

    // This part up to the while loop is just to compute the eddy viscosity field at time = startTime if time is not equal to zero
    if (time != 0)
    {
        volScalarField nutTmp("nutRec", problem->nutModes[0] * 0);
        Eigen::VectorXd gNut0;
        gNut0.resize(nphiNut);
        std::vector<double> tv0;
        tv0.resize(1);
        tv0[0] = time;

        for (label i = 0; i < nphiNut; i++)
        {
            gNut0(i) = problem->rbfSplines[i]->eval(tv0);
        }

        for (label j = 0; j < nphiNut; j++)
        {
            nutTmp += problem->nutModes[j] * gNut0(j);
        }

        nutRec.append(nutTmp);
    }

    // Start the time loop
    while (time < finalTime + dt)
    {
        time = time + dt;
        std::vector<double> tv;
        tv.resize(1);
        tv[0] = time;

        for (label i = 0; i < nphiNut; i++)
        {
            newtonObjectPPE.gNut(i) = problem->rbfSplines[i]->eval(tv);
        }

        volScalarField nutTmp("nutRec", problem->nutModes[0] * 0);

        for (label j = 0; j < nphiNut; j++)
        {
            nutTmp += problem->nutModes[j] * newtonObjectPPE.gNut(j);
        }

        nutRec.append(nutTmp);
        Eigen::VectorXd res(y);
        res.setZero();
        hnls.solve(y);

        for (label j = 0; j < N_BC; j++)
        {
            y(j) = velNow(j, 0);
        }

        Info << "before the operator" << endl;
        newtonObjectPPE.operator()(y, res);
        newtonObjectPPE.y_old = y;
        std::cout << "################## Online solve N° " << count_online_solve <<
                  " ##################" << std::endl;
        Info << "Time = " << time << endl;
        std::cout << "Solving for the parameter: " << velNow << std::endl;

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
    ITHACAstream::exportMatrix(online_solution, "red_coeff", "python",
                               "./ITHACAoutput/red_coeff");
    ITHACAstream::exportMatrix(online_solution, "red_coeff", "matlab",
                               "./ITHACAoutput/red_coeff");
    count_online_solve += 1;
}

void ReducedUnsteadyNSTurb::reconstructPPE(fileName folder, int printEvery)
{
    mkDir(folder);
    ITHACAutilities::createSymLink(folder);
    int counter = 0;
    int nextWrite = 0;
    int counter2 = 1;

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
            ITHACAstream::exportSolution(nutRec[nextWrite], name(counter2), folder);
            nextWrite += printEvery;
            double timeNow = online_solution[i](0, 0);
            std::ofstream of(folder + name(counter2) + "/" + name(timeNow));
            counter2 ++;
        }

        counter++;
    }
}

void ReducedUnsteadyNSTurb::reconstructSUP(fileName folder, int printEvery)
{
    mkDir(folder);
    ITHACAutilities::createSymLink(folder);
    int counter = 0;
    int nextWrite = 0;
    int counter2 = 1;

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
            ITHACAstream::exportSolution(nutRec[nextWrite], name(counter2), folder);
            nextWrite += printEvery;
            double timeNow = online_solution[i](0, 0);
            std::ofstream of(folder + name(counter2) + "/" + name(timeNow));
            counter2 ++;
        }

        counter++;
    }
}
// ************************************************************************* //

