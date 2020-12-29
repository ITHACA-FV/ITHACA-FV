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
    dimA = problem->velRBF.cols();
    interChoice = problem->interChoice;

    // Create locally the velocity modes
    for (int k = 0; k < problem->liftfield.size(); k++)
    {
        Umodes.append(tmp<volVectorField>(problem->liftfield[k]));
    }

    for (int k = 0; k < problem->NUmodes; k++)
    {
        Umodes.append(tmp<volVectorField>(problem->Umodes[k]));
    }

    for (int k = 0; k < problem->NSUPmodes; k++)
    {
        Umodes.append(tmp<volVectorField>(problem->supmodes[k]));
    }

    // Create locally the pressure modes
    for (int k = 0; k < problem->NPmodes; k++)
    {
        Pmodes.append(tmp<volScalarField>(problem->Pmodes[k]));
    }

    // Create locally the eddy viscosity modes
    for (int k = 0; k < problem->nNutModes; k++)
    {
        nutModes.append(tmp<volScalarField>(problem->nutModes[k]));
    }

    newtonObjectSUP = newtonUnsteadyNSTurbSUP(Nphi_u + Nphi_p, Nphi_u + Nphi_p,
                      fomProblem);
    newtonObjectPPE = newtonUnsteadyNSTurbPPE(Nphi_u + Nphi_p, Nphi_u + Nphi_p,
                      fomProblem);
    newtonObjectSUPAve = newtonUnsteadyNSTurbSUPAve(Nphi_u + Nphi_p,
                         Nphi_u + Nphi_p,
                         fomProblem);
    newtonObjectPPEAve = newtonUnsteadyNSTurbPPEAve(Nphi_u + Nphi_p,
                         Nphi_u + Nphi_p,
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
    // Penalty term
    Eigen::MatrixXd penaltyU = Eigen::MatrixXd::Zero(Nphi_u, N_BC);

    // Term for penalty method
    if (problem->bcMethod == "penalty")
    {
        for (int l = 0; l < N_BC; l++)
        {
            penaltyU.col(l) = bc(l) * problem->bcVelVec[l] - problem->bcVelMat[l] *
                              aTmp;
        }
    }

    for (int i = 0; i < Nphi_u; i++)
    {
        cc = aTmp.transpose() * Eigen::SliceFromTensor(problem->C_tensor, 0,
                i) * aTmp - gNut.transpose() *
             Eigen::SliceFromTensor(problem->cTotalTensor, 0, i) * aTmp;
        fvec(i) = - m5(i) + m1(i) - cc(0, 0) - m2(i);

        if (problem->bcMethod == "penalty")
        {
            fvec(i) += ((penaltyU * tauU)(i, 0));
        }
    }

    for (int j = 0; j < Nphi_p; j++)
    {
        int k = j + Nphi_u;
        fvec(k) = m3(j);
    }

    if (problem->bcMethod == "lift")
    {
        for (int j = 0; j < N_BC; j++)
        {
            fvec(j) = x(j) - bc(j);
        }
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

// Operator to evaluate the residual for the supremizer approach
int newtonUnsteadyNSTurbSUPAve::operator()(const Eigen::VectorXd& x,
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
    // Penalty term
    Eigen::MatrixXd penaltyU = Eigen::MatrixXd::Zero(Nphi_u, N_BC);

    // Term for penalty method
    if (problem->bcMethod == "penalty")
    {
        for (int l = 0; l < N_BC; l++)
        {
            penaltyU.col(l) = bc(l) * problem->bcVelVec[l] - problem->bcVelMat[l] *
                              aTmp;
        }
    }

    for (int i = 0; i < Nphi_u; i++)
    {
        cc = aTmp.transpose() * Eigen::SliceFromTensor(problem->C_tensor, 0,
                i) * aTmp - gNut.transpose() *
             Eigen::SliceFromTensor(problem->cTotalTensor, 0,
                                    i) * aTmp - gNutAve.transpose() *
             Eigen::SliceFromTensor(problem->cTotalAveTensor, 0, i) * aTmp;
        fvec(i) = - m5(i) + m1(i) - cc(0, 0) - m2(i);

        if (problem->bcMethod == "penalty")
        {
            fvec(i) += ((penaltyU * tauU)(i, 0));
        }
    }

    for (int j = 0; j < Nphi_p; j++)
    {
        int k = j + Nphi_u;
        fvec(k) = m3(j);
    }

    if (problem->bcMethod == "lift")
    {
        for (int j = 0; j < N_BC; j++)
        {
            fvec(j) = x(j) - bc(j);
        }
    }

    return 0;
}

// Operator to evaluate the Jacobian for the supremizer approach
int newtonUnsteadyNSTurbSUPAve::df(const Eigen::VectorXd& x,
                                   Eigen::MatrixXd& fjac) const
{
    Eigen::NumericalDiff<newtonUnsteadyNSTurbSUPAve> numDiff(*this);
    numDiff.df(x, fjac);
    return 0;
}

// * * * * * * * * * * * * * * * Operators PPE * * * * * * * * * * * * * * * //

// Operator to evaluate the residual for the supremizer approach
int newtonUnsteadyNSTurbPPE::operator()(const Eigen::VectorXd& x,
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
    Eigen::VectorXd m1 = problem->bTotalMatrix * aTmp * nu;
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
    // Penalty term
    Eigen::MatrixXd penaltyU = Eigen::MatrixXd::Zero(Nphi_u, N_BC);

    // Term for penalty method
    if (problem->bcMethod == "penalty")
    {
        for (int l = 0; l < N_BC; l++)
        {
            penaltyU.col(l) = bc(l) * problem->bcVelVec[l] - problem->bcVelMat[l] *
                              aTmp;
        }
    }

    for (int i = 0; i < Nphi_u; i++)
    {
        cc = aTmp.transpose() * Eigen::SliceFromTensor(problem->C_tensor, 0,
                i) * aTmp - gNut.transpose() *
             Eigen::SliceFromTensor(problem->cTotalTensor, 0, i) * aTmp;
        fvec(i) = - m5(i) + m1(i) - cc(0, 0) - m2(i);

        if (problem->bcMethod == "penalty")
        {
            fvec(i) += ((penaltyU * tauU)(i, 0));
        }
    }

    for (int j = 0; j < Nphi_p; j++)
    {
        int k = j + Nphi_u;
        gg = aTmp.transpose() * Eigen::SliceFromTensor(problem->gTensor, 0,
                j) * aTmp;
        bb = aTmp.transpose() * Eigen::SliceFromTensor(problem->bc2Tensor, 0,
                j) * aTmp;
        //fvec(k) = m3(j, 0) - gg(0, 0) - m6(j, 0) + bb(0, 0);
        fvec(k) = m3(j, 0) + gg(0, 0) - m7(j, 0);
    }

    if (problem->bcMethod == "lift")
    {
        for (int j = 0; j < N_BC; j++)
        {
            fvec(j) = x(j) - bc(j);
        }
    }

    return 0;
}

// Operator to evaluate the Jacobian for the supremizer approach
int newtonUnsteadyNSTurbPPE::df(const Eigen::VectorXd& x,
                                Eigen::MatrixXd& fjac) const
{
    Eigen::NumericalDiff<newtonUnsteadyNSTurbPPE> numDiff(*this);
    numDiff.df(x, fjac);
    return 0;
}

int newtonUnsteadyNSTurbPPEAve::operator()(const Eigen::VectorXd& x,
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
    Eigen::MatrixXd nn(1, 1);
    // Mom Term
    Eigen::VectorXd m1 = problem->bTotalMatrix * aTmp * nu;
    // Gradient of pressure
    Eigen::VectorXd m2 = problem->K_matrix * bTmp;
    // Mass Term
    Eigen::VectorXd m5 = problem->M_matrix * a_dot;
    // Time-derivative of the divergence Term
    //Eigen::VectorXd m5 = problem->M_matrix * a_dot;
    // Pressure Term
    Eigen::VectorXd m3 = problem->D_matrix * bTmp;
    // BC PPE
    Eigen::VectorXd m6 = problem->BC1_matrix * aTmp * nu;
    // BC PPE
    Eigen::VectorXd m7 = problem->BC3_matrix * aTmp * nu;
    // Penalty term
    Eigen::MatrixXd penaltyU = Eigen::MatrixXd::Zero(Nphi_u, N_BC);

    // Term for penalty method
    if (problem->bcMethod == "penalty")
    {
        for (int l = 0; l < N_BC; l++)
        {
            penaltyU.col(l) = bc(l) * problem->bcVelVec[l] - problem->bcVelMat[l] *
                              aTmp;
        }
    }

    for (int i = 0; i < Nphi_u; i++)
    {
        cc = aTmp.transpose() * Eigen::SliceFromTensor(problem->C_tensor, 0,
                i) * aTmp - gNut.transpose() *
             Eigen::SliceFromTensor(problem->cTotalTensor, 0,
                                    i) * aTmp - gNutAve.transpose() *
             Eigen::SliceFromTensor(problem->cTotalAveTensor, 0, i) * aTmp;
        fvec(i) = - m5(i) + m1(i) - cc(0, 0) - m2(i);

        if (problem->bcMethod == "penalty")
        {
            fvec(i) += ((penaltyU * tauU)(i, 0));
        }
    }

    for (int j = 0; j < Nphi_p; j++)
    {
        int k = j + Nphi_u;
        gg = aTmp.transpose() * Eigen::SliceFromTensor(problem->gTensor, 0,
                j) * aTmp;
        bb = aTmp.transpose() * Eigen::SliceFromTensor(problem->bc2Tensor, 0,
                j) * aTmp;
        nn = gNut.transpose() *
             Eigen::SliceFromTensor(problem->cTotalPPETensor, 0,
                                    j) * aTmp + gNutAve.transpose() *
             Eigen::SliceFromTensor(problem->cTotalPPEAveTensor, 0, j) * aTmp;
        //fvec(k) = m3(j, 0) + gg(0, 0) - m6(j, 0) + bb(0, 0) - nn(0, 0);
        fvec(k) = m3(j, 0) + gg(0, 0) - m7(j, 0) - nn(0, 0);
        //fvec(k) = m3(j, 0) + gg(0, 0) - m7(j, 0);
    }

    if (problem->bcMethod == "lift")
    {
        for (int j = 0; j < N_BC; j++)
        {
            fvec(j) = x(j) - bc(j);
        }
    }

    return 0;
}

// Operator to evaluate the Jacobian for the PPE approach
int newtonUnsteadyNSTurbPPEAve::df(const Eigen::VectorXd& x,
                                   Eigen::MatrixXd& fjac) const
{
    Eigen::NumericalDiff<newtonUnsteadyNSTurbPPEAve> numDiff(*this);
    numDiff.df(x, fjac);
    return 0;
}


// * * * * * * * * * * * * * * * Solve Functions  * * * * * * * * * * * * * //
void ReducedUnsteadyNSTurb::solveOnlineSUP(Eigen::MatrixXd vel)
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

    if (problem->bcMethod == "lift")
    {
        vel_now = setOnlineVelocity(vel);
    }
    else if (problem->bcMethod == "penalty")
    {
        vel_now = vel;
    }

    // Create and resize the solution vector
    y.resize(Nphi_u + Nphi_p, 1);
    y.setZero();
    // Set Initial Conditions
    y.head(Nphi_u) = initCond.col(0).head(Nphi_u);
    y.tail(Nphi_p) = initCond.col(0).tail(Nphi_p);
    int nextStore = 0;
    int counter2 = 0;

    // Change initial condition for the lifting function
    if (problem->bcMethod == "lift")
    {
        for (int j = 0; j < N_BC; j++)
        {
            y(j) = vel_now(j, 0);
        }
    }

    int firstRBFInd;

    if (skipLift == true && problem->bcMethod == "lift")
    {
        firstRBFInd = N_BC;
    }
    else
    {
        firstRBFInd = 0;
    }

    // Set some properties of the newton object
    newtonObjectSUP.nu = nu;
    newtonObjectSUP.y_old = y;
    newtonObjectSUP.yOldOld = newtonObjectSUP.y_old;
    newtonObjectSUP.dt = dt;
    newtonObjectSUP.bc.resize(N_BC);
    newtonObjectSUP.tauU = tauU;
    newtonObjectSUP.gNut = nut0;

    for (int j = 0; j < N_BC; j++)
    {
        newtonObjectSUP.bc(j) = vel_now(j, 0);
    }

    // Set number of online solutions
    int Ntsteps = static_cast<int>((finalTime - tstart) / dt);
    int onlineSize = static_cast<int>(Ntsteps / numberOfStores);
    online_solution.resize(onlineSize);
    rbfCoeffMat.resize(nphiNut + 1, onlineSize + 3);
    // Set the initial time
    time = tstart;
    // Counting variable
    int counter = 0;
    // Create vector to store temporal solution and save initial condition as first solution
    Eigen::MatrixXd tmp_sol(Nphi_u + Nphi_p + 1, 1);
    tmp_sol(0) = time;
    tmp_sol.col(0).tail(y.rows()) = y;

    // This part up to the while loop is just to compute the eddy viscosity field at time = startTime if time is not equal to zero
    if ((time != 0) || (startFromZero == true))
    {
        online_solution[counter] = tmp_sol;
        counter ++;
        rbfCoeffMat(0, counter2) = time;
        rbfCoeffMat.block(1, counter2, nphiNut, 1) = nut0;
        counter2++;
        nextStore += numberOfStores;
    }

    // Create nonlinear solver object
    Eigen::HybridNonLinearSolver<newtonUnsteadyNSTurbSUP> hnls(newtonObjectSUP);
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
        newtonObjectSUP.operator()(y, res);
        Eigen::VectorXd tv;
        Eigen::VectorXd aDer;
        aDer = (y.head(Nphi_u) - newtonObjectSUP.y_old.head(Nphi_u)) / dt;
        tv.resize(dimA);

        switch (interChoice)
        {
            case 1:
                tv << y.segment(firstRBFInd, dimA);
                break;

            case 2:
                tv << muStar, y.segment(firstRBFInd, dimA - muStar.size());
                break;

            case 3:
                tv << y.segment(firstRBFInd, dimA / 2), aDer.segment(firstRBFInd, dimA / 2);
                break;

            case 4:
                tv << muStar, y.segment(firstRBFInd, (dimA - muStar.size()) / 2),
                aDer.segment(firstRBFInd, (dimA - muStar.size()) / 2);
                break;

            default:
                tv << y.segment(firstRBFInd, dimA);
                break;
        }

        for (int i = 0; i < nphiNut; i++)
        {
            newtonObjectSUP.gNut(i) = problem->rbfSplines[i]->eval(tv);
        }

        // Change initial condition for the lifting function
        if (problem->bcMethod == "lift")
        {
            for (int j = 0; j < N_BC; j++)
            {
                y(j) = vel_now(j, 0);
            }
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

            rbfCoeffMat(0, counter2) = time;
            rbfCoeffMat.block(1, counter2, nphiNut, 1) = newtonObjectSUP.gNut;
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

void ReducedUnsteadyNSTurb::solveOnlineSUPAve(Eigen::MatrixXd vel)
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

    if (problem->bcMethod == "lift")
    {
        vel_now = setOnlineVelocity(vel);
    }
    else if (problem->bcMethod == "penalty")
    {
        vel_now = vel;
    }

    // Create and resize the solution vector
    y.resize(Nphi_u + Nphi_p, 1);
    y.setZero();
    // Set Initial Conditions
    y.head(Nphi_u) = initCond.col(0).head(Nphi_u);
    y.tail(Nphi_p) = initCond.col(0).tail(Nphi_p);
    int nextStore = 0;
    int counter2 = 0;

    // Change initial condition for the lifting function
    if (problem->bcMethod == "lift")
    {
        for (int j = 0; j < N_BC; j++)
        {
            y(j) = vel_now(j, 0);
        }
    }

    int firstRBFInd;

    if (skipLift == true && problem->bcMethod == "lift")
    {
        firstRBFInd = N_BC;
    }
    else
    {
        firstRBFInd = 0;
    }

    // Set some properties of the newton object
    newtonObjectSUPAve.nu = nu;
    newtonObjectSUPAve.y_old = y;
    newtonObjectSUPAve.yOldOld = newtonObjectSUPAve.y_old;
    newtonObjectSUPAve.dt = dt;
    newtonObjectSUPAve.bc.resize(N_BC);
    newtonObjectSUPAve.tauU = tauU;
    newtonObjectSUPAve.gNut = nut0;
    newtonObjectSUPAve.gNutAve = gNutAve;

    for (int j = 0; j < N_BC; j++)
    {
        newtonObjectSUPAve.bc(j) = vel_now(j, 0);
    }

    // Set number of online solutions
    int Ntsteps = static_cast<int>((finalTime - tstart) / dt);
    int onlineSize = static_cast<int>(Ntsteps / numberOfStores);
    online_solution.resize(onlineSize);
    rbfCoeffMat.resize(nphiNut + 1, onlineSize + 3);
    // Set the initial time
    time = tstart;
    // Counting variable
    int counter = 0;
    // Create vector to store temporal solution and save initial condition as first solution
    Eigen::MatrixXd tmp_sol(Nphi_u + Nphi_p + 1, 1);
    tmp_sol(0) = time;
    tmp_sol.col(0).tail(y.rows()) = y;

    // This part up to the while loop is just to compute the eddy viscosity field at time = startTime if time is not equal to zero
    if ((time != 0) || (startFromZero == true))
    {
        online_solution[counter] = tmp_sol;
        counter ++;
        rbfCoeffMat(0, counter2) = time;
        rbfCoeffMat.block(1, counter2, nphiNut, 1) = nut0;
        counter2++;
        nextStore += numberOfStores;
    }

    // Create nonlinear solver object
    Eigen::HybridNonLinearSolver<newtonUnsteadyNSTurbSUPAve> hnls(
        newtonObjectSUPAve);
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
        newtonObjectSUPAve.operator()(y, res);
        Eigen::VectorXd tv;
        Eigen::VectorXd aDer;
        aDer = (y.head(Nphi_u) - newtonObjectSUPAve.y_old.head(Nphi_u)) / dt;
        tv.resize(dimA);

        switch (interChoice)
        {
            case 1:
                tv << y.segment(firstRBFInd, dimA);
                break;

            case 2:
                tv << muStar, y.segment(firstRBFInd, dimA - muStar.size());
                break;

            case 3:
                tv << y.segment(firstRBFInd, dimA / 2), aDer.segment(firstRBFInd, dimA / 2);
                break;

            case 4:
                tv << muStar, y.segment(firstRBFInd, (dimA - muStar.size()) / 2),
                aDer.segment(firstRBFInd, (dimA - muStar.size()) / 2);
                break;

            default:
                tv << y.segment(firstRBFInd, dimA);
                break;
        }

        for (int i = 0; i < nphiNut; i++)
        {
            newtonObjectSUPAve.gNut(i) = problem->rbfSplines[i]->eval(tv);
        }

        // Change initial condition for the lifting function
        if (problem->bcMethod == "lift")
        {
            for (int j = 0; j < N_BC; j++)
            {
                y(j) = vel_now(j, 0);
            }
        }

        newtonObjectSUPAve.yOldOld = newtonObjectSUPAve.y_old;
        newtonObjectSUPAve.y_old = y;
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

            rbfCoeffMat(0, counter2) = time;
            rbfCoeffMat.block(1, counter2, nphiNut, 1) = newtonObjectSUPAve.gNut;
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
void ReducedUnsteadyNSTurb::solveOnlinePPE(Eigen::MatrixXd vel)
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

    if (problem->bcMethod == "lift")
    {
        vel_now = setOnlineVelocity(vel);
    }
    else if (problem->bcMethod == "penalty")
    {
        vel_now = vel;
    }

    // Create and resize the solution vector
    y.resize(Nphi_u + Nphi_p, 1);
    y.setZero();
    // Set Initial Conditions
    y.head(Nphi_u) = initCond.col(0).head(Nphi_u);
    y.tail(Nphi_p) = initCond.col(0).tail(Nphi_p);
    int nextStore = 0;
    int counter2 = 0;

    // Change initial condition for the lifting function
    if (problem->bcMethod == "lift")
    {
        for (int j = 0; j < N_BC; j++)
        {
            y(j) = vel_now(j, 0);
        }
    }

    int firstRBFInd;

    if (skipLift == true && problem->bcMethod == "lift")
    {
        firstRBFInd = N_BC;
    }
    else
    {
        firstRBFInd = 0;
    }

    // Set some properties of the newton object
    newtonObjectPPE.nu = nu;
    newtonObjectPPE.y_old = y;
    newtonObjectPPE.yOldOld = newtonObjectPPE.y_old;
    newtonObjectPPE.dt = dt;
    newtonObjectPPE.bc.resize(N_BC);
    newtonObjectPPE.tauU = tauU;
    newtonObjectPPE.gNut = nut0;

    for (int j = 0; j < N_BC; j++)
    {
        newtonObjectPPE.bc(j) = vel_now(j, 0);
    }

    // Set number of online solutions
    int Ntsteps = static_cast<int>((finalTime - tstart) / dt);
    int onlineSize = static_cast<int>(Ntsteps / numberOfStores);
    online_solution.resize(onlineSize);
    rbfCoeffMat.resize(nphiNut + 1, onlineSize + 3);
    // Set the initial time
    time = tstart;
    // Counting variable
    int counter = 0;
    // Create vectpr to store temporal solution and save initial condition as first solution
    Eigen::MatrixXd tmp_sol(Nphi_u + Nphi_p + 1, 1);
    tmp_sol(0) = time;
    tmp_sol.col(0).tail(y.rows()) = y;

    // This part up to the while loop is just to compute the eddy viscosity field at time = startTime if time is not equal to zero
    if ((time != 0) || (startFromZero == true))
    {
        online_solution[counter] = tmp_sol;
        counter ++;
        rbfCoeffMat(0, counter2) = time;
        rbfCoeffMat.block(1, counter2, nphiNut, 1) = nut0;
        counter2++;
        nextStore += numberOfStores;
    }

    // Create nonlinear solver object
    Eigen::HybridNonLinearSolver<newtonUnsteadyNSTurbPPE> hnls(newtonObjectPPE);
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
        newtonObjectPPE.operator()(y, res);
        Eigen::VectorXd tv;
        Eigen::VectorXd aDer;
        aDer = (y.head(Nphi_u) - newtonObjectPPE.y_old.head(Nphi_u)) / dt;
        tv.resize(dimA);

        switch (interChoice)
        {
            case 1:
                tv << y.segment(firstRBFInd, dimA);
                break;

            case 2:
                tv << muStar, y.segment(firstRBFInd, dimA - muStar.size());
                break;

            case 3:
                tv << y.segment(firstRBFInd, dimA / 2), aDer.segment(firstRBFInd, dimA / 2);
                break;

            case 4:
                tv << muStar, y.segment(firstRBFInd, (dimA - muStar.size()) / 2),
                aDer.segment(firstRBFInd, (dimA - muStar.size()) / 2);
                break;

            default:
                tv << y.segment(firstRBFInd, dimA);
                break;
        }

        for (int i = 0; i < nphiNut; i++)
        {
            newtonObjectPPE.gNut(i) = problem->rbfSplines[i]->eval(tv);
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

            rbfCoeffMat(0, counter2) = time;
            rbfCoeffMat.block(1, counter2, nphiNut, 1) = newtonObjectPPE.gNut;
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

void ReducedUnsteadyNSTurb::solveOnlinePPEAve(Eigen::MatrixXd vel)
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

    if (problem->bcMethod == "lift")
    {
        vel_now = setOnlineVelocity(vel);
    }
    else if (problem->bcMethod == "penalty")
    {
        vel_now = vel;
    }

    // Create and resize the solution vector
    y.resize(Nphi_u + Nphi_p, 1);
    y.setZero();
    // Set Initial Conditions
    y.head(Nphi_u) = initCond.col(0).head(Nphi_u);
    y.tail(Nphi_p) = initCond.col(0).tail(Nphi_p);
    int nextStore = 0;
    int counter2 = 0;

    // Change initial condition for the lifting function
    if (problem->bcMethod == "lift")
    {
        for (int j = 0; j < N_BC; j++)
        {
            y(j) = vel_now(j, 0);
        }
    }

    int firstRBFInd;

    if (skipLift == true && problem->bcMethod == "lift")
    {
        firstRBFInd = N_BC;
    }
    else
    {
        firstRBFInd = 0;
    }

    // Set some properties of the newton object
    newtonObjectPPEAve.nu = nu;
    newtonObjectPPEAve.y_old = y;
    newtonObjectPPEAve.yOldOld = newtonObjectPPEAve.y_old;
    newtonObjectPPEAve.dt = dt;
    newtonObjectPPEAve.bc.resize(N_BC);
    newtonObjectPPEAve.tauU = tauU;
    newtonObjectPPEAve.gNut = nut0;
    newtonObjectPPEAve.gNutAve = gNutAve;

    for (int j = 0; j < N_BC; j++)
    {
        newtonObjectPPEAve.bc(j) = vel_now(j, 0);
    }

    // Set number of online solutions
    int Ntsteps = static_cast<int>((finalTime - tstart) / dt);
    int onlineSize = static_cast<int>(Ntsteps / numberOfStores);
    online_solution.resize(onlineSize);
    rbfCoeffMat.resize(nphiNut + 1, onlineSize + 3);
    // Set the initial time
    time = tstart;
    // Counting variable
    int counter = 0;
    // Create vectpr to store temporal solution and save initial condition as first solution
    Eigen::MatrixXd tmp_sol(Nphi_u + Nphi_p + 1, 1);
    tmp_sol(0) = time;
    tmp_sol.col(0).tail(y.rows()) = y;

    // This part up to the while loop is just to compute the eddy viscosity field at time = startTime if time is not equal to zero
    if ((time != 0) || (startFromZero == true))
    {
        online_solution[counter] = tmp_sol;
        counter ++;
        rbfCoeffMat(0, counter2) = time;
        rbfCoeffMat.block(1, counter2, nphiNut, 1) = nut0;
        counter2++;
        nextStore += numberOfStores;
    }

    // Create nonlinear solver object
    Eigen::HybridNonLinearSolver<newtonUnsteadyNSTurbPPEAve> hnls(
        newtonObjectPPEAve);
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
        newtonObjectPPEAve.operator()(y, res);
        Eigen::VectorXd tv;
        Eigen::VectorXd aDer;
        aDer = (y.head(Nphi_u) - newtonObjectPPEAve.y_old.head(Nphi_u)) / dt;
        tv.resize(dimA);

        switch (interChoice)
        {
            case 1:
                tv << y.segment(firstRBFInd, dimA);
                break;

            case 2:
                tv << muStar, y.segment(firstRBFInd, dimA - muStar.size());
                break;

            case 3:
                tv << y.segment(firstRBFInd, dimA / 2), aDer.segment(firstRBFInd, dimA / 2);
                break;

            case 4:
                tv << muStar, y.segment(firstRBFInd, (dimA - muStar.size()) / 2),
                aDer.segment(firstRBFInd, (dimA - muStar.size()) / 2);
                break;

            default:
                tv << y.segment(firstRBFInd, dimA);
                break;
        }

        for (int i = 0; i < nphiNut; i++)
        {
            newtonObjectPPEAve.gNut(i) = problem->rbfSplines[i]->eval(tv);
        }

        // Change initial condition for the lifting function
        if (problem->bcMethod == "lift")
        {
            for (int j = 0; j < N_BC; j++)
            {
                y(j) = vel_now(j, 0);
            }
        }

        newtonObjectPPEAve.yOldOld = newtonObjectPPEAve.y_old;
        newtonObjectPPEAve.y_old = y;
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

            rbfCoeffMat(0, counter2) = time;
            rbfCoeffMat.block(1, counter2, nphiNut, 1) = newtonObjectPPEAve.gNut;
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

void ReducedUnsteadyNSTurb::reconstruct(bool exportFields, fileName folder)
{
    if (exportFields)
    {
        mkDir(folder);
        ITHACAutilities::createSymLink(folder);
    }

    int counter = 0;
    int nextWrite = 0;
    int counter2 = 1;
    int exportEveryIndex = round(exportEvery / storeEvery);
    volScalarField nutAveNow("nutAveNow", nutModes[0] * 0);
    List < Eigen::MatrixXd> CoeffU;
    List < Eigen::MatrixXd> CoeffP;
    List < Eigen::MatrixXd> CoeffNut;
    CoeffU.resize(0);
    CoeffP.resize(0);
    CoeffNut.resize(0);

    for (int k = 0; k < problem->nutAve.size(); k++)
    {
        nutAveNow += gNutAve(k) * problem->nutAve[k];
    }

    for (int i = 0; i < online_solution.size(); i++)
    {
        if (counter == nextWrite)
        {
            Eigen::MatrixXd currentUCoeff;
            Eigen::MatrixXd currentPCoeff;
            Eigen::MatrixXd currentNutCoeff;
            currentUCoeff = online_solution[i].block(1, 0, Nphi_u, 1);
            currentPCoeff = online_solution[i].bottomRows(Nphi_p);
            currentNutCoeff = rbfCoeffMat.block(1, i, nphiNut, 1);
            CoeffU.append(currentUCoeff);
            CoeffP.append(currentPCoeff);
            CoeffNut.append(currentNutCoeff);
            nextWrite += exportEveryIndex;
        }

        counter++;
    }

    volVectorField uRec("uRec", Umodes[0]);
    volScalarField pRec("pRec", problem->Pmodes[0]);
    volScalarField nutRec("nutRec", problem->nutModes[0]);
    uRecFields = problem->L_U_SUPmodes.reconstruct(uRec, CoeffU, "uRec");
    pRecFields = problem->Pmodes.reconstruct(pRec, CoeffP, "pRec");
    nutRecFields = problem->nutModes.reconstruct(nutRec, CoeffNut, "nutRec");

    for (int k = 0; k < nutRecFields.size(); k++)
    {
        nutRecFields[k] += nutAveNow;
    }

    if (exportFields)
    {
        ITHACAstream::exportFields(uRecFields, folder,
                                   "uRec");
        ITHACAstream::exportFields(pRecFields, folder,
                                   "pRec");
        ITHACAstream::exportFields(nutRecFields, folder,
                                   "nutRec");
    }
}

Eigen::MatrixXd ReducedUnsteadyNSTurb::setOnlineVelocity(Eigen::MatrixXd vel)
{
    assert(problem->inletIndex.rows() == vel.rows()
           && "Imposed boundary conditions dimensions do not match given values matrix dimensions");
    Eigen::MatrixXd vel_scal;
    vel_scal.resize(vel.rows(), vel.cols());

    for (int k = 0; k < problem->inletIndex.rows(); k++)
    {
        int p = problem->inletIndex(k, 0);
        int l = problem->inletIndex(k, 1);
        scalar area = gSum(problem->liftfield[0].mesh().magSf().boundaryField()[p]);
        scalar u_lf = gSum(problem->liftfield[k].mesh().magSf().boundaryField()[p] *
                           problem->liftfield[k].boundaryField()[p]).component(l) / area;
        vel_scal(k, 0) = vel(k, 0) / u_lf;
    }

    return vel_scal;
}
// ************************************************************************* //

