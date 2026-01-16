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
    N_BC    = problem->inletIndex.rows();
    Nphi_u  = problem->B_matrix.rows();     // DO NOT use NUmodes; B_matrix.rows() is correct
    Nphi_p  = problem->K_matrix.cols();
    interChoice = problem->interChoice;

    nphiNutAvg   = problem->nutAvgModes.size();
    nphiNutFluct = problem->nutFluctModes.size();
    nphiNut      = nphiNutAvg + nphiNutFluct;

    // Local velocity modes (lift + physical + SUP, for full y)
    Umodes.clear();
    for (int k = 0; k < problem->liftfield.size(); ++k)
        Umodes.append(problem->liftfield[k].clone());
    for (int k = 0; k < problem->NUmodes; ++k)
        Umodes.append(problem->Umodes[k].clone());
    for (int k = 0; k < problem->NSUPmodes; ++k)
        Umodes.append(problem->supmodes[k].clone());

    Pmodes.clear();
    for (int k = 0; k < problem->NPmodes; ++k)
        Pmodes.append(problem->Pmodes[k].clone());

    nutAvgModes.clear();
    for (int k = 0; k < problem->nutAvgModes.size(); ++k)
        nutAvgModes.append(problem->nutAvgModes[k].clone());
    nutFluctModes.clear();
    for (int k = 0; k < problem->nutFluctModes.size(); ++k)
        nutFluctModes.append(problem->nutFluctModes[k].clone());

    // Copy RBF splines
    rbfSplinesNutAvg   = problem->rbfSplinesNutAvg;
    rbfSplinesNutFluct = problem->rbfSplinesNutFluct;
    samplesNutAvg      = problem->samplesNutAvg;
    samplesNutFluct    = problem->samplesNutFluct;

    coeffNutAvg.resize(nphiNutAvg);
    coeffNutFluct.resize(nphiNutFluct);

    // Set dimA = input size for fluctuation RBF (typically 2*physical modes)
    // If original was using dimA = problem->velRBF.cols(), you should set accordingly
    // For dual: aNow (Nphys), aDotNow (Nphys)
    dimA = 2* Nphi_u;

    // Newton objects (as in original)
    newtonObjectSUP     = newtonUnsteadyNSTurbSUP(Nphi_u + Nphi_p, Nphi_u + Nphi_p, fomProblem);
    newtonObjectPPE     = newtonUnsteadyNSTurbPPE(Nphi_u + Nphi_p, Nphi_u + Nphi_p, fomProblem);
    newtonObjectSUPAve  = newtonUnsteadyNSTurbSUPAve(Nphi_u + Nphi_p, Nphi_u + Nphi_p, fomProblem);
    newtonObjectPPEAve  = newtonUnsteadyNSTurbPPEAve(Nphi_u + Nphi_p, Nphi_u + Nphi_p, fomProblem);

    std::cout << "[DUAL-ROM] Initialized. Nphi_u (with BCs) = " << Nphi_u
              << ", Nphi_p = " << Nphi_p
              << ", nphiNutAvg = " << nphiNutAvg
              << ", nphiNutFluct = " << nphiNutFluct
              << ", dimA = " << dimA << std::endl;
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

// ===== SUP-AVG residual (use gNutFluct with cTotalFluctTensor) =====
int newtonUnsteadyNSTurbSUPAve::operator()(const Eigen::VectorXd& x,
                                           Eigen::VectorXd& fvec) const
{
    // unpack
    Eigen::VectorXd aTmp = x.head(Nphi_u);
    Eigen::VectorXd bTmp = x.tail(Nphi_p);

    // time derivative
    Eigen::VectorXd a_dot(Nphi_u);
    if (problem->timeDerivativeSchemeOrder == "first")
        a_dot = (x.head(Nphi_u) - y_old.head(Nphi_u)) / dt;
    else
        a_dot = (1.5 * x.head(Nphi_u) - 2.0 * y_old.head(Nphi_u) + 0.5 * yOldOld.head(Nphi_u)) / dt;

    // core operators
    Eigen::VectorXd m1 = problem->bTotalMatrix * aTmp * nu;  // ν B a (incl. bt)
    Eigen::VectorXd m2 = problem->K_matrix     * bTmp;       // grad p
    Eigen::VectorXd m5 = problem->M_matrix     * a_dot;      // mass
    Eigen::VectorXd m3 = problem->P_matrix     * aTmp;       // SUP pressure block

    // penalty (unchanged)
    Eigen::MatrixXd penaltyU = Eigen::MatrixXd::Zero(Nphi_u, N_BC);
    if (problem->bcMethod == "penalty")
        for (int l = 0; l < N_BC; ++l)
            penaltyU.col(l) = bc(l) * problem->bcVelVec[l] - problem->bcVelMat[l] * aTmp;

    // momentum residual
    Eigen::MatrixXd cc(1,1);
    for (int i = 0; i < Nphi_u; ++i)
    {
        cc =  aTmp.transpose() * Eigen::SliceFromTensor(problem->C_tensor, 0, i) * aTmp
            // use FLUCT here (replaces old gNut · cTotalTensor)
            - gNutFluct.transpose()
              * Eigen::SliceFromTensor(problem->cTotalFluctTensor, 0, i) * aTmp
            // keep AVG term unchanged
            - gNutAve.transpose()
              * Eigen::SliceFromTensor(problem->cTotalAveTensor,   0, i) * aTmp;

        fvec(i) = -m5(i) + m1(i) - cc(0,0) - m2(i);

        if (problem->bcMethod == "penalty")
            fvec(i) += ((penaltyU * tauU)(i, 0));
    }

    // SUP “pressure” equations (unchanged)
    for (int j = 0; j < Nphi_p; ++j)
    {
        const int k = j + Nphi_u;
        fvec(k) = m3(j);
    }

    // lift overwrite (unchanged)
    if (problem->bcMethod == "lift")
        for (int j = 0; j < N_BC; ++j)
            fvec(j) = x(j) - bc(j);

    return 0;
}

// Jacobian (unchanged)
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

    if (problem->timeDerivativeSchemeOrder == "first")
    {
        a_dot = (x.head(Nphi_u) - y_old.head(Nphi_u)) / dt;
    }
    else
    {
        a_dot = (1.5 * x.head(Nphi_u) - 2.0 * y_old.head(Nphi_u) + 0.5 * yOldOld.head(Nphi_u)) / dt;
    }

    // Temporary 1x1 holders for quadratic/bilinear forms
    Eigen::MatrixXd cc(1, 1);
    Eigen::MatrixXd gg(1, 1);
    Eigen::MatrixXd bb(1, 1);
    Eigen::MatrixXd nn(1, 1);

    // Momentum-side vectors
    Eigen::VectorXd m1 = problem->bTotalMatrix * aTmp * nu;   // ν B a
    Eigen::VectorXd m2 = problem->K_matrix     * bTmp;        // K b   (grad p)
    Eigen::VectorXd m5 = problem->M_matrix     * a_dot;       // M ȧ

    // PPE-side vectors/matrices
    Eigen::VectorXd m3 = problem->D_matrix     * bTmp;        // D b   (Δp term)
    Eigen::VectorXd m6 = problem->BC1_matrix   * aTmp * nu;   // (not used here)
    Eigen::VectorXd m7 = problem->BC3_matrix   * aTmp * nu;   // ν N a

    // Penalty term (if used)
    Eigen::MatrixXd penaltyU = Eigen::MatrixXd::Zero(Nphi_u, N_BC);
    if (problem->bcMethod == "penalty")
    {
        for (int l = 0; l < N_BC; ++l)
        {
            penaltyU.col(l) = bc(l) * problem->bcVelVec[l] - problem->bcVelMat[l] * aTmp;
        }
    }

    // ===== Momentum equations residual =====
    for (int i = 0; i < Nphi_u; ++i)
    {
        cc = aTmp.transpose() * Eigen::SliceFromTensor(problem->C_tensor, 0, i) * aTmp
           - gNutAve.transpose()   * Eigen::SliceFromTensor(problem->cTotalAveTensor,   0, i) * aTmp
           - gNutFluct.transpose() * Eigen::SliceFromTensor(problem->cTotalFluctTensor, 0, i) * aTmp;

        // f_u = -M ȧ + ν B a - (nonlinear conv + ν_t terms) - K b
        fvec(i) = - m5(i) + m1(i) - cc(0, 0) - m2(i);

        if (problem->bcMethod == "penalty")
        {
            fvec(i) += ((penaltyU * tauU)(i, 0));
        }
    }

    // ===== Pressure (PPE) equations residual =====
    // Thesis: D b + aᵀ G a - ν N a - L = 0
    for (int j = 0; j < Nphi_p; ++j)
    {
        const int k = j + Nphi_u;

        gg = aTmp.transpose() * Eigen::SliceFromTensor(problem->gTensor, 0, j) * aTmp; // aᵀ G_j a
        bb = aTmp.transpose() * Eigen::SliceFromTensor(problem->bc2Tensor, 0, j) * aTmp; // (kept computed; not used per thesis)
        nn = gNutAve.transpose()   * Eigen::SliceFromTensor(problem->cTotalPPEAveTensor,   0, j) * aTmp
           + gNutFluct.transpose() * Eigen::SliceFromTensor(problem->cTotalPPEFluctTensor, 0, j) * aTmp;

        // Safe access to L(j): if not present or wrong-sized, use 0.0
        double Lj = 0.0;
        if (problem->L_vector.rows() == Nphi_p && problem->L_vector.cols() >= 1)
        {
            Lj = problem->L_vector(j, 0);
        }

        // f_p = D b + aᵀ G a - ν N a - (ν_t PPE terms) - L
        fvec(k) = m3(j, 0) + gg(0, 0) - m7(j, 0) - nn(0, 0) - Lj;
    }

    // Lift BC overwrite on velocity DOFs, as in your current implementation
    if (problem->bcMethod == "lift")
    {
        for (int j = 0; j < N_BC; ++j)
        {
            fvec(j) = x(j) - bc(j);
        }
    }

    return 0;
}

// Jacobian via numerical differentiation (unchanged)
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
    std::cout << "\n=== ENTERING solveOnlineSUP ===" << std::endl;

    M_Assert(exportEvery >= dt, "dt must be smaller than exportEvery");
    M_Assert(storeEvery >= dt, "dt must be smaller than storeEvery");
    M_Assert(ITHACAutilities::isInteger(storeEvery / dt), "storeEvery must be a multiple of dt");
    M_Assert(ITHACAutilities::isInteger(exportEvery / dt), "exportEvery must be a multiple of dt");
    M_Assert(ITHACAutilities::isInteger(exportEvery / storeEvery), "exportEvery must be a multiple of storeEvery");

    int numberOfStores = round(storeEvery / dt);
    std::cout << "[DEBUG] storeEvery = " << storeEvery << ", dt = " << dt
              << ", numberOfStores = " << numberOfStores << std::endl;

    std::cout << "[DEBUG] vel shape: (" << vel.rows() << ", " << vel.cols() << ")\n";
    std::cout << "[DEBUG] vel input:\n" << vel << std::endl;

    vel_now = vel;  // Direct assignment
    std::cout << "[DEBUG] vel_now set directly:\n" << vel_now << std::endl;

    if (problem->bcMethod == "lift")
        std::cout << "[DEBUG] Using lifting function for BCs\n";
    else if (problem->bcMethod == "penalty")
        std::cout << "[DEBUG] Using penalty method\n";

    // Initial y
    y.resize(Nphi_u + Nphi_p, 1);
    y.setZero();
    y.head(Nphi_u) = ITHACAutilities::getCoeffs(problem->Ufield[0], Umodes);
    y.tail(Nphi_p) = ITHACAutilities::getCoeffs(problem->Pfield[0], Pmodes);
    nut0 = ITHACAutilities::getCoeffs(problem->nutFields[0], nutModes);

    std::cout << "[DEBUG] Nphi_u = " << Nphi_u << ", Nphi_p = " << Nphi_p << std::endl;
    std::cout << "[DEBUG] y.size() = " << y.rows() << " x " << y.cols() << std::endl;
    std::cout << "[DEBUG] Initial y.head(Nphi_u): " << y.head(Nphi_u).transpose() << std::endl;

    // Update with vel_now for lifting
    if (problem->bcMethod == "lift")
    {
        for (int j = 0; j < N_BC; ++j)
            y(j) = vel_now(j, 0);
        std::cout << "[DEBUG] Updated y(0:N_BC) with vel_now\n";
    }

    const int firstRBFInd = (skipLift && problem->bcMethod == "lift") ? N_BC : 0;
    newtonObjectSUP.nu = nu;
    newtonObjectSUP.y_old = y;
    newtonObjectSUP.yOldOld = y;
    newtonObjectSUP.dt = dt;
    newtonObjectSUP.bc.resize(N_BC);
    newtonObjectSUP.tauU = tauU;
    newtonObjectSUP.gNut = nut0;

    for (int j = 0; j < N_BC; ++j)
        newtonObjectSUP.bc(j) = vel_now(j, 0);

    const int Ntsteps = static_cast<int>((finalTime - tstart) / dt);
    const int onlineSize = static_cast<int>(Ntsteps / numberOfStores);
    std::cout << "[DEBUG] Ntsteps = " << Ntsteps << ", onlineSize = " << onlineSize << std::endl;

    online_solution.resize(onlineSize);
    rbfCoeffMat.resize(nphiNut + 1, onlineSize + 3);
    std::cout << "[DEBUG] rbfCoeffMat shape = " << rbfCoeffMat.rows() << " x " << rbfCoeffMat.cols() << "\n";

    time = tstart;
    int counter = 0, counter2 = 0, nextStore = 0;
    Eigen::MatrixXd tmp_sol(Nphi_u + Nphi_p + 1, 1);

    tmp_sol(0) = time;
    tmp_sol.col(0).tail(y.rows()) = y;

    if ((time != 0) || (startFromZero == true))
    {
        online_solution[counter] = tmp_sol;
        rbfCoeffMat(0, counter2) = time;
        rbfCoeffMat.block(1, counter2, nphiNut, 1) = nut0;
        ++counter;
        ++counter2;
        nextStore += numberOfStores;
    }

    Eigen::HybridNonLinearSolver<newtonUnsteadyNSTurbSUP> hnls(newtonObjectSUP);
    Color::Modifier red(Color::FG_RED), green(Color::FG_GREEN), def(Color::FG_DEFAULT);

    Eigen::VectorXd muStar(1);
    muStar(0) = vel_now(0, 0);

    while (time < finalTime)
    {
        time += dt;
        std::cout << "\n[DEBUG] ==== TIME STEP @ t = " << time << " ====" << std::endl;

        Eigen::VectorXd res = Eigen::VectorXd::Zero(y.rows());
        hnls.solve(y);
        newtonObjectSUP.operator()(y, res);

        Eigen::VectorXd aDer = (y.head(Nphi_u) - newtonObjectSUP.y_old.head(Nphi_u)) / dt;
        Eigen::VectorXd tv(dimA);

        std::cout << "[DEBUG] Constructing tv with interChoice = " << interChoice << std::endl;

        switch (interChoice)
        {
            case 1:
                tv = y.segment(firstRBFInd, dimA);
                break;
            case 2:
                tv.head(1) = muStar;
                tv.tail(dimA - 1) = y.segment(firstRBFInd, dimA - 1);
                break;
            case 3:
                tv.head(dimA / 2) = y.segment(firstRBFInd, dimA / 2);
                tv.tail(dimA / 2) = aDer.segment(firstRBFInd, dimA / 2);
                break;
            case 4:
                {
                    const int half = (dimA - 1) / 2;
                    tv(0) = muStar(0);
                    tv.segment(1, half) = y.segment(firstRBFInd, half);
                    tv.segment(1 + half, half) = aDer.segment(firstRBFInd, half);
                }
                break;
            default:
                FatalErrorInFunction << "Invalid interChoice: " << interChoice << abort(FatalError);
        }

        for (int i = 0; i < nphiNut; ++i)
            newtonObjectSUP.gNut(i) = problem->rbfSplines[i]->eval(tv);

        if (problem->bcMethod == "lift")
        {
            for (int j = 0; j < N_BC; ++j)
                y(j) = vel_now(j, 0);
        }

        newtonObjectSUP.operator()(y, res);
        newtonObjectSUP.yOldOld = newtonObjectSUP.y_old;
        newtonObjectSUP.y_old = y;

        std::cout << "################## Online solve N° " << count_online_solve << " ##################\n";
        std::cout << "Time = " << time << ", Parameter: " << vel_now(0, 0) << std::endl;

        if (res.norm() < 1e-5)
            std::cout << green << "|F(x)| = " << res.norm() << " - Converged in " << hnls.iter << " iterations" << def << std::endl;
        else
            std::cout << red << "|F(x)| = " << res.norm() << " - Not converged in " << hnls.iter << " iterations" << def << std::endl;

        ++count_online_solve;

        tmp_sol(0) = time;
        tmp_sol.col(0).tail(y.rows()) = y;

        if (counter == nextStore)
        {
            if (counter2 >= online_solution.size())
                online_solution.append(tmp_sol);
            else
                online_solution[counter2] = tmp_sol;

            rbfCoeffMat(0, counter2) = time;
            rbfCoeffMat.block(1, counter2, nphiNut, 1) = newtonObjectSUP.gNut;

            nextStore += numberOfStores;
            ++counter2;
        }

        ++counter;
    }

    ITHACAstream::exportMatrix(online_solution, "red_coeff", "python", "./ITHACAoutput/red_coeff");
    ITHACAstream::exportMatrix(online_solution, "red_coeff", "matlab", "./ITHACAoutput/red_coeff");

    ++count_online_solve;
    std::cout << "=== Exiting solveOnlineSUP ===\n" << std::endl;
}

// ================== ONLINE SOLVER: SUP + (AVG linear-μ, FLUCT RBF) ==================
void ReducedUnsteadyNSTurb::solveOnlineSUPAve(Eigen::MatrixXd vel)
{
    using std::cout; using std::endl;
    #include <fstream>
    #include <iomanip>

    cout << "\n\n[DEBUG] === Entering solveOnlineSUPAve() (AVG linear-μ + FLUCT RBF) ===\n";

    // ---- open debug logs (overwrite each run) ----
    std::ofstream avgLog("./ITHACAoutput/Coefficients/avg_coeffs_online_log_SUP.txt");
    std::ofstream flLog ("./ITHACAoutput/Coefficients/fluct_coeffs_online_log_SUP.txt");
    avgLog << std::setprecision(15);
    flLog  << std::setprecision(15);
    avgLog << "# step time mu*  gNutAvg[0..nphiNutAvg-1]\n";
    flLog  << "# step time      gNutFluct[0..nphiNutFluct-1]\n";

    // --- Setup/checks ---
    cout << "[DEBUG] dt=" << dt << ", storeEvery=" << storeEvery << ", exportEvery=" << exportEvery << endl;
    M_Assert(exportEvery >= dt && storeEvery >= dt, "dt must be <= exportEvery/storeEvery");
    M_Assert(ITHACAutilities::isInteger(storeEvery / dt), "storeEvery must be integer multiple of dt.");
    M_Assert(ITHACAutilities::isInteger(exportEvery / dt), "exportEvery must be integer multiple of dt.");
    M_Assert(ITHACAutilities::isInteger(exportEvery / storeEvery), "exportEvery must be integer multiple of storeEvery.");
    const int numberOfStores = round(storeEvery / dt);
    cout << "[DEBUG] numberOfStores = " << numberOfStores << endl;

    // --- Velocity input/BCs ---
    if (problem->bcMethod == "lift") {
        vel_now = setOnlineVelocity(vel);
        cout << "[DEBUG] vel_now set via setOnlineVelocity(vel), shape = (" << vel_now.rows() << "," << vel_now.cols() << ")\n";
    } else {
        vel_now = vel;
        cout << "[DEBUG] vel_now = vel, shape = (" << vel_now.rows() << "," << vel_now.cols() << ")\n";
    }

    // --- Initial projection (t=0, FOM) ---
    cout << "[DEBUG] Projecting Ufield[0] and Pfield[0] onto modes (t=0)\n";
    y.resize(Nphi_u + Nphi_p, 1);
    y.setZero();
    y.head(Nphi_u) = ITHACAutilities::getCoeffs(problem->Uomfield[0], Umodes);
    cout << "[DEBUG] y.head(Nphi_u) = " << y.head(Nphi_u).transpose() << endl;
    y.tail(Nphi_p) = ITHACAutilities::getCoeffs(problem->Pfield[0], Pmodes);
    cout << "[DEBUG] y.tail(Nphi_p) = " << y.tail(Nphi_p).transpose() << endl;

    if (problem->bcMethod == "lift") {
        for (int j = 0; j < N_BC; j++) {
            y(j) = vel_now(j, 0);
            cout << "[DEBUG] Applied lift BC: y(" << j << ") = " << y(j) << endl;
        }
    }

    // --- PHYSICAL DOFs (exclude BC, exclude SUP) for RBF fluct input ---
    const bool useLift = (problem->bcMethod == "lift");
    const int firstRBFInd = useLift ? N_BC : 0;
    const int Nphys       = problem->NUmodes;           // only physical velocity modes
    cout << "[DEBUG] firstRBFInd = " << firstRBFInd
         << ", Nphys = " << Nphys
         << ", Nphi_u (lift+phys+SUP) = " << Nphi_u << "\n";
    M_Assert(firstRBFInd + Nphys <= Nphi_u, "RBF input slice exceeds velocity DOF range.");

    // --- RBF inputs ---
    Eigen::VectorXd muStar(1); muStar(0) = vel(0, 0);
    const double muStarVal = muStar(0);
    cout << "[DEBUG] muStar = " << muStar.transpose() << endl;

    // === Load AVG interpolation tables (once) ===
    static bool avgTblLoaded = false;
    static Eigen::VectorXd muAvg;        // [M]
    static Eigen::MatrixXd Cavg;         // [nphiNutAvg x M]
    if (!avgTblLoaded) {
        const word coeffDir = "./ITHACAoutput/Coefficients/";
        Eigen::MatrixXd muAvgMat = ITHACAstream::readMatrix(coeffDir + "NutAvg_mu_unique_mat.txt");
        Cavg                    = ITHACAstream::readMatrix(coeffDir + "NutAvg_coeffs_by_mu_mat.txt");

        // Convert muAvgMat to a column VectorXd (handle 1xM or Mx1)
        if (muAvgMat.rows() == 1 && muAvgMat.cols() > 1)
            muAvg = muAvgMat.transpose();
        else
            muAvg = muAvgMat.col(0);

        M_Assert(Cavg.rows() == nphiNutAvg, "NutAvg coeff table row mismatch vs nphiNutAvg.");
        M_Assert(Cavg.cols() == muAvg.size(), "NutAvg table columns must match μ-grid length.");

        avgTblLoaded = true;
        cout << "[DEBUG] Loaded NutAvg tables: M=" << muAvg.size()
             << "  Cavg=" << Cavg.rows() << "x" << Cavg.cols() << endl;
    }

    // === Helper: linear interpolation in μ (for AVG only) ===
    auto interpNutAvgLinear = [&](double mu)->Eigen::VectorXd {
        const int M = static_cast<int>(muAvg.size());
        if (mu <= muAvg(0))   return Cavg.col(0);
        if (mu >= muAvg(M-1)) return Cavg.col(M-1);
        int j = 1; while (j < M && muAvg(j) < mu) ++j;    // first μ_j >= mu
        const double muL = muAvg(j-1), muR = muAvg(j);
        const double t = (mu - muL) / (muR - muL);
        return (1.0 - t)*Cavg.col(j-1) + t*Cavg.col(j);
    };

    // aOld/aNow/aDot for PHYSICAL part only (exclude lift & SUP)
    Eigen::VectorXd aOld    = y.segment(firstRBFInd, Nphys);
    Eigen::VectorXd aNow    = aOld;
    Eigen::VectorXd aDotNow = Eigen::VectorXd::Zero(Nphys);

    // Prepare 2*Nphys input vector for fluctuation RBF
    Eigen::VectorXd rbfInputFluct(2 * Nphys);
    rbfInputFluct.head(Nphys) = aNow;
    rbfInputFluct.tail(Nphys) = aDotNow;

    auto printHead = [](const Eigen::VectorXd& v, int k)
    {
        std::stringstream ss; ss << std::setprecision(6);
        k = std::min<int>(k, v.size());
        for (int i = 0; i < k; ++i) { ss << v(i) << (i+1<k ? ", " : ""); }
        return ss.str();
    };
    auto vecMin = [](const Eigen::VectorXd& v){ return v.minCoeff(); };
    auto vecMax = [](const Eigen::VectorXd& v){ return v.maxCoeff(); };

    // --- Dual evaluations at t=0 ---
    Eigen::VectorXd coeffNutAvg(nphiNutAvg), coeffNutFluct(nphiNutFluct);

    // AVG via linear μ-interpolation
    coeffNutAvg = interpNutAvgLinear(muStarVal);

    // FLUCT at t=0: zero (consistent with many offline setups)
    coeffNutFluct.setZero();

    cout << "[DEBUG] t=0: ||gAvg||2="   << coeffNutAvg.norm()
         << "  min/max=(" << vecMin(coeffNutAvg)   << "," << vecMax(coeffNutAvg)   << ")  head=["
         << printHead(coeffNutAvg, 5) << "]\n";
    cout << "[DEBUG] t=0: ||gFluct||2=" << coeffNutFluct.norm()
         << "  min/max=(" << vecMin(coeffNutFluct) << "," << vecMax(coeffNutFluct) << ")  head=["
         << printHead(coeffNutFluct, 5) << "]\n";

    avgLog << 0 << " " << tstart << " " << muStarVal << " ";
    for (int k = 0; k < coeffNutAvg.size(); ++k) avgLog << coeffNutAvg(k) << (k+1<coeffNutAvg.size() ? ' ' : '\n');
    flLog  << 0 << " " << tstart << " ";
    for (int k = 0; k < coeffNutFluct.size(); ++k) flLog  << coeffNutFluct(k) << (k+1<coeffNutFluct.size() ? ' ' : '\n');

    // --- Setup Newton object, BCs, outputs ---
    newtonObjectSUPAve.nu        = nu;
    newtonObjectSUPAve.y_old     = y;
    newtonObjectSUPAve.yOldOld   = y;
    newtonObjectSUPAve.dt        = dt;
    newtonObjectSUPAve.bc.resize(N_BC);
    newtonObjectSUPAve.tauU      = tauU;
    newtonObjectSUPAve.gNutAve   = coeffNutAvg;
    newtonObjectSUPAve.gNutFluct = coeffNutFluct;   // <<< CHANGED: use gNutFluct
    for (int j = 0; j < N_BC; j++)
        newtonObjectSUPAve.bc(j) = vel_now(j, 0);

    const int Ntsteps    = static_cast<int>((finalTime - tstart) / dt);
    const int onlineSize = static_cast<int>(Ntsteps / numberOfStores) + 3; // slack
    online_solution.resize(onlineSize);
    rbfCoeffMat.setZero(1 + nphiNutAvg + nphiNutFluct, onlineSize);

    // --- Export t=0 solution ---
    double time = tstart;
    int step = 0, nextStore = 0, counter2 = 0, counter = 0;
    Eigen::VectorXd tmp_sol(Nphi_u + Nphi_p + 1);
    tmp_sol(0) = time;
    tmp_sol.tail(y.rows()) = y;
    online_solution[counter] = tmp_sol;
    rbfCoeffMat(0, counter2) = time;
    rbfCoeffMat.block(1, counter2, nphiNutAvg, 1) = coeffNutAvg;
    rbfCoeffMat.block(1 + nphiNutAvg, counter2, nphiNutFluct, 1) = coeffNutFluct;
    cout << "[DEBUG] t=0: Stored solution and RBF coeffs at counter2=" << counter2 << endl;
    counter2++; counter++; nextStore += numberOfStores;

    // --- Create nonlinear solver object
    Eigen::HybridNonLinearSolver<newtonUnsteadyNSTurbSUPAve> hnls(newtonObjectSUPAve);

    // === Main time loop ===
    cout << "\n[DEBUG] == ENTERING MAIN TIME LOOP t>0 ==\n";
    while (time < finalTime)
    {
        time += dt;
        step++;
        cout << "\n---------------------------------------------\n";
        cout << "[DEBUG] === Time step: t=" << time << ", step=" << step << " ===\n";

        problem->timeDerivativeSchemeOrder = "first";
        cout << "[DEBUG] y (before Newton) = " << y.transpose() << endl;

        // Nonlinear solve
        Eigen::VectorXd res(y); res.setZero();
        hnls.solve(y);
        cout << "[DEBUG] y (after  Newton) = " << y.transpose() << endl;
        newtonObjectSUPAve.operator()(y, res);

        // --- PHYSICAL coefficients only for RBF fluctuation input ---
        Eigen::VectorXd aNowFull = y.head(Nphi_u);
        Eigen::VectorXd aOldFull = newtonObjectSUPAve.y_old.head(Nphi_u);
        Eigen::VectorXd aNowPhys = aNowFull.segment(firstRBFInd, Nphys);
        Eigen::VectorXd aDotPhys = (aNowFull.segment(firstRBFInd, Nphys) - aOldFull.segment(firstRBFInd, Nphys)) / dt;

        rbfInputFluct.head(Nphys) = aNowPhys;
        rbfInputFluct.tail(Nphys) = aDotPhys;

        cout << "[DEBUG] aNow_phys = " << aNowPhys.transpose() << endl;
        cout << "[DEBUG] aDot_phys = " << aDotPhys.transpose() << endl;

        // --- Dual evaluation ---
        // AVG via linear μ-interpolation
        coeffNutAvg = interpNutAvgLinear(muStarVal);

        // FLUCT via RBF([a,aDot])
        for (int i = 0; i < nphiNutFluct; i++) {
            coeffNutFluct(i) = rbfSplinesNutFluct[i]->eval(rbfInputFluct);
        }

        // <<< CHANGED: assign to gNutFluct (not gNut)
        newtonObjectSUPAve.gNutAve   = coeffNutAvg;
        newtonObjectSUPAve.gNutFluct = coeffNutFluct;

        // ---- per-step debug summary to console ----
        cout << std::setprecision(6)
             << "[DEBUG] coeffNutAvg:  ||·||2=" << coeffNutAvg.norm()
             << "  min/max=(" << vecMin(coeffNutAvg) << "," << vecMax(coeffNutAvg)
             << ")  head=[" << printHead(coeffNutAvg, 5) << "]\n";
        cout << "[DEBUG] coeffNutFluct:||·||2=" << coeffNutFluct.norm()
             << "  min/max=(" << vecMin(coeffNutFluct) << "," << vecMax(coeffNutFluct)
             << ")  head=[" << printHead(coeffNutFluct, 5) << "]\n";

        // ---- full vectors to log files every step ----
        avgLog << step << " " << time << " " << muStarVal << " ";
        for (int k = 0; k < coeffNutAvg.size(); ++k)
            avgLog << coeffNutAvg(k) << (k+1<coeffNutAvg.size() ? ' ' : '\n');

        flLog << step << " " << time << " ";
        for (int k = 0; k < coeffNutFluct.size(); ++k)
            flLog << coeffNutFluct(k) << (k+1<coeffNutFluct.size() ? ' ' : '\n');

        // --- Re-apply lift BCs in y (keeps stored solution consistent)
        if (useLift) {
            for (int j = 0; j < N_BC; j++) y(j) = vel_now(j, 0);
        }

        // Update old states
        newtonObjectSUPAve.yOldOld = newtonObjectSUPAve.y_old;
        newtonObjectSUPAve.y_old   = y;

        // --- Store/Export schedule ---
        tmp_sol(0) = time;
        tmp_sol.tail(y.rows()) = y;
        if (counter == nextStore)
        {
            if (counter2 >= online_solution.size()) {
                online_solution.append(tmp_sol);
            } else {
                online_solution[counter2] = tmp_sol;
            }

            if (counter2 < rbfCoeffMat.cols()) {
                rbfCoeffMat(0, counter2) = time;
                rbfCoeffMat.block(1, counter2, nphiNutAvg, 1) = coeffNutAvg;
                rbfCoeffMat.block(1 + nphiNutAvg, counter2, nphiNutFluct, 1) = coeffNutFluct;
            }

            cout << "[DEBUG] Stored solution and RBF coeffs at counter2=" << counter2 << endl;
            nextStore += numberOfStores;
            counter2++;
        }
        counter++;
    }

    // --- Trim to actually stored columns BEFORE export
    const int storedCols = counter2;
    cout << "[DEBUG] Trimming outputs: storedCols=" << storedCols
         << ", rbfCoeffMat orig cols=" << rbfCoeffMat.cols()
         << ", online_solution orig size=" << online_solution.size() << endl;

    if (storedCols > 0 && storedCols <= rbfCoeffMat.cols()) {
        rbfCoeffMat.conservativeResize(rbfCoeffMat.rows(), storedCols);
    }
    if (storedCols >= 0 && storedCols <= online_solution.size()) {
        online_solution.resize(storedCols);
    }

    // --- Export
    cout << "\n[DEBUG] Exporting online_solution and rbfCoeffMat...\n";
    ITHACAstream::exportMatrix(online_solution, "red_coeff_SUP",        "python", "./ITHACAoutput/red_coeff");
    ITHACAstream::exportMatrix(online_solution, "red_coeff_SUP",        "matlab", "./ITHACAoutput/red_coeff");
    ITHACAstream::exportMatrix(rbfCoeffMat,     "rbf_coeffs_online_SUP","python", "./ITHACAoutput/red_coeff");
    ITHACAstream::exportMatrix(rbfCoeffMat,     "rbf_coeffs_online_SUP","matlab", "./ITHACAoutput/red_coeff");

    // close logs
    avgLog.close();
    flLog.close();

    cout << "[DEBUG] solveOnlineSUPAve() completed.\n\n";
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
    using std::cout; using std::endl;
    #include <fstream>
    #include <iomanip>

    cout << "\n\n[DEBUG] === Entering solveOnlinePPEAve() (AVG linear-μ + FLUCT RBF) ===\n";

    // ---- open debug logs (overwrite each run) ----
    std::ofstream avgLog("./ITHACAoutput/Coefficients/avg_coeffs_online_log.txt");
    std::ofstream flLog ("./ITHACAoutput/Coefficients/fluct_coeffs_online_log.txt");
    avgLog << std::setprecision(15);
    flLog  << std::setprecision(15);
    avgLog << "# step time mu*  gNutAvg[0..nphiNutAvg-1]\n";
    flLog  << "# step time      gNutFluct[0..nphiNutFluct-1]\n";

    // --- Setup/checks ---
    cout << "[DEBUG] dt=" << dt << ", storeEvery=" << storeEvery << ", exportEvery=" << exportEvery << endl;
    M_Assert(exportEvery >= dt && storeEvery >= dt, "dt must be <= exportEvery/storeEvery");
    M_Assert(ITHACAutilities::isInteger(storeEvery / dt), "storeEvery must be integer multiple of dt.");
    M_Assert(ITHACAutilities::isInteger(exportEvery / dt), "exportEvery must be integer multiple of dt.");
    M_Assert(ITHACAutilities::isInteger(exportEvery / storeEvery), "exportEvery must be integer multiple of storeEvery.");
    const int numberOfStores = round(storeEvery / dt);
    cout << "[DEBUG] numberOfStores = " << numberOfStores << endl;

    // --- Velocity input/BCs ---
    if (problem->bcMethod == "lift") {
        vel_now = setOnlineVelocity(vel);
        cout << "[DEBUG] vel_now set via setOnlineVelocity(vel), shape = (" << vel_now.rows() << "," << vel_now.cols() << ")\n";
    } else {
        vel_now = vel;
        cout << "[DEBUG] vel_now = vel, shape = (" << vel_now.rows() << "," << vel_now.cols() << ")\n";
    }

    // --- Initial projection (t=0, FOM) ---
    cout << "[DEBUG] Projecting Ufield[0] and Pfield[0] onto modes (t=0)\n";
    y.resize(Nphi_u + Nphi_p, 1);
    y.setZero();
    y.head(Nphi_u) = ITHACAutilities::getCoeffs(problem->Uomfield[0], Umodes);
    cout << "[DEBUG] y.head(Nphi_u) = " << y.head(Nphi_u).transpose() << endl;
    y.tail(Nphi_p) = ITHACAutilities::getCoeffs(problem->Pfield[0], Pmodes);
    cout << "[DEBUG] y.tail(Nphi_p) = " << y.tail(Nphi_p).transpose() << endl;
    cout << "[DEBUG] y (full) = " << y.transpose() << endl;

    if (problem->bcMethod == "lift") {
        for (int j = 0; j < N_BC; j++) {
            y(j) = vel_now(j, 0);
            cout << "[DEBUG] Applied lift BC: y(" << j << ") = " << y(j) << endl;
        }
    }

    // --- PHYSICAL DOFs (exclude BC and SUP for RBF fluct input) ---
    const bool useLift = (problem->bcMethod == "lift");
    const int firstRBFInd = useLift ? N_BC : 0;
    const int Nphys      = problem->NUmodes;
    cout << "[DEBUG] firstRBFInd = " << firstRBFInd
         << ", Nphys = " << Nphys
         << ", Nphi_u (lift+phys+SUP) = " << Nphi_u << "\n";

    M_Assert(firstRBFInd + Nphys <= Nphi_u,
             "RBF input slice exceeds velocity DOF range. Check N_BC/NUmodes/Nphi_u.");

    // --- RBF inputs ---
    Eigen::VectorXd muStar(1); muStar(0) = vel(0, 0);
    const double muStarVal = muStar(0);
    cout << "[DEBUG] muStar = " << muStar.transpose() << endl;

    // === Load AVG interpolation tables (once) ===
    //     (keeps FLUCT path unchanged)
    static bool avgTblLoaded = false;
    static Eigen::VectorXd muAvg;        // [M]
    static Eigen::MatrixXd Cavg;         // [nphiNutAvg x M]
    if (!avgTblLoaded) {
        const word coeffDir = "./ITHACAoutput/Coefficients/";
        Eigen::MatrixXd muAvgMat = ITHACAstream::readMatrix(coeffDir + "NutAvg_mu_unique_mat.txt");
        Cavg                    = ITHACAstream::readMatrix(coeffDir + "NutAvg_coeffs_by_mu_mat.txt");

        // Convert muAvgMat to a column VectorXd (handle 1xM or Mx1)
        if (muAvgMat.rows() == 1 && muAvgMat.cols() > 1)
            muAvg = muAvgMat.transpose();
        else
            muAvg = muAvgMat.col(0);

        M_Assert(Cavg.rows() == nphiNutAvg, "NutAvg coeff table row mismatch vs nphiNutAvg.");
        M_Assert(Cavg.cols() == muAvg.size(), "NutAvg table columns must match mu grid length.");

        avgTblLoaded = true;
        cout << "[DEBUG] Loaded NutAvg tables: M=" << muAvg.size()
             << "  Cavg=" << Cavg.rows() << "x" << Cavg.cols() << endl;
    }

    // === Helper: linear interpolation in μ (for AVG only) ===
    auto interpNutAvgLinear = [&](double mu)->Eigen::VectorXd {
        const int M = static_cast<int>(muAvg.size());
        if (mu <= muAvg(0))   return Cavg.col(0);
        if (mu >= muAvg(M-1)) return Cavg.col(M-1);
        int j = 1; while (j < M && muAvg(j) < mu) ++j;    // first μ_j >= mu
        const double muL = muAvg(j-1), muR = muAvg(j);
        const double t = (mu - muL) / (muR - muL);
        // (Optional tiny debug — comment out if noisy)
        cout << "[DEBUG] AVG μ-bracket: jL=" << (j-1) << " jR=" << j
             << "  muL=" << muL << " muR=" << muR << "  t=" << t << endl;
        return (1.0 - t)*Cavg.col(j-1) + t*Cavg.col(j);
    };

    // aOld/aNow/aDot for PHYSICAL part only
    Eigen::VectorXd aOld    = y.segment(firstRBFInd, Nphys);
    Eigen::VectorXd aNow    = aOld;
    Eigen::VectorXd aDotNow = Eigen::VectorXd::Zero(Nphys);

    // Prepare 2*Nphys input vector for fluctuation RBF
    Eigen::VectorXd rbfInputFluct(2 * Nphys);
    rbfInputFluct.head(Nphys) = aNow;
    rbfInputFluct.tail(Nphys) = aDotNow;

    cout << "[DEBUG] t=0: aOld_phys = " << aOld.transpose() << endl;
    cout << "[DEBUG] t=0: aDot_phys = " << aDotNow.transpose() << endl;

    // --- Dual evaluations at t=0 ---
    Eigen::VectorXd coeffNutAvg(nphiNutAvg), coeffNutFluct(nphiNutFluct);

    // *** AVG via linear μ-interpolation (replaces RBF call) ***
    coeffNutAvg = interpNutAvgLinear(muStarVal);

    // (FLUCT unchanged)
    coeffNutFluct.setZero(); // consistent with offline phase choice

    auto printHead = [](const Eigen::VectorXd& v, int k)
    {
        std::stringstream ss; ss << std::setprecision(6);
        k = std::min<int>(k, v.size());
        for (int i = 0; i < k; ++i) { ss << v(i) << (i+1<k ? ", " : ""); }
        return ss.str();
    };
    auto vecMin = [](const Eigen::VectorXd& v){ return v.minCoeff(); };
    auto vecMax = [](const Eigen::VectorXd& v){ return v.maxCoeff(); };

    cout << "[DEBUG] t=0: ||gAvg||2="   << coeffNutAvg.norm()
         << "  min/max=(" << vecMin(coeffNutAvg)   << "," << vecMax(coeffNutAvg)   << ")  head=["
         << printHead(coeffNutAvg, 5) << "]\n";
    cout << "[DEBUG] t=0: ||gFluct||2=" << coeffNutFluct.norm()
         << "  min/max=(" << vecMin(coeffNutFluct) << "," << vecMax(coeffNutFluct) << ")  head=["
         << printHead(coeffNutFluct, 5) << "]\n";

    avgLog << 0 << " " << tstart << " " << muStarVal << " ";
    for (int k = 0; k < coeffNutAvg.size(); ++k) avgLog << coeffNutAvg(k) << (k+1<coeffNutAvg.size() ? ' ' : '\n');
    flLog  << 0 << " " << tstart << " ";
    for (int k = 0; k < coeffNutFluct.size(); ++k) flLog  << coeffNutFluct(k) << (k+1<coeffNutFluct.size() ? ' ' : '\n');

    newtonObjectPPEAve.gNutAve   = coeffNutAvg;
    newtonObjectPPEAve.gNutFluct = coeffNutFluct;

    // --- Setup Newton object, BCs, outputs ---
    newtonObjectPPEAve.nu      = nu;
    newtonObjectPPEAve.y_old   = y;
    newtonObjectPPEAve.yOldOld = y;
    newtonObjectPPEAve.dt      = dt;
    newtonObjectPPEAve.bc.resize(N_BC);
    newtonObjectPPEAve.tauU    = tauU;
    newtonObjectPPEAve.gNutAve.resize(nphiNutAvg);
    newtonObjectPPEAve.gNutFluct.resize(nphiNutFluct);
    for (int j = 0; j < N_BC; j++)
        newtonObjectPPEAve.bc(j) = vel_now(j, 0);

    const int Ntsteps    = static_cast<int>((finalTime - tstart) / dt);
    const int onlineSize = static_cast<int>(Ntsteps / numberOfStores) + 3; // slack
    online_solution.resize(onlineSize);
    rbfCoeffMat.setZero(1 + nphiNutAvg + nphiNutFluct, onlineSize);

    // --- Export t=0 solution ---
    double time = tstart;
    int step = 0, nextStore = 0, counter2 = 0, counter = 0;
    Eigen::VectorXd tmp_sol(Nphi_u + Nphi_p + 1);
    tmp_sol(0) = time;
    tmp_sol.tail(y.rows()) = y;
    online_solution[counter] = tmp_sol;
    rbfCoeffMat(0, counter2) = time;
    rbfCoeffMat.block(1, counter2, nphiNutAvg, 1) = coeffNutAvg;
    rbfCoeffMat.block(1 + nphiNutAvg, counter2, nphiNutFluct, 1) = coeffNutFluct;
    cout << "[DEBUG] t=0: Stored solution and RBF coeffs at counter2=" << counter2 << endl;
    counter2++; counter++; nextStore += numberOfStores;

    // --- Create nonlinear solver object
    Eigen::HybridNonLinearSolver<newtonUnsteadyNSTurbPPEAve> hnls(newtonObjectPPEAve);

    // === Main time loop ===
    cout << "\n[DEBUG] == ENTERING MAIN TIME LOOP t>0 ==\n";
    while (time < finalTime)
    {
        time += dt;
        step++;
        cout << "\n---------------------------------------------\n";
        cout << "[DEBUG] === Time step: t=" << time << ", step=" << step << " ===\n";

        problem->timeDerivativeSchemeOrder = "first";
        cout << "[DEBUG] y (before Newton) = " << y.transpose() << endl;

        // Nonlinear solve
        Eigen::VectorXd res(y); res.setZero();
        hnls.solve(y);
        cout << "[DEBUG] y (after  Newton) = " << y.transpose() << endl;
        newtonObjectPPEAve.operator()(y, res);

        // --- PHYSICAL coefficients only for RBF fluctuation input ---
        aNow    = y.segment(firstRBFInd, Nphys);
        aDotNow = (aNow - aOld) / dt;

        rbfInputFluct.head(Nphys) = aNow;
        rbfInputFluct.tail(Nphys) = aDotNow;

        cout << "[DEBUG] aNow_phys = " << aNow.transpose() << endl;
        cout << "[DEBUG] aDot_phys = " << aDotNow.transpose() << endl;

        // --- Dual evaluation ---
        // *** AVG via linear μ-interpolation (no RBF) ***
        coeffNutAvg = interpNutAvgLinear(muStarVal);

        // (FLUCT unchanged)
        for (int i = 0; i < nphiNutFluct; i++) {
            coeffNutFluct(i) = rbfSplinesNutFluct[i]->eval(rbfInputFluct);
        }

        newtonObjectPPEAve.gNutAve   = coeffNutAvg;
        newtonObjectPPEAve.gNutFluct = coeffNutFluct;

        // ---- per-step debug summary to console ----
        cout << std::setprecision(6)
             << "[DEBUG] coeffNutAvg:  ||·||2=" << coeffNutAvg.norm()
             << "  min/max=(" << vecMin(coeffNutAvg) << "," << vecMax(coeffNutAvg)
             << ")  head=[" << printHead(coeffNutAvg, 5) << "]\n";
        cout << "[DEBUG] coeffNutFluct:||·||2=" << coeffNutFluct.norm()
             << "  min/max=(" << vecMin(coeffNutFluct) << "," << vecMax(coeffNutFluct)
             << ")  head=[" << printHead(coeffNutFluct, 5) << "]\n";

        // ---- full vectors to log files every step ----
        avgLog << step << " " << time << " " << muStarVal << " ";
        for (int k = 0; k < coeffNutAvg.size(); ++k)
            avgLog << coeffNutAvg(k) << (k+1<coeffNutAvg.size() ? ' ' : '\n');

        flLog << step << " " << time << " ";
        for (int k = 0; k < coeffNutFluct.size(); ++k)
            flLog << coeffNutFluct(k) << (k+1<coeffNutFluct.size() ? ' ' : '\n');

        // --- Re-apply lift BCs in y (keeps stored solution consistent)
        if (useLift) {
            for (int j = 0; j < N_BC; j++) y(j) = vel_now(j, 0);
        }

        // Update old states
        newtonObjectPPEAve.yOldOld = newtonObjectPPEAve.y_old;
        newtonObjectPPEAve.y_old   = y;
        aOld = aNow;

        // --- Store/Export schedule ---
        tmp_sol(0) = time;
        tmp_sol.tail(y.rows()) = y;
        if (counter == nextStore)
        {
            if (counter2 >= online_solution.size()) {
                online_solution.append(tmp_sol);
            } else {
                online_solution[counter2] = tmp_sol;
            }

            if (counter2 < rbfCoeffMat.cols()) {
                rbfCoeffMat(0, counter2) = time;
                rbfCoeffMat.block(1, counter2, nphiNutAvg, 1) = coeffNutAvg;
                rbfCoeffMat.block(1 + nphiNutAvg, counter2, nphiNutFluct, 1) = coeffNutFluct;
            } else {
                cout << "[WARN] rbfCoeffMat is full; skipping store at t=" << time << endl;
            }

            cout << "[DEBUG] Stored solution and RBF coeffs at counter2=" << counter2 << endl;
            nextStore += numberOfStores;
            counter2++;
        }
        counter++;
    }

    // --- Trim to actually stored columns BEFORE export
    const int storedCols = counter2;
    cout << "[DEBUG] Trimming outputs: storedCols=" << storedCols
         << ", rbfCoeffMat orig cols=" << rbfCoeffMat.cols()
         << ", online_solution orig size=" << online_solution.size() << endl;

    if (storedCols > 0 && storedCols <= rbfCoeffMat.cols()) {
        rbfCoeffMat.conservativeResize(rbfCoeffMat.rows(), storedCols);
    }
    if (storedCols >= 0 && storedCols <= online_solution.size()) {
        online_solution.resize(storedCols);
    }

    // --- Export
    cout << "\n[DEBUG] Exporting online_solution and rbfCoeffMat...\n";
    ITHACAstream::exportMatrix(online_solution, "red_coeff",          "python", "./ITHACAoutput/red_coeff");
    ITHACAstream::exportMatrix(online_solution, "red_coeff",          "matlab", "./ITHACAoutput/red_coeff");
    ITHACAstream::exportMatrix(rbfCoeffMat,     "rbf_coeffs_online",  "python", "./ITHACAoutput/red_coeff");
    ITHACAstream::exportMatrix(rbfCoeffMat,     "rbf_coeffs_online",  "matlab", "./ITHACAoutput/red_coeff");

    // close logs
    avgLog.close();
    flLog.close();

    cout << "[DEBUG] solveOnlinePPEAve() completed.\n\n";
}


void ReducedUnsteadyNSTurb::reconstruct(bool exportFields, fileName folder)
{
    if (exportFields) { mkDir(folder); ITHACAutilities::createSymLink(folder); }

    const int nphiNutAvg   = nutAvgModes.size();
    const int nphiNutFluct = nutFluctModes.size();
    const int nphiU        = Nphi_u;
    const int nphiP        = Nphi_p;

    const int nSteps = std::min<int>(online_solution.size(), rbfCoeffMat.cols());

    // FILL CLASS MEMBERS
    uRecFields.clear();  pRecFields.clear();  nutRecFields.clear();
    uRecFields.setSize(nSteps);
    pRecFields.setSize(nSteps);
    nutRecFields.setSize(nSteps);

    for (int i = 0; i < nSteps; ++i)
    {
        const Eigen::VectorXd &yi = online_solution[i];
        Eigen::VectorXd currentUCoeff = yi.block(1, 0, nphiU, 1);
        Eigen::VectorXd currentPCoeff = yi.block(1 + nphiU, 0, nphiP, 1);

        Eigen::VectorXd nut_avg_coeffs   = rbfCoeffMat.block(1,              i, nphiNutAvg,   1);
        Eigen::VectorXd nut_fluct_coeffs = rbfCoeffMat.block(1 + nphiNutAvg, i, nphiNutFluct, 1);

        // U
        auto* uPtr = new volVectorField("uRec", Umodes[0] * 0.0);
        for (int j = 0; j < nphiU; ++j) { *uPtr += currentUCoeff(j) * Umodes[j]; }
        uRecFields.set(i, uPtr);

        // P
        auto* pPtr = new volScalarField("pRec", Pmodes[0] * 0.0);
        for (int j = 0; j < nphiP; ++j) { *pPtr += currentPCoeff(j) * Pmodes[j]; }
        pRecFields.set(i, pPtr);

        // nut = avg + fluct
        auto* nutAvgPtr   = new volScalarField("nutRecAvg",   nutAvgModes[0]   * 0.0);
        for (int j = 0; j < nphiNutAvg; ++j)   { *nutAvgPtr   += nut_avg_coeffs(j)   * nutAvgModes[j]; }
        auto* nutFluctPtr = new volScalarField("nutRecFluct", nutFluctModes[0] * 0.0);
        for (int j = 0; j < nphiNutFluct; ++j) { *nutFluctPtr += nut_fluct_coeffs(j) * nutFluctModes[j]; }

        auto* nutPtr = new volScalarField("nutRec", (*nutAvgPtr) + (*nutFluctPtr));
        delete nutAvgPtr; delete nutFluctPtr; // keep only final field
        nutRecFields.set(i, nutPtr);
    }

    if (exportFields) {
        ITHACAstream::exportFields(uRecFields,   folder, "uRec");
        ITHACAstream::exportFields(pRecFields,   folder, "pRec");
        ITHACAstream::exportFields(nutRecFields, folder, "nutRec");
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
