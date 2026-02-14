/*---------------------------------------------------------------------------*\
     ██╗████████╗██╗  ██╗ █████╗  ██████╗ █████╗       ███████╗██╗   ██╗
     ██║╚══██╔══╝██║  ██║██╔══██╗██╔════╝██╔══██╗      ██╔════╝██║   ██║
     ██║   ██║   ███████║███████║██║     ███████║█████╗█████╗  ██║   ██║
     ██║   ██║   ██╔══██║██╔══██║██║     ██╔══██║╚════╝██╔══╝  ╚██╗ ██╔╝
     ██║   ██║   ██║  ██║██║  ██║╚██████╗██║  ██║      ██║      ╚████╔╝
     ╚═╝   ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝      ╚═╝       ╚═══╝
 *
 * In real Time Highly Advanced Computational Applications for Finite Volumes
 * Copyright (C) 2017 by the ITHACA-FV authors
 *-------------------------------------------------------------------------------
 *
 * License
 *     This file is part of ITHACA-FV
 *
 *     ITHACA-FV is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     ITHACA-FV is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *     GNU Lesser General Public License for more details.
 *
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with ITHACA-FV. If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

/// \file
/// Source file of the reducedUnsteadyNS class

#include "ReducedUnsteadyNSTurb.H"

#include <cassert>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

ReducedUnsteadyNSTurb::ReducedUnsteadyNSTurb() = default;

ReducedUnsteadyNSTurb::ReducedUnsteadyNSTurb(UnsteadyNSTurb& fomProblem)
{
    problem      = &fomProblem;
    N_BC         = problem->inletIndex.rows();
    Nphi_u       = problem->B_matrix.rows();  // DO NOT use NUmodes; B_matrix.rows() is correct
    Nphi_p       = problem->K_matrix.cols();
    interChoice  = problem->interChoice;

    nphiNutAvg   = problem->nutAvgModes.size();
    nphiNutFluct = problem->nutFluctModes.size();
    nphiNut      = nphiNutAvg + nphiNutFluct;

    // Local velocity modes (lift + physical + SUP, for full y)
    Umodes.clear();
    for (int k = 0; k < problem->liftfield.size(); ++k)
    {
        Umodes.append(problem->liftfield[k].clone());
    }
    for (int k = 0; k < problem->NUmodes; ++k)
    {
        Umodes.append(problem->Umodes[k].clone());
    }
    for (int k = 0; k < problem->NSUPmodes; ++k)
    {
        Umodes.append(problem->supmodes[k].clone());
    }

    Pmodes.clear();
    for (int k = 0; k < problem->NPmodes; ++k)
    {
        Pmodes.append(problem->Pmodes[k].clone());
    }

    nutAvgModes.clear();
    for (int k = 0; k < problem->nutAvgModes.size(); ++k)
    {
        nutAvgModes.append(problem->nutAvgModes[k].clone());
    }

    nutFluctModes.clear();
    for (int k = 0; k < problem->nutFluctModes.size(); ++k)
    {
        nutFluctModes.append(problem->nutFluctModes[k].clone());
    }

    // Copy RBF splines
    rbfSplinesNutAvg   = problem->rbfSplinesNutAvg;
    rbfSplinesNutFluct = problem->rbfSplinesNutFluct;
    samplesNutAvg      = problem->samplesNutAvg;
    samplesNutFluct    = problem->samplesNutFluct;

    coeffNutAvg.resize(nphiNutAvg);
    coeffNutFluct.resize(nphiNutFluct);

    // Set dimA = input size for fluctuation RBF (typically 2*physical modes)
    // For dual: aNow (Nphys), aDotNow (Nphys)
    dimA = 2 * Nphi_u;

    // Newton objects (as in original)
    newtonObjectSUP    = newtonUnsteadyNSTurbSUP(Nphi_u + Nphi_p, Nphi_u + Nphi_p, fomProblem);
    newtonObjectPPE    = newtonUnsteadyNSTurbPPE(Nphi_u + Nphi_p, Nphi_u + Nphi_p, fomProblem);
    newtonObjectSUPAve = newtonUnsteadyNSTurbSUPAve(Nphi_u + Nphi_p, Nphi_u + Nphi_p, fomProblem);
    newtonObjectPPEAve = newtonUnsteadyNSTurbPPEAve(Nphi_u + Nphi_p, Nphi_u + Nphi_p, fomProblem);

    std::cout << "[DUAL-ROM] Initialized. Nphi_u (with BCs) = " << Nphi_u
              << ", Nphi_p = " << Nphi_p
              << ", nphiNutAvg = " << nphiNutAvg
              << ", nphiNutFluct = " << nphiNutFluct
              << ", dimA = " << dimA << std::endl;
}

// * * * * * * * * * * * * * * Operators supremizer * * * * * * * * * * * * * //

int newtonUnsteadyNSTurbSUP::operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const
{
    Eigen::VectorXd aTmp = x.head(Nphi_u);
    Eigen::VectorXd bTmp = x.tail(Nphi_p);

    Eigen::VectorXd a_dot(Nphi_u);
    if (problem->timeDerivativeSchemeOrder == "first")
    {
        a_dot = (x.head(Nphi_u) - y_old.head(Nphi_u)) / dt;
    }
    else
    {
        a_dot = (1.5 * x.head(Nphi_u) - 2.0 * y_old.head(Nphi_u) + 0.5 * yOldOld.head(Nphi_u)) / dt;
    }

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

    if (problem->bcMethod == "penalty")
    {
        for (int l = 0; l < N_BC; ++l)
        {
            penaltyU.col(l) = bc(l) * problem->bcVelVec[l] - problem->bcVelMat[l] * aTmp;
        }
    }

    for (int i = 0; i < Nphi_u; ++i)
    {
        cc = aTmp.transpose() * Eigen::SliceFromTensor(problem->C_tensor, 0, i) * aTmp
           - gNut.transpose() * Eigen::SliceFromTensor(problem->cTotalTensor, 0, i) * aTmp;

        fvec(i) = -m5(i) + m1(i) - cc(0, 0) - m2(i);

        if (problem->bcMethod == "penalty")
        {
            fvec(i) += ((penaltyU * tauU)(i, 0));
        }
    }

    for (int j = 0; j < Nphi_p; ++j)
    {
        const int k = j + Nphi_u;
        fvec(k) = m3(j);
    }

    if (problem->bcMethod == "lift")
    {
        for (int j = 0; j < N_BC; ++j)
        {
            fvec(j) = x(j) - bc(j);
        }
    }

    return 0;
}

int newtonUnsteadyNSTurbSUP::df(const Eigen::VectorXd& x, Eigen::MatrixXd& fjac) const
{
    Eigen::NumericalDiff<newtonUnsteadyNSTurbSUP> numDiff(*this);
    numDiff.df(x, fjac);
    return 0;
}

// ===== SUP-AVG residual (use gNutFluct with cTotalFluctTensor) =====
int newtonUnsteadyNSTurbSUPAve::operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const
{
    Eigen::VectorXd aTmp = x.head(Nphi_u);
    Eigen::VectorXd bTmp = x.tail(Nphi_p);

    Eigen::VectorXd a_dot(Nphi_u);
    if (problem->timeDerivativeSchemeOrder == "first")
    {
        a_dot = (x.head(Nphi_u) - y_old.head(Nphi_u)) / dt;
    }
    else
    {
        a_dot = (1.5 * x.head(Nphi_u) - 2.0 * y_old.head(Nphi_u) + 0.5 * yOldOld.head(Nphi_u)) / dt;
    }

    Eigen::VectorXd m1 = problem->bTotalMatrix * aTmp * nu;  // ν B a (incl. bt)
    Eigen::VectorXd m2 = problem->K_matrix * bTmp;           // grad p
    Eigen::VectorXd m5 = problem->M_matrix * a_dot;          // mass
    Eigen::VectorXd m3 = problem->P_matrix * aTmp;           // SUP pressure block

    Eigen::MatrixXd penaltyU = Eigen::MatrixXd::Zero(Nphi_u, N_BC);
    if (problem->bcMethod == "penalty")
    {
        for (int l = 0; l < N_BC; ++l)
        {
            penaltyU.col(l) = bc(l) * problem->bcVelVec[l] - problem->bcVelMat[l] * aTmp;
        }
    }

    Eigen::MatrixXd cc(1, 1);
    for (int i = 0; i < Nphi_u; ++i)
    {
        cc = aTmp.transpose() * Eigen::SliceFromTensor(problem->C_tensor, 0, i) * aTmp
           - gNutFluct.transpose()
                 * Eigen::SliceFromTensor(problem->cTotalFluctTensor, 0, i) * aTmp
           - gNutAve.transpose()
                 * Eigen::SliceFromTensor(problem->cTotalAveTensor, 0, i) * aTmp;

        fvec(i) = -m5(i) + m1(i) - cc(0, 0) - m2(i);

        if (problem->bcMethod == "penalty")
        {
            fvec(i) += ((penaltyU * tauU)(i, 0));
        }
    }

    for (int j = 0; j < Nphi_p; ++j)
    {
        const int k = j + Nphi_u;
        fvec(k) = m3(j);
    }

    if (problem->bcMethod == "lift")
    {
        for (int j = 0; j < N_BC; ++j)
        {
            fvec(j) = x(j) - bc(j);
        }
    }

    return 0;
}

int newtonUnsteadyNSTurbSUPAve::df(const Eigen::VectorXd& x, Eigen::MatrixXd& fjac) const
{
    Eigen::NumericalDiff<newtonUnsteadyNSTurbSUPAve> numDiff(*this);
    numDiff.df(x, fjac);
    return 0;
}

// * * * * * * * * * * * * * * * Operators PPE * * * * * * * * * * * * * * * //

int newtonUnsteadyNSTurbPPE::operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const
{
    Eigen::VectorXd aTmp = x.head(Nphi_u);
    Eigen::VectorXd bTmp = x.tail(Nphi_p);

    Eigen::VectorXd a_dot(Nphi_u);
    if (problem->timeDerivativeSchemeOrder == "first")
    {
        a_dot = (x.head(Nphi_u) - y_old.head(Nphi_u)) / dt;
    }
    else
    {
        a_dot = (1.5 * x.head(Nphi_u) - 2.0 * y_old.head(Nphi_u) + 0.5 * yOldOld.head(Nphi_u)) / dt;
    }

    Eigen::MatrixXd cc(1, 1);
    Eigen::MatrixXd gg(1, 1);
    Eigen::MatrixXd bb(1, 1);

    Eigen::VectorXd m1 = problem->bTotalMatrix * aTmp * nu;
    Eigen::VectorXd m2 = problem->K_matrix * bTmp;
    Eigen::VectorXd m5 = problem->M_matrix * a_dot;

    Eigen::VectorXd m3 = problem->D_matrix * bTmp;
    Eigen::VectorXd m6 = problem->BC1_matrix * aTmp * nu;
    Eigen::VectorXd m7 = problem->BC3_matrix * aTmp * nu;

    Eigen::MatrixXd penaltyU = Eigen::MatrixXd::Zero(Nphi_u, N_BC);
    if (problem->bcMethod == "penalty")
    {
        for (int l = 0; l < N_BC; ++l)
        {
            penaltyU.col(l) = bc(l) * problem->bcVelVec[l] - problem->bcVelMat[l] * aTmp;
        }
    }

    for (int i = 0; i < Nphi_u; ++i)
    {
        cc = aTmp.transpose() * Eigen::SliceFromTensor(problem->C_tensor, 0, i) * aTmp
           - gNut.transpose() * Eigen::SliceFromTensor(problem->cTotalTensor, 0, i) * aTmp;

        fvec(i) = -m5(i) + m1(i) - cc(0, 0) - m2(i);

        if (problem->bcMethod == "penalty")
        {
            fvec(i) += ((penaltyU * tauU)(i, 0));
        }
    }

    for (int j = 0; j < Nphi_p; ++j)
    {
        const int k = j + Nphi_u;

        gg = aTmp.transpose() * Eigen::SliceFromTensor(problem->gTensor, 0, j) * aTmp;
        bb = aTmp.transpose() * Eigen::SliceFromTensor(problem->bc2Tensor, 0, j) * aTmp;

        // fvec(k) = m3(j, 0) - gg(0, 0) - m6(j, 0) + bb(0, 0);
        fvec(k) = m3(j, 0) + gg(0, 0) - m7(j, 0);
    }

    if (problem->bcMethod == "lift")
    {
        for (int j = 0; j < N_BC; ++j)
        {
            fvec(j) = x(j) - bc(j);
        }
    }

    return 0;
}

int newtonUnsteadyNSTurbPPE::df(const Eigen::VectorXd& x, Eigen::MatrixXd& fjac) const
{
    Eigen::NumericalDiff<newtonUnsteadyNSTurbPPE> numDiff(*this);
    numDiff.df(x, fjac);
    return 0;
}

int newtonUnsteadyNSTurbPPEAve::operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const
{
    Eigen::VectorXd aTmp = x.head(Nphi_u);
    Eigen::VectorXd bTmp = x.tail(Nphi_p);

    Eigen::VectorXd a_dot(Nphi_u);
    if (problem->timeDerivativeSchemeOrder == "first")
    {
        a_dot = (x.head(Nphi_u) - y_old.head(Nphi_u)) / dt;
    }
    else
    {
        a_dot = (1.5 * x.head(Nphi_u) - 2.0 * y_old.head(Nphi_u) + 0.5 * yOldOld.head(Nphi_u)) / dt;
    }

    Eigen::MatrixXd cc(1, 1);
    Eigen::MatrixXd gg(1, 1);
    Eigen::MatrixXd bb(1, 1);
    Eigen::MatrixXd nn(1, 1);

    Eigen::VectorXd m1 = problem->bTotalMatrix * aTmp * nu;  // ν B a
    Eigen::VectorXd m2 = problem->K_matrix * bTmp;           // K b   (grad p)
    Eigen::VectorXd m5 = problem->M_matrix * a_dot;          // M ȧ

    Eigen::VectorXd m3 = problem->D_matrix * bTmp;           // D b   (Δp term)
    Eigen::VectorXd m6 = problem->BC1_matrix * aTmp * nu;    // (not used here)
    Eigen::VectorXd m7 = problem->BC3_matrix * aTmp * nu;    // ν N a

    Eigen::MatrixXd penaltyU = Eigen::MatrixXd::Zero(Nphi_u, N_BC);
    if (problem->bcMethod == "penalty")
    {
        for (int l = 0; l < N_BC; ++l)
        {
            penaltyU.col(l) = bc(l) * problem->bcVelVec[l] - problem->bcVelMat[l] * aTmp;
        }
    }

    for (int i = 0; i < Nphi_u; ++i)
    {
        cc = aTmp.transpose() * Eigen::SliceFromTensor(problem->C_tensor, 0, i) * aTmp
           - gNutAve.transpose()
                 * Eigen::SliceFromTensor(problem->cTotalAveTensor, 0, i) * aTmp
           - gNutFluct.transpose()
                 * Eigen::SliceFromTensor(problem->cTotalFluctTensor, 0, i) * aTmp;

        fvec(i) = -m5(i) + m1(i) - cc(0, 0) - m2(i);

        if (problem->bcMethod == "penalty")
        {
            fvec(i) += ((penaltyU * tauU)(i, 0));
        }
    }

    // Thesis: D b + aᵀ G a - ν N a - L = 0
    for (int j = 0; j < Nphi_p; ++j)
    {
        const int k = j + Nphi_u;

        gg = aTmp.transpose() * Eigen::SliceFromTensor(problem->gTensor, 0, j) * aTmp;
        bb = aTmp.transpose() * Eigen::SliceFromTensor(problem->bc2Tensor, 0, j) * aTmp;

        nn = gNutAve.transpose()
                 * Eigen::SliceFromTensor(problem->cTotalPPEAveTensor, 0, j) * aTmp
           + gNutFluct.transpose()
                 * Eigen::SliceFromTensor(problem->cTotalPPEFluctTensor, 0, j) * aTmp;

        double Lj = 0.0;
        if (problem->L_vector.rows() == Nphi_p && problem->L_vector.cols() >= 1)
        {
            Lj = problem->L_vector(j, 0);
        }

        fvec(k) = m3(j, 0) + gg(0, 0) - m7(j, 0) - nn(0, 0) - Lj;
    }

    if (problem->bcMethod == "lift")
    {
        for (int j = 0; j < N_BC; ++j)
        {
            fvec(j) = x(j) - bc(j);
        }
    }

    return 0;
}

int newtonUnsteadyNSTurbPPEAve::df(const Eigen::VectorXd& x, Eigen::MatrixXd& fjac) const
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

    const int numberOfStores = round(storeEvery / dt);
    std::cout << "[DEBUG] storeEvery = " << storeEvery << ", dt = " << dt
              << ", numberOfStores = " << numberOfStores << std::endl;

    std::cout << "[DEBUG] vel shape: (" << vel.rows() << ", " << vel.cols() << ")\n";
    std::cout << "[DEBUG] vel input:\n" << vel << std::endl;

    vel_now = vel;
    std::cout << "[DEBUG] vel_now set directly:\n" << vel_now << std::endl;

    if (problem->bcMethod == "lift")
    {
        std::cout << "[DEBUG] Using lifting function for BCs\n";
    }
    else if (problem->bcMethod == "penalty")
    {
        std::cout << "[DEBUG] Using penalty method\n";
    }

    y.resize(Nphi_u + Nphi_p, 1);
    y.setZero();
    y.head(Nphi_u) = ITHACAutilities::getCoeffs(problem->Ufield[0], Umodes);
    y.tail(Nphi_p) = ITHACAutilities::getCoeffs(problem->Pfield[0], Pmodes);
    nut0           = ITHACAutilities::getCoeffs(problem->nutFields[0], nutModes);

    std::cout << "[DEBUG] Nphi_u = " << Nphi_u << ", Nphi_p = " << Nphi_p << std::endl;
    std::cout << "[DEBUG] y.size() = " << y.rows() << " x " << y.cols() << std::endl;
    std::cout << "[DEBUG] Initial y.head(Nphi_u): " << y.head(Nphi_u).transpose() << std::endl;

    if (problem->bcMethod == "lift")
    {
        for (int j = 0; j < N_BC; ++j)
        {
            y(j) = vel_now(j, 0);
        }
        std::cout << "[DEBUG] Updated y(0:N_BC) with vel_now\n";
    }

    const int firstRBFInd = (skipLift && problem->bcMethod == "lift") ? N_BC : 0;

    newtonObjectSUP.nu      = nu;
    newtonObjectSUP.y_old   = y;
    newtonObjectSUP.yOldOld = y;
    newtonObjectSUP.dt      = dt;
    newtonObjectSUP.bc.resize(N_BC);
    newtonObjectSUP.tauU = tauU;
    newtonObjectSUP.gNut = nut0;

    for (int j = 0; j < N_BC; ++j)
    {
        newtonObjectSUP.bc(j) = vel_now(j, 0);
    }

    const int Ntsteps    = static_cast<int>((finalTime - tstart) / dt);
    const int onlineSize = static_cast<int>(Ntsteps / numberOfStores);

    std::cout << "[DEBUG] Ntsteps = " << Ntsteps << ", onlineSize = " << onlineSize << std::endl;

    online_solution.resize(onlineSize);
    rbfCoeffMat.resize(nphiNut + 1, onlineSize + 3);

    std::cout << "[DEBUG] rbfCoeffMat shape = " << rbfCoeffMat.rows() << " x " << rbfCoeffMat.cols() << "\n";

    time = tstart;
    int counter = 0;
    int counter2 = 0;
    int nextStore = 0;

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

    Color::Modifier red(Color::FG_RED);
    Color::Modifier green(Color::FG_GREEN);
    Color::Modifier def(Color::FG_DEFAULT);

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
                break;
            }

            default:
                FatalErrorInFunction << "Invalid interChoice: " << interChoice << abort(FatalError);
        }

        for (int i = 0; i < nphiNut; ++i)
        {
            newtonObjectSUP.gNut(i) = problem->rbfSplines[i]->eval(tv);
        }

        if (problem->bcMethod == "lift")
        {
            for (int j = 0; j < N_BC; ++j)
            {
                y(j) = vel_now(j, 0);
            }
        }

        newtonObjectSUP.operator()(y, res);

        newtonObjectSUP.yOldOld = newtonObjectSUP.y_old;
        newtonObjectSUP.y_old   = y;

        std::cout << "################## Online solve N° " << count_online_solve
                  << " ##################\n";
        std::cout << "Time = " << time << ", Parameter: " << vel_now(0, 0) << std::endl;

        if (res.norm() < 1e-5)
        {
            std::cout << green << "|F(x)| = " << res.norm()
                      << " - Converged in " << hnls.iter << " iterations" << def << std::endl;
        }
        else
        {
            std::cout << red << "|F(x)| = " << res.norm()
                      << " - Not converged in " << hnls.iter << " iterations" << def << std::endl;
        }

        ++count_online_solve;

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
    using std::cout;
    using std::endl;

    cout << "\n\n[DEBUG] === Entering solveOnlineSUPAve() (AVG linear-μ + FLUCT RBF) ===\n";

    std::ofstream avgLog("./ITHACAoutput/Coefficients/avg_coeffs_online_log_SUP.txt");
    std::ofstream flLog("./ITHACAoutput/Coefficients/fluct_coeffs_online_log_SUP.txt");
    avgLog << std::setprecision(15);
    flLog << std::setprecision(15);
    avgLog << "# step time mu*  gNutAvg[0..nphiNutAvg-1]\n";
    flLog << "# step time      gNutFluct[0..nphiNutFluct-1]\n";

    cout << "[DEBUG] dt=" << dt << ", storeEvery=" << storeEvery << ", exportEvery=" << exportEvery << endl;
    M_Assert(exportEvery >= dt && storeEvery >= dt, "dt must be <= exportEvery/storeEvery");
    M_Assert(ITHACAutilities::isInteger(storeEvery / dt), "storeEvery must be integer multiple of dt.");
    M_Assert(ITHACAutilities::isInteger(exportEvery / dt), "exportEvery must be integer multiple of dt.");
    M_Assert(
        ITHACAutilities::isInteger(exportEvery / storeEvery),
        "exportEvery must be integer multiple of storeEvery.");

    const int numberOfStores = round(storeEvery / dt);
    cout << "[DEBUG] numberOfStores = " << numberOfStores << endl;

    if (problem->bcMethod == "lift")
    {
        vel_now = setOnlineVelocity(vel);
        cout << "[DEBUG] vel_now set via setOnlineVelocity(vel), shape = (" << vel_now.rows() << ","
             << vel_now.cols() << ")\n";
    }
    else
    {
        vel_now = vel;
        cout << "[DEBUG] vel_now = vel, shape = (" << vel_now.rows() << "," << vel_now.cols() << ")\n";
    }

    cout << "[DEBUG] Projecting Ufield[0] and Pfield[0] onto modes (t=0)\n";
    y.resize(Nphi_u + Nphi_p, 1);
    y.setZero();

    y.head(Nphi_u) = ITHACAutilities::getCoeffs(problem->Uomfield[0], Umodes);
    cout << "[DEBUG] y.head(Nphi_u) = " << y.head(Nphi_u).transpose() << endl;

    y.tail(Nphi_p) = ITHACAutilities::getCoeffs(problem->Pfield[0], Pmodes);
    cout << "[DEBUG] y.tail(Nphi_p) = " << y.tail(Nphi_p).transpose() << endl;

    const bool useLift = (problem->bcMethod == "lift");
    if (useLift)
    {
        for (int j = 0; j < N_BC; ++j)
        {
            y(j) = vel_now(j, 0);
            cout << "[DEBUG] Applied lift BC: y(" << j << ") = " << y(j) << endl;
        }
    }

    const int firstRBFInd = useLift ? N_BC : 0;
    const int Nphys       = problem->NUmodes;  // only physical velocity modes

    cout << "[DEBUG] firstRBFInd = " << firstRBFInd
         << ", Nphys = " << Nphys
         << ", Nphi_u (lift+phys+SUP) = " << Nphi_u << "\n";

    M_Assert(firstRBFInd + Nphys <= Nphi_u, "RBF input slice exceeds velocity DOF range.");

    Eigen::VectorXd muStar(1);
    muStar(0) = vel(0, 0);
    const double muStarVal = muStar(0);
    cout << "[DEBUG] muStar = " << muStar.transpose() << endl;

    static bool avgTblLoaded = false;
    static Eigen::VectorXd muAvg;
    static Eigen::MatrixXd Cavg;

    if (!avgTblLoaded)
    {
        const word coeffDir = "./ITHACAoutput/Coefficients/";
        const Eigen::MatrixXd muAvgMat = ITHACAstream::readMatrix(coeffDir + "NutAvg_mu_unique_mat.txt");
        Cavg = ITHACAstream::readMatrix(coeffDir + "NutAvg_coeffs_by_mu_mat.txt");

        if (muAvgMat.rows() == 1 && muAvgMat.cols() > 1)
        {
            muAvg = muAvgMat.transpose();
        }
        else
        {
            muAvg = muAvgMat.col(0);
        }

        M_Assert(Cavg.rows() == nphiNutAvg, "NutAvg coeff table row mismatch vs nphiNutAvg.");
        M_Assert(Cavg.cols() == muAvg.size(), "NutAvg table columns must match μ-grid length.");

        avgTblLoaded = true;
        cout << "[DEBUG] Loaded NutAvg tables: M=" << muAvg.size()
             << "  Cavg=" << Cavg.rows() << "x" << Cavg.cols() << endl;
    }

    auto interpNutAvgLinear = [&](double mu) -> Eigen::VectorXd
    {
        const int M = static_cast<int>(muAvg.size());
        if (mu <= muAvg(0))   return Cavg.col(0);
        if (mu >= muAvg(M - 1)) return Cavg.col(M - 1);

        int j = 1;
        while (j < M && muAvg(j) < mu) ++j;

        const double muL = muAvg(j - 1);
        const double muR = muAvg(j);
        const double t = (mu - muL) / (muR - muL);

        return (1.0 - t) * Cavg.col(j - 1) + t * Cavg.col(j);
    };

    Eigen::VectorXd aOld    = y.segment(firstRBFInd, Nphys);
    Eigen::VectorXd aNow    = aOld;
    Eigen::VectorXd aDotNow = Eigen::VectorXd::Zero(Nphys);

    Eigen::VectorXd rbfInputFluct(2 * Nphys);
    rbfInputFluct.head(Nphys) = aNow;
    rbfInputFluct.tail(Nphys) = aDotNow;

    auto printHead = [](const Eigen::VectorXd& v, int k)
    {
        std::stringstream ss;
        ss << std::setprecision(6);
        k = std::min<int>(k, v.size());
        for (int i = 0; i < k; ++i)
        {
            ss << v(i);
            if (i + 1 < k) ss << ", ";
        }
        return ss.str();
    };

    auto vecMin = [](const Eigen::VectorXd& v) { return v.minCoeff(); };
    auto vecMax = [](const Eigen::VectorXd& v) { return v.maxCoeff(); };

    Eigen::VectorXd coeffNutAvgLocal(nphiNutAvg);
    Eigen::VectorXd coeffNutFluctLocal(nphiNutFluct);

    coeffNutAvgLocal = interpNutAvgLinear(muStarVal);
    coeffNutFluctLocal.setZero();

    cout << "[DEBUG] t=0: ||gAvg||2=" << coeffNutAvgLocal.norm()
         << "  min/max=(" << vecMin(coeffNutAvgLocal) << "," << vecMax(coeffNutAvgLocal)
         << ")  head=[" << printHead(coeffNutAvgLocal, 5) << "]\n";
    cout << "[DEBUG] t=0: ||gFluct||2=" << coeffNutFluctLocal.norm()
         << "  min/max=(" << vecMin(coeffNutFluctLocal) << "," << vecMax(coeffNutFluctLocal)
         << ")  head=[" << printHead(coeffNutFluctLocal, 5) << "]\n";

    avgLog << 0 << " " << tstart << " " << muStarVal << " ";
    for (int k = 0; k < coeffNutAvgLocal.size(); ++k)
    {
        avgLog << coeffNutAvgLocal(k) << (k + 1 < coeffNutAvgLocal.size() ? ' ' : '\n');
    }

    flLog << 0 << " " << tstart << " ";
    for (int k = 0; k < coeffNutFluctLocal.size(); ++k)
    {
        flLog << coeffNutFluctLocal(k) << (k + 1 < coeffNutFluctLocal.size() ? ' ' : '\n');
    }

    newtonObjectSUPAve.nu        = nu;
    newtonObjectSUPAve.y_old     = y;
    newtonObjectSUPAve.yOldOld   = y;
    newtonObjectSUPAve.dt        = dt;
    newtonObjectSUPAve.bc.resize(N_BC);
    newtonObjectSUPAve.tauU      = tauU;
    newtonObjectSUPAve.gNutAve   = coeffNutAvgLocal;
    newtonObjectSUPAve.gNutFluct = coeffNutFluctLocal;

    for (int j = 0; j < N_BC; ++j)
    {
        newtonObjectSUPAve.bc(j) = vel_now(j, 0);
    }

    const int Ntsteps    = static_cast<int>((finalTime - tstart) / dt);
    const int onlineSize = static_cast<int>(Ntsteps / numberOfStores) + 3;

    online_solution.resize(onlineSize);
    rbfCoeffMat.setZero(1 + nphiNutAvg + nphiNutFluct, onlineSize);

    double timeLocal = tstart;
    int step = 0;
    int nextStore = 0;
    int counter2 = 0;
    int counter = 0;

    Eigen::VectorXd tmp_sol(Nphi_u + Nphi_p + 1);
    tmp_sol(0) = timeLocal;
    tmp_sol.tail(y.rows()) = y;

    online_solution[counter] = tmp_sol;
    rbfCoeffMat(0, counter2) = timeLocal;
    rbfCoeffMat.block(1, counter2, nphiNutAvg, 1) = coeffNutAvgLocal;
    rbfCoeffMat.block(1 + nphiNutAvg, counter2, nphiNutFluct, 1) = coeffNutFluctLocal;

    cout << "[DEBUG] t=0: Stored solution and RBF coeffs at counter2=" << counter2 << endl;

    ++counter2;
    ++counter;
    nextStore += numberOfStores;

    Eigen::HybridNonLinearSolver<newtonUnsteadyNSTurbSUPAve> hnls(newtonObjectSUPAve);

    cout << "\n[DEBUG] == ENTERING MAIN TIME LOOP t>0 ==\n";
    while (timeLocal < finalTime)
    {
        timeLocal += dt;
        ++step;

        cout << "\n---------------------------------------------\n";
        cout << "[DEBUG] === Time step: t=" << timeLocal << ", step=" << step << " ===\n";

        problem->timeDerivativeSchemeOrder = "first";
        cout << "[DEBUG] y (before Newton) = " << y.transpose() << endl;

        Eigen::VectorXd res(y);
        res.setZero();

        hnls.solve(y);
        cout << "[DEBUG] y (after  Newton) = " << y.transpose() << endl;

        newtonObjectSUPAve.operator()(y, res);

        const Eigen::VectorXd aNowFull = y.head(Nphi_u);
        const Eigen::VectorXd aOldFull = newtonObjectSUPAve.y_old.head(Nphi_u);

        const Eigen::VectorXd aNowPhys = aNowFull.segment(firstRBFInd, Nphys);
        const Eigen::VectorXd aDotPhys =
            (aNowFull.segment(firstRBFInd, Nphys) - aOldFull.segment(firstRBFInd, Nphys)) / dt;

        rbfInputFluct.head(Nphys) = aNowPhys;
        rbfInputFluct.tail(Nphys) = aDotPhys;

        cout << "[DEBUG] aNow_phys = " << aNowPhys.transpose() << endl;
        cout << "[DEBUG] aDot_phys = " << aDotPhys.transpose() << endl;

        coeffNutAvgLocal = interpNutAvgLinear(muStarVal);

        for (int i = 0; i < nphiNutFluct; ++i)
        {
            coeffNutFluctLocal(i) = rbfSplinesNutFluct[i]->eval(rbfInputFluct);
        }

        newtonObjectSUPAve.gNutAve   = coeffNutAvgLocal;
        newtonObjectSUPAve.gNutFluct = coeffNutFluctLocal;

        cout << std::setprecision(6)
             << "[DEBUG] coeffNutAvg:  ||·||2=" << coeffNutAvgLocal.norm()
             << "  min/max=(" << vecMin(coeffNutAvgLocal) << "," << vecMax(coeffNutAvgLocal)
             << ")  head=[" << printHead(coeffNutAvgLocal, 5) << "]\n";
        cout << "[DEBUG] coeffNutFluct:||·||2=" << coeffNutFluctLocal.norm()
             << "  min/max=(" << vecMin(coeffNutFluctLocal) << "," << vecMax(coeffNutFluctLocal)
             << ")  head=[" << printHead(coeffNutFluctLocal, 5) << "]\n";

        avgLog << step << " " << timeLocal << " " << muStarVal << " ";
        for (int k = 0; k < coeffNutAvgLocal.size(); ++k)
        {
            avgLog << coeffNutAvgLocal(k) << (k + 1 < coeffNutAvgLocal.size() ? ' ' : '\n');
        }

        flLog << step << " " << timeLocal << " ";
        for (int k = 0; k < coeffNutFluctLocal.size(); ++k)
        {
            flLog << coeffNutFluctLocal(k) << (k + 1 < coeffNutFluctLocal.size() ? ' ' : '\n');
        }

        if (useLift)
        {
            for (int j = 0; j < N_BC; ++j)
            {
                y(j) = vel_now(j, 0);
            }
        }

        newtonObjectSUPAve.yOldOld = newtonObjectSUPAve.y_old;
        newtonObjectSUPAve.y_old   = y;

        tmp_sol(0) = timeLocal;
        tmp_sol.tail(y.rows()) = y;

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

            if (counter2 < rbfCoeffMat.cols())
            {
                rbfCoeffMat(0, counter2) = timeLocal;
                rbfCoeffMat.block(1, counter2, nphiNutAvg, 1) = coeffNutAvgLocal;
                rbfCoeffMat.block(1 + nphiNutAvg, counter2, nphiNutFluct, 1) = coeffNutFluctLocal;
            }

            cout << "[DEBUG] Stored solution and RBF coeffs at counter2=" << counter2 << endl;

            nextStore += numberOfStores;
            ++counter2;
        }

        ++counter;
    }

    const int storedCols = counter2;

    cout << "[DEBUG] Trimming outputs: storedCols=" << storedCols
         << ", rbfCoeffMat orig cols=" << rbfCoeffMat.cols()
         << ", online_solution orig size=" << online_solution.size() << endl;

    if (storedCols > 0 && storedCols <= rbfCoeffMat.cols())
    {
        rbfCoeffMat.conservativeResize(rbfCoeffMat.rows(), storedCols);
    }
    if (storedCols >= 0 && storedCols <= online_solution.size())
    {
        online_solution.resize(storedCols);
    }

    cout << "\n[DEBUG] Exporting online_solution and rbfCoeffMat...\n";
    ITHACAstream::exportMatrix(online_solution, "red_coeff_SUP", "python", "./ITHACAoutput/red_coeff");
    ITHACAstream::exportMatrix(online_solution, "red_coeff_SUP", "matlab", "./ITHACAoutput/red_coeff");
    ITHACAstream::exportMatrix(rbfCoeffMat, "rbf_coeffs_online_SUP", "python", "./ITHACAoutput/red_coeff");
    ITHACAstream::exportMatrix(rbfCoeffMat, "rbf_coeffs_online_SUP", "matlab", "./ITHACAoutput/red_coeff");

    avgLog.close();
    flLog.close();

    cout << "[DEBUG] solveOnlineSUPAve() completed.\n\n";
}

void ReducedUnsteadyNSTurb::solveOnlinePPE(Eigen::MatrixXd vel)
{
    M_Assert(exportEvery >= dt, "The time step dt must be smaller than exportEvery.");
    M_Assert(storeEvery >= dt, "The time step dt must be smaller than storeEvery.");
    M_Assert(
        ITHACAutilities::isInteger(storeEvery / dt) == true,
        "The variable storeEvery must be an integer multiple of the time step dt.");
    M_Assert(
        ITHACAutilities::isInteger(exportEvery / dt) == true,
        "The variable exportEvery must be an integer multiple of the time step dt.");
    M_Assert(
        ITHACAutilities::isInteger(exportEvery / storeEvery) == true,
        "The variable exportEvery must be an integer multiple of the variable storeEvery.");

    const int numberOfStores = round(storeEvery / dt);

    if (problem->bcMethod == "lift")
    {
        vel_now = setOnlineVelocity(vel);
    }
    else if (problem->bcMethod == "penalty")
    {
        vel_now = vel;
    }

    y.resize(Nphi_u + Nphi_p, 1);
    y.setZero();

    y.head(Nphi_u) = initCond.col(0).head(Nphi_u);
    y.tail(Nphi_p) = initCond.col(0).tail(Nphi_p);

    int nextStore = 0;
    int counter2 = 0;

    if (problem->bcMethod == "lift")
    {
        for (int j = 0; j < N_BC; ++j)
        {
            y(j) = vel_now(j, 0);
        }
    }

    const int firstRBFInd = (skipLift == true && problem->bcMethod == "lift") ? N_BC : 0;

    newtonObjectPPE.nu      = nu;
    newtonObjectPPE.y_old   = y;
    newtonObjectPPE.yOldOld = newtonObjectPPE.y_old;
    newtonObjectPPE.dt      = dt;
    newtonObjectPPE.bc.resize(N_BC);
    newtonObjectPPE.tauU = tauU;
    newtonObjectPPE.gNut = nut0;

    for (int j = 0; j < N_BC; ++j)
    {
        newtonObjectPPE.bc(j) = vel_now(j, 0);
    }

    const int Ntsteps    = static_cast<int>((finalTime - tstart) / dt);
    const int onlineSize = static_cast<int>(Ntsteps / numberOfStores);

    online_solution.resize(onlineSize);
    rbfCoeffMat.resize(nphiNut + 1, onlineSize + 3);

    time = tstart;

    int counter = 0;

    Eigen::MatrixXd tmp_sol(Nphi_u + Nphi_p + 1, 1);
    tmp_sol(0) = time;
    tmp_sol.col(0).tail(y.rows()) = y;

    if ((time != 0) || (startFromZero == true))
    {
        online_solution[counter] = tmp_sol;
        ++counter;

        rbfCoeffMat(0, counter2) = time;
        rbfCoeffMat.block(1, counter2, nphiNut, 1) = nut0;

        ++counter2;
        nextStore += numberOfStores;
    }

    Eigen::HybridNonLinearSolver<newtonUnsteadyNSTurbPPE> hnls(newtonObjectPPE);

    Color::Modifier red(Color::FG_RED);
    Color::Modifier green(Color::FG_GREEN);
    Color::Modifier def(Color::FG_DEFAULT);

    while (time < finalTime)
    {
        time += dt;

        Eigen::VectorXd res(y);
        res.setZero();

        hnls.solve(y);
        newtonObjectPPE.operator()(y, res);

        Eigen::VectorXd tv(dimA);
        const Eigen::VectorXd aDer = (y.head(Nphi_u) - newtonObjectPPE.y_old.head(Nphi_u)) / dt;

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
                tv << muStar,
                    y.segment(firstRBFInd, (dimA - muStar.size()) / 2),
                    aDer.segment(firstRBFInd, (dimA - muStar.size()) / 2);
                break;

            default:
                tv << y.segment(firstRBFInd, dimA);
                break;
        }

        for (int i = 0; i < nphiNut; ++i)
        {
            newtonObjectPPE.gNut(i) = problem->rbfSplines[i]->eval(tv);
        }

        newtonObjectPPE.operator()(y, res);

        newtonObjectPPE.yOldOld = newtonObjectPPE.y_old;
        newtonObjectPPE.y_old   = y;

        std::cout << "################## Online solve N° " << count_online_solve
                  << " ##################" << std::endl;

        Info << "Time = " << time << endl;
        std::cout << "Solving for the parameter: " << vel_now << std::endl;

        if (res.norm() < 1e-5)
        {
            std::cout << green << "|F(x)| = " << res.norm() << " - Minimun reached in " << hnls.iter
                      << " iterations " << def << std::endl << std::endl;
        }
        else
        {
            std::cout << red << "|F(x)| = " << res.norm() << " - Minimun reached in " << hnls.iter
                      << " iterations " << def << std::endl << std::endl;
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
            ++counter2;
        }

        ++counter;
    }

    ITHACAstream::exportMatrix(online_solution, "red_coeff", "python", "./ITHACAoutput/red_coeff");
    ITHACAstream::exportMatrix(online_solution, "red_coeff", "matlab", "./ITHACAoutput/red_coeff");

    count_online_solve += 1;
}

void ReducedUnsteadyNSTurb::solveOnlinePPEAve(Eigen::MatrixXd vel)
{
    using std::cout;
    using std::endl;

    cout << "\n\n[DEBUG] === Entering solveOnlinePPEAve() (AVG linear-μ + FLUCT RBF) ===\n";

    std::ofstream avgLog("./ITHACAoutput/Coefficients/avg_coeffs_online_log.txt");
    std::ofstream flLog("./ITHACAoutput/Coefficients/fluct_coeffs_online_log.txt");
    avgLog << std::setprecision(15);
    flLog << std::setprecision(15);
    avgLog << "# step time mu*  gNutAvg[0..nphiNutAvg-1]\n";
    flLog << "# step time      gNutFluct[0..nphiNutFluct-1]\n";

    cout << "[DEBUG] dt=" << dt << ", storeEvery=" << storeEvery << ", exportEvery=" << exportEvery << endl;
    M_Assert(exportEvery >= dt && storeEvery >= dt, "dt must be <= exportEvery/storeEvery");
    M_Assert(ITHACAutilities::isInteger(storeEvery / dt), "storeEvery must be integer multiple of dt.");
    M_Assert(ITHACAutilities::isInteger(exportEvery / dt), "exportEvery must be integer multiple of dt.");
    M_Assert(
        ITHACAutilities::isInteger(exportEvery / storeEvery),
        "exportEvery must be integer multiple of storeEvery.");

    const int numberOfStores = round(storeEvery / dt);
    cout << "[DEBUG] numberOfStores = " << numberOfStores << endl;

    if (problem->bcMethod == "lift")
    {
        vel_now = setOnlineVelocity(vel);
        cout << "[DEBUG] vel_now set via setOnlineVelocity(vel), shape = (" << vel_now.rows() << ","
             << vel_now.cols() << ")\n";
    }
    else
    {
        vel_now = vel;
        cout << "[DEBUG] vel_now = vel, shape = (" << vel_now.rows() << "," << vel_now.cols() << ")\n";
    }

    cout << "[DEBUG] Projecting Ufield[0] and Pfield[0] onto modes (t=0)\n";
    y.resize(Nphi_u + Nphi_p, 1);
    y.setZero();

    y.head(Nphi_u) = ITHACAutilities::getCoeffs(problem->Uomfield[0], Umodes);
    cout << "[DEBUG] y.head(Nphi_u) = " << y.head(Nphi_u).transpose() << endl;

    y.tail(Nphi_p) = ITHACAutilities::getCoeffs(problem->Pfield[0], Pmodes);
    cout << "[DEBUG] y.tail(Nphi_p) = " << y.tail(Nphi_p).transpose() << endl;
    cout << "[DEBUG] y (full) = " << y.transpose() << endl;

    const bool useLift = (problem->bcMethod == "lift");
    if (useLift)
    {
        for (int j = 0; j < N_BC; ++j)
        {
            y(j) = vel_now(j, 0);
            cout << "[DEBUG] Applied lift BC: y(" << j << ") = " << y(j) << endl;
        }
    }

    const int firstRBFInd = useLift ? N_BC : 0;
    const int Nphys = problem->NUmodes;

    cout << "[DEBUG] firstRBFInd = " << firstRBFInd
         << ", Nphys = " << Nphys
         << ", Nphi_u (lift+phys+SUP) = " << Nphi_u << "\n";

    M_Assert(
        firstRBFInd + Nphys <= Nphi_u,
        "RBF input slice exceeds velocity DOF range. Check N_BC/NUmodes/Nphi_u.");

    Eigen::VectorXd muStar(1);
    muStar(0) = vel(0, 0);
    const double muStarVal = muStar(0);
    cout << "[DEBUG] muStar = " << muStar.transpose() << endl;

    static bool avgTblLoaded = false;
    static Eigen::VectorXd muAvg;
    static Eigen::MatrixXd Cavg;

    if (!avgTblLoaded)
    {
        const word coeffDir = "./ITHACAoutput/Coefficients/";
        const Eigen::MatrixXd muAvgMat = ITHACAstream::readMatrix(coeffDir + "NutAvg_mu_unique_mat.txt");
        Cavg = ITHACAstream::readMatrix(coeffDir + "NutAvg_coeffs_by_mu_mat.txt");

        if (muAvgMat.rows() == 1 && muAvgMat.cols() > 1)
        {
            muAvg = muAvgMat.transpose();
        }
        else
        {
            muAvg = muAvgMat.col(0);
        }

        M_Assert(Cavg.rows() == nphiNutAvg, "NutAvg coeff table row mismatch vs nphiNutAvg.");
        M_Assert(Cavg.cols() == muAvg.size(), "NutAvg table columns must match mu grid length.");

        avgTblLoaded = true;
        cout << "[DEBUG] Loaded NutAvg tables: M=" << muAvg.size()
             << "  Cavg=" << Cavg.rows() << "x" << Cavg.cols() << endl;
    }

    auto interpNutAvgLinear = [&](double mu) -> Eigen::VectorXd
    {
        const int M = static_cast<int>(muAvg.size());
        if (mu <= muAvg(0))     return Cavg.col(0);
        if (mu >= muAvg(M - 1)) return Cavg.col(M - 1);

        int j = 1;
        while (j < M && muAvg(j) < mu) ++j;

        const double muL = muAvg(j - 1);
        const double muR = muAvg(j);
        const double t = (mu - muL) / (muR - muL);

        cout << "[DEBUG] AVG μ-bracket: jL=" << (j - 1) << " jR=" << j
             << "  muL=" << muL << " muR=" << muR << "  t=" << t << endl;

        return (1.0 - t) * Cavg.col(j - 1) + t * Cavg.col(j);
    };

    Eigen::VectorXd aOld = y.segment(firstRBFInd, Nphys);
    Eigen::VectorXd aNow = aOld;
    Eigen::VectorXd aDotNow = Eigen::VectorXd::Zero(Nphys);

    Eigen::VectorXd rbfInputFluct(2 * Nphys);
    rbfInputFluct.head(Nphys) = aNow;
    rbfInputFluct.tail(Nphys) = aDotNow;

    cout << "[DEBUG] t=0: aOld_phys = " << aOld.transpose() << endl;
    cout << "[DEBUG] t=0: aDot_phys = " << aDotNow.transpose() << endl;

    Eigen::VectorXd coeffNutAvgLocal(nphiNutAvg);
    Eigen::VectorXd coeffNutFluctLocal(nphiNutFluct);

    coeffNutAvgLocal = interpNutAvgLinear(muStarVal);
    coeffNutFluctLocal.setZero();

    auto printHead = [](const Eigen::VectorXd& v, int k)
    {
        std::stringstream ss;
        ss << std::setprecision(6);
        k = std::min<int>(k, v.size());
        for (int i = 0; i < k; ++i)
        {
            ss << v(i);
            if (i + 1 < k) ss << ", ";
        }
        return ss.str();
    };

    auto vecMin = [](const Eigen::VectorXd& v) { return v.minCoeff(); };
    auto vecMax = [](const Eigen::VectorXd& v) { return v.maxCoeff(); };

    cout << "[DEBUG] t=0: ||gAvg||2=" << coeffNutAvgLocal.norm()
         << "  min/max=(" << vecMin(coeffNutAvgLocal) << "," << vecMax(coeffNutAvgLocal)
         << ")  head=[" << printHead(coeffNutAvgLocal, 5) << "]\n";
    cout << "[DEBUG] t=0: ||gFluct||2=" << coeffNutFluctLocal.norm()
         << "  min/max=(" << vecMin(coeffNutFluctLocal) << "," << vecMax(coeffNutFluctLocal)
         << ")  head=[" << printHead(coeffNutFluctLocal, 5) << "]\n";

    avgLog << 0 << " " << tstart << " " << muStarVal << " ";
    for (int k = 0; k < coeffNutAvgLocal.size(); ++k)
    {
        avgLog << coeffNutAvgLocal(k) << (k + 1 < coeffNutAvgLocal.size() ? ' ' : '\n');
    }

    flLog << 0 << " " << tstart << " ";
    for (int k = 0; k < coeffNutFluctLocal.size(); ++k)
    {
        flLog << coeffNutFluctLocal(k) << (k + 1 < coeffNutFluctLocal.size() ? ' ' : '\n');
    }

    newtonObjectPPEAve.gNutAve   = coeffNutAvgLocal;
    newtonObjectPPEAve.gNutFluct = coeffNutFluctLocal;

    newtonObjectPPEAve.nu      = nu;
    newtonObjectPPEAve.y_old   = y;
    newtonObjectPPEAve.yOldOld = y;
    newtonObjectPPEAve.dt      = dt;
    newtonObjectPPEAve.bc.resize(N_BC);
    newtonObjectPPEAve.tauU = tauU;
    newtonObjectPPEAve.gNutAve.resize(nphiNutAvg);
    newtonObjectPPEAve.gNutFluct.resize(nphiNutFluct);

    for (int j = 0; j < N_BC; ++j)
    {
        newtonObjectPPEAve.bc(j) = vel_now(j, 0);
    }

    const int Ntsteps = static_cast<int>((finalTime - tstart) / dt);
    const int onlineSize = static_cast<int>(Ntsteps / numberOfStores) + 3;

    online_solution.resize(onlineSize);
    rbfCoeffMat.setZero(1 + nphiNutAvg + nphiNutFluct, onlineSize);

    double timeLocal = tstart;
    int step = 0;
    int nextStore = 0;
    int counter2 = 0;
    int counter = 0;

    Eigen::VectorXd tmp_sol(Nphi_u + Nphi_p + 1);
    tmp_sol(0) = timeLocal;
    tmp_sol.tail(y.rows()) = y;

    online_solution[counter] = tmp_sol;
    rbfCoeffMat(0, counter2) = timeLocal;
    rbfCoeffMat.block(1, counter2, nphiNutAvg, 1) = coeffNutAvgLocal;
    rbfCoeffMat.block(1 + nphiNutAvg, counter2, nphiNutFluct, 1) = coeffNutFluctLocal;

    cout << "[DEBUG] t=0: Stored solution and RBF coeffs at counter2=" << counter2 << endl;

    ++counter2;
    ++counter;
    nextStore += numberOfStores;

    Eigen::HybridNonLinearSolver<newtonUnsteadyNSTurbPPEAve> hnls(newtonObjectPPEAve);

    cout << "\n[DEBUG] == ENTERING MAIN TIME LOOP t>0 ==\n";
    while (timeLocal < finalTime)
    {
        timeLocal += dt;
        ++step;

        cout << "\n---------------------------------------------\n";
        cout << "[DEBUG] === Time step: t=" << timeLocal << ", step=" << step << " ===\n";

        problem->timeDerivativeSchemeOrder = "first";
        cout << "[DEBUG] y (before Newton) = " << y.transpose() << endl;

        Eigen::VectorXd res(y);
        res.setZero();

        hnls.solve(y);
        cout << "[DEBUG] y (after  Newton) = " << y.transpose() << endl;

        newtonObjectPPEAve.operator()(y, res);

        aNow = y.segment(firstRBFInd, Nphys);
        aDotNow = (aNow - aOld) / dt;

        rbfInputFluct.head(Nphys) = aNow;
        rbfInputFluct.tail(Nphys) = aDotNow;

        cout << "[DEBUG] aNow_phys = " << aNow.transpose() << endl;
        cout << "[DEBUG] aDot_phys = " << aDotNow.transpose() << endl;

        coeffNutAvgLocal = interpNutAvgLinear(muStarVal);

        for (int i = 0; i < nphiNutFluct; ++i)
        {
            coeffNutFluctLocal(i) = rbfSplinesNutFluct[i]->eval(rbfInputFluct);
        }

        newtonObjectPPEAve.gNutAve   = coeffNutAvgLocal;
        newtonObjectPPEAve.gNutFluct = coeffNutFluctLocal;

        cout << std::setprecision(6)
             << "[DEBUG] coeffNutAvg:  ||·||2=" << coeffNutAvgLocal.norm()
             << "  min/max=(" << vecMin(coeffNutAvgLocal) << "," << vecMax(coeffNutAvgLocal)
             << ")  head=[" << printHead(coeffNutAvgLocal, 5) << "]\n";
        cout << "[DEBUG] coeffNutFluct:||·||2=" << coeffNutFluctLocal.norm()
             << "  min/max=(" << vecMin(coeffNutFluctLocal) << "," << vecMax(coeffNutFluctLocal)
             << ")  head=[" << printHead(coeffNutFluctLocal, 5) << "]\n";

        avgLog << step << " " << timeLocal << " " << muStarVal << " ";
        for (int k = 0; k < coeffNutAvgLocal.size(); ++k)
        {
            avgLog << coeffNutAvgLocal(k) << (k + 1 < coeffNutAvgLocal.size() ? ' ' : '\n');
        }

        flLog << step << " " << timeLocal << " ";
        for (int k = 0; k < coeffNutFluctLocal.size(); ++k)
        {
            flLog << coeffNutFluctLocal(k) << (k + 1 < coeffNutFluctLocal.size() ? ' ' : '\n');
        }

        if (useLift)
        {
            for (int j = 0; j < N_BC; ++j)
            {
                y(j) = vel_now(j, 0);
            }
        }

        newtonObjectPPEAve.yOldOld = newtonObjectPPEAve.y_old;
        newtonObjectPPEAve.y_old   = y;
        aOld = aNow;

        tmp_sol(0) = timeLocal;
        tmp_sol.tail(y.rows()) = y;

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

            if (counter2 < rbfCoeffMat.cols())
            {
                rbfCoeffMat(0, counter2) = timeLocal;
                rbfCoeffMat.block(1, counter2, nphiNutAvg, 1) = coeffNutAvgLocal;
                rbfCoeffMat.block(1 + nphiNutAvg, counter2, nphiNutFluct, 1) = coeffNutFluctLocal;
            }
            else
            {
                cout << "[WARN] rbfCoeffMat is full; skipping store at t=" << timeLocal << endl;
            }

            cout << "[DEBUG] Stored solution and RBF coeffs at counter2=" << counter2 << endl;

            nextStore += numberOfStores;
            ++counter2;
        }

        ++counter;
    }

    const int storedCols = counter2;

    cout << "[DEBUG] Trimming outputs: storedCols=" << storedCols
         << ", rbfCoeffMat orig cols=" << rbfCoeffMat.cols()
         << ", online_solution orig size=" << online_solution.size() << endl;

    if (storedCols > 0 && storedCols <= rbfCoeffMat.cols())
    {
        rbfCoeffMat.conservativeResize(rbfCoeffMat.rows(), storedCols);
    }
    if (storedCols >= 0 && storedCols <= online_solution.size())
    {
        online_solution.resize(storedCols);
    }

    cout << "\n[DEBUG] Exporting online_solution and rbfCoeffMat...\n";
    ITHACAstream::exportMatrix(online_solution, "red_coeff", "python", "./ITHACAoutput/red_coeff");
    ITHACAstream::exportMatrix(online_solution, "red_coeff", "matlab", "./ITHACAoutput/red_coeff");
    ITHACAstream::exportMatrix(rbfCoeffMat, "rbf_coeffs_online", "python", "./ITHACAoutput/red_coeff");
    ITHACAstream::exportMatrix(rbfCoeffMat, "rbf_coeffs_online", "matlab", "./ITHACAoutput/red_coeff");

    avgLog.close();
    flLog.close();

    cout << "[DEBUG] solveOnlinePPEAve() completed.\n\n";
}

void ReducedUnsteadyNSTurb::reconstruct(bool exportFields, fileName folder)
{
    if (exportFields)
    {
        mkDir(folder);
        ITHACAutilities::createSymLink(folder);
    }

    const int nphiNutAvgLocal   = nutAvgModes.size();
    const int nphiNutFluctLocal = nutFluctModes.size();
    const int nphiULocal        = Nphi_u;
    const int nphiPLocal        = Nphi_p;

    const int nSteps = std::min<int>(online_solution.size(), rbfCoeffMat.cols());

    uRecFields.clear();
    pRecFields.clear();
    nutRecFields.clear();

    uRecFields.setSize(nSteps);
    pRecFields.setSize(nSteps);
    nutRecFields.setSize(nSteps);

    for (int i = 0; i < nSteps; ++i)
    {
        const Eigen::VectorXd& yi = online_solution[i];

        const Eigen::VectorXd currentUCoeff = yi.block(1, 0, nphiULocal, 1);
        const Eigen::VectorXd currentPCoeff = yi.block(1 + nphiULocal, 0, nphiPLocal, 1);

        const Eigen::VectorXd nut_avg_coeffs =
            rbfCoeffMat.block(1, i, nphiNutAvgLocal, 1);
        const Eigen::VectorXd nut_fluct_coeffs =
            rbfCoeffMat.block(1 + nphiNutAvgLocal, i, nphiNutFluctLocal, 1);

        // U
        auto* uPtr = new volVectorField("uRec", Umodes[0] * 0.0);
        for (int j = 0; j < nphiULocal; ++j)
        {
            *uPtr += currentUCoeff(j) * Umodes[j];
        }
        uRecFields.set(i, uPtr);

        // P
        auto* pPtr = new volScalarField("pRec", Pmodes[0] * 0.0);
        for (int j = 0; j < nphiPLocal; ++j)
        {
            *pPtr += currentPCoeff(j) * Pmodes[j];
        }
        pRecFields.set(i, pPtr);

        // nut = avg + fluct
        auto* nutAvgPtr = new volScalarField("nutRecAvg", nutAvgModes[0] * 0.0);
        for (int j = 0; j < nphiNutAvgLocal; ++j)
        {
            *nutAvgPtr += nut_avg_coeffs(j) * nutAvgModes[j];
        }

        auto* nutFluctPtr = new volScalarField("nutRecFluct", nutFluctModes[0] * 0.0);
        for (int j = 0; j < nphiNutFluctLocal; ++j)
        {
            *nutFluctPtr += nut_fluct_coeffs(j) * nutFluctModes[j];
        }

        auto* nutPtr = new volScalarField("nutRec", (*nutAvgPtr) + (*nutFluctPtr));
        delete nutAvgPtr;
        delete nutFluctPtr;

        nutRecFields.set(i, nutPtr);
    }

    if (exportFields)
    {
        ITHACAstream::exportFields(uRecFields, folder, "uRec");
        ITHACAstream::exportFields(pRecFields, folder, "pRec");
        ITHACAstream::exportFields(nutRecFields, folder, "nutRec");
    }
}

Eigen::MatrixXd ReducedUnsteadyNSTurb::setOnlineVelocity(Eigen::MatrixXd vel)
{
    assert(problem->inletIndex.rows() == vel.rows()
           && "Imposed boundary conditions dimensions do not match given values matrix dimensions");

    Eigen::MatrixXd vel_scal(vel.rows(), vel.cols());

    for (int k = 0; k < problem->inletIndex.rows(); ++k)
    {
        const int p = problem->inletIndex(k, 0);
        const int l = problem->inletIndex(k, 1);

        const scalar area =
            gSum(problem->liftfield[0].mesh().magSf().boundaryField()[p]);

        const scalar u_lf =
            gSum(problem->liftfield[k].mesh().magSf().boundaryField()[p]
                 * problem->liftfield[k].boundaryField()[p]).component(l) / area;

        vel_scal(k, 0) = vel(k, 0) / u_lf;
    }

    return vel_scal;
}

// ************************************************************************* //
