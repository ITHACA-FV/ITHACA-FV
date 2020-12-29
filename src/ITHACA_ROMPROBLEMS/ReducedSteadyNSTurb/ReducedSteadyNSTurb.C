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

#include "ReducedSteadyNSTurb.H"

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructor
ReducedSteadyNSTurb::ReducedSteadyNSTurb()
{
}

ReducedSteadyNSTurb::ReducedSteadyNSTurb(SteadyNSTurb& fomProblem)
    :
    problem(&fomProblem)
{
    N_BC = problem->inletIndex.rows();
    Nphi_u = problem->B_matrix.rows();
    Nphi_p = problem->K_matrix.cols();
    nphiNut = problem->cTotalTensor.dimension(1);

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

    newtonObject = newtonSteadyNSTurb(Nphi_u + Nphi_p, Nphi_u + Nphi_p,
                                      fomProblem);
}

int newtonSteadyNSTurb::operator()(const Eigen::VectorXd& x,
                                   Eigen::VectorXd& fvec) const
{
    Eigen::VectorXd aTmp(Nphi_u);
    Eigen::VectorXd bTmp(Nphi_p);
    aTmp = x.head(Nphi_u);
    bTmp = x.tail(Nphi_p);
    // Convective term
    Eigen::MatrixXd cc(1, 1);
    // Mom Term
    Eigen::VectorXd m1 = problem->bTotalMatrix * aTmp * nu;
    // Gradient of pressure
    Eigen::VectorXd m2 = problem->K_matrix * bTmp;
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
        fvec(i) = m1(i) - cc(0, 0) - m2(i);

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

int newtonSteadyNSTurb::df(const Eigen::VectorXd& x,
                           Eigen::MatrixXd& fjac) const
{
    Eigen::NumericalDiff<newtonSteadyNSTurb> numDiff(*this);
    numDiff.df(x, fjac);
    return 0;
}


// * * * * * * * * * * * * * * * Solve Functions  * * * * * * * * * * * * * //


void ReducedSteadyNSTurb::solveOnlineSUP(Eigen::MatrixXd vel)
{
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

    // Change initial condition for the lifting function
    if (problem->bcMethod == "lift")
    {
        for (int j = 0; j < N_BC; j++)
        {
            y(j) = vel_now(j, 0);
        }
    }

    Color::Modifier red(Color::FG_RED);
    Color::Modifier green(Color::FG_GREEN);
    Color::Modifier def(Color::FG_DEFAULT);
    Eigen::HybridNonLinearSolver<newtonSteadyNSTurb> hnls(newtonObject);
    newtonObject.bc.resize(N_BC);
    newtonObject.tauU = tauU;

    for (int j = 0; j < N_BC; j++)
    {
        newtonObject.bc(j) = vel_now(j, 0);
    }

    if (problem->viscCoeff == "L2")
    {
        for (int i = 0; i < nphiNut; i++)
        {
            newtonObject.gNut = problem->nutCoeff;
        }
    }
    else if (problem->viscCoeff == "RBF")
    {
        for (int i = 0; i < nphiNut; i++)
        {
            newtonObject.gNut(i) = problem->rbfSplines[i]->eval(vel_now);
            rbfCoeff = newtonObject.gNut;
        }
    }
    else
    {
        Info << "The way to compute the eddy viscosity coefficients has to be either L2 or RBF"
             << endl;
        exit(0);
    }

    newtonObject.nu = nu;
    hnls.solve(y);
    Eigen::VectorXd res(y);
    newtonObject.operator()(y, res);
    std::cout << "################## Online solve N° " << count_online_solve <<
              " ##################" << std::endl;
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
}


void ReducedSteadyNSTurb::reconstruct(bool exportFields, fileName folder,
                                      int printevery)
{
    if (exportFields)
    {
        mkDir(folder);
        ITHACAutilities::createSymLink(folder);
    }

    int counter = 0;
    int nextWrite = 0;
    List <Eigen::MatrixXd> CoeffU;
    List <Eigen::MatrixXd> CoeffP;
    List <Eigen::MatrixXd> CoeffNut;
    CoeffU.resize(0);
    CoeffP.resize(0);
    CoeffNut.resize(0);

    for (int i = 0; i < online_solution.size(); i++)
    {
        if (counter == nextWrite)
        {
            Eigen::MatrixXd currentUCoeff;
            Eigen::MatrixXd currentPCoeff;
            Eigen::MatrixXd currentNutCoeff;
            currentUCoeff = online_solution[i].block(1, 0, Nphi_u, 1);
            currentPCoeff = online_solution[i].bottomRows(Nphi_p);
            currentNutCoeff = rbfCoeffMat.block(0, i, nphiNut, 1);
            CoeffU.append(currentUCoeff);
            CoeffP.append(currentPCoeff);
            CoeffNut.append(currentNutCoeff);
            nextWrite += printevery;
        }

        counter++;
    }

    volVectorField uRec("uRec", Umodes[0]);
    volScalarField pRec("pRec", problem->Pmodes[0]);
    volScalarField nutRec("nutRec", problem->nutModes[0]);
    uRecFields = problem->L_U_SUPmodes.reconstruct(uRec, CoeffU, "uRec");
    pRecFields = problem->Pmodes.reconstruct(pRec, CoeffP, "pRec");
    nutRecFields = problem->nutModes.reconstruct(nutRec, CoeffNut, "nutRec");

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

Eigen::MatrixXd ReducedSteadyNSTurb::setOnlineVelocity(Eigen::MatrixXd vel)
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

