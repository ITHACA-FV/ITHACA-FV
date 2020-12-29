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

#include "ReducedSteadyNSTurbIntrusive.H"

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructor
ReducedSteadyNSTurbIntrusive::ReducedSteadyNSTurbIntrusive()
{
}

ReducedSteadyNSTurbIntrusive::ReducedSteadyNSTurbIntrusive(
    SteadyNSTurbIntrusive& fomProblem)
    :
    problem(&fomProblem)
{
    N_BC = problem->inletIndex.rows();
    Nphi_u = problem->bMatrix.rows();

    for (int k = 0; k < Nphi_u; k++)
    {
        Umodes.append(tmp<volVectorField>(problem->Umodes[k]));
    }

    newtonObject = newtonSteadyNSTurbIntrusive(Nphi_u, Nphi_u,
                   fomProblem);
}

int newtonSteadyNSTurbIntrusive::operator()(const Eigen::VectorXd& x,
        Eigen::VectorXd& fvec) const
{
    Eigen::VectorXd aTmp(Nphi_u);
    aTmp = x;
    // Convective term
    Eigen::MatrixXd cc(1, 1);
    // Mom Term
    Eigen::VectorXd m1 = problem->bTotalMatrix * aTmp * nu;
    // Gradient of pressure
    Eigen::VectorXd m2 = problem->kMatrix * aTmp;
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
        cc = aTmp.transpose() * Eigen::SliceFromTensor(problem->cTotalTensor, 0,
                i) * aTmp;
        fvec(i) = m1(i) - cc(0, 0) - m2(i);

        if (problem->bcMethod == "penalty")
        {
            fvec(i) += ((penaltyU * tauU)(i, 0));
        }
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

int newtonSteadyNSTurbIntrusive::df(const Eigen::VectorXd& x,
                                    Eigen::MatrixXd& fjac) const
{
    Eigen::NumericalDiff<newtonSteadyNSTurbIntrusive> numDiff(*this);
    numDiff.df(x, fjac);
    return 0;
}


// * * * * * * * * * * * * * * * Solve Functions  * * * * * * * * * * * * * //


void ReducedSteadyNSTurbIntrusive::solveOnline(Eigen::MatrixXd vel)
{
    if (problem->bcMethod == "lift")
    {
        vel_now = setOnlineVelocity(vel);
    }
    else if (problem->bcMethod == "penalty")
    {
        vel_now = vel;
    }

    y.resize(Nphi_u, 1);
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
    Eigen::HybridNonLinearSolver<newtonSteadyNSTurbIntrusive> hnls(newtonObject);
    newtonObject.bc.resize(N_BC);
    newtonObject.tauU = tauU;

    for (int j = 0; j < N_BC; j++)
    {
        newtonObject.bc(j) = vel_now(j, 0);
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


void ReducedSteadyNSTurbIntrusive::reconstruct(fileName folder,
        int printEvery)
{
    mkDir(folder);
    ITHACAutilities::createSymLink(folder);
    int counter = 0;
    int nextWrite = 0;

    for (int i = 0; i < online_solution.size(); i++)
    {
        if (counter == nextWrite)
        {
            volVectorField uRec("uRec", Umodes[0] * 0);
            volScalarField pRec("pRec", problem->Pmodes[0] * 0);
            volScalarField nutTemp("nutTemp", problem->nutModes[0] * 0);

            for (int j = 0; j < Nphi_u; j++)
            {
                uRec += Umodes[j] * online_solution[i](j + 1, 0);
                pRec += problem->Pmodes[j] * online_solution[i](j + 1, 0);
                nutTemp += problem->nutModes[j] * online_solution[i](j + 1, 0);
            }

            ITHACAstream::exportSolution(uRec, name(online_solution[i](0, 0)), folder);
            ITHACAstream::exportSolution(pRec, name(online_solution[i](0, 0)), folder);
            nextWrite += printEvery;
            UREC.append(tmp<volVectorField>(uRec));
            PREC.append(tmp<volScalarField>(pRec));
            nutRec.append(tmp<volScalarField>(nutTemp));
        }

        counter++;
    }
}

Eigen::MatrixXd ReducedSteadyNSTurbIntrusive::setOnlineVelocity(
    Eigen::MatrixXd vel)
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

void ReducedSteadyNSTurbIntrusive::reconstructLiftAndDrag(
    SteadyNSTurbIntrusive& problem,
    fileName folder)
{
    mkDir(folder);
    system("ln -s ../../constant " + folder + "/constant");
    system("ln -s ../../0 " + folder + "/0");
    system("ln -s ../../system " + folder + "/system");
    Eigen::VectorXd cl(online_solution.size());
    Eigen::VectorXd cd(online_solution.size());
    //Read FORCESdict
    IOdictionary FORCESdict
    (
        IOobject
        (
            "FORCESdict",
            "./system",
            Umodes[0].mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    fTau.setZero(online_solution.size(), 3);
    fN.setZero(online_solution.size(), 3);

    for (int i = 0; i < online_solution.size(); i++)
    {
        for (int j = 0; j < Nphi_u; j++)
        {
            fTau.row(i) += problem.tauMatrix.row(j) * online_solution[i](j + 1, 0);
        }

        for (int j = 0; j < Nphi_u; j++)
        {
            fN.row(i) += problem.nMatrix.row(j) * online_solution[i](j + 1, 0);
        }
    }

    // Export the matrices
    if (para->exportPython)
    {
        ITHACAstream::exportMatrix(fTau, "fTau", "python", folder);
        ITHACAstream::exportMatrix(fN, "fN", "python", folder);
    }

    if (para->exportMatlab)
    {
        ITHACAstream::exportMatrix(fTau, "fTau", "matlab", folder);
        ITHACAstream::exportMatrix(fN, "fN", "matlab", folder);
    }

    if (para->exportTxt)
    {
        ITHACAstream::exportMatrix(fTau, "fTau", "eigen", folder);
        ITHACAstream::exportMatrix(fN, "fN", "eigen", folder);
    }
}
// ************************************************************************* //


// ************************************************************************* //

