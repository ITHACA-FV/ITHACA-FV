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
/// Source file of the reducedSteadyNS class

#include "ReducedSimpleSteadyNS.H"

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructor
reducedSimpleSteadyNS::reducedSimpleSteadyNS()
{
}

reducedSimpleSteadyNS::reducedSimpleSteadyNS(SteadyNSSimple& FOMproblem)
    :
    problem(&FOMproblem)
{
    // Create a new Umodes set where the first ones are the lift functions
    for (int i = 0; i < problem->inletIndex.rows(); i++)
    {
        ULmodes.append(problem->liftfield[i]);
    }

    for (int i = 0; i < problem->Umodes.size(); i++)
    {
        ULmodes.append(problem->Umodes.toPtrList()[i]);
    }
}

// * * * * * * * * * * * * * * * Solve Functions  * * * * * * * * * * * * * //

void reducedSimpleSteadyNS::solveOnline_Simple(scalar mu_now,
        scalar NmodesUproj, scalar NmodesPproj)
{
    counter++;
    scalar UprojN;
    scalar PprojN;

    if (NmodesUproj == 0)
    {
        UprojN = ULmodes.size();
    }
    else
    {
        UprojN = NmodesUproj;
    }

    if (NmodesPproj == 0)
    {
        PprojN = problem->Pmodes.size();
    }
    else
    {
        PprojN = NmodesPproj;
    }

    Eigen::VectorXd uresidualOld = Eigen::VectorXd::Zero(UprojN);
    Eigen::VectorXd presidualOld = Eigen::VectorXd::Zero(UprojN);
    Eigen::VectorXd uresidual;
    Eigen::VectorXd presidual;
    scalar U_norm_res(1);
    scalar P_norm_res(1);
    Eigen::MatrixXd a = Eigen::VectorXd::Zero(UprojN);
    Eigen::MatrixXd b = Eigen::VectorXd::Zero(PprojN);
    ITHACAparameters para;
    float residualJumpLim =
        para.ITHACAdict->lookupOrDefault<float>("residualJumpLim", 1e-5);
    float normalizedResidualLim =
        para.ITHACAdict->lookupOrDefault<float>("normalizedResidualLim", 1e-5);
    scalar residual_jump(1 + residualJumpLim);
    volVectorField Uaux("Uaux", problem->_U());
    volScalarField Paux("Paux", problem->_p());
    int iter = 0;

    while (residual_jump > residualJumpLim
            || std::max(U_norm_res, P_norm_res) > normalizedResidualLim)
    {
        iter++;
        Uaux = ULmodes.reconstruct(a, "Uaux");
        Paux = problem->Pmodes.reconstruct(b, "Paux");
        simpleControl& simple = problem->_simple();
        setRefCell(Paux, simple.dict(), problem->pRefCell, problem->pRefValue);
        problem->_phi() = linearInterpolate(Uaux) & problem->_U().mesh().Sf();
        fvVectorMatrix Au(get_Umatrix_Online(Uaux, Paux));
        List<Eigen::MatrixXd> RedLinSysU = ULmodes.project(Au, UprojN);
        a = reducedProblem::solveLinearSys(RedLinSysU, a, uresidual, vel_now);
        //Info << uresidual.norm() << endl;
        Uaux = ULmodes.reconstruct(a, "Uaux");
        problem->_phi() = linearInterpolate(Uaux) & problem->_U().mesh().Sf();
        fvScalarMatrix Ap(get_Pmatrix_Online(Uaux, Paux));
        List<Eigen::MatrixXd> RedLinSysP = problem->Pmodes.project(Ap, PprojN);
        b = reducedProblem::solveLinearSys(RedLinSysP, b, presidual);
        //Info << presidual.norm() << endl;
        uresidualOld = uresidualOld - uresidual;
        presidualOld = presidualOld - presidual;
        uresidualOld = uresidualOld.cwiseAbs();
        presidualOld = presidualOld.cwiseAbs();
        //std::cout << uresidualOld.sum() << std::endl;
        //std::cout << presidualOld.sum() << std::endl;
        residual_jump = std::max(uresidualOld.sum(), presidualOld.sum());
        uresidualOld = uresidual;
        presidualOld = presidual;
        uresidual = uresidual.cwiseAbs();
        presidual = presidual.cwiseAbs();
        U_norm_res = uresidual.sum() / (RedLinSysU[1].cwiseAbs()).sum();
        P_norm_res = presidual.sum() / (RedLinSysP[1].cwiseAbs()).sum();
        // std::cout << "Residual jump = " << residual_jump << std::endl;
        // std::cout << "Normalized residual = " << std::max(U_norm_res,P_norm_res) << std::endl;
    }

    std::cout << "Solution " << counter << " converged in " << iter <<
              " iterations." << std::endl;
    std::cout << "Final normalized residual for velocity: " << U_norm_res <<
              std::endl;
    std::cout << "Final normalized residual for pressure: " << P_norm_res <<
              std::endl;
    Uaux = ULmodes.reconstruct(a, "Uaux");
    Paux = problem->Pmodes.reconstruct(b, "Paux");
    ITHACAstream::exportSolution(Uaux, name(counter),
                                 "./ITHACAoutput/Reconstruct/");
    ITHACAstream::exportSolution(Paux, name(counter),
                                 "./ITHACAoutput/Reconstruct/");
}

fvVectorMatrix reducedSimpleSteadyNS::get_Umatrix_Online(volVectorField& U,
        volScalarField& p)
{
    IOMRFZoneList& MRF = problem->_MRF();
    surfaceScalarField& phi = problem->_phi();
    fv::options& fvOptions = problem->_fvOptions();
    MRF.correctBoundaryVelocity(U);
    fvVectorMatrix Ueqn
    (
        fvm::div(phi, U)
        + MRF.DDt(U)
        + problem->turbulence->divDevReff(U)
        ==
        fvOptions(U)
    );
    problem->Ueqn_global = &Ueqn;
    return Ueqn;
}

fvScalarMatrix reducedSimpleSteadyNS::get_Pmatrix_Online(volVectorField& U,
        volScalarField& p)
{
    IOMRFZoneList& MRF = problem->_MRF();
    MRF.correctBoundaryVelocity(U);
    volScalarField rAU(1.0 / problem->Ueqn_global->A());
    volVectorField HbyA(constrainHbyA(rAU * problem->Ueqn_global->H(), U, p));
    surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
    MRF.makeRelative(phiHbyA);
    adjustPhi(phiHbyA, U, p);
    tmp<volScalarField> rAtU(rAU);
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA)
    );
    return pEqn;
}


void reducedSimpleSteadyNS::setOnlineVelocity(Eigen::MatrixXd vel)
{
    M_Assert(problem->inletIndex.rows() == vel.size(),
             "Imposed boundary conditions dimensions do not match given values matrix dimensions");
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
    vel_now = vel_scal;
}
