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
/// Source file of the HyperReducedCompressibleUnSteadyNS class

#include "HyperReducedCompressibleUnSteadyNS.H"

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructor
HyperReducedCompressibleUnSteadyNS::HyperReducedCompressibleUnSteadyNS()
{
}

HyperReducedCompressibleUnSteadyNS::HyperReducedCompressibleUnSteadyNS(
    CompressibleUnSteadyRhoPimple& FOMproblem)
    :
    problem(&FOMproblem)
{
    for (int i = 0; i < problem->Umodes.size(); i++)
    {
        Umodes.append((problem->Umodes.toPtrList()[i]).clone());
    }

    for (int i = 0; i < problem->Pmodes.size(); i++)
    {
        Pmodes.append((problem->Pmodes.toPtrList()[i]).clone());
    }

    for (int i = 0; i < problem->Emodes.size(); i++)
    {
        Emodes.append((problem->Emodes.toPtrList()[i]).clone());
    }

    //problem->restart();
    std::cout << "################ HyperReduced Ctor called ##################" <<
              std::endl;
}



// * * * * * * * * * * * * * * * Solve Functions  * * * * * * * * * * * * * //

void HyperReducedCompressibleUnSteadyNS::SolveHyperReducedSys(int NmodesUproj,
        int NmodesPproj,
        int NmodesEproj,
        fileName folder)
{
    Time& runTime = problem->_runTime();
    volVectorField& U = problem->_U();
    volScalarField& p = problem->pThermo().p();
    volScalarField& E = problem->pThermo().he();
    pimpleControl& pimple = problem->_pimple();
    volScalarField _nut(problem->turbulence->nut());
    instantList Times = runTime.times();
    runTime.setEndTime(finalTime);
    runTime.setTime(Times[1], 1);
    runTime.setDeltaT(timeStep);
    nextWrite = startTime; // timeStep initialization
    int counter = 1;
    bool  correctPhi = problem->correctPhi;
    bool  checkMeshCourantNo = problem->checkMeshCourantNo;
    bool  moveMeshOuterCorrectors = problem->moveMeshOuterCorrectors;
    scalar  cumulativeContErr = problem->cumulativeContErr;
#include "HyperRedSolvers.H"
}
std::tuple<Eigen::MatrixXd, Eigen::VectorXd>
HyperReducedCompressibleUnSteadyNS::HyperReducedSys(Eigen::SparseMatrix<double>&
        S,
        Eigen::VectorXd& se,
        List<label>& uniqueMagicPoints,
        Eigen::MatrixXd& Modes)
{
    dynamicFvMesh& mesh = problem->meshPtr();
    int s = uniqueMagicPoints.size();
    Eigen::MatrixXd B;
    int m = Modes.cols(); // number of modes
    B.setZero(s, m);
    Eigen::VectorXd b;
    b.setZero(s);

    for (int k = 0; k < uniqueMagicPoints.size(); ++k)
    {
        int magicPoint = uniqueMagicPoints[k];
        Eigen::SparseVector<double> vec = S.row(magicPoint);
        Eigen::RowVectorXd tempB = Eigen::RowVectorXd::Zero(m);

        for (Eigen::SparseVector<double>::InnerIterator it(vec); it; ++it)
        {
            tempB += it.value() * Modes.row(it.index());
        }

        // Apply scaling by volume
        double volume = 1.0 / mesh.V()[magicPoint];
        B.row(k) = tempB * volume;
        b(k) = volume * se(magicPoint);
    }

    std::tuple<Eigen::MatrixXd, Eigen::VectorXd> HRSys;
    HRSys = std::make_tuple(B, b);
    return HRSys;
}
bool HyperReducedCompressibleUnSteadyNS::checkWrite(Time& timeObject)
{
    scalar diffnow = mag(nextWrite - atof(timeObject.timeName().c_str()));
    scalar diffnext = mag(nextWrite - atof(timeObject.timeName().c_str()) -
                          timeObject.deltaTValue());

    if ( diffnow < diffnext)
    {
        return true;
    }
    else
    {
        return false;
    }
}


// void HyperReducedCompressibleUnSteadyNS::setOnlineVelocity(Eigen::MatrixXd vel)
// {
//     M_Assert(problem->inletIndex.rows() == vel.size(),
//              "Imposed boundary conditions dimensions do not match given values matrix dimensions into setOnlineVelocity");
//     Eigen::MatrixXd vel_scal;
//     vel_scal.resize(vel.rows(), vel.cols());

//     for (int k = 0; k < problem->inletIndex.rows(); k++)
//     {
//         label p = problem->inletIndex(k, 0);
//         label l = problem->inletIndex(k, 1);
//         scalar area = gSum(problem->liftfield[0].mesh().magSf().boundaryField()[p]);
//         scalar u_lf = gSum(problem->liftfield[k].mesh().magSf().boundaryField()[p] *
//                            problem->liftfield[k].boundaryField()[p]).component(l) / area;
//         vel_scal(k, 0) = vel(k, 0) / u_lf;
//     }

//     vel_now = vel_scal;
// }


// void HyperReducedCompressibleUnSteadyNS::projectReducedOperators(int NmodesUproj, int NmodesPproj, int NmodesEproj)
// {
//     PtrList<volVectorField> gradModP;
//     for (label i = 0; i < NmodesPproj; i++)
//     {
//         gradModP.append(fvc::grad(problem->Pmodes[i]));
//     }
//     projGradModP = problem->Umodes.project(gradModP, NmodesUproj);
// }
