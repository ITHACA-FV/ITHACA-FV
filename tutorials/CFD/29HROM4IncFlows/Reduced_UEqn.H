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
// Solve the Momentum equation
//MRF.correctBoundaryVelocity(U);
fvVectorMatrix UEqn
(
    fvm::ddt(U) 
    + fvm::div(phi, U)
    //+ fvc::div(phi, U)
    //+ MRF.DDt(U)
    + turbulence->divDevReff(U)
    == fvOptions(U) //-fvc::div(phi, U) + fvOptions(U)
);
//fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();
fvOptions.constrain(UEqn);

std::tuple<Eigen::MatrixXd,Eigen::VectorXd> HRSysU;
//List<Eigen::MatrixXd> RedLinSysU;
if (pimple.momentumPredictor())
{
    // solve(UEqn == -fvc::grad(p));
    fvVectorMatrix UEqnWithPressure(UEqn == -fvc::grad(p));    
    Foam2Eigen::fvMat2Eigen(UEqnWithPressure, Su, su);  // 8.16779--333.  
    HRSysU = HyperReducedSys(Su, su, UMagicPoints, DeimU);
    //auto t2 = std::chrono::high_resolution_clock::now();
    //auto time_span_u = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    //time_rom_u +=time_span_u.count(); 
    // Hyper-reduced system
    Eigen::MatrixXd A_u = std::get<0>(HRSysU);
    Eigen::VectorXd a_u = std::get<1>(HRSysU);
    //Eigen::MatrixXd I_u = Eigen::MatrixXd::Identity(A_u.cols(), A_u.cols());
    Eigen::VectorXd c_u = A_u.colPivHouseholderQr().solve(a_u);
    //Eigen::VectorXd c_u =  (A_u.transpose()*A_u + I_u).inverse()*A_u.transpose()*a_u;
    //Eigen::VectorXd c_u =  (A_u.transpose()*A_u).inverse()*A_u.transpose()*a_u;
    //c_u = (c_u + lambda_t * u_ref) / (1.0 + lambda_t);
    Eigen::VectorXd z_u = problem->DeimU->U*c_u;
    U = Foam2Eigen::Eigen2field(U, z_u);
    U.correctBoundaryConditions();
    fvOptions.correct(U);
}
