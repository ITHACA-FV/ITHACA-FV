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

// Solve pressure equation
volScalarField rAU(1.0 / UEqn.A());
volVectorField HbyA(constrainHbyA(rAU* UEqn.H(), U, p));
surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::flux(HbyA)
    + fvc::interpolate(rAU) * fvc::ddtCorr(U, phi)
);

//MRF.makeRelative(phiHbyA);

adjustPhi(phiHbyA, U, p);

tmp<volScalarField> rAtU(rAU);
/*
if (pimple.consistent())
{
    rAtU = 1.0 / max(1.0 / rAU - UEqn.H1(), 0.1 / rAU);
    phiHbyA +=
        fvc::interpolate(rAtU() - rAU) * fvc::snGrad(p) * mesh.magSf();
    HbyA -= (rAU - rAtU()) * fvc::grad(p);
}

if (pimple.nCorrPISO() <= 1)
{
    tUEqn.clear();
}
*/
// Update the pressure BCs to ensure flux consistency
//constrainPressure(p, U, phiHbyA, rAtU(), MRF);
std::tuple<Eigen::MatrixXd,Eigen::VectorXd> HRSysP;
// Non-orthogonal pressure corrector loop
while (pimple.correctNonOrthogonal())
{
    // Pressure corrector
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAtU(), p) 
        == 
        fvc::div(phiHbyA)
    );
    pEqn.setReference(pRefCell, pRefValue);
    //pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));
    
    Foam2Eigen::fvMat2Eigen(pEqn, Sp, sp);
    HRSysP = HyperReducedSys(Sp, sp, PMagicPoints, DeimP);
    // Hyper-reduced system
    Eigen::MatrixXd A_p = std::get<0>(HRSysP);
    Eigen::VectorXd a_p = std::get<1>(HRSysP);
    /// Solve the Hyper-reduced system for p
    //Eigen::VectorXd c_p = A_p.completeOrthogonalDecomposition().solve(a_p);
    Eigen::VectorXd c_p = A_p.colPivHouseholderQr().solve(a_p);
    //Eigen::MatrixXd I_p = Eigen::MatrixXd::Identity(A_p.cols(), A_p.cols());
    //Eigen::VectorXd c_p = (A_p.transpose()*A_p + I_p).inverse()*A_p.transpose()*a_p;
    //Eigen::VectorXd c_p = (A_p.transpose()*A_p).inverse()*A_p.transpose()*a_p;
    //c_p = (c_p + lambda_t * p_ref) / (1.0 + lambda_t);
    //Eigen::VectorXd c_p = A_p.ldlt().solve(a_p);
    //Eigen::VectorXd c_p = A_p.ldlt().solve(a_p);
    Eigen::VectorXd z_p = problem->DeimP->U*c_p;
    /// Eigen2field
    p = Foam2Eigen::Eigen2field(p, z_p); //Okay
    //std::copy(z_p.data(), z_p.data() + z_p.size(), p.ref().begin()); // Okay
    p.correctBoundaryConditions();
    if (pimple.finalNonOrthogonalIter())
    {
        phi = phiHbyA - pEqn.flux();
    }
}

//#include "continuityErrs.H"
// Explicitly relax pressure for momentum corrector
p.relax();
U = HbyA - rAtU() * fvc::grad(p);
U.correctBoundaryConditions();
fvOptions.correct(U);
