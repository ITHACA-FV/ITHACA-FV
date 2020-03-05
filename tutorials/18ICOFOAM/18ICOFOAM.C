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
Description
    Example of steady NS Reduction Problem
SourceFiles
    03steadyNS.C
\*---------------------------------------------------------------------------*/

#include "unsteadyNS.H"
#include "ITHACAPOD.H"
#include "ReducedUnsteadyNS.H"
#include "ITHACAstream.H"
#include <chrono>
#include<math.h>
#include<iomanip>


class tutorial18 : public unsteadyNS
{
    public:
        /// Constructor
        explicit tutorial18(int argc, char* argv[])
        {
            _args = autoPtr<argList>
                    (
                        new argList(argc, argv)
                    );

            if (!_args->checkRootCase())
            {
                Foam::FatalError.exit();
            }

            argList& args = _args();
#include "createTime.H"
#include "createMesh.H"
            std::cerr << "debug point 1" << std::endl;
#include "createFields.H"
            std::cerr << "debug point 2" << std::endl;
            _piso = autoPtr<pisoControl>
                    (
                        new pisoControl
                        (
                            mesh
                        )
                    );
            para = new ITHACAparameters;
            offline = ITHACAutilities::check_off();
            podex = ITHACAutilities::check_pod();
            std::cerr << "debug point 1" << std::endl;
        }

        autoPtr<IOdictionary> _transportProperties;
        autoPtr<dimensionedScalar> _nu;
        autoPtr<pisoControl> _piso;





        /// Perform an Offline solve
        void offlineSolve()
        {
            volVectorField& U = _U();
            volScalarField& p = _p();
            Vector<double> inl(1, 0, 0);
            List<scalar> mu_now(1);

            if (offline)
            {
                ITHACAstream::read_fields(Ufield, U, "./ITHACAoutput/Offline/");
                ITHACAstream::read_fields(Pfield, p, "./ITHACAoutput/Offline/");
            }

            else
            {
                for (label i = 0; i < mu.cols(); i++)
                {
                    mu_now[0] = mu(0, i);
                    truthSolve(mu_now, "./ITHACAoutput/Offline/");
                }
            }
        }

        void truthSolve(List<scalar> mu_now, word folder)
        {
            Time& runTime = _runTime();
            surfaceScalarField& phi = _phi();
            fvMesh& mesh = _mesh();
#include "initContinuityErrs.H"
            pisoControl piso = _piso();
            volVectorField& U = _U();
            volScalarField& p = _p();
            dimensionedScalar nu = _nu();
            instantList Times = runTime.times();
            runTime.setEndTime(finalTime);
            runTime.setTime(Times[1], 1);
            runTime.setDeltaT(timeStep);
            nextWrite = startTime;
            // Export and store the initial conditions for velocity and pressure
            ITHACAstream::exportSolution(U, name(counter), folder);
            ITHACAstream::exportSolution(p, name(counter), folder);
            std::ofstream of(folder + name(counter) + "/" +
                             runTime.timeName());
            Ufield.append(U);
            Pfield.append(p);
            counter = 0;
            counter++;
            nextWrite += writeEvery;

            while (runTime.run())
            {
                runTime.setEndTime(finalTime);
                runTime++;
                Info << "Time = " << runTime.timeName() << nl << endl;
                // Momentum predictor
                fvVectorMatrix UEqn
                (
                    fvm::ddt(U)
                    + fvm::div(phi, U)
                    - fvm::laplacian(nu, U)
                );

                if (piso.momentumPredictor())
                {
                    solve(UEqn == -fvc::grad(p));
                }

                // --- PISO loop
                while (piso.correct())
                {
                    volScalarField rAU(1.0 / UEqn.A());
                    volVectorField HbyA(constrainHbyA(rAU * UEqn.H(), U, p));
                    surfaceScalarField phiHbyA
                    (
                        "phiHbyA",
                        fvc::flux(HbyA)
                        + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
                    );
                    adjustPhi(phiHbyA, U, p);
                    // Update the pressure BCs to ensure flux consistency
                    constrainPressure(p, U, phiHbyA, rAU);

                    // Non-orthogonal pressure corrector loop
                    while (piso.correctNonOrthogonal())
                    {
                        // Pressure corrector
                        fvScalarMatrix pEqn
                        (
                            fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
                        );
                        pEqn.setReference(0, 0);
                        pEqn.solve();

                        if (piso.finalNonOrthogonalIter())
                        {
                            phi = phiHbyA - pEqn.flux();
                        }
                    }

                    U = HbyA - rAU * fvc::grad(p);
                    U.correctBoundaryConditions();
                }

                if (checkWrite(runTime))
                {
                    ITHACAstream::exportSolution(U, name(counter), folder);
                    ITHACAstream::exportSolution(p, name(counter), folder);
                    std::ofstream of(folder + name(counter) + "/" +
                                     runTime.timeName());
                    Ufield.append(U);
                    Pfield.append(p);
                    counter++;
                    nextWrite += writeEvery;
                }
            }

            Info << "End\n" << endl;
        }
        void restart()
        {
            volScalarField& p = _p();
            volScalarField& p0 = _p0();
            volVectorField& U = _U();
            volVectorField& U0 = _U0();
            surfaceScalarField& phi = _phi();
            surfaceScalarField& phi0 = _phi0();
            p = p0;
            U = U0;
            phi = phi0;
        }
};

class tutorial18red : public reducedUnsteadyNS
{
    public:

        tutorial18red(tutorial18& FOMproblem)
            :
            problem(&FOMproblem)
        {}

        int counter;

        Modes<vector> ULmodes;
        /// Imposed boundary conditions.
        Eigen::MatrixXd vel_now;
        Eigen::MatrixXd redGradP;
        tutorial18* problem;

        void solveOnlineICO(int NmodesUproj, int NmodesPproj, word folder)
        {
            problem->restart();
            ULmodes.resize(0);

            for (int i = 0; i < problem->inletIndex.rows(); i++)
            {
                ULmodes.append(problem->liftfield[i]);
            }

            for (int i = 0; i < NmodesUproj; i++)
            {
                ULmodes.append(problem->Umodes.toPtrList()[i]);
            }

            int UprojN = ULmodes.size();
            int PprojN = NmodesPproj;
            Eigen::VectorXd uresidualOld = Eigen::VectorXd::Zero(UprojN);
            Eigen::VectorXd presidualOld = Eigen::VectorXd::Zero(PprojN);
            Eigen::VectorXd uresidual = Eigen::VectorXd::Zero(UprojN);
            Eigen::VectorXd presidual = Eigen::VectorXd::Zero(PprojN);
            scalar U_norm_res(1);
            scalar P_norm_res(1);
            Eigen::MatrixXd a = Eigen::VectorXd::Zero(UprojN);
            Eigen::MatrixXd b = Eigen::VectorXd::Zero(PprojN);
            a(0) = vel_now(0, 0);
            /// Velocity field
            volVectorField& U = problem->_U();
            /// Pressure field
            volScalarField& p = problem->_p();
            ///
            surfaceScalarField& phi = problem->_phi();
            Time& runTime = problem->_runTime();
            fvMesh& mesh = problem->_mesh();
            pisoControl piso = problem->_piso();
            dimensionedScalar nu = problem->_nu();
            instantList Times = runTime.times();
            runTime.setEndTime(problem->finalTime);
            runTime.setTime(Times[1], 1);
            runTime.timeName();
            runTime.setDeltaT(problem->timeStep);
            problem->nextWrite = problem->startTime;
            U.rename("Ured");
            p.rename("Pred");
            // Export and store the initial conditions for velocity and pressure
            counter = 0;
            ITHACAstream::exportSolution(U, name(counter), folder);
            ITHACAstream::exportSolution(p, name(counter), folder);
            U.rename("U");
            p.rename("p");
            counter++;
            problem->nextWrite += problem->writeEvery;
            std::cerr << "debug point 1" << std::endl;

            while (runTime.run())
            {
                runTime++;
                Info << "Time = " << runTime.timeName() << nl << endl;
                // Momentum predictor
                fvVectorMatrix UEqn
                (
                    fvm::ddt(U)
                    + fvm::div(phi, U)
                    - fvm::laplacian(nu, U)
                );
                UEqn.relax();
                List<Eigen::MatrixXd> RedLinSysU = ULmodes.project(UEqn, UprojN);
                volVectorField H = UEqn.H();
                volScalarField A = UEqn.A();
                RedLinSysU[1] += redGradP * b;
                a = reducedProblem::solveLinearSys(RedLinSysU, a, uresidual, vel_now);
                ULmodes.reconstruct(U, a, "U");
                List<Eigen::MatrixXd> RedLinSysP;

                // if (piso.momentumPredictor())
                // {
                //     solve(UEqn == -fvc::grad(p));
                // }

                // --- PISO loop
                while (piso.correct())
                {
                    volScalarField rAU(1.0 / UEqn.A());
                    volVectorField HbyA(constrainHbyA(rAU * UEqn.H(), U, p));
                    surfaceScalarField phiHbyA
                    (
                        "phiHbyA",
                        fvc::flux(HbyA)
                        + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
                    );
                    adjustPhi(phiHbyA, U, p);
                    // Update the pressure BCs to ensure flux consistency
                    constrainPressure(p, U, phiHbyA, rAU);

                    // Non-orthogonal pressure corrector loop
                    while (piso.correctNonOrthogonal())
                    {
                        // Pressure corrector
                        fvScalarMatrix pEqn
                        (
                            fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
                        );
                        pEqn.setReference(0, 0);
                        RedLinSysP = problem->Pmodes.project(pEqn, PprojN);
                        b = reducedProblem::solveLinearSys(RedLinSysP, b, presidual);
                        problem->Pmodes.reconstruct(p, b, "p");

                        //pEqn.solve();

                        if (piso.finalNonOrthogonalIter())
                        {
                            phi = phiHbyA - pEqn.flux();
                        }
                    }

                    U = HbyA - rAU * fvc::grad(p);
                    U.correctBoundaryConditions();
                }

                if (problem->checkWrite(runTime))
                {
                    U.rename("Ured");
                    p.rename("Pred");
                    ITHACAstream::exportSolution(U, name(counter), folder);
                    ITHACAstream::exportSolution(p, name(counter), folder);
                    std::ofstream of(folder + name(counter) + "/" +
                                     runTime.timeName());
                    counter++;
                    problem->nextWrite += problem->writeEvery;
                    U.rename("U");
                    p.rename("p");
                }
            }
        }

        void project(int nModesU, int nModesP)
        {
            ULmodes.resize(0);

            for (int i = 0; i < problem->inletIndex.rows(); i++)
            {
                ULmodes.append(problem->liftfield[i]);
            }

            for (int i = 0; i < nModesU; i++)
            {
                ULmodes.append(problem->Umodes.toPtrList()[i]);
            }

            ULmodes.toEigen();
            volScalarField& P = problem->_p();
            volVectorField& U = problem->_U();
            PtrList<volVectorField> gradP;
            gradP.resize(0);

            for (int i = 0; i < nModesP; i++)
            {
                gradP.append(-fvc::grad(problem->Pmodes[i]));
            }

            redGradP = ULmodes.project(gradP);
        }

        Eigen::MatrixXd setOnlineVelocity(Eigen::MatrixXd vel)
        {
            assert(problem->inletIndex.rows() == vel.rows()
                   && "Imposed boundary conditions dimensions do not match given values matrix dimensions");
            Eigen::MatrixXd vel_scal;
            vel_scal.resize(vel.rows(), vel.cols());

            for (int k = 0; k < problem->inletIndex.rows(); k++)
            {
                label p = problem->inletIndex(k, 0);
                label l = problem->inletIndex(k, 1);
                scalar area = gSum(problem->liftfield[0].mesh().magSf().boundaryField()[p]);
                scalar u_lf = gSum(problem->liftfield[k].mesh().magSf().boundaryField()[p] *
                                   problem->liftfield[k].boundaryField()[p]).component(l) / area;

                for (int i = 0; i < vel.cols(); i++)
                {
                    vel_scal(k, i) = vel(k, i) / u_lf;
                }
            }

            return vel_scal;
        }


};

int main(int argc, char* argv[])
{
    // Construct the tutorial object
    tutorial18 example(argc, argv);
    // Read some parameters from file
    ITHACAparameters para;
    int NmodesUout = para.ITHACAdict->lookupOrDefault<int>("NmodesUout", 15);
    int NmodesPout = para.ITHACAdict->lookupOrDefault<int>("NmodesPout", 15);
    int NmodesUproj = para.ITHACAdict->lookupOrDefault<int>("NmodesUproj", 10);
    int NmodesPproj = para.ITHACAdict->lookupOrDefault<int>("NmodesPproj", 10);
    /// Set the number of parameters
    example.Pnumber = 1;
    /// Set samples
    example.Tnumber = 1;
    /// Set the parameters infos
    example.setParameters();
    // Set the parameter ranges
    example.mu_range(0, 0) = 0.001;
    example.mu_range(0, 1) = 0.001;
    // Generate equispaced samples inside the parameter range
    example.genEquiPar();
    // Set the inlet boundaries where we have non homogeneous boundary conditions
    example.inletIndex.resize(1, 2);
    example.inletIndex(0, 0) = 0;
    example.inletIndex(0, 1) = 0;
    // Time parameters
    example.startTime = 0;
    example.finalTime = 1;
    example.timeStep = 0.001;
    example.writeEvery = 0.01;
    // Perform The Offline Solve;
    example.offlineSolve();
    example.liftfield.append(ITHACAutilities::computeAverage(example.Ufield));
    ITHACAutilities::normalizeFields(example.liftfield);
    example.liftfield[0].rename("Ulift0");
    //ITHACAutilities::changeBCtype(example.liftfield[0],"fixedValue",0);
    ITHACAstream::exportSolution(example.liftfield[0], "0", "./");
    // Homogenize the snapshots
    example.computeLift(example.Ufield, example.liftfield, example.Uomfield);
    // Perform POD on velocity and pressure and store the first 10 modes
    ITHACAPOD::getModes(example.Uomfield, example.Umodes, example.podex, 0, 0,
                        NmodesUout);
    ITHACAPOD::getModes(example.Pfield, example.Pmodes, example.podex, 0, 0,
                        NmodesPout);
    // Create the reduced object
    Eigen::MatrixXd vel(1, 1);
    vel(0, 0) = 1.0;
    tutorial18red reduced(example);
    reduced.vel_now = reduced.setOnlineVelocity(vel);
    reduced.project(10, 10);
    reduced.solveOnlineICO(10, 10, "./ITHACAoutput/Offline/");
    //example.liftSolve();
    exit(0);
}
