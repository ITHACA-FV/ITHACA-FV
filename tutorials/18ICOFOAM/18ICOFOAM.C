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
            :
            unsteadyNS(argc, argv),
            U(_U()),
            p(_p()),
            phi(_phi())
        {}

        /// Velocity field
        volVectorField& U;
        /// Pressure field
        volScalarField& p;
        ///
        surfaceScalarField& phi;



        /// Perform an Offline solve
        void offlineSolve()
        {
            Vector<double> inl(1, 0, 0);
            List<scalar> mu_now(1);

            if (offline)
            {
                ITHACAstream::read_fields(Ufield, U, "./ITHACAoutput/Offline/");
                ITHACAstream::read_fields(Pfield, p, "./ITHACAoutput/Offline/");
                mu_samples =
                    ITHACAstream::readMatrix("./ITHACAoutput/Offline/mu_samples_mat.txt");
            }

            else
            {
                for (label i = 0; i < mu.cols(); i++)
                {
                    //inl[0] = mu(0, i);
                    mu_now[0] = mu(0, i);
                    //assignBC(U, BCind, inl);
                    assignIF(U, inl);
                    change_viscosity( mu(0, i));
                    truthSolve(mu_now);
                }
            }
        }

        void truthSolve(List<scalar> mu_now)
        {
            Time& runTime = _runTime();
            surfaceScalarField& phi = _phi();
            fvMesh& mesh = _mesh();
#include "initContinuityErrs.H"
            fv::options& fvOptions = _fvOptions();
            pisoControl piso(mesh);
            singlePhaseTransportModel& laminarTransport = _laminarTransport();
            volScalarField nu = laminarTransport.nu();

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
                        pEqn.setReference(pRefCell, pRefValue);
                        pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));

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
                    ITHACAstream::exportSolution(U, name(counter), "ITHACAoutput/Offline/");
                    ITHACAstream::exportSolution(p, name(counter), "ITHACAoutput/Offline/");
                    std::ofstream of("ITHACAoutput/Offline/" + name(counter) + "/" +
                                     runTime.timeName());
                    Ufield.append(U);
                    Pfield.append(p);
                    counter++;
                    nextWrite += writeEvery;
                    writeMu(mu_now);
                    // --- Fill in the mu_samples with parameters (time, mu) to be used for the PODI sample points
                    mu_samples.conservativeResize(mu_samples.rows() + 1, mu_now.size() + 1);
                    mu_samples(mu_samples.rows() - 1, 0) = atof(runTime.timeName().c_str());

                    for (int i = 0; i < mu_now.size(); i++)
                    {
                        mu_samples(mu_samples.rows() - 1, i + 1) = mu_now[i];
                    }
                }
                // Resize to Unitary if not initialized by user (i.e. non-parametric problem)
                if (mu.cols() == 0)
                {
                    mu.resize(1, 1);
                }

                if (mu_samples.rows() == counter * mu.cols())
                {
                    ITHACAstream::exportMatrix(mu_samples, "mu_samples", "eigen",
                                               "ITHACAoutput/Offline/");
                }
            }

            Info << "End\n" << endl;
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
    //example.liftSolve();
    exit(0);
}
