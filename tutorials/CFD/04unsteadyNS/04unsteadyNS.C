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
    Example of an unsteady NS Reduction Problem
SourceFiles
    04unsteadyNS.C
\*---------------------------------------------------------------------------*/

#include "unsteadyNS.H"
#include "ITHACAPOD.H"
#include "ReducedUnsteadyNS.H"
#include "ITHACAstream.H"
#include <chrono>
#include <math.h>
#include <iomanip>
#include <iostream>

class tutorial04 : public unsteadyNS
{
    public:
        explicit tutorial04(int argc, char* argv[])
            : unsteadyNS(argc, argv), U(_U()), p(_p()) {}

        // Fields To Perform
        volVectorField& U;
        volScalarField& p;

        void offlineSolve()
        {
            Vector<double> inl(1, 0, 0);
            List<scalar> mu_now(1);
            label BCind = 0;

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
                    assignIF(U, inl);
                    change_viscosity(mu(0, i));
                    truthSolve(mu_now);
                }
            }
        }
};

/*---------------------------------------------------------------------------*\
                               Starting the MAIN
\*---------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
    // Construct the tutorial04 object
    tutorial04 example(argc, argv);
    // Read parameters from ITHACAdict file
    ITHACAparameters* para = ITHACAparameters::getInstance(example._mesh(),
                             example._runTime());
    int NmodesUout = para->ITHACAdict->lookupOrDefault<int>("NmodesUout", 15);
    int NmodesPout = para->ITHACAdict->lookupOrDefault<int>("NmodesPout", 15);
    int NmodesSUPout = para->ITHACAdict->lookupOrDefault<int>("NmodesSUPout", 15);
    int NmodesUproj = para->ITHACAdict->lookupOrDefault<int>("NmodesUproj", 10);
    int NmodesPproj = para->ITHACAdict->lookupOrDefault<int>("NmodesPproj", 10);
    int NmodesSUPproj = para->ITHACAdict->lookupOrDefault<int>("NmodesSUPproj", 10);
    // Set the number of parameters
    example.Pnumber = 1;
    // Set samples
    example.Tnumber = 1;
    // Set the parameters infos
    example.setParameters();
    // Set the parameter ranges
    example.mu_range(0, 0) = 0.005;
    example.mu_range(0, 1) = 0.005;
    // Generate equispaced samples inside the parameter range
    example.genEquiPar();
    // Set the inlet boundaries where we have non homogeneous boundary conditions
    example.inletIndex.resize(1, 2);
    example.inletIndex(0, 0) = 0;
    example.inletIndex(0, 1) = 0;
    // Time parameters
    example.startTime = 60;
    example.finalTime = 70;
    example.timeStep = 0.01;
    example.writeEvery = 0.1;
    // Perform The Offline Solve
    example.offlineSolve();

    // Check if lift or penalty method should be used
    if (example.bcMethod == "lift")
    {
        // Search the lift function
        example.liftSolve();
        // Normalize the lifting function
        ITHACAutilities::normalizeFields(example.liftfield);
        // Create homogeneous basis functions for velocity
        example.computeLift(example.Ufield, example.liftfield, example.Uomfield);
        // Perform a POD decomposition for velocity and pressure
        ITHACAPOD::getModes(example.Uomfield, example.Umodes, example._U().name(),
                            example.podex, 0, 0, NmodesUout);
        ITHACAPOD::getModes(example.Pfield, example.Pmodes, example._p().name(),
                            example.podex, 0, 0, NmodesPout);
    }
    else
    {
        // Perform a POD decomposition for velocity and pressure without lifting function
        ITHACAPOD::getModes(example.Ufield, example.Umodes, example._U().name(),
                            example.podex, 0, 0, NmodesUout);
        ITHACAPOD::getModes(example.Pfield, example.Pmodes, example._p().name(),
                            example.podex, 0, 0, NmodesPout);
    }

    // Check the method and perform the appropriate projection
    if (example.method == "supremizer")
    {
        example.solvesupremizer();
        ITHACAPOD::getModes(example.supfield, example.supmodes, example._U().name(),
                            example.podex, example.supex, 1, NmodesSUPout);
        example.projectSUP("./Matrices", NmodesUproj, NmodesPproj, NmodesSUPproj);
    }
    else if (example.method == "PPE")
    {
        example.projectPPE("./Matrices", NmodesUproj, NmodesPproj);
    }

    reducedUnsteadyNS reduced(example);
    // Set values of the reduced stuff
    reduced.nu = 0.005;
    reduced.tstart = 60;
    reduced.finalTime = 70;
    reduced.dt = 0.005;
    reduced.storeEvery = 0.005;
    reduced.exportEvery = 0.1;
    // Set the online velocity
    Eigen::MatrixXd vel_now(1, 1);
    vel_now(0, 0) = 1;

    // If using the penalty method, set tauU
    if (example.bcMethod == "penalty")
    {
        reduced.tauU = Eigen::MatrixXd::Zero(1, 1);
        reduced.tauU(0, 0) = 1e-1;  // Example penalty coefficient
    }

    // Perform the online solve based on the method
    if (example.method == "supremizer")
    {
        reduced.solveOnline_sup(vel_now, 1);
    }
    else if (example.method == "PPE")
    {
        reduced.solveOnline_PPE(vel_now, 1);
    }

    // Reconstruct the solution and export it
    reduced.reconstruct(true,
                        "./ITHACAoutput/Reconstruction" + example.method + "/");
    return 0;
}

