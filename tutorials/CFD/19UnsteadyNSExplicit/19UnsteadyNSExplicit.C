/*---------------------------------------------------------------------------*\
     ██╗████████╗██╗  ██╗ █████╗  ██████╗ █████╗       ███████╗██╗   ██╗
     ██║╚══██╔══╝██║  ██║██╔══██╗██╔════╝██╔══██╗      ██╔════╝██║   ██║
     ██║   ██║   ███████║███████║██║     ███████║█████╗█████╗  ██║   ██║
     ██║   ██║   ██╔══██║██╔══██║██║     ██╔══██║╚════╝██╔══╝  ╚██╗ ██╔╝
     ██║   ██║   ██║  ██║██║  ██║╚██████╗██║  ██║      ██║      ╚████╔╝
     ╚═╝   ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝      ╚═╝       ╚═══╝

 * In real Time Highly Advanced Computational Applications for Finite Volumes
 * Copyright (C) 2020 by the ITHACA-FV authors
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
    Example of an unsteady NS Reduction Problem with the Discretize-Then-
    Project approach for the lid driven cavity problem
SourceFiles
    19UnsteadyNSExplicit.C
\*---------------------------------------------------------------------------*/

#include "UnsteadyNSExplicit.H"
#include "ITHACAPOD.H"
#include "ReducedUnsteadyNSExplicit.H"
#include "ITHACAstream.H"
#include <chrono>
#include<math.h>
#include<iomanip>

class tutorial19: public UnsteadyNSExplicit
{
    public:
        explicit tutorial19(int argc, char* argv[])
            :
            UnsteadyNSExplicit(argc, argv),
            U(_U()),
            p(_p())
        {}

        // Fields To Perform
        volVectorField& U;
        volScalarField& p;

        void offlineSolve()
        {
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
    // Construct the tutorial19 object
    tutorial19 example(argc, argv);
    // Read parameters from ITHACAdict file
    ITHACAparameters* para = ITHACAparameters::getInstance(example._mesh(),
                             example._runTime());
    int NmodesUout = para->ITHACAdict->lookupOrDefault<int>("NmodesUout", 15);
    int NmodesPout = para->ITHACAdict->lookupOrDefault<int>("NmodesPout", 15);
    int NmodesSUPout = 0;
    int NmodesUproj = para->ITHACAdict->lookupOrDefault<int>("NmodesUproj", 10);
    int NmodesPproj = para->ITHACAdict->lookupOrDefault<int>("NmodesPproj", 10);
    int NmodesSUPproj = 0;
    /// Set the number of parameters
    example.Pnumber = 1;
    /// Set samples
    example.Tnumber = 1;
    /// Set the parameters infos
    example.setParameters();
    // Set the parameter ranges
    example.mu_range(0, 0) = 0.01;
    example.mu_range(0, 1) = 0.01;
    // Generate equispaced samples inside the parameter range
    example.genEquiPar();
    // Set the inlet boundaries where we have non homogeneous boundary conditions
    example.inletIndex.resize(1, 2);
    example.inletIndex(0, 0) = 0;
    example.inletIndex(0, 1) = 0;
    // Time parameters
    example.startTime = 0.0;
    example.finalTime = 1.0;
    example.timeStep = 0.005;
    example.writeEvery = 0.005;
    // Perform The Offline Solve;
    example.offlineSolve();
    // Perform a POD decomposition for velocity and pressure
    ITHACAPOD::getModes(example.Ufield, example.Umodes, example._U().name(),
                        example.podex, 0, 0,
                        NmodesUout);
    ITHACAPOD::getModes(example.Pfield, example.Pmodes, example._p().name(),
                        example.podex, 0, 0,
                        NmodesPout);
    // Galerkin Projection
    example.discretizeThenProject("./Matrices", NmodesUproj, NmodesPproj,
                                  NmodesSUPproj);
    ReducedUnsteadyNSExplicit reduced(example);
    // Set values of the reduced order model
    reduced.nu = 0.01;
    reduced.tstart = 0.0;
    reduced.finalTime = 1.0;
    reduced.dt = 0.005;
    reduced.storeEvery = 0.005;
    reduced.exportEvery = 0.005;
    // Set the online velocity
    Eigen::MatrixXd vel_now(1, 1);
    vel_now(0, 0) = 1;
    reduced.solveOnline(vel_now, 1);
    // Reconstruct the solution and export it
    reduced.reconstruct(false, "./ITHACAoutput/Reconstruction/");
    exit(0);
}




