/*---------------------------------------------------------------------------*\
Copyright (C) 2017 by the ITHACA-FV authors

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
    Example of NS-Stokes and heat transport equation Reduction Problem
SourceFiles
    06NonIsothermalMixingElbow.C


\*---------------------------------------------------------------------------*/

#include "unsteadyNST.H"
#include "ITHACAPOD.H"
#include "ReducedUnsteadyNS.H"
#include "ReducedUnsteadyNST.H"
#include "ITHACAstream.H"
#include <chrono>
#include<math.h>
#include<iomanip>

class tutorial07: public unsteadyNST
{
    public:
        explicit tutorial07(int argc, char* argv[])
            :
            unsteadyNST(argc, argv),
            U(_U()),
            p(_p()),
            T(_T())
        {}

        // Fields To Perform
        volVectorField& U;
        volScalarField& p;
        volScalarField& T;

        void offlineSolve()
        {
            Vector<double> inl(0, 0, 0);
            List<scalar> mu_now(1);
            Info << "here" << endl;

            if (offline)
            {
                ITHACAstream::read_fields(Ufield, U, "./ITHACAoutput/Offline/");
                ITHACAstream::read_fields(Pfield, p, "./ITHACAoutput/Offline/");
                ITHACAstream::read_fields(Tfield, T, "./ITHACAoutput/Offline/");
                mu_samples =
                    ITHACAstream::readMatrix("./ITHACAoutput/Offline/mu_samples_mat.txt");
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
    // Construct the tutorial object
    tutorial07 example(argc, argv);
    // Read some parameters from file
    ITHACAparameters* para = ITHACAparameters::getInstance(example._mesh(),
                             example._runTime());
    int NmodesUout = para->ITHACAdict->lookupOrDefault<int>("NmodesUout", 5);
    int NmodesPout = para->ITHACAdict->lookupOrDefault<int>("NmodesPout", 5);
    int NmodesTout = para->ITHACAdict->lookupOrDefault<int>("NmodesTout", 5);
    int NmodesSUPout = para->ITHACAdict->lookupOrDefault<int>("NmodesSUPout", 5);
    int NmodesUproj = para->ITHACAdict->lookupOrDefault<int>("NmodesUproj", 5);
    int NmodesPproj = para->ITHACAdict->lookupOrDefault<int>("NmodesPproj", 15);
    int NmodesTproj = para->ITHACAdict->lookupOrDefault<int>("NmodesTproj", 5);
    int NmodesSUPproj = para->ITHACAdict->lookupOrDefault<int>("NmodesSUPproj", 5);
    /// Set the number of parameters
    example.Pnumber = 1;
    /// Set samples
    example.Tnumber = 1;
    /// Set the parameters infos
    example.setParameters();
    // Set the parameter ranges
    example.mu_range(0, 0) = 0.1;
    example.mu_range(0, 1) = 0.1;
    //Generate equispaced samples inside the parameter range
    example.genEquiPar();
    // Set the inlet boundaries where you have non homogeneous boundary conditions
    example.inletIndex.resize(2, 2);
    example.inletIndex << 3, 0, 2, 1;
    example.inletIndexT.resize(3, 1);
    example.inletIndexT << 3, 2, 0;
    // Time parameters
    example.startTime = 0;
    example.finalTime = 50;
    example.timeStep = 0.05;
    example.writeEvery = 0.1;
    // Perform The Offline Solve;
    example.offlineSolve();
    // Solve the supremizer problem
    example.solvesupremizer();
    // Search the lift function for the velocity
    example.liftSolve();
    // Search the lift function for the temperature
    example.liftSolveT();
    // Create homogeneous basis functions for velocity
    example.computeLift(example.Ufield, example.liftfield, example.Uomfield);
    // Create homogeneous basis functions for temperature
    example.computeLiftT(example.Tfield, example.liftfieldT, example.Tomfield);
    // Perform a POD decomposition for velocity temperature and pressure fields
    ITHACAPOD::getModes(example.Uomfield, example.Umodes, example._U().name(),
                        example.podex, 0, 0,
                        NmodesUout);
    ITHACAPOD::getModes(example.Pfield, example.Pmodes, example._p().name(),
                        example.podex, 0, 0,
                        NmodesPout);
    ITHACAPOD::getModes(example.Tomfield, example.Tmodes, example._T().name(),
                        example.podex, 0, 0,
                        NmodesTout);
    ITHACAPOD::getModes(example.supfield, example.supmodes, example._U().name(),
                        example.podex,
                        example.supex, 1, NmodesSUPout);
    example.projectSUP("./Matrices", NmodesUproj, NmodesPproj, NmodesTproj,
                       NmodesSUPproj);
    reducedUnsteadyNST reduced(example);
    // Set values of the ridotto stuff
    reduced.nu = 0.1;
    reduced.tstart = 0;
    reduced.finalTime = 50;
    reduced.dt = 0.05;
    reduced.DT = 1e-06;
    // Set the online velocity
    Eigen::MatrixXd vel_now;
    vel_now.resize(2, 1);
    vel_now << 0.6, 1.2;
    // Set the online temperature and the value of the internal field
    Eigen::MatrixXd temp_now;
    temp_now.resize(3, 1);
    temp_now << 60, 70, 60;
    reduced.solveOnline_sup(vel_now, temp_now);
    // Reconstruct the solution and export it
    reduced.reconstruct_sup("./ITHACAoutput/ReconstructionSUP/", 2);
    reduced.reconstruct_supt("./ITHACAoutput/ReconstructionSUP/", 2);
    exit(0);
}

