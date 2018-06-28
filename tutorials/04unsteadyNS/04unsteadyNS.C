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
#include "reducedUnsteadyNS.H"
#include "ITHACAstream.H"
#include <chrono>
#include<math.h>
#include<iomanip>

class tutorial04: public unsteadyNS
{
public:
    explicit tutorial04(int argc, char *argv[])
        :
        unsteadyNS(argc, argv),
        U(_U()),
        p(_p())
    {}

    // Fields To Perform
    volVectorField& U;
    volScalarField& p;

        void offlineSolve()
    {
        Vector<double> inl(0, 0, 0);
        label BCind = 1;
        List<scalar> mu_now(1);
        Info << "here" << endl;
        if (offline)
        {
            ITHACAstream::read_fields(Ufield, U, "./ITHACAoutput/Offline/");
            ITHACAstream::read_fields(Pfield, p, "./ITHACAoutput/Offline/");
        }
        else
        {
            for (label i = 0; i < mu.cols(); i++)
            {
                inl[0] = mu(0, i);
                mu_now[0] = mu(0, i);
                //assignBC(U, BCind, inl);
                assignIF(U, inl);
                change_viscosity( mu(0, i));
                truthSolve(mu_now);
            }
        }
    }
};

/*---------------------------------------------------------------------------*\
                               Starting the MAIN
\*---------------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
    // Construct the tutorial object
    tutorial04 example(argc, argv);

    // Read some parameters from file
    ITHACAparameters para;
    int NmodesUout = para.ITHACAdict->lookupOrDefault<int>("NmodesUout", 15);
    int NmodesPout = para.ITHACAdict->lookupOrDefault<int>("NmodesPout", 15);
    int NmodesSUPout = para.ITHACAdict->lookupOrDefault<int>("NmodesSUPout", 15);
    int NmodesUproj = para.ITHACAdict->lookupOrDefault<int>("NmodesUproj", 10);
    int NmodesPproj = para.ITHACAdict->lookupOrDefault<int>("NmodesPproj", 10);
    int NmodesSUPproj = para.ITHACAdict->lookupOrDefault<int>("NmodesSUPproj", 10);


    /// Set the number of parameters
    example.Pnumber = 1;
    /// Set samples
    example.Tnumber = 5;
    /// Set the parameters infos
    example.setParameters();
    // Set the parameter ranges
    example.mu_range(0, 0) = 0.005;
    example.mu_range(0, 1) = 0.01;
    // Generate equispaced samples inside the parameter range
    example.genEquiPar();

    // Set the inlet boundaries where you have non homogeneous boundary conditions
    example.inletIndex.resize(1, 2);
    example.inletIndex(0, 0) = 0;
    example.inletIndex(0, 1) = 0;

    // Time parameters
    example.startTime = 180;
    example.finalTime = 190;
    example.timeStep = 0.01;
    example.writeEvery = 0.05;

    // Perform The Offline Solve;
    example.offlineSolve();

    // Solve the supremizer problem
    example.solvesupremizer();

    // Search the lift function
    example.liftSolve();

    // Create homogeneous basis functions for velocity
    example.computeLift(example.Ufield, example.liftfield, example.Uomfield);
    
    // Perform a POD decomposition for velocity and pressure
    ITHACAPOD::getModes(example.Uomfield, example.Umodes, example.podex, 0, 0, NmodesUout);
    ITHACAPOD::getModes(example.Pfield, example.Pmodes, example.podex, 0, 0, NmodesPout);
    ITHACAPOD::getModes(example.supfield, example.supmodes, example.podex, example.supex, 1, NmodesSUPout);

    example.projectSUP("./Matrices", NmodesUproj, NmodesPproj, NmodesSUPproj);
    
    reducedUnsteadyNS ridotto(example, "SUP");
    //unsteadyNSreduced ridotto(example, "PPE");

    // Set values of the ridotto stuff
    ridotto.nu = 0.005;
    ridotto.tstart = 0;
    ridotto.finalTime = 10;
    ridotto.dt = 0.01;

    // Set the online velocity
    Eigen::MatrixXd vel_now(1, 1);
    vel_now(0, 0) = 1;
    ridotto.solveOnline_sup(vel_now);
    // Reconstruct the solution and export it
    ridotto.reconstruct_sup(example, "./ITHACAoutput/ReconstructionSUP/", 5);
    //ridotto.reconstruct_PPE(example,"./ITHACAoutput/Reconstruction/",4);
    exit(0);
}


