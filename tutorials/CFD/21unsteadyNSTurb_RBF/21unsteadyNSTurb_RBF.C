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
#include "UnsteadyNSTurb.H"
#include "ITHACAPOD.H"
#include "ReducedUnsteadyNS.H"
#include "ReducedUnsteadyNSTurb.H"
#include "ITHACAstream.H"
#include <chrono>
#include<math.h>
#include<iomanip>

#include "fvCFD.H"
#include "fvOptions.H"
#include "IOmanip.H"
#include "ITHACAutilities.H"
#include "Foam2Eigen.H"
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <Eigen/SparseLU>
#include "EigenFunctions.H"
#include "reductionProblem.H"
#include "steadyNS.H"

class tutorial21: public UnsteadyNSTurb
{
    public:
        explicit tutorial21(int argc, char* argv[])
            :
            UnsteadyNSTurb(argc, argv),
            U(_U()),
            p(_p()),
            nut(_nut())
        {}

        // Fields To Perform
        volVectorField& U;
        volScalarField& p;
        volScalarField& nut;

        void offlineSolve(std::string offlinepath)
        {
            Vector<double> inl(1, 0, 0);
            List<scalar> mu_now(1);
            label BCind = 0;

            if ((offline) && (ITHACAutilities::check_folder(offlinepath) == true))
            {
                ITHACAstream::read_fields(Ufield, U, offlinepath);
                ITHACAstream::read_fields(Pfield, p, offlinepath);
                ITHACAstream::read_fields(nutFields, nut, offlinepath);
            }
            else
            {
                for (label i = 0; i < mu.cols(); i++)
                {
                    inl[0] = mu(0, i);
                    mu_now[0] = mu(0, i);
                    assignBC(U, BCind, inl);
                    truthSolve(mu_now,offlinepath);
                }
            }
        }
};

/*---------------------------------------------------------------------------*\
                               Starting the MAIN
\*---------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
    // Construct the tutorial21 object
    tutorial21 example(argc, argv);
    // Read parameters from ITHACAdict file
    ITHACAparameters* para = ITHACAparameters::getInstance(example._mesh(),
                             example._runTime());
    int NmodesU = para->ITHACAdict->lookupOrDefault<int>("NmodesU", 10);
    int NmodesP = para->ITHACAdict->lookupOrDefault<int>("NmodesP", 10);
    int NmodesSUP = para->ITHACAdict->lookupOrDefault<int>("NmodesSUP", 10);
    int NmodesNUT = para->ITHACAdict->lookupOrDefault<int>("NmodesNUT", 10);
    int NmodesProject = para->ITHACAdict->lookupOrDefault<int>("NmodesProject", 10);

    /// Set the number of parameters
    example.Pnumber = 1;
    /// Set samples
    example.Tnumber = 2;
    /// Set the parameters infos
    example.setParameters();
    // Set the parameter ranges
    example.mu_range(0, 0) = 1.0;
    example.mu_range(0, 1) = 1.1;
    // Generate equispaced samples inside the parameter range
    example.genEquiPar();
    // Set the inlet boundaries where we have non homogeneous boundary conditions
    example.inletIndex.resize(1, 2);
    example.inletIndex(0, 0) = 0;
    example.inletIndex(0, 1) = 0;
    // Time parameters
    example.startTime = 0;
    example.finalTime = 5;
    example.timeStep = 0.001;
    example.writeEvery = 0.1;
    // Perform The Offline Solve;
    example.offlineSolve("./ITHACAoutput/Offline/");

    // Define velRBF
    example.velRBF = Eigen::MatrixXd(102,1);
    for (int i=0; i<51; i++)
        example.velRBF(i,0) = 1.0;
    
    for (int i=51; i<102; i++)
        example.velRBF(i,0) = 1.1;
    
    // Solve the supremizer problem
    example.solvesupremizer();
    // Search the lift function
    example.liftSolve();
    // Normalize the lifting function
    ITHACAutilities::normalizeFields(example.liftfield);
    // Create homogeneous basis functions for velocity
    example.computeLift(example.Ufield, example.liftfield, example.Uomfield);
    // Export the homogeneous velocity snapshots
    ITHACAstream::exportFields(example.Uomfield, "./ITHACAoutput/Offline",
                               "Uofield");
    // Perform a POD decomposition for the velocity, the pressure and the eddy viscosity
    ITHACAPOD::getModes(example.Uomfield, example.Umodes, example._U().name(),
                        example.podex,
                        example.supex, 0, NmodesProject);
    ITHACAPOD::getModes(example.Pfield, example.Pmodes, example.p().name(),
                        example.podex,
                        example.supex, 0, NmodesProject);
    ITHACAPOD::getModes(example.nutFields, example.nutModes, example._nut().name(),
                        example.podex,
                        example.supex, 0, NmodesProject);
    // Solve the supremizer problem based on the pressure modes
    example.solvesupremizer("modes");
    // Compute the reduced order matrices
    example.projectSUP("./Matrices", NmodesU, NmodesP, NmodesSUP,
                      NmodesNUT,true);

    // Create an object of the turbulent class
    ReducedUnsteadyNSTurb pod_rbf(
        example);

    // Set value of the reduced viscosity and the penalty factor
    pod_rbf.nu = 1e-05;
    pod_rbf.tauU.resize(1, 1);
    pod_rbf.tstart = 0;
    pod_rbf.finalTime = 10.1;
    pod_rbf.dt = 0.005;
    pod_rbf.storeEvery = 0.005;
    pod_rbf.exportEvery = 0.1;

    // We create the matrix rbfCoeff which will store the values of the interpolation results for the eddy viscosity field
    Eigen::MatrixXd rbfCoeff;
    rbfCoeff.resize(NmodesNUT+1, 1);

    // Perform an online solve for the new values of inlet velocities
    for (label k = 0; k < 1; k++)
    {
        Eigen::MatrixXd velNow(1, 1);
        velNow(0, 0) = 1.05;
        pod_rbf.tauU(0, 0) = 0;
        pod_rbf.solveOnlineSUP(velNow);
        rbfCoeff.col(k) = pod_rbf.rbfCoeffMat.col(k);
        Eigen::MatrixXd tmp_sol(pod_rbf.y.rows() + 1, 1);
        tmp_sol(0) = k + 1;
        tmp_sol.col(0).tail(pod_rbf.y.rows()) = pod_rbf.y;
        pod_rbf.online_solution.append(tmp_sol);
    }

    // Save the matrix of interpolated eddy viscosity coefficients
    ITHACAstream::exportMatrix(rbfCoeff, "rbfCoeff", "python",
                               "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(rbfCoeff, "rbfCoeff", "matlab",
                               "./ITHACAoutput/Matrices/");
    // Save the online solution
    ITHACAstream::exportMatrix(pod_rbf.online_solution, "red_coeff", "python",
                               "./ITHACAoutput/red_coeff");
    ITHACAstream::exportMatrix(pod_rbf.online_solution, "red_coeff", "matlab",
                               "./ITHACAoutput/red_coeff");
    ITHACAstream::exportMatrix(pod_rbf.online_solution, "red_coeff", "eigen",
                               "./ITHACAoutput/red_coeff");

    // Reconstruct and export the solution
    pod_rbf.reconstruct(true, "./ITHACAoutput/Reconstruction/");
    
    // Solve the full order problem for the online velocity values for the purpose of comparison
    tutorial21 example2(argc, argv);
    /// Set the number of parameters
    example2.Pnumber = 1;
    /// Set samples
    example2.Tnumber = 2;
    /// Set the parameters infos
    example2.setParameters();
    // Set the parameter ranges
    example2.mu_range(0, 0) = 1.05;
    example2.mu_range(0, 1) = 1.05;
    // Generate equispaced samples inside the parameter range
    example2.genEquiPar();
    // Set the inlet boundaries where we have non homogeneous boundary conditions
    example2.inletIndex.resize(1, 2);
    example2.inletIndex(0, 0) = 0;
    example2.inletIndex(0, 1) = 0;
    // Time parameters
    example2.startTime = 0;
    example2.finalTime = 5;
    example2.timeStep = 0.001;
    example2.writeEvery = 0.1;
    example2.offlineSolve("./ITHACAoutput/Offline_check/");

    Eigen::MatrixXd errFrobU = ITHACAutilities::errorFrobRel(example2.Ufield,
                               pod_rbf.uRecFields);
    Eigen::MatrixXd errFrobP =  ITHACAutilities::errorFrobRel(example2.Pfield,
                                pod_rbf.pRecFields);
    Eigen::MatrixXd errFrobNut =  ITHACAutilities::errorFrobRel(example2.nutFields,
                                  pod_rbf.nutRecFields);
    ITHACAstream::exportMatrix(errFrobU, "errFrobU", "matlab",
                               "./ITHACAoutput/ErrorsFrob/");
    ITHACAstream::exportMatrix(errFrobP, "errFrobP", "matlab",
                               "./ITHACAoutput/ErrorsFrob/");
    ITHACAstream::exportMatrix(errFrobNut, "errFrobNut", "matlab",
                               "./ITHACAoutput/ErrorsFrob/");
    Eigen::MatrixXd errL2U = ITHACAutilities::errorL2Rel(example2.Ufield,
                             pod_rbf.uRecFields);
    Eigen::MatrixXd errL2P =  ITHACAutilities::errorL2Rel(example2.Pfield,
                              pod_rbf.pRecFields);
    Eigen::MatrixXd errL2Nut =  ITHACAutilities::errorL2Rel(example2.nutFields,
                                pod_rbf.nutRecFields);
    ITHACAstream::exportMatrix(errL2U, "errL2U", "matlab",
                               "./ITHACAoutput/ErrorsL2/");
    ITHACAstream::exportMatrix(errL2P, "errL2P", "matlab",
                               "./ITHACAoutput/ErrorsL2/");
    ITHACAstream::exportMatrix(errL2Nut, "errL2Nut", "matlab",
                               "./ITHACAoutput/ErrorsL2/");
    exit(0);
}