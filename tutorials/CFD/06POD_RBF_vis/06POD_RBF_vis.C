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
    Example of NS-Stokes Reduction Problem for Turbulent Flow Case
SourceFiles
    06POD_RBF.C
\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "fvOptions.H"
#include "IOmanip.H"
#include "simpleControl.H"
#include "ITHACAutilities.H"
#include "Foam2Eigen.H"
#include "ITHACAstream.H"
#include "ITHACAPOD.H"
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <Eigen/SparseLU>
#include "EigenFunctions.H"
#include <chrono>
#include "reductionProblem.H"
#include "steadyNS.H"
#include "SteadyNSTurb.H"
#include "ReducedSteadyNS.H"
#include "ReducedSteadyNSTurb.H"
#include "IFstream.H"

/// \brief Class where the tutorial number 6 is implemented.
/// \details It is a child of the laplacianProblem class and some of its
/// functions are overridden to be adapted to the specific case.
class tutorial06 : public SteadyNSTurb
{
    public:
        explicit tutorial06(int argc, char* argv[])
            :
            SteadyNSTurb(argc, argv),
            U(_U()),
            p(_p()),
            nut(_nut())
        {}
        //! [tutorial06]
        // Relevant Fields
        volVectorField& U;
        volScalarField& p;
        volScalarField& nut;
        /// Perform an Offline solve
        void offlineSolve()
        {
            Vector<double> inl(0, 0, 0);
            List<scalar> mu_now(1);
            label BCind = 1;

            // if the offline solution is already performed read the fields
            if (offline)
            {
                ITHACAstream::read_fields(Ufield, U, "./ITHACAoutput/Offline/");
                ITHACAstream::read_fields(Pfield, p, "./ITHACAoutput/Offline/");
                ITHACAstream::read_fields(nutFields, nut, "./ITHACAoutput/Offline/");

                // if the offline stage is not completed then resume it
                if (Ufield.size() < mu.rows())
                {                    
                    for (label i = Ufield.size(); i < mu.rows(); i++)
                    {
                        mu_now[0] = mu(i, 0);
                        change_viscosity(mu(i, 0));
                        counter = Ufield.size() + 1;
                        Info<< "The viscosity of Offline case No. " << counter 
                            << " is changed to " << mu(i, 0) << endl;
                        truthSolve("./ITHACAoutput/Offline/");
                        restart();
                    }
                }
            }
            else
            {
                for (label i = 0; i < mu.rows(); i++)
                {
                    mu_now[0] = mu(i, 0);
                    change_viscosity(mu(i, 0));
                    Info<< "The viscosity of Offline case No. " << counter 
                        << " is changed to " << mu(i, 0) << endl;
                    truthSolve("./ITHACAoutput/Offline/");
                    restart();
                }
            }
        }

        /// Perform an Offline solve for a special set of parameter samples called par
        void offlineSolve(Eigen::MatrixXd par, fileName folder)
        {
            Vector<double> inl(0, 0, 0);
            List<scalar> mu_now(1);

            for (label i = 0; i < par.rows(); i++)
            {
                mu_now[0] = mu(i, 0);
                change_viscosity(mu(i, 0));
                Info<< "The viscosity of Offline case No. " << counter 
                    << " is changed to " << mu(i, 0) << endl;
                truthSolve(folder);
                restart();
            }
        }

        /// Perform an Online solve for a special set of parameter samples called par
        void onlineSolve(fileName folder)
        {
            Vector<double> inl(0, 0, 0);
            List<scalar> mu_now(1);
            label BCind = 1;

            fileName filePath = folder + "1/U";
            IFstream exFileOff(filePath);
            filePath = folder + "processor0/1/U";
            IFstream exFileOff2(filePath);
            
            // if the offline solution is already performed read the fields
            if (exFileOff.good() || exFileOff2.good())
            {
                ITHACAstream::read_fields(Ufield, U, folder);
                ITHACAstream::read_fields(Pfield, p, folder);
                ITHACAstream::read_fields(nutFields, nut, folder);

                // if the offline stage is not completed then resume it
                if (Ufield.size() < mu.rows())
                {                    
                    for (label i = Ufield.size(); i < mu.rows(); i++)
                    {
                        mu_now[0] = mu(i, 0);
                        change_viscosity(mu(i, 0));
                        counter = Ufield.size() + 1;
                        Info<< "The viscosity of Online case No. " << counter 
                            << " is changed to " << mu(i, 0) << endl;
                        truthSolve(folder);
                        restart();
                    }
                }
            }
            else
            {
                for (label i = 0; i < mu.rows(); i++)
                {
                    mu_now[0] = mu(i, 0);
                    change_viscosity(mu(i, 0));                    
                    Info<< "The viscosity of Online case No. " << counter 
                        << " is changed to " << mu(i, 0) << endl;
                    truthSolve(folder);
                    restart();
                }
            }
        }

        void truthSolve(fileName folder)
        {
            Time& runTime = _runTime();
            fvMesh& mesh = _mesh();
            volScalarField& p = _p();
            volVectorField& U = _U();
            surfaceScalarField& phi = _phi();
            fv::options& fvOptions = _fvOptions();
            simpleControl& simple = _simple();
            IOMRFZoneList& MRF = _MRF();
            singlePhaseTransportModel& laminarTransport = _laminarTransport();
#include "NLsolveSteadyNSTurb.H"
            ITHACAstream::exportSolution(U, name(counter), folder);
            ITHACAstream::exportSolution(p, name(counter), folder);
            volScalarField _nut(turbulence->nut());
            ITHACAstream::exportSolution(_nut, name(counter), folder);
            Ufield.append((U).clone());
            Pfield.append((p).clone());
            nutFields.append((_nut).clone());
            counter++;
        }
};

int main(int argc, char* argv[])
{
    // Construct the tutorial06 object
    tutorial06 example(argc, argv);
    // Read parameters samples for the offline stage and the ones for the online stages
    word par_offline("./par_offline");
    word par_new("./par_online");
    Eigen::MatrixXd par_online = ITHACAstream::readMatrix(par_new);
    example.mu = ITHACAstream::readMatrix(par_offline);

    // Set the inlet boundaries where we have parameterized boundary conditions
    example.inletIndex.resize(1, 2);
    example.inletIndex(0, 0) = 2;
    example.inletIndex(0, 1) = 2;

    // create a list to store the time of different steps
    List<scalar> timeList;
    // The list for name of the steps
    List<word> nameList;

    ITHACAparameters* para = ITHACAparameters::getInstance(example._mesh(),
                             example._runTime());
    // Read parameters from ITHACAdict file
    int NmodesU = para->ITHACAdict->lookupOrDefault<int>("NmodesU", 5);
    int NmodesP = para->ITHACAdict->lookupOrDefault<int>("NmodesP", 5);
    int NmodesSUP = para->ITHACAdict->lookupOrDefault<int>("NmodesSUP", 5);
    int NmodesNUT = para->ITHACAdict->lookupOrDefault<int>("NmodesNUT", 5);
    int NmodesProject = para->ITHACAdict->lookupOrDefault<int>("NmodesProject", 5);
    word stabilization = para->ITHACAdict->lookupOrDefault<word>("Stabilization",
                         "supremizer");
    bool supremizerConsistent = para->ITHACAdict->lookupOrDefault<bool>("supremizerConsistent",
                        false);
    bool exportrecField = para->ITHACAdict->lookupOrDefault<bool>("exportrecField",
                            false);
    bool exportErrorField = para->ITHACAdict->lookupOrDefault<bool>("exportErrorField",
                            false);
    // Perform The Offline Solve;
    example.offlineSolve();
    // Read the lift functions
    ITHACAstream::read_fields(example.liftfield, example.U, "./lift/");
    timeList.append(example._runTime().elapsedCpuTime());
    nameList.append("OfflineSolve");

    // Create the homogeneous set of snapshots for the velocity field
    // ITHACAutilities::normalizeFields(example.liftfield);
    // Homogenize the snapshots
    if (example.nonUniformbc)
    {
        // The shape of mu is different from the old code, so transpose it
        Eigen::MatrixXd muTranspose (1, example.mu.rows());
        // set all value of mu to 1
        muTranspose.setConstant(1.0);
        example.computeLift(example.Ufield, example.liftfield, example.Uomfield, muTranspose);
    }
    else
    {
        example.computeLift(example.Ufield, example.liftfield, example.Uomfield);
    }

    IFstream exUomfield1(fileName ("./ITHACAoutput/Offline/1/Uofield"));
    IFstream exUomfield2(fileName ("./ITHACAoutput/Offline/processor0/1/Uofield"));
    if (!exUomfield1.good() && !exUomfield2.good())
    {
        // Export the homogeneous velocity snapshots
        ITHACAstream::exportFields(example.Uomfield, "./ITHACAoutput/Offline",
            "Uofield");
    }
    timeList.append(example._runTime().elapsedCpuTime());
    nameList.append("HomogenizeUofield");
    
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
    timeList.append(example._runTime().elapsedCpuTime());
    nameList.append("POD");

    example.Ufield.clear();
    example.Phifield.clear();
    example.Uomfield.clear();

    // Solve the supremizer problem based on the pressure modes
    if (stabilization == "supremizer")
    {
        if (supremizerConsistent)
        {
            example.solvesupremizerConsistent("modes");
        }
        else
        {
            example.solvesupremizer("modes");
        }        
    }
    timeList.append(example._runTime().elapsedCpuTime());
    nameList.append("SolveSupremizer");

    example.Pfield.clear();

    // Compute the reduced order matrices
    // Get reduced matrices
    if (stabilization == "supremizer")
    {
        example.projectSUP("./Matrices", NmodesU, NmodesP, NmodesSUP,
                           NmodesNUT);
    }
    else if (stabilization == "PPE")
    {
        example.projectPPE("./Matrices", NmodesU, NmodesP, NmodesSUP,
                           NmodesNUT);
    }
    timeList.append(example._runTime().elapsedCpuTime());
    nameList.append("Project");

    // Create an object of the turbulent class
    ReducedSteadyNSTurb pod_rbf(
        example);
    // We create the matrix rbfCoeff which will store the values of the interpolation results for the eddy viscosity field
    Eigen::MatrixXd rbfCoeff;
    rbfCoeff.resize(NmodesNUT, par_online.rows());
    pod_rbf.onlineMu = par_online;

    Eigen::MatrixXd velNow(1, 1);
    velNow(0, 0) = 1.0;

    // Perform an online solve for the new values of inlet velocities
    for (label k = 0; k < par_online.rows(); k++)
    {
        pod_rbf.tauU.resize(2, 1);
        pod_rbf.tauU(0, 0) = 0;
        pod_rbf.tauU(1, 0) = 0;

        // Set value of the reduced viscosity and the penalty factor
        pod_rbf.nu = par_online(k, 0);

        if (stabilization == "supremizer")
        {
            pod_rbf.solveOnlineSUP(velNow);
        }
        else if (stabilization == "PPE")
        {
            pod_rbf.solveOnlinePPE(velNow);
        }

        rbfCoeff.col(k) = pod_rbf.rbfCoeff;
        Eigen::MatrixXd tmp_sol(pod_rbf.y.rows() + 1, 1);
        tmp_sol(0) = k + 1;
        tmp_sol.col(0).tail(pod_rbf.y.rows()) = pod_rbf.y;
        pod_rbf.online_solution.append(tmp_sol);

    }
    timeList.append(example._runTime().elapsedCpuTime());
    nameList.append("OnlineSolve");

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
    pod_rbf.rbfCoeffMat = rbfCoeff;
    // Reconstruct and export the solution
    pod_rbf.reconstruct(exportrecField, "./ITHACAoutput/Online/");
    timeList.append(example._runTime().elapsedCpuTime());
    nameList.append("Reconstruct");

    example.Umodes.clear();
    example.Pmodes.clear();
    example.nutModes.clear();
    example.supmodes.clear();
    example.liftfield.clear();
    example.nutFields.clear();
    example.L_U_SUPmodes.clear();
    example.L_PHImodes.clear();

    example.mu = par_online;
    example.onlineSolve("./ITHACAoutput/Online/");

    // Write error of online solutions
    Eigen::MatrixXd errL2U = ITHACAutilities::errorL2Rel(example.Ufield,
                                                            pod_rbf.uRecFields);
    Eigen::MatrixXd errL2P = ITHACAutilities::errorL2Rel(example.Pfield,
                                                            pod_rbf.pRecFields);
    ITHACAstream::exportMatrix(errL2U, "errL2U", "python",
                                "./ITHACAoutput/Online/ErrorsL2/");
    ITHACAstream::exportMatrix(errL2P, "errL2P", "python",
                                "./ITHACAoutput/Online/ErrorsL2/");

    if (exportErrorField)
    {
        // Export errorfields
        for (label k = 0; k < example.mu.rows(); k++)
        {
            volVectorField Uerror("Uerror", example.Ufield[k] - pod_rbf.uRecFields[k]);
            volScalarField perror("perror", example.Pfield[k] - pod_rbf.pRecFields[k]);
            ITHACAstream::exportSolution(Uerror,
                                        name(k+1),
                                        "./ITHACAoutput/Online/");
            ITHACAstream::exportSolution(perror,
                                        name(k+1),
                                        "./ITHACAoutput/Online/");
        }
    }
    timeList.append(example._runTime().elapsedCpuTime());
    nameList.append("done");

    List<scalar> globalTimeList(timeList.size(), 0.0);
    forAll(timeList, i)
    {
        scalar localVal = timeList[i];
        reduce(localVal, maxOp<scalar>());
        globalTimeList[i] = localVal;
    }

    Info<< "The elapsed time for the different steps is:\n"
        << "-----------------------------------------------------\n";
    for (label i = 0; i < globalTimeList.size(); i++)
    {
        Info<< nameList[i] << " = " << globalTimeList[i] << endl;
    }
    Info<< "-----------------------------------------------------\n";

    exit(0);
}

//--------
/// \dir 06POD_RBF Folder of the turorial 6
/// \file
/// \brief Implementation of a tutorial of turbulent blabla

/// \example 06POD_RBF.C
/// \section intro_podRBF Introduction to turorial 6
/// The problems consists of steady Navier-Stokes problem with a parameterized velocity at the inlet.
/// The physical problem is the pitz-daily depicted in the following image
/// \image html pitzdaily.png
/// The velocity at the inlet is parameterized in both x and y directions.
/// In other words, the parameters in this setting are the magnitude of the velocity at the inlet and
/// the inclination of the velocity with respect to the inlet.
/// \section CODE A detailed look into the code
///
/// This section explains the main steps necessary to construct the tutorial N°6.
///
/// \subsection header The necessary header files
///
/// First of all let's have a look into the header files which have to be included, indicating what they are responsible for:
///
/// \dontinclude 06POD_RBF.C
/// The OpenFOAM header files:
/// \skip fvCFD
/// \until simpleControl.H
///
/// ITHACA-FV header files, <ITHACAstream.H> is responsible for reading and exporting the fields and other sorts of data.
/// <ITHACAutilities.H> has the utilities which compute the mass matrices, fields norms, fields error,...etc.
/// <Foam2Eigen.H> is for converting fields and other data from OpenFOAM to Eigen format.
/// <ITHACAPOD.H> is for the computation of the POD modes.
/// \until #include "ITHACAPOD.H"
/// The Eigen library for matrix manipulation and linear and non-linear algebra operations:
/// \until EigenFunctions.H
/// chrono to compute the speedup
/// \until chrono
/// The header files of ITHACA-FV necessary for this tutorial are: <reductionProblem.H> A general class for
/// the implementation of a full order parameterized problem. <steadyNS.H> is
/// for the full order steady NS problem , <SteadyNSTurb.H> is the child of <steadyNS.H>
/// and it is the class for the full order steady NS turbulent problem.
/// Finally <reducedSteadyNS.H> and <ReducedSteadyNSTurb.H> are for the construction of the reduced order problems
/// \until ReducedSteadyNSTurb
///
/// \subsection classtutorial06 Definition of the tutorial06 class
/// We define the tutorial06 class as a child of the SteadyNSTurb class.
/// The constructor is defined with members that are the fields required to be manipulated
/// during the resolution of the full order problem using simpleFoam. Such fields are
/// also initialized with the same initial conditions in the solver.
/// \skipline tutorial06
/// \until {}
///
/// Inside the tutorial06 class we define the offlineSolve method according to the
/// specific parameterized problem that needs to be solved. If the offline solve has
/// been previously performed then the method just reads the existing snapshots from the Offline directory,
/// and if the offline solve has been started but not completed then it continues the offline stage
/// from the last snapshot computed. If the procedure has not started at all, the method
/// loops over all the parameters samples, changes the inlet velocity components at the inlet
/// with the iterable parameter sample for both components of the velocity
/// then it performs the offline solve.
/// \skipline offlineSolve
/// \until }
/// \skipline else
/// \until }
/// \skipline }
/// \skipline }
///
/// This offlineSolve method is just designed to compute the full order solutions
/// for a set of parameters samples called par, the goal is to compute the full order solution for
/// a cross validation test samples which will be used to test the reduced order model
/// in the online stage.
///
/// \skipline offlineSolve
/// \until };
/// \subsection main Definition of the main function
///
/// In this section we address the definition of the main function.
/// First we construct the object "example" of type tutorial06:
///
/// \skipline example
///
///
/// Then we read the parameters samples for the offline stage and the ones for the online stages
/// using method readMatrix in the ITHACAstream class
///
/// \skipline word
/// \until (par_offline);
///
/// Then we set the inlet boundaries where we have parameterized boundary conditions
/// and we define in which directions we have the parameterization.
/// \skipline example.inletIndex
/// \until  example.inletIndex(1, 1) = 1;
///
///
/// Then we parse the ITHACAdict file to determine the number of modes to be written out
/// and also the ones to be used for the projection of the velocity, pressure, supremizer and the eddy viscosity:
///
/// \skipline ITHACAparameters
/// \until NmodesProject
///
/// Now we are ready to perform the offline stage:
/// \skipline   example.offlineSolve();
///
/// The next step is to read the lifting functions
/// \skipline example.liftfield
///
/// Then we compute the homogeneous velocity field snapshots, and then we output them
/// \skipline Uomfield
/// \until ITHACAstream::exportFields(example.Uomfield
///
/// After that, the modes for velocity, pressure and the eddy viscosity are obtained:
/// \skipline ITHACAPOD
/// \until nutModes
///
/// Then we compute the supremizer modes on the basis of the POD pressure modes obtained from the last step
///
/// \skipline solvesup
///
/// then we perform the projection onto the POD modes
/// \skipline projectSUP
///
/// Now we proceed to the ROM part of the tutorial, at first we construct an object of the class <ReducedSteadyNSTurb.H>.
///
/// \skipline reduced
///
/// Now we set the value of the reduced viscosity and we initialize the penalty factor which will be zero
///
/// \skipline nu
/// \until tauU
///
/// We create the matrix rbfCoeff in order to store the values of the interpolated eddy viscosity
/// coefficient using the RBF in the online stage
///
/// \skipline rbfCoeff
/// \until rbfCoeff
///
/// Now we solve the online reduced system for the parameter values stored in par_online which is a
/// different set of samples for the parameters than the one used in the offline stage.
/// This represent a cross validation test for assessing the reduced order model in a better way.
///
/// \skipline for
/// \until }
///
/// Now we output the matrix rbfCoeff and the online solution
///
/// \skipline rbfCoeff
/// \until eigen
///
/// The last step is to reconstruct the velocity and pressure fields using the reduced solution obtained in
/// the online stage and with the POD modes computed earlier on
///
/// \skipline reconstruct
/// \until nuRec
///
/// After we carried out the online stage using an object of the turbulent class <ReducedSteadyNSTurb.H>, now we will repeat
/// same procedure but for the class that does not take into consideration turbulence at the reduced order level which is
/// <reducedSteadyNS.H>. This allows us at the end of comparing the reduced results obtained from both classes.
///
///
/// \skipline pod_normal
/// \until pod_normal.reconstruct_sup
///
/// We finally solve the offline stage for the checking parameters set that was used for validating the reduced order
/// model
/// \skipline offlineSolve
/// \until createSymLink
///




