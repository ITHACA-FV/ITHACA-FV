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

#include "IOmanip.H"
#include "ITHACAPOD.H"
#include "ITHACAstream.H"
#include "ReducedSteadyNS.H"
#include "forces.H"
#include "steadyNS.H"

class tutorial03 : public steadyNS
{
    public:
        /// Constructor
        explicit tutorial03(int argc, char* argv[])
            : steadyNS(argc, argv), U(_U()), p(_p()), args(_args()) {}

        /// Velocity field
        volVectorField& U;
        /// Pressure field
        volScalarField& p;
        /// Arg List
        argList& args;
        /// Inlet BC lable
        label inlet_idx = 0;

        /// Perform an Offline solve
        void offlineSolve()
        {
            List<scalar> mu_now(mu.rows());
            vector vel_new_ini;
            List<vector> vel_new_list;

            // if the offline solution is already performed read the fields
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
                    for(label k=0; k < mu.rows(); k++)
                    {
                        mu_now[k] = mu(k, i);
                    }
                    
                    vel_new_ini = vector(mu(0, i), 0, 0);
                    Info << "The velocity coeffs of case " << i << ": " << mu_now << nl;
                    List<scalar> inlet_yAxis = U.mesh().Cf().boundaryField()[inlet_idx].component(1);
                    vel_new_list.resize(inlet_yAxis.size());
                    scalar midInlet = gAverage(inlet_yAxis);
                    scalar minInlet = gMin(inlet_yAxis);
                    scalar maxInlet = gMax(inlet_yAxis);
                    forAll(vel_new_list, faceI)
                    {
                        scalar Ux_tmp = mu(0, i) 
                                      + mu(1, i)*inlet_yAxis[faceI]
                                      + mu(2, i)*pow(midInlet - inlet_yAxis[faceI], 2)
                                      + mu(3, i)*pow(minInlet - inlet_yAxis[faceI], 2)
                                      + mu(4, i)*pow(maxInlet - inlet_yAxis[faceI], 2);
                        vel_new_list[faceI] = vector(Ux_tmp, 0, 0);
                    }
                    ITHACAutilities::assignBC(_U(), inlet_idx, vel_new_list);
                    assignIF(_U(), vel_new_ini);
                    truthSolve(mu_now);
                    restart();
                }
            }
        }

        /// Perform truthSolve for online parameters
        void onlineSolve()
        {
            List<scalar> mu_now(mu.rows());
            vector vel_new_ini;
            List<vector> vel_new_list;

            std::ifstream exFileOff("./ITHACAoutput/Online/");

            if (exFileOff)
            {
                ITHACAstream::read_fields(Ufield, U, "./ITHACAoutput/Online/");
                ITHACAstream::read_fields(Pfield, p, "./ITHACAoutput/Online/");
            }
            else
            {
                for (label i = 0; i < mu.cols(); i++)
                {
                    for(label k=0; k < mu.rows(); k++)
                    {
                        mu_now[k] = mu(k, i);
                    }
                    
                    vel_new_ini = vector(mu(0, i), 0, 0);
                    Info << "The velocity coeffs of case " << i << ": " << mu_now << nl;
                    List<scalar> inlet_yAxis = U.mesh().Cf().boundaryField()[inlet_idx].component(1);
                    vel_new_list.resize(inlet_yAxis.size());
                    scalar midInlet = gAverage(inlet_yAxis);
                    scalar minInlet = gMin(inlet_yAxis);
                    scalar maxInlet = gMax(inlet_yAxis);
                    forAll(vel_new_list, faceI)
                    {
                        scalar Ux_tmp = mu(0, i) 
                                      + mu(1, i)*inlet_yAxis[faceI]
                                      + mu(2, i)*pow(midInlet - inlet_yAxis[faceI], 2)
                                      + mu(3, i)*pow(minInlet - inlet_yAxis[faceI], 2)
                                      + mu(4, i)*pow(maxInlet - inlet_yAxis[faceI], 2);
                        vel_new_list[faceI] = vector(Ux_tmp, 0, 0);
                    }
                    ITHACAutilities::assignBC(_U(), inlet_idx, vel_new_list);
                    assignIF(_U(), vel_new_ini);
                    truthSolve(mu_now, "Online");
                    restart();
                }
            }            
        }
};

void offline_stage(tutorial03& example);
void online_stage(tutorial03& example, tutorial03& onlineExample);

int main(int argc, char* argv[])
{
#if OPENFOAM > 1812
    // load stage to perform
    argList::addOption("stage", "offline", "Perform offline stage");
    argList::addOption("stage", "online", "Perform online stage");
    // add options for parallel run
    HashTable<string> validParOptions;
    validParOptions.set
    (
        "stage",
        "offline"
    );
    validParOptions.set
    (
        "stage",
        "online"
    );
    Pstream::addValidParOptions(validParOptions);
    // Construct the tutorial object
    tutorial03 example(argc, argv);

    if (example._args().get("stage").match("offline"))
    {
        // perform the offline stage, extracting the modes from the snapshots'
        // dataset corresponding to parOffline
        offline_stage(example);
    }
    else if (example._args().get("stage").match("online"))
    {
        // load precomputed modes and reduced matrices
        offline_stage(example);
        // Construct the tutorial object for online truthSolve
        tutorial03 onlineExample(argc, argv);
        online_stage(example, onlineExample);
    }
    else
    {
        Info << "Pass '-stage offline', '-stage online'" << endl;
    }

#else

    if (argc == 1)
    {
        Info << "Pass 'offline' or 'online' as first arguments."
             << endl;
        exit(0);
    }

    // process arguments removing "offline" or "online" keywords
    int argc_proc = argc - 1;
    char* argv_proc[argc_proc];
    argv_proc[0] = argv[0];

    if (argc > 2)
    {
        std::copy(argv + 2, argv + argc, argv_proc + 1);
    }

    argc--;
    // Construct the tutorial object
    tutorial03 example(argc, argv);

    if (std::strcmp(argv[1], "offline") == 0)
    {
        // perform the offline stage, extracting the modes from the snapshots' dataset corresponding to parOffline
        offline_stage(example);
    }
    else if (std::strcmp(argv[1], "online") == 0)
    {
        // load precomputed modes and reduced matrices
        offline_stage(example);
        // Construct the tutorial object for online truthSolve
        tutorial03 onlineExample(argc, argv);
        online_stage(example, onlineExample);
    }
    else
    {
        Info << "Pass offline, online" << endl;
    }

#endif
    exit(0);
}

void offline_stage(tutorial03& example)
{
    // Read some parameters from file
    ITHACAparameters* para =
        ITHACAparameters::getInstance(example._mesh(), example._runTime());
    int NmodesUout = para->ITHACAdict->lookupOrDefault<int>("NmodesUout", 15);
    int NmodesPout = para->ITHACAdict->lookupOrDefault<int>("NmodesPout", 15);
    int NmodesSUPout = para->ITHACAdict->lookupOrDefault<int>("NmodesSUPout", 15);
    int NmodesUproj = para->ITHACAdict->lookupOrDefault<int>("NmodesUproj", 10);
    int NmodesPproj = para->ITHACAdict->lookupOrDefault<int>("NmodesPproj", 10);
    int NmodesSUPproj =
        para->ITHACAdict->lookupOrDefault<int>("NmodesSUPproj", 10);
    // Read the par file where the training parameters are stored
    word filename("./parOffline");
    example.mu = ITHACAstream::readMatrix(filename);
    // Set the inlet boundaries patch 0 directions x and y
    example.inletIndex.resize(5, 2);
    example.inletIndex.setZero();
    // example.inletIndex(0, 0) = 0;
    // example.inletIndex(0, 1) = 0;
    // example.inletIndex(1, 0) = 0;
    // example.inletIndex(1, 1) = 0;
    // Perform the offline solve
    example.offlineSolve();
    // Solve the supremizer problem
    example.solvesupremizer();

    if (example.bcMethod == "lift")
    {
        ITHACAstream::read_fields(example.liftfield, example._U(), "./lift/");
        // ITHACAutilities::normalizeFields(example.liftfield);
        // Homogenize the snapshots
        if (example.nonUniformbc)
        {
            example.computeLift(example.Ufield, example.liftfield, example.Uomfield, example.mu);
        }
        else
        {
            example.computeLift(example.Ufield, example.liftfield, example.Uomfield);
        }
        // Perform POD on the velocity snapshots
        ITHACAPOD::getModes(example.Uomfield, example.Umodes, example._U().name(),
                            example.podex, 0, 0, NmodesUout);
    }
    else
    {
        ITHACAPOD::getModes(example.Ufield, example.Umodes, example._U().name(),
                            example.podex, 0, 0, NmodesUout);
    }

    // Perform POD on velocity pressure and supremizers and store the first 10
    // modes
    ITHACAPOD::getModes(example.Pfield, example.Pmodes, example._p().name(),
                        example.podex, 0, 0, NmodesPout);
    ITHACAPOD::getModes(example.supfield, example.supmodes, example._U().name(),
                        example.podex, example.supex, 1, NmodesSUPout);
    // Perform the Galerkin Projection
    example.projectSUP("./Matrices", NmodesUproj, NmodesPproj, NmodesSUPproj);
}

void online_stage(tutorial03 &example, tutorial03 &onlineExample)
{
    // Create the reduced object
    reducedSteadyNS reduced(example);
    // Read the par file where the test parameters are stored
    word filename("./parOnline");
    example.mu = ITHACAstream::readMatrix(filename);

    // Set the inlet velocity
    Eigen::MatrixXd vel_now(5, 1);

    // Run truthsolve for online parameters
    onlineExample.mu = ITHACAstream::readMatrix(filename);
    onlineExample.onlineSolve();

    // used only for penalty approach
    reduced.tauU = Eigen::MatrixXd::Zero(1, 1);
    reduced.tauU(0, 0) = 1e-1;

    // Perform an online solve for the new values of inlet velocities
    for (label k = 0; k < example.mu.cols(); k++)
    {
        for(label n=0; n < example.mu.rows(); n++)
        {
            vel_now(n, 0) = example.mu(n, k);
        } 
        Info << "Evaluation of the reduced order model on the test set" << endl;
        Info << "Inlet Ux = " << vel_now << " nu = " << 1
             << endl;
        // Set the reduced viscosity
        reduced.nu = 1.0;
        reduced.solveOnline_sup(vel_now);
        Eigen::MatrixXd tmp_sol(reduced.y.rows() + 1, 1);
        tmp_sol(0) = k + 1;
        tmp_sol.col(0).tail(reduced.y.rows()) = reduced.y;
        reduced.online_solution.append(tmp_sol);
    }

    // Save the online solution
    ITHACAstream::exportMatrix(reduced.online_solution, "red_coeff", "python",
                               "./ITHACAoutput/Online/red_coeff");
    ITHACAstream::exportMatrix(reduced.online_solution, "red_coeff", "matlab",
                               "./ITHACAoutput/Online/red_coeff");
    ITHACAstream::exportMatrix(reduced.online_solution, "red_coeff", "eigen",
                               "./ITHACAoutput/Online/red_coeff");
    // Reconstruct and export the solution
    reduced.reconstruct(true, "./ITHACAoutput/Online/");

    // Write error of online solutions
    Eigen::MatrixXd errL2U = ITHACAutilities::errorL2Rel(onlineExample.Ufield,
                                                         reduced.uRecFields);
    Eigen::MatrixXd errL2P = ITHACAutilities::errorL2Rel(onlineExample.Pfield,
                                                         reduced.pRecFields);
    ITHACAstream::exportMatrix(errL2U, "errL2U", "python",
                               "./ITHACAoutput/Online/ErrorsL2/");
    ITHACAstream::exportMatrix(errL2P, "errL2P", "python",
                               "./ITHACAoutput/Online/ErrorsL2/");

    // Export errorfields
    // Perform an online solve for the new values of inlet velocities
    for (label k = 0; k < example.mu.cols(); k++)
    {
        volVectorField Uerror("Uerror", onlineExample.Ufield[k] - reduced.uRecFields[k]);
        volScalarField perror("perror", onlineExample.Pfield[k] - reduced.pRecFields[k]);
        ITHACAstream::exportSolution(Uerror,
                                     name(k+1),
                                     "./ITHACAoutput/Online/");
        ITHACAstream::exportSolution(perror,
                                     name(k+1),
                                     "./ITHACAoutput/Online/");
    }

    if (Pstream::parRun())
    {
        bool endedOnline = true;
        reduce(endedOnline, sumOp<label>());
    }
}

//--------
/// \dir 03steadyNS Folder of the turorial 3
/// \file
/// \brief Implementation of a tutorial of a steady Navier-Stokes problem

/// \example 03steadyNS.C
/// \section intro_sreadyNS Introduction to tutorial 3
/// The problems consists of steady Navier-Stokes problem with parametrized
/// viscosity. The physical problem is the backward facing step depicted in the
/// following image: \image html step.png At the inlet a uniform and constant
/// velocity equal to 1 m/s is prescribed.
///
/// \section code03 A detailed look into the code
///
/// In this section are explained the main steps necessary to construct the
/// tutorial N°3
///
/// \subsection header The necessary header files
///
/// First of all let's have a look to the header files that needs to be included
/// and what they are responsible for:
///
/// The header file of ITHACA-FV necessary for this tutorial
///
/// \dontinclude 03steadyNS.C
/// \skip steadyNS
/// \until reducedSteady
///
/// \subsection classtuto03 Implementation of the tutorial03 class
///
/// Then we can define the tutorial03 class as a child of the steadyNS class
/// \skipline tutorial03
/// \until {}
///
/// The members of the class are the fields that needs to be manipulated during
/// the resolution of the problem
///
/// Inside the class it is defined the offlineSolve method according to the
/// specific parametrized problem that needs to be solved.
///
/// \skipline offlineSolve
/// \until {
///
///
/// If the offline solve has already been performed than read the existing
/// snapshots
///
/// \skipline offline
/// \until }
///
/// else perform the offline solve where a loop over all the parameters is
/// performed:
///
/// \skipline else
/// \until }
/// \skipline }
///
/// See also the steadyNS class for the definition of the methods.
///
/// \subsection main Definition of the main function
///
/// Once the tutorial03 class is defined the main function is defined,
/// an example of type tutorial03 is constructed:
///
/// \skipline tutorial03
///
/// In this case the vector of parameter is read from a txt file
///
/// \skipline word
/// \until example.mu
///
/// The inlet boundary is set:
///
/// \skipline example.inlet
/// \until example.inletIndex(0, 1) = 0;
///
/// and the offline stage is performed:
///
/// \skipline Solve()
///
/// and the supremizer problem is solved:
///
/// \skipline supremizer()
///
/// In order to show the functionality of reading fields in this case the
/// lifting function is read from a precomputed simulation with a unitary inlet
/// velocity:
///
/// \skipline stream
///
/// Then the snapshots matrix is homogenized:
///
/// \skipline computeLift
///
/// and the modes for velocity, pressure and supremizers are obtained:
///
/// \skipline getModes
/// \until supfield
///
/// then the projection onto the POD modes is performed with:
///
/// \skipline projectSUP
///
/// the reduced object is constructed:
///
/// \skipline reducedSteady
///
/// and the online solve is performed for some values of the viscosity:
///
/// \skipline Eigen::
/// \until }
///
/// The vel_now matrix in this case is not used since there are no parametrized
/// boundary conditions.
///
/// The viscosity is set with the command:
///
/// \code
/// reduced.nu = example.mu(k,0)
/// \endcode
///
/// finally the online solution stored during the online solve is exported to
/// file in three different formats with the lines:
///
/// \skipline exportMatrix
/// \until "eigen"
///
/// and the online solution is reconstructed and exported to file
///
/// \skipline reconstruct
///
///
///
///
///
/// \section plaincode The plain program
/// Here there's the plain code
///
