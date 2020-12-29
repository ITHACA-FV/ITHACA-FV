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
    Example of Boussinesq approximation for two way coupling NS-momentum equation
    and heat transport equation for enclosed flows.
SourceFiles
    10UnsteadyBBEnclosed.C
\*---------------------------------------------------------------------------*/

#include "UnsteadyBB.H"
#include "ITHACAPOD.H"
#include "ReducedUnsteadyBB.H"
#include "ITHACAstream.H"
#include <chrono>
#include <math.h>
#include <iomanip>

class tutorial10: public UnsteadyBB
{
    public:
        explicit tutorial10(int argc, char* argv[])
            :
            UnsteadyBB(argc, argv),
            U(_U()),
            p(_p()),
            p_rgh(_p_rgh()),
            T(_T())
        {}

        // Fields To Perform
        volVectorField& U;
        volScalarField& p;
        volScalarField& p_rgh;
        volScalarField& T;

        void offlineSolve(Eigen::MatrixXd par_BC)
        {
            List<scalar> mu_now(1);

            if (offline)
            {
                ITHACAstream::read_fields(Ufield, U, "./ITHACAoutput/Offline/");
                ITHACAstream::read_fields(Tfield, T, "./ITHACAoutput/Offline/");
            }
            else
            {
                for (label k = 0; k < par_BC.rows(); k++)
                {
                    for (label j = 0; j < par_BC.cols(); j++)
                    {
                        for (label i = 0; i < mu.cols(); i++)
                        {
                            mu_now[0] = mu(0, i);
                        }

                        assignBC(T, inletIndexT(j, 0), par_BC(k, j));
                    }

                    truthSolve(mu_now);
                }
            }
        }


        void onlineSolveFull(Eigen::MatrixXd par_BC, label para_set_BC, fileName folder)
        {
            if (ITHACAutilities::check_folder(folder))
            {
            }
            else
            {
                mkDir(folder);
                ITHACAutilities::createSymLink(folder);
                label i = para_set_BC;

                for (label j = 0; j < par_BC.cols(); j++)
                {
                    assignBC(T, inletIndexT(j, 0), par_BC(i, j));
                }

                truthSolve(folder);
            }
        }

        void onlineSolveRead(fileName folder)
        {
            if (ITHACAutilities::check_folder(folder))
            {
                ITHACAstream::read_fields(Ufield_on, U, folder);
                ITHACAstream::read_fields(Tfield_on, T, folder);
            }
            else
            {
            }
        }


        // Method to compute the lifting function for temperature
        void liftSolveT()
        {
            for (label k = 0; k < inletIndexT.rows(); k++)
            {
                Time& runTime = _runTime();
                fvMesh& mesh = _mesh();
                volScalarField T = _T();
                simpleControl simple(mesh);
                volScalarField& alphat = _alphat();
                dimensionedScalar& Pr = _Pr();
                dimensionedScalar& Prt = _Prt();
                label BCind = inletIndexT(k, 0);
                volScalarField Tlift("Tlift" + name(k), T);
                instantList Times = runTime.times();
                runTime.setTime(Times[1], 1);
                Info << "Solving a lifting Problem" << endl;
                scalar t1 = 1;
                scalar t0 = 0;

                for (label j = 0; j < T.boundaryField().size(); j++)
                {
                    assignIF(Tlift, t0);

                    if (j == BCind)
                    {
                        assignBC(Tlift, j, t1);
                    }
                    else if (T.boundaryField()[BCind].type() == "fixedValue")
                    {
                        assignBC(Tlift, j, t0);
                    }
                    else
                    {
                    }
                }

                while (simple.correctNonOrthogonal())
                {
                    alphat = turbulence->nut() / Prt;
                    alphat.correctBoundaryConditions();
                    volScalarField alphaEff("alphaEff", turbulence->nu() / Pr + alphat);
                    fvScalarMatrix TEqn
                    (
                        fvm::laplacian(alphaEff, Tlift)
                    );
                    TEqn.solve();
                    Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                         << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                         << nl << endl;
                }

                Tlift.write();
                liftfieldT.append(tmp<volScalarField>(Tlift));
            }
        }
};


/*---------------------------------------------------------------------------*\
                               Starting the MAIN
\*---------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
    // Construct the tutorial object
    tutorial10 example(argc, argv);
    // the offline samples for the boundary conditions
    word par_offline_BC("./par_offline_BC");
    Eigen::MatrixXd par_off_BC = ITHACAstream::readMatrix(par_offline_BC);
    // the samples which will be used for setting the boundary condition in the online stage
    word par_online_BC("./par_online_BC");
    Eigen::MatrixXd par_on_BC = ITHACAstream::readMatrix(par_online_BC);
    // Read some parameters from file
    ITHACAparameters* para = ITHACAparameters::getInstance(example._mesh(),
                             example._runTime());
    int NmodesUproj   = para->ITHACAdict->lookupOrDefault<int>("NmodesUproj", 5);
    int NmodesPproj   = para->ITHACAdict->lookupOrDefault<int>("NmodesPproj", 5);
    int NmodesTproj   = para->ITHACAdict->lookupOrDefault<int>("NmodesTproj", 5);
    int NmodesSUPproj = para->ITHACAdict->lookupOrDefault<int>("NmodesSUPproj", 5);
    int NmodesOut     = para->ITHACAdict->lookupOrDefault<int>("NmodesOut", 15);
    /// Set the number of parameters
    example.Pnumber = 1;
    /// Set samples
    example.Tnumber = 1;
    /// Set the parameters infos
    example.setParameters();
    /// Set the parameter ranges
    example.mu_range(0, 0) = 0.00001;
    example.mu_range(0, 1) = 0.00001;
    // Generate equispaced samples inside the parameter range
    example.genEquiPar();
    /// Set the inlet Temperature boundaries where there are non homogeneous boundary conditions
    example.inletIndexT.resize(2, 1);
    example.inletIndexT << 1, 2;
    /// Time parameters
    example.startTime = 0.0;
    example.finalTime = 10.0;
    example.timeStep = 0.005;
    example.writeEvery = 0.01;
    // Perform the Offline Solve;
    example.offlineSolve(par_off_BC);
    // Search the lift function for the temperature
    example.liftSolveT();
    // Create homogeneous basis functions for temperature
    example.computeLiftT(example.Tfield, example.liftfieldT, example.Tomfield);
    // Perform a POD decomposition for velocity temperature and pressure fields
    ITHACAPOD::getModes(example.Ufield, example.Umodes, example._U().name(),
                        example.podex, 0, 0,
                        NmodesOut);
    ITHACAPOD::getModes(example.Tomfield, example.Tmodes, example._T().name(),
                        example.podex, 0, 0,
                        NmodesOut);
    // Create a list with number of modes for which the projection needs to be performed
    Eigen::MatrixXd List_of_modes(NmodesOut - 5, 1);

    for (int i = 0; i < List_of_modes.rows(); i++)
    {
        List_of_modes(i, 0) = i + 1;
    }

    // Export with number of modes for which the projection needs to be performed
    ITHACAstream::exportMatrix(List_of_modes, "List_of_modes", "eigen",
                               "./ITHACAoutput/l2error");
    // Create locally the temperature modes
    PtrList<volScalarField> TLmodes;

    for (label k = 0; k < example.liftfieldT.size(); k++)
    {
        TLmodes.append(tmp<volScalarField>(example.liftfieldT[k]));
    }

    for (label k = 0; k < List_of_modes.size(); k++)
    {
        TLmodes.append(tmp<volScalarField>(example.Tmodes[k]));
    }

    // Perform the projection for all number of modes in list List_of_modes
    Eigen::MatrixXd L2errorProjMatrixU(example.Ufield.size(), List_of_modes.rows());
    Eigen::MatrixXd L2errorProjMatrixT(example.Tfield.size(), List_of_modes.rows());

    // Calculate the coefficients and L2 error and store the error in a matrix for each number of modes
    for (int i = 0; i < List_of_modes.rows(); i++)
    {
        Eigen::MatrixXd coeffU = ITHACAutilities::getCoeffs(example.Ufield,
                                 example.Umodes,
                                 List_of_modes(i, 0) + example.liftfield.size() + NmodesSUPproj);
        Eigen::MatrixXd coeffT = ITHACAutilities::getCoeffs(example.Tfield, TLmodes,
                                 List_of_modes(i, 0) + example.liftfieldT.size());
        PtrList<volVectorField> rec_fieldU = ITHACAutilities::reconstructFromCoeff(
                example.Umodes, coeffU, List_of_modes(i, 0));
        PtrList<volScalarField> rec_fieldT = ITHACAutilities::reconstructFromCoeff(
                TLmodes, coeffT, List_of_modes(i, 0) + example.liftfieldT.size());
        Eigen::MatrixXd L2errorProjU = ITHACAutilities::errorL2Rel(example.Ufield,
                                       rec_fieldU);
        Eigen::MatrixXd L2errorProjT = ITHACAutilities::errorL2Rel(example.Tfield,
                                       rec_fieldT);
        L2errorProjMatrixU.col(i) = L2errorProjU;
        L2errorProjMatrixT.col(i) = L2errorProjT;
    }

    // Export the matrix containing the error
    ITHACAstream::exportMatrix(L2errorProjMatrixU, "L2errorProjMatrixU", "eigen",
                               "./ITHACAoutput/l2error");
    ITHACAstream::exportMatrix(L2errorProjMatrixT, "L2errorProjMatrixT", "eigen",
                               "./ITHACAoutput/l2error");
    // Get reduced matrices
    example.projectSUP("./Matrices", NmodesUproj, NmodesPproj, NmodesTproj,
                       NmodesSUPproj);
    // Resize the modes for projection
    example.Tmodes.resize(NmodesTproj);
    example.Umodes.resize(NmodesUproj);
    // Online part
    ReducedUnsteadyBB reduced(example);
    // Set values of the online solve
    reduced.nu = 0.00001;
    reduced.Pr = 0.71;
    reduced.tstart = 0.0;
    reduced.finalTime = 10;
    reduced.dt = 0.005;
    // No parametrization of velocity on boundary
    Eigen::MatrixXd vel_now_BC(0, 0);

    // Set the online temperature BC and solve reduced model
    for (label k = 0; k < (par_on_BC.rows()); k++)
    {
        Eigen::MatrixXd temp_now_BC(2, 1);
        temp_now_BC(0, 0) = par_on_BC(k, 0);
        temp_now_BC(1, 0) = par_on_BC(k, 1);
        reduced.solveOnline_sup(temp_now_BC, vel_now_BC, k, par_on_BC.rows());
        reduced.reconstruct_sup("./ITHACAoutput/ReconstructionSUP", 2);
    }

    // Performing full order simulation for second parameter set - temp_BC
    tutorial10 HFonline2(argc, argv);
    HFonline2.Pnumber = 1;
    HFonline2.Tnumber = 1;
    HFonline2.setParameters();
    HFonline2.mu_range(0, 0) = 0.00001;
    HFonline2.mu_range(0, 1) = 0.00001;
    HFonline2.genEquiPar();
    HFonline2.inletIndexT.resize(2, 1);
    HFonline2.inletIndexT << 1, 2;
    HFonline2.startTime = 0.0;
    HFonline2.finalTime = 10;
    HFonline2.timeStep = 0.005;
    HFonline2.writeEvery = 0.01;
    // Reconstruct the online solution
    HFonline2.onlineSolveFull(par_on_BC, 1,
                              "./ITHACAoutput/HFonline2");
    // Performing full order simulation for third parameter set - temp_BC
    tutorial10 HFonline3(argc, argv);
    HFonline3.Pnumber = 1;
    HFonline3.Tnumber = 1;
    HFonline3.setParameters();
    HFonline3.mu_range(0, 0) = 0.00001;
    HFonline3.mu_range(0, 1) = 0.00001;
    HFonline3.genEquiPar();
    HFonline3.inletIndexT.resize(2, 1);
    HFonline3.inletIndexT << 1, 2;
    HFonline3.startTime = 0.0;
    HFonline3.finalTime = 10.0;
    HFonline3.timeStep = 0.005;
    HFonline3.writeEvery = 0.01;
    // Reconstruct the online solution
    HFonline3.onlineSolveFull(par_on_BC, 2,
                              "./ITHACAoutput/HFonline3");
    // Reading in the high-fidelity solutions for the parameter set
    // for which the offline solve has been performed
    example.onlineSolveRead("./ITHACAoutput/Offline/");
    // Reading in the high-fidelity solutions for the second parameter set
    example.onlineSolveRead("./ITHACAoutput/HFonline2/");
    // Reading in the high-fidelity solutions for the second parameter set
    example.onlineSolveRead("./ITHACAoutput/HFonline3/");
    // Calculate error between online- and corresponding full order solution
    Eigen::MatrixXd L2errorMatrixU = ITHACAutilities::errorL2Rel(
                                         example.Ufield_on, reduced.UREC);
    Eigen::MatrixXd L2errorMatrixT = ITHACAutilities::errorL2Rel(
                                         example.Tfield_on, reduced.TREC);
    //Export the matrix containing the error
    ITHACAstream::exportMatrix(L2errorMatrixU, "L2errorMatrixU", "eigen",
                               "./ITHACAoutput/l2error");
    ITHACAstream::exportMatrix(L2errorMatrixT, "L2errorMatrixT", "eigen",
                               "./ITHACAoutput/l2error");
    exit(0);
}

/// \dir 10unsteadyBB_open Folder of the tutorial 10
/// \file
/// \brief Implementation of tutorial 10 for an unsteady Buoyant Boussinesq problem for an enclosed cavity
///
/// \example 10UnsteadyBBEnclosed.C
/// \section intro_unsteadyBB Introduction to tutorial 10
/// In this tutorial an unsteady Buoyant Boussinesq (BB) 2D problem with paramerized temperature boundary conditions is implemented.
/// The physical problem represents a differentially heated cavity. A uniform temperature is
/// set to the left (hot) and the right (cold) sides of the cavity while the other sides are set to adiabatic. The cavity aspect ratio is 1.0.
/// The flow is laminar and the working fluid is air with Pr = 0.7. The ambient temperature is 300 K. The hot wall, Th, has a temperature of 301.5 K,
/// while the cold wall, Tc, is set to 298.5 K. The initial condition for the velocity is (0.0001, 0, 0) m/s.
///
/// The following image illustrates the geometry with the boundary conditions.
/// \image html setup_diff_heated_cavity.png
///
/// \section code04 A detailed look into the code
///
/// In this section we explain the main steps necessary to construct the tutorial N°10
///
/// \subsection header ITHACA-FV header files
///
/// First of all let's have a look at the header files that need to be included and what they are responsible for.
///
/// The header files of ITHACA-FV necessary for this tutorial are: <UnsteadyBB.H> for the full order unsteady BB problem,
/// <ITHACAPOD.H> for the POD decomposition, <ReducedUnsteadyBB.H> for the construction of the reduced order problem,
/// and finally <ITHACAstream.H> for some ITHACA input-output operations.
///
/// \dontinclude 10UnsteadyBBEnclosed.C
/// \skip UnsteadyBB
/// \until ITHACAstream
///
/// \subsection classtutorial10 Definition of the tutorial10 class
///
/// We define the tutorial10 class as a child of the UnsteadyBB class.
/// The constructor is defined with members that are the fields need to be manipulated
/// during the resolution of the full order problem using the BuoyancyBoussinesqPimpleFoam solver. Such fields are
/// also initialized with the same initial conditions in the solver.
/// \skipline tutorial10
/// \until {}
///
/// Inside the tutorial10 class we define the offlineSolve method according to the
/// specific parametrized problem that needs to be solved. If the offline solve has
/// been previously performed then the method simply reads in the existing snapshots from the Offline directory.
/// If not, the offlineSolve changes the values of the temperature boundary conditions before performing the offline solve.
///
/// \skipline offlineSolve
/// \until }
/// \skipline else
/// \until }
/// \skipline assignBC
/// \skipline }
/// \skipline truthSolve
/// \skipline }
/// \skipline }
/// \skipline }
///
/// In order to calculate the L2 error between the online solved reduced order model and the corresponding full order solution
/// at the end of this tutorial, the full order solve can be performed with the onlineSolveFull for a parameter set for which
/// the reduced order model has been solved online.
///
/// \skipline onlineSolveFull
/// \until }
/// \skipline else
/// \until }
/// \skipline truthSolve
/// \skipline }
/// \skipline }
///
/// If the full order solve has been performed previously then the existing snapshots can be read in from the
/// specified directory.
///
/// \skipline onlineSolveRead
/// \until }
/// \skipline else
/// \until }
/// \skipline }
///
/// The liftingfunctions for temperature in this problem are determined by solving a steady state laplacian function
///
/// \skipline liftSolveT
/// \until liftfieldT
/// \skipline }
/// \until }
///
/// \subsection main Definition of the main function
///
/// In this section we show the definition of the main function.
/// First we construct the object "example" of type tutorial10:
///
/// \skipline example
///
/// The parameter sets for the temperature BC are read in from file, both for the offline and online phase.
/// Each row corresponds to a parameter set and each column to a specific boundary.
///
/// \skipline par_offline_BC
/// \until Eigen
/// \skipline par_online_BC
/// \until Eigen
///
/// Then we parse the ITHACAdict file to determine the number of modes
/// to be written out and also the ones to be used for projection of
/// the velocity, pressure, temperature and the supremizer:
/// \skipline ITHACAparameters
/// \until NmodesOut
///
/// We note that a default value can be assigned in case the parser did
/// not find the corresponding string in the ITHACAdict file.
///
/// Now the kinematic viscosity is set to 0.00001 m^2/s. It is possible to parametrize the viscosity. For more
/// info take a look at tutorial04.
///
/// \skipline example.Pnumber
/// \until example.genEquiPar()
///
/// After that we set the boundaries where we have a non homogeneous BC for temperature. Patch 1 corresponds to the
/// hot boundary and Patch 2 to the cold boundary.
///
/// \skipline example.inletIndexT
/// \until example.inletIndexT << 1, 2;
///
/// Furthermore, we set the parameters for the time integration. In this case the simulation time is
/// 10 seconds, with a step size = 0.005 seconds, and the data is written every 0.01 seconds, i.e.
///
/// \skipline example.startTime
/// \until example.writeEvery
///
/// Now we are ready to perform the offline stage:
///
/// \skipline example.offlineSolve
///
/// No lifting function has to be computed for velocity as the velocity is zero everywhere on the boundary.
/// For the temperature, the lifting functions are determined by:
///
/// \skipline liftSolveT()
///
/// Then we create homogenuous basis functions for the temperature:
///
/// \skipline computeLiftT(example.Tfield, example.liftfieldT, example.Tomfield)
///
/// And after that, we obtain the modes for velocity and temperature:
///
/// \skipline getModes
/// \skipline getModes
///
/// Before continuiting, the projection error is calculated for certain number of modes for both
/// temperature and velocity.
///
/// \skipline List_of_modes
/// \until }
/// \skipline exportMatrix
///
/// The temperature modes are created locally
///
/// \skipline TLmodes
/// \until }
/// \skipline for
/// \until }
///
/// and the projection onto the POD modes is performed with:
///
/// \skipline L2errorProjMatrixU
/// \skipline L2errorProjMatrixT
///
/// Finally the L2 error between full order solution and projection of the basis are calculated
/// \skipline for
/// \until }
///
/// and exported
///
/// \skipline L2errorProjMatrixU
/// \until L2errorProjMatrixT
/// \skipline "./ITHACAoutput/l2error");
///
/// Then the projection onto the POD modes is performed to get the reduced matrices
///
/// \skipline projectSUP
/// \until NmodesSUPproj
///
/// and the modes are resized to the number for which the projection has been performed
///
/// \skipline Tmodes
/// \until Umodes
///
/// Now that we obtained all the necessary information from the POD decomposition and the reduced matrices,
/// we are ready to construct the dynamical system for the reduced order model (ROM). We proceed
/// by constructing the object "reduced" of type ReducedUnsteadyBB:
///
/// \skipline ReducedUnsteadyBB
///
/// We can use the new constructed ROM to perform the online procedure, from which we can simulate the
/// problem at new set of parameters. For instance, we solve the problem for 10 seconds of physical time:
///
/// \skipline reduced.nu
/// \until reduced.dt
///
/// And then we can use the new constructed ROM to perform the online procedure, from which we can simulate the
/// problem at new set of parameters.
/// As the velocity is homogenous on the boundary, it is not parametrized and therefore an empty matrix is created:
///
/// \skipline Eigen::MatrixXd
///
/// Lastly the online procedure is performed for all temperature boundary sets as defined in the par_online_BC
/// file. And he ROM solution is reconstructed and exported:
///
/// \skipline for
/// \until }
///
/// For the second and third parameter set the full order solution is calculated
///
/// \skipline // Performing full order simulation for second parameter set - temp_BC
/// \until "./ITHACAoutput/HFonline3/");
///
/// and the solutions are read in:
/// \skipline example.onlineSolveRead("./ITHACAoutput/Offline/");
/// \until example.onlineSolveRead("./ITHACAoutput/HFonline3/");
///
/// Finally the L2 error between full and reduced order solutions is calculated
///
/// \skipline Ufield_on
/// \until Tfield_on
/// \skipline L2errorMatrixU
/// \until L2errorMatrixT
/// \skipline "./ITHACAoutput/l2error");
///
/// \section plaincode The plain program
/// Here's the plain code
///
