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
    and heat transport equation for open flows where pressure connot be neglected.
SourceFiles
    11UnsteadyBBOpen.C
\*---------------------------------------------------------------------------*/

#include "UnsteadyBB.H"
#include "ITHACAPOD.H"
#include "ReducedUnsteadyBB.H"
#include "ITHACAstream.H"
#include <chrono>
#include <math.h>
#include <iomanip>

class tutorial11: public UnsteadyBB
{
    public:
        explicit tutorial11(int argc, char* argv[])
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
                ITHACAstream::read_fields(Pfield, p, "./ITHACAoutput/Offline/");
                ITHACAstream::read_fields(Prghfield, p_rgh, "./ITHACAoutput/Offline/");
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
                volScalarField& T = _T();
                volVectorField& U = _U();
                surfaceScalarField& phi = _phi();
                phi = linearInterpolate(U) & mesh.Sf();
                simpleControl simple(mesh);
                // IOMRFZoneList& MRF = _MRF();
                // singlePhaseTransportModel& laminarTransport = _laminarTransport();
                // volScalarField& nut = _nut();
                volScalarField& alphat = _alphat();
                // dimensionedScalar& nu = _nu();
                dimensionedScalar& Pr = _Pr();
                dimensionedScalar& Prt = _Prt();
                label BCind = inletIndexT(k, 0);
                volScalarField Tlift("Tlift" + name(k), T);
                instantList Times = runTime.times();
                runTime.setTime(Times[1], 1);
                Info << "Solving a lifting Problem" << endl;
                scalar t1 = 1;
                scalar t0 = 0;
                alphat = turbulence->nut() / Prt;
                alphat.correctBoundaryConditions();
                volScalarField alphaEff("alphaEff", turbulence->nu() / Pr + alphat);

                for (label j = 0; j < T.boundaryField().size(); j++)
                {
                    if (j == BCind)
                    {
                        assignBC(Tlift, j, t1);
                        assignIF(Tlift, t0);
                    }
                    else if (T.boundaryField()[BCind].type() == "fixedValue")
                    {
                        assignBC(Tlift, j, t0);
                        assignIF(Tlift, t0);
                    }
                    else
                    {
                    }
                }

                while (simple.correctNonOrthogonal())
                {
                    fvScalarMatrix TEqn
                    (
                        fvm::div(phi, Tlift)
                        - fvm::laplacian(alphaEff, Tlift)
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
    tutorial11 example(argc, argv);
    // the offline samples for the boundary conditions
    word par_offline_BC("./par_offline_BC");
    // the samples which will be used for setting the boundary condition in the online stage
    word par_online_BC("./par_online_BC");
    Eigen::MatrixXd par_off_BC = ITHACAstream::readMatrix(par_offline_BC);
    Eigen::MatrixXd par_on_BC = ITHACAstream::readMatrix(par_online_BC);
    // Read some parameters from file
    ITHACAparameters* para = ITHACAparameters::getInstance(example._mesh(),
                             example._runTime());
    int NmodesUproj   = para->ITHACAdict->lookupOrDefault<int>("NmodesUproj", 5);
    int NmodesPproj   = para->ITHACAdict->lookupOrDefault<int>("NmodesPproj", 5);
    int NmodesPrghproj = para->ITHACAdict->lookupOrDefault<int>("NmodesPrghproj",
                         5);
    int NmodesTproj   = para->ITHACAdict->lookupOrDefault<int>("NmodesTproj", 5);
    int NmodesSUPproj = para->ITHACAdict->lookupOrDefault<int>("NmodesSUPproj", 5);
    int NmodesOut     = para->ITHACAdict->lookupOrDefault<int>("NmodesOut", 15);
    // Set the number of parameters
    example.Pnumber = 1;
    // Set samples
    example.Tnumber = 1;
    // Set the parameters infos
    example.setParameters();
    // Set the parameter ranges
    example.mu_range(0, 0) = 0.00001;
    example.mu_range(0, 1) = 0.00001;
    // Generate equispaced samples inside the parameter range
    example.genEquiPar();
    // Set the inlet Temperature boundaries where we have non homogeneous boundary conditions
    example.inletIndexT.resize(3, 1);
    example.inletIndexT << 1, 2, 3;
    // Set the inlet Velocity boundaries where we have non homogeneous boundary conditions
    example.inletIndex.resize(1, 2);
    example.inletIndex << 3, 0;
    // Time parameters
    example.startTime = 0.0;
    example.finalTime = 5.0;
    example.timeStep = 0.002;
    example.writeEvery = 0.01;
    // Perform the Offline Solve;
    example.offlineSolve(par_off_BC);
    // Search the lift function for the velocity
    example.liftSolve();
    // Search the lift function for the temperature
    example.liftSolveT();
    // Create homogeneous basis functions for velocity
    example.computeLift(example.Ufield, example.liftfield, example.Uomfield);
    // Create homogeneous basis functions for temperature
    example.computeLiftT(example.Tfield, example.liftfieldT, example.Tomfield);
    // Perform a POD decomposition for velocity ,temperature and pressure fields
    ITHACAPOD::getModes(example.Uomfield, example.Umodes, example._U().name(),
                        example.podex, 0, 0,
                        NmodesOut);
    ITHACAPOD::getModes(example.Pfield, example.Pmodes, example._p().name(),
                        example.podex, 0, 0,
                        NmodesOut);
    ITHACAPOD::getModes(example.Prghfield, example.Prghmodes,
                        example._p_rgh().name(), example.podex, 0, 0,
                        NmodesOut);
    ITHACAPOD::getModes(example.Tomfield, example.Tmodes, example._T().name(),
                        example.podex, 0, 0,
                        NmodesOut);
    // Solve the supremizer problem
    example.solvesupremizer("modes");
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

    // Create locally the velocity modes
    PtrList<volVectorField> ULmodes;

    for (label k = 0; k < example.liftfield.size(); k++)
    {
        ULmodes.append(tmp<volVectorField>(example.liftfield[k]));
    }

    for (label k = 0; k < List_of_modes.size(); k++)
    {
        ULmodes.append(tmp<volVectorField>(example.Umodes[k]));
    }

    for (label k = 0; k < NmodesSUPproj; k++)
    {
        ULmodes.append(tmp<volVectorField>(example.supmodes[k]));
    }

    // Perform the projection for all number of modes in List_of_modes for temperature and velocity
    Eigen::MatrixXd L2errorProjMatrixU(example.Ufield.size(), List_of_modes.rows());
    Eigen::MatrixXd L2errorProjMatrixT(example.Tfield.size(), List_of_modes.rows());

    // Calculate the coefficients and L2 error and store the error in a matrix for each number of modes
    for (int i = 0; i < List_of_modes.rows(); i++)
    {
        Eigen::MatrixXd coeffU = ITHACAutilities::getCoeffs(example.Ufield, ULmodes,
                                 List_of_modes(i, 0) + example.liftfield.size() + NmodesSUPproj);
        Eigen::MatrixXd coeffT = ITHACAutilities::getCoeffs(example.Tfield, TLmodes,
                                 List_of_modes(i, 0) + example.liftfieldT.size());
        PtrList<volVectorField> rec_fieldU = ITHACAutilities::reconstructFromCoeff(
                ULmodes, coeffU, List_of_modes(i,
                                               0) + example.liftfield.size() + NmodesSUPproj);
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
    example.projectSUP("./Matrices", NmodesUproj, NmodesPrghproj, NmodesTproj,
                       NmodesSUPproj);
    // Resize the modes for projection
    example.Tmodes.resize(NmodesTproj);
    example.Umodes.resize(NmodesUproj);
    example.Pmodes.resize(NmodesPproj);
    example.Prghmodes.resize(NmodesPrghproj);
    // Online part
    ReducedUnsteadyBB reduced(example);
    // Set values of the online solve
    reduced.nu = 0.00001;
    reduced.Pr = 0.71;
    reduced.tstart = 0.0;
    reduced.finalTime = 5.0;
    reduced.dt = 0.002;
    // Set the online velocity
    Eigen::MatrixXd vel_now_BC(1, 1);
    vel_now_BC(0, 0) = 0.0157;

    // Set the online temperature BC and solve reduced model
    for (label k = 0; k < par_on_BC.rows(); k++)
    {
        Eigen::MatrixXd temp_now_BC(3, 1);
        temp_now_BC(0, 0) = par_on_BC(k, 0);
        temp_now_BC(1, 0) = par_on_BC(k, 1);
        temp_now_BC(2, 0) = par_on_BC(k, 2);
        reduced.solveOnline_sup(temp_now_BC, vel_now_BC, k, par_on_BC.rows());
        reduced.reconstruct_sup("./ITHACAoutput/ReconstructionSUP/", 5);
    }

    // Performing full order simulation for second parameter set of temp_now_BC
    tutorial11 HFonline2(argc, argv);
    HFonline2.Pnumber = 1;
    HFonline2.Tnumber = 1;
    HFonline2.setParameters();
    HFonline2.mu_range(0, 0) = 0.00001;
    HFonline2.mu_range(0, 1) = 0.00001;
    HFonline2.genEquiPar();
    HFonline2.inletIndexT.resize(3, 1);
    HFonline2.inletIndexT << 1, 2, 3;
    HFonline2.inletIndex.resize(1, 2);
    HFonline2.inletIndex << 3, 0;
    HFonline2.startTime = 0.0;
    HFonline2.finalTime = 5.0;
    HFonline2.timeStep = 0.002;
    HFonline2.writeEvery = 0.01;
    // Reconstruct the online solution
    HFonline2.onlineSolveFull(par_on_BC, 1,
                              "./ITHACAoutput/high_fidelity_online2");
    // Reading in the high-fidelity solutions for the parameter set
    // for which the offline solve has been performed
    example.onlineSolveRead("./ITHACAoutput/Offline/");
    // Reading in the high-fidelity solutions for the second parameter set
    example.onlineSolveRead("./ITHACAoutput/high_fidelity_online2/");
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

/// \dir 11UnsteadyBB Folder of the turorial 11
/// \file
/// \brief Implementation of tutorial 11 for an unsteady Buoyant Boussinesq problem for an open ended cavity

/// \example 11UnsteadyBBOpen.C
/// \section intro_unsteadyBB Introduction to tutorial 11
/// In this tutorial an unsteady Buoyant Boussinesq (BB) 2D problem with paramerized temperature boundary conditions is implemented.
/// The physical problem represents an differentially heated cavity with an inlet and outlet. A uniform temperature is
/// set to the left (hot) and the right (cold) sides of the cavity while the other surfaces are set to adiabatic.
/// The cavity aspect ratio is 1.0 and the ratio between the height of the inlet/outlet and the cavity is 0.1.
/// The flow is laminar with Re = 10 and Ri = 1.2*10^4. The working fluid is air with Pr = 0.72.
/// The ambient temperature is 300 K. The hot wall, Th, has a temperature of 310 K, while the cold wall, Tc,
/// and inlet BC's are set to 300 K. The inlet velocity is 0.0157 m/s.
///
/// The following image illustrates the simulated system at time = 5 s:
/// \image html open_ended_cavity.png
///
/// \section code04 A detailed look into the code
///
/// In this section we explain the main steps necessary to construct the tutorial N°11
///
/// \subsection header ITHACA-FV header files
///
/// First of all let's have a look at the header files that need to be included and what they are responsible for.
///
/// The header files of ITHACA-FV necessary for this tutorial are: <UnsteadyBB.H> for the full order unsteady BB problem,
/// <ITHACAPOD.H> for the POD decomposition, <ReducedUnsteadyBB.H> for the construction of the reduced order problem,
/// and finally <ITHACAstream.H> for some ITHACA input-output operations.
///
/// \dontinclude 11UnsteadyBBOpen.C
/// \skip UnsteadyBB
/// \until ITHACAstream
///
/// \subsection classtutorial11 Definition of the tutorial11 class
///
/// We define the tutorial11 class as a child of the UnsteadyBB class.
/// The constructor is defined with members that are the fields need to be manipulated
/// during the resolution of the full order problem using pimpleFoam. Such fields are
/// also initialized with the same initial conditions in the solver.
/// \skipline tutorial11
/// \until {}
///
/// Inside the tutorial11 class we define the offlineSolve method according to the
/// specific parametrized problem that needs to be solved. If the offline solve has
/// been previously performed then the method just reads the existing snapshots from the Offline directory.
/// Otherwise it loops over all the parameters, changes the values of the temperature boundary conditions and
/// then performs the offline solve.
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
/// First we construct the object "example" of type tutorial11:
///
/// \skipline example
///
/// The parameter sets for the temperature BC are read in from file, both for the offline and online phase.
/// Each row corresponds to a parameter set and each column to a specific boundary.
///
/// \skipline par_offline_BC
/// \until Eigen
/// \skipline Eigen
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
/// info take a look at tutorial04
///
/// \skipline example.Pnumber
/// \until example.genEquiPar()
///
/// After that we set the boundaries where we have a non homogeneous BC for temperature. Patch 1 corresponds to the
/// hot boundary and Patch 2 to the cold boundary and Patch 3 to the inlet boundary.
///
/// \skipline example.inletIndexT
/// \until example.inletIndexT << 1, 2, 3;
///
/// At inlet boundaries (patch 3) where we have the non homogeneous BC (in x-direction, 0):
///
/// \skipline example.inlet
/// \until example.inletIndex << 3, 0;
///
/// And we set the parameters for the time integration, so as to simulate 5 seconds for each
/// simulation, with a step size = 0.002 seconds, nd the data is written every 0.01 seconds, i.e.
///
/// \skipline example.startTime
/// \until example.writeEvery
///
/// Now we are ready to perform the offline stage:
///
/// \skipline example.offlineSolve
///
/// In order to search and compute the velocity lifting function (which should be a step function of value
/// equals to the unitary inlet velocity), we perform the following:
///
/// \skipline liftSolve()
///
/// And for the temperature lifting function:
///
/// \skipline liftSolveT()
///
/// Then we create homogeneous basis functions for the velocity:
///
/// \skipline computeLift(example.Ufield, example.liftfield, example.Uomfield)
///
/// And temperature:
///
/// \skipline computeLiftT(example.Tfield, example.liftfieldT, example.Tomfield)
///
/// And after that, we obtain the modes for velocity, pressure, shifted pressure and temperature:
///
/// \skipline getModes
/// \skipline getModes
/// \skipline getModes
/// \skipline getModes
///
/// Then the supremizer problem is solved;
///
/// \skipline example.solvesupremizer("modes")
///
/// Before continuiting, the projection error is calculated for a range of number of modes for
/// temperature and velocity.
///
/// \skipline List_of_modes
/// \until }
/// \skipline exportMatrix
///
/// The temperature and velcoity modes are created locally
///
/// \skipline TLmodes
/// \until }
/// \skipline for
/// \until }
/// \skipline ULmodes
/// \until }
/// \skipline for
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
/// \until Prghmodes
///
/// Now that we obtained all the necessary information from the POD decomposition and the reduced matrices,
/// we are now ready to construct the dynamical system for the reduced order model (ROM). We proceed
/// by constructing the object "reduced" of type ReducedUnsteadyBB:
///
/// \skipline ReducedUnsteadyBB
///
/// We can use the new constructed ROM to perform the online procedure, from which we can simulate the
/// problem at new set of parameters. For instance, we solve the problem for 5 seconds of physical time:
///
/// \skipline reduced.nu
/// \until reduced.dt
///
/// And then we can use the new constructed ROM to perform the online procedure, from which we can simulate the
/// problem at new set of parameters. First the online inlet velocity is set:
///
/// \skipline Eigen::MatrixXd
/// \until    vel_now_BC
///
/// Then the online procedure is performed for all temperature boundary sets as defined in the par_online_BC
/// file. And he ROM solution is reconstructed and exported:
///
/// \skipline for
/// \until }
///
/// For the second parameter set the full order solution is calculated
///
/// \skipline // Performing full order simulation for second parameter set of temp_now_BC
/// \until "./ITHACAoutput/HFonline2/");
///
/// and the solutions are read in:
/// \skipline example.onlineSolveRead("./ITHACAoutput/Offline/");
/// \until example.onlineSolveRead("./ITHACAoutput/HFonline2/");
///
/// Finally the L2 error between full and reduced order solutions is calculated
///
/// \skipline Ufield_on
/// \until Tfield_on
/// \skipline L2errorMatrixU
/// \until L2errorMatrixT
/// \skipline "./ITHACAoutput/l2error");
///
///
/// \section plaincode The plain program
/// Here there's the plain code
///
