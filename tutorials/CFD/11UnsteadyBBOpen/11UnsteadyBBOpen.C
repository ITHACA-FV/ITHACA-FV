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
                liftfieldT.append((Tlift).clone());
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
    word stabilization = para->ITHACAdict->lookupOrDefault<word>("Stabilization",
                         "supremizer");
    int NmodesUproj   = para->ITHACAdict->lookupOrDefault<int>("NmodesUproj", 5);
    int NmodesPproj   = para->ITHACAdict->lookupOrDefault<int>("NmodesPproj", 5);
    int NmodesPrghproj = para->ITHACAdict->lookupOrDefault<int>("NmodesPrghproj",
                         5);
    int NmodesTproj   = para->ITHACAdict->lookupOrDefault<int>("NmodesTproj", 5);
    int NmodesSUPproj = 0;

    if (stabilization == "supremizer")
    {
        NmodesSUPproj = para->ITHACAdict->lookupOrDefault<int>("NmodesSUPproj", 5);
    }

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
    ITHACAPOD::getModes(
        example.Uomfield, example.Umodes, example._U().name(),
        example.podex, 0, 0, NmodesOut, false);
    ITHACAPOD::getModes(
        example.Pfield, example.Pmodes, example._p().name(),
        example.podex, 0, 0, NmodesOut, false);
    ITHACAPOD::getModes(
        example.Prghfield, example.Prghmodes, example._p_rgh().name(),
        example.podex, 0, 0, NmodesOut, false);
    ITHACAPOD::getModes(
        example.Tomfield, example.Tmodes, example._T().name(),
        example.podex, 0, 0, NmodesOut, false);

    // Solve the supremizer problem
    if (stabilization == "supremizer")
    {
        example.solvesupremizer("modes");
    }

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
        TLmodes.append((example.liftfieldT[k]).clone());
    }

    for (label k = 0; k < List_of_modes.size(); k++)
    {
        TLmodes.append((example.Tmodes[k]).clone());
    }

    // Create locally the velocity modes
    PtrList<volVectorField> ULmodes;

    for (label k = 0; k < example.liftfield.size(); k++)
    {
        ULmodes.append((example.liftfield[k]).clone());
    }

    for (label k = 0; k < List_of_modes.size(); k++)
    {
        ULmodes.append((example.Umodes[k]).clone());
    }

    if (stabilization == "supremizer")
    {
        for (label k = 0; k < NmodesSUPproj; k++)
        {
            ULmodes.append((example.supmodes[k]).clone());
        }
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
    if (stabilization == "supremizer")
    {
        example.projectSUP("./Matrices", NmodesUproj, NmodesPrghproj, NmodesTproj,
                           NmodesSUPproj);
    }
    else if (stabilization == "PPE")
    {
        example.projectPPE("./Matrices", NmodesUproj, NmodesPrghproj, NmodesTproj,
                           NmodesSUPproj);
    }
    else
    {
        // TODO: warning message
    }

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

        if (stabilization == "supremizer")
        {
            reduced.solveOnline_sup(temp_now_BC, vel_now_BC, k, par_on_BC.rows());
            reduced.reconstruct_sup("./ITHACAoutput/ReconstructionSUP/", 5);
        }
        else if (stabilization == "PPE")
        {
            reduced.solveOnline_PPE(temp_now_BC, vel_now_BC, k, par_on_BC.rows());
            reduced.reconstruct_PPE("./ITHACAoutput/ReconstructionPPE/", 5);
        }
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
    return 0;
}

