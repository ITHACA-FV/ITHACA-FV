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
    Example of univariate PODI Problem
SourceFiles
    05PODI_univ.C
\*---------------------------------------------------------------------------*/

#include "unsteadyNS.H"
#include "ITHACAPOD.H"
#include "ITHACAstream.H"
#include "forces.H"
#include "IOmanip.H"

class tutorial05 : public unsteadyNS
{
    public:
        /// Constructor
        explicit tutorial05(int argc, char* argv[])
            :
            unsteadyNS(argc, argv),
            U(_U()),
            p(_p())
        {}

        /// Velocity field
        volVectorField& U;
        /// Pressure field
        volScalarField& p;

        /// Perform an Offline solve
        void offlineSolve()
        {
            Vector<double> Uinl(0, 0, 0);
            List<scalar> mu_now(1);
            mu_now[0] = 1;

            // if the offline solution is already performed read the fields
            if (offline)
            {
                ITHACAstream::read_fields(Ufield, U, "./ITHACAoutput/Offline/");
                ITHACAstream::read_fields(Pfield, p, "./ITHACAoutput/Offline/");
            }
            else
            {
                assignIF(U, Uinl);
                truthSolve(mu_now);
            }
        }
};

int main(int argc, char* argv[])
{
    // Construct the tutorial object
    tutorial05 example(argc, argv);
    // Read some parameters from file
    ITHACAparameters para;
    int NmodesUout = para.ITHACAdict->lookupOrDefault<int>("NmodesUout", 15);
    int NmodesPout = para.ITHACAdict->lookupOrDefault<int>("NmodesPout", 15);
    int NmodesUproj = para.ITHACAdict->lookupOrDefault<int>("NmodesUproj", 10);
    int NmodesPproj = para.ITHACAdict->lookupOrDefault<int>("NmodesPproj", 10);
    example.startTime = 1.0;
    example.finalTime = 20;
    example.timeStep = 0.005;
    example.writeEvery = 0.5;
    // Perform the offline solve
    example.offlineSolve();
    // original fields to exploit for interpolation
    PtrList <volVectorField> U_1;
    PtrList <volScalarField> P_1;

    for (int k = 0; k < example.Ufield.size(); k += 2)
    {
        U_1.append(example.Ufield[k]);
        P_1.append(example.Pfield[k]);
    }

    // expected fields to check against
    PtrList <volVectorField> U_2;
    PtrList <volScalarField> P_2;

    for (int k = 1; k < example.Ufield.size(); k += 2)
    {
        U_2.append(example.Ufield[k]);
        P_2.append(example.Pfield[k]);
    }

    // Perform POD decomposition on given fields (i.e. U_1 & P_1, assuming we do not know U_2 & P_2)
    // then store the first 15 modes
    PtrList <volVectorField> Umodes_1;
    PtrList <volScalarField> Pmodes_1;
    ITHACAPOD::getModes(U_1, Umodes_1, example.podex, 0, 0, NmodesUproj);
    ITHACAPOD::getModes(P_1, Pmodes_1, example.podex, 0, 0, NmodesPproj);
    /////////////////////
    // ITHACAUTILITIES //
    /////////////////////
    // Get coefficients (Orthonormal) -- given
    Eigen::MatrixXd coeff_U = ITHACAutilities::get_coeffs_ortho(U_1, Umodes_1);
    Eigen::MatrixXd coeff_P = ITHACAutilities::get_coeffs_ortho(P_1, Pmodes_1);
    // Get coefficients (Orthonormal) -- ref
    Eigen::MatrixXd coeff_U_ref = ITHACAutilities::get_coeffs_ortho(U_2, Umodes_1);
    Eigen::MatrixXd coeff_P_ref = ITHACAutilities::get_coeffs_ortho(P_2, Pmodes_1);
    // Vector for the given time snapshots (20 snapshots)
    Eigen::VectorXd time(20);

    for (int i = 0; i < 20; i++)
    {
        time(i) = i + 1;
    }

    // Vector for the expected (interpolated) snapshots (19 snapshots)
    Eigen::VectorXd time_interp(19);

    for (int i = 0; i < 19; i++)
    {
        time_interp(i) = i + 1.5;
    }

    Eigen::MatrixXd coeff_U_spl = ITHACAutilities::bsplineInterp_univar(time,
                                  coeff_U, time_interp, 3);
    Eigen::MatrixXd coeff_P_spl = ITHACAutilities::bsplineInterp_univar(time,
                                  coeff_P, time_interp, 3);
    // Export U & P coeffs of snapshots for the modes 1 to 15
    ITHACAstream::exportMatrix(coeff_U, "coeff_U_given", "eigen",
                               "ITHACAoutput/coefficients");
    ITHACAstream::exportMatrix(coeff_P, "coeff_P_given", "eigen",
                               "ITHACAoutput/coefficients");
    ITHACAstream::exportMatrix(coeff_U_ref, "coeff_U_ref", "eigen",
                               "ITHACAoutput/coefficients");
    ITHACAstream::exportMatrix(coeff_P_ref, "coeff_P_ref", "eigen",
                               "ITHACAoutput/coefficients");
    ITHACAstream::exportMatrix(coeff_U_spl, "coeff_U_spl", "eigen",
                               "ITHACAoutput/coefficients");
    ITHACAstream::exportMatrix(coeff_P_spl, "coeff_P_spl", "eigen",
                               "ITHACAoutput/coefficients");
    // Compute relative interpolation error in the L2 norm
    double interp_error_U = (coeff_U_ref - coeff_U_spl).norm() / coeff_U_ref.norm()
                            * 100.;
    double interp_error_P = (coeff_P_ref - coeff_P_spl).norm() / coeff_P_ref.norm()
                            * 100.;
    // Write out interpolation error
    std::ofstream file("ITHACAoutput/interp_error.txt");
    file << "# Relative error of the Velocity coefficients in L2 norm = " <<
         interp_error_U << "\n";
    file << "# Relative error of the Pressure coefficients in L2 norm = " <<
         interp_error_P << "\n";
    file.close();
    // Compute interpolated Velocity Snapshots using the first 15 POD modes
    PtrList < volVectorField > U_recon;
    example.reconstruct_from_matrix(U_recon, Umodes_1, 15, coeff_U_spl);
    // Compute interpolated Pressure Snapshots using the first 15 POD modes
    PtrList < volScalarField > P_recon;
    example.reconstruct_from_matrix(P_recon, Pmodes_1, 15, coeff_P_spl);
    // Export unknown fields (true)
    ITHACAstream::exportFields(U_2, "ITHACAoutput/interp_fields", "U");
    ITHACAstream::exportFields(P_2, "ITHACAoutput/interp_fields", "P");
    // Export unknown fields (reconstructed using interpolated coefficients)
    ITHACAstream::exportFields(U_recon, "ITHACAoutput/interp_fields", "U_recon");
    ITHACAstream::exportFields(P_recon, "ITHACAoutput/interp_fields", "P_recon");
    // Error Fields
    PtrList <volVectorField> U_error;
    PtrList <volScalarField> P_error;

    for (int i = 0; i < U_recon.size(); i++)
    {
        U_error.append(U_2[i] - U_recon[i]);
        P_error.append(P_2[i] - P_recon[i]);
    }

    ITHACAstream::exportFields(U_error, "ITHACAoutput/interp_fields", "U_error");
    ITHACAstream::exportFields(P_error, "ITHACAoutput/interp_fields", "P_error");
    // Compute fields relative error in the L2 norm
    Eigen::MatrixXd error_norm_U = ITHACAutilities::error_listfields(U_2, U_recon);
    Eigen::MatrixXd error_norm_P = ITHACAutilities::error_listfields(P_2, P_recon);
    // Write out the L2 relative error for the Velocity Fields
    std::ofstream file2("ITHACAoutput/Ufield_L2_error.txt");
    file2 << error_norm_U << "\n";
    file2.close();
    // \Write out the L2 relative error for the Pressure Fields
    std::ofstream file3("ITHACAoutput/Pfield_L2_error.txt");
    file3 << error_norm_P << "\n";
    file3.close();
    exit(0);
}



/// \dir 05PODI_univ Folder of the turorial 5
/// \file
/// \brief Implementation of tutorial 5 for POD with univariate interpolation for an unsteady Navier-Stokes problem

/// \example 05PODI_univ.C
/// \section intro_PODI_univ Introduction to tutorial 5
/// In this tutorial we implement a cubic BSpline *univariate* interpolation of the POD temporal coefficients for an unsteady Navier-Stokes 2D problem.
/// Similar to tutorial 4, the physical problem represents an incompressible flow passing around a very long cylinder. The simulation domain is rectangular
/// with spatial bounds of [-4, 30], and [-5, 5] in the X and Y directions, respectively. The cylinder has a radius of
/// 0.5 unit length and is located at the origin. The system has a prescribed uniform inlet velocity of 1 m/s which is constant through the whole simulation.
///
/// \section code05 A detailed look into the code
///
/// In this section we explain the main steps necessary to construct the tutorial N°5
///
/// \subsection header ITHACA-FV header files
///
/// First of all let's have a look at the header files that need to be included and what they are responsible for.
///
/// The header files of ITHACA-FV necessary for this tutorial are: <unsteadyNS.H> for the full order unsteady NS problem,
/// <ITHACAPOD.H> for the POD decomposition, <reducedUnsteadyNS.H> for the construction of the reduced order problem,
/// and finally <ITHACAstream.H> for some ITHACA input-output operations.
///
/// \dontinclude 05PODI_univ.C
/// \skip unsteadyNS
/// \until ITHACAstream
///
/// \subsection classtutorial05 Definition of the tutorial05 class
///
/// We define the tutorial05 class as a child of the unsteadyNS class.
/// The constructor is defined with members that are the fields need to be manipulated
/// during the resolution of the full order problem using pimpleFoam. Such fields are
/// also initialized with the same initial conditions in the solver.
/// \skipline tutorial05
/// \until {}
///
/// Inside the tutorial05 class we define the offlineSolve method according to the
/// specific parametrized problem that needs to be solved. If the offline solve has
/// been previously performed then the method just reads the existing snapshots from the Offline directory.
/// Otherwise it performs the offline solve.
///
/// \skipline offlineSolve
/// \until }
/// \skipline else
/// \until }
/// \skipline }
/// \skipline }
///
/// \subsection main Definition of the main function
///
/// In this section we show the definition of the main function.
/// First we construct the object "example" of type tutorial05:
///
/// \skipline example
///
/// Then we parse the ITHACAdict file to determine the number of modes
/// to be written out and also the ones to be used for projection of
/// the velocity and pressure:
/// \skipline ITHACAparameters
/// \until NmodesPproj
///
/// we note that a default value can be assigned in case the parser did
/// not find the corresponding string in the ITHACAdict file.
///
/// After that we set the parameters for the time integration, so as to simulate 20 seconds
/// with a step size = 0.005 seconds, and the data are dumped every 0.5 seconds, i.e.
///
/// \skipline example.startTime
/// \until example.writeEvery
///
/// and we are ready to perform the offline stage:
///
/// \skipline Solve()
///
/// Now let us assume we are given only the field variables at time snapshots t={1, 2, 3, ..., 20}
/// and we need to estimate the fields at t={1.5, 2.5, 3.5, ..., 19.5}. This can be achieved by interpolating
/// the POD temporal coefficients of the given snapshots at the time instances of the required snapshots.
/// Before we proceed, we have first to separate the offline solutions (39 snapshots) into a given set (20 snapshots)
/// and a required or unknown set (19 snapshots) by using PtrList:
///
/// \skipline PtrList <volVectorField> U_1;
/// \until }
/// \skipline PtrList <volVectorField> U_2;
/// \until }
///
/// After that, we perform POD decomposition on the set of the given fields (with the 20 snapshots)
/// and we obtain the energetic spatial modes of the Velocity and the Pressure:
///
/// \skipline Umodes_1
/// \until getModes(P_1
///
/// Having the POD modes of the given dataset, we can compute the temporal coefficients (amplitudes)
/// for the given snapshots. Also we compute the temporal coefficients of the required snapshots so as
/// to use as reference to compare against the interpolated coefficients.
///
/// \skipline coeff_U
/// \until coeff_P
///
/// \skipline coeff_U_ref
/// \until coeff_P_ref
///
/// We note that due to the orthonormality of the spatial modes, the temporal coefficients
/// are computed just by performing an inner product in the Hilbert space between the POD modes and the field snapshots.
///
/// Now we proceed to the coefficients interpolation. First we create a vector for t={1, 2, 3, ..., 20} to be used as the original input,
/// and a vector for t={1.5, 2.5, 3.5, ..., 19.5} to be used as the interpolation input. Then we perform the BSpline interpolation with degree=3:
///
/// \skipline VectorXd time
/// \until }
/// \skipline VectorXd time_interp
/// \until }
/// \skipline coeff_U_spl
/// \until coeff_P_spl
///
/// The plots of the coefficients for the first 6 POD modes are shown in the following picture
/// \image html coeff.png
///
/// The coefficients matrices that correspond to the given data, reference data, and the interpolated data are then exported:
///
/// \skipline exportMatrix
/// \until coeff_P_spl
///
/// Then we compute the interpolation relative error in the L2-norm and we write out the data, i.e.
///
/// \skipline interp_error_U
/// \until interp_error_P
/// \skipline ofstream
/// \until close()
///
/// Now we reconstruct the Velocity and Pressure fields for the snapshots t={1.5, 2.5, ..., 19.5}
/// by projecting the interpolated coefficients onto the first 15 modes, then we export such fields
/// into a new directory <interp_fields>, i.e.
///
/// \skipline U_recon
/// \until example
/// \skipline P_recon
/// \until example
/// \skipline U_2
/// \until P_2
/// \skipline U_recon
/// \until P_recon
///
/// We can also compute the error fields between the reconstructed and the reference ones:
///
/// \skipline U_error
/// \until "P_error"
///
/// Finally we compute the error fields in the relative L2 norm, then we write out the results, i.e.
///
/// \skipline error_norm_U
/// \until error_norm_P
/// \skipline file2
/// \until close()
/// \skipline file3
/// \until close()
///
/// \section plaincode The plain program
/// Here there's the plain code
///
