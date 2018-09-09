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
    Example of bivariate PODI Problem
SourceFiles
    06PODI_biv.C
\*---------------------------------------------------------------------------*/

#include <vector>

#include "unsteadyNS.H"
#include "ITHACAPOD.H"
#include "ITHACAstream.H"
#include "forces.H"
#include "IOmanip.H"

class tutorial06 : public unsteadyNS
{
    public:
        /// Constructor
        explicit tutorial06(int argc, char* argv[])
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
            Vector<double> inl(5, 0, 0);
            // Inlet BC: batch index = 0
            label BCind = 0;
            assignBC(U, BCind, inl);
            List<scalar> mu_now(1);

            // If the offline solution is already performed read the fields
            if (offline)
            {
                ITHACAstream::read_fields(Ufield, U, "./ITHACAoutput/Offline/");
                ITHACAstream::read_fields(Pfield, p, "./ITHACAoutput/Offline/");
            }
            else
            {
                Vector<double> Uinit(0, 0, 0);

                for (label i = 0; i < mu.rows(); i++)
                {
                    change_viscosity(mu(i, 0));
                    assignIF(U, Uinit);
                    mu_now[0] = mu(i, 0);
                    truthSolve(mu_now);
                }
            }
        }
};

int main(int argc, char* argv[])
{
    // Construct the tutorial object
    tutorial06 example(argc, argv);
    // Read some parameters from file
    ITHACAparameters para;
    int NmodesUout = para.ITHACAdict->lookupOrDefault<int>("NmodesUout", 15);
    int NmodesPout = para.ITHACAdict->lookupOrDefault<int>("NmodesPout", 15);
    int NmodesUproj = para.ITHACAdict->lookupOrDefault<int>("NmodesUproj", 10);
    int NmodesPproj = para.ITHACAdict->lookupOrDefault<int>("NmodesPproj", 10);
    example.startTime = 10.0;
    example.finalTime = 20.0;
    example.timeStep = 0.0025;
    example.writeEvery = 0.5;
    // Read the par file where the parameters are stored
    word filename("./par");
    example.mu = ITHACAstream::readMatrix(filename);
    // Perform the offline solve
    example.offlineSolve();
    ////////////////////////////////////////////////////////////////////////////////
    // From the fields extract 2 sets of indices: givenIndices and unknownIndices //
    ////////////////////////////////////////////////////////////////////////////////
    // Nsnapshots for one parametric simulation
    int Nsnapshots = (example.finalTime - example.startTime) / example.writeEvery +
                     1;
    int i;
    // Set of indices of given fields -> given mu given T (25*11 snapshots)
    std::vector< int > givenIndices;
    i = 0;

    while (i < example.Ufield.size())
    {
        for (int j = 0; j < Nsnapshots; j++)
        {
            // Select only time snapshots with local even indices
            if (j % 2 == 0)
            {
                givenIndices.push_back(i);
            }

            // Iterate snapshots indices
            i += 1;
        }

        // Skip all time snapshots corresponding to next mu,
        // so that we account only to every other parametric simulation
        i += Nsnapshots;
    }

    // Set of indices of unknown fields -> unknown mu unknown T (25*10 snapshots)
    std::vector< int > unknownIndices;
    i = 0;
    i += Nsnapshots;

    while (i < example.Ufield.size())
    {
        for (int j = 0; j < Nsnapshots; j++)
        {
            // Select only time snapshots with local odd indices
            if (j % 2 == 1)
            {
                unknownIndices.push_back(i);
            }

            // Iterate snapshots indices
            i += 1;
        }

        // Skip all time snapshots corresponding to next mu,
        // so that we account only to every other parametric simulation
        i += Nsnapshots;
    }

    // Fields to be regarded as the given data (given mu given T)
    PtrList <volVectorField> U_1;
    PtrList <volScalarField> P_1;
    // Fields to be regarded as the unknown data (unknown mu unknown T)
    PtrList <volVectorField> U_2;
    PtrList <volScalarField> P_2;
    bool check_1, check_2;

    for (int k = 0; k < example.Ufield.size(); k++)
    {
        // True if the current iterate (i.e. snapshot index) is contained in vector "givenIndices"
        check_1 = std::find(givenIndices.begin(), givenIndices.end(),
                            k) != givenIndices.end();
        // True if the current iterate (i.e. snapshot index) is contained in vector "unknownIndices"
        check_2 = std::find(unknownIndices.begin(), unknownIndices.end(),
                            k) != unknownIndices.end();

        if (check_1 == true)
        {
            U_1.append(example.Ufield[k]);
            P_1.append(example.Pfield[k]);
        }

        if (check_2 == true)
        {
            U_2.append(example.Ufield[k]);
            P_2.append(example.Pfield[k]);
        }
    }

    // Perform POD decomposition on given fields (i.e. U_1 & P_1, assuming we do not know all other fields)
    // then store the first 15 modes
    PtrList <volVectorField> Umodes_1;
    PtrList <volScalarField> Pmodes_1;
    ITHACAPOD::getModes(U_1, Umodes_1, example.podex, 0, 0, NmodesUproj);
    ITHACAPOD::getModes(P_1, Pmodes_1, example.podex, 0, 0, NmodesPproj);
    /////////////////////
    // ITHACAUTILITIES //
    /////////////////////
    // Get coefficients (Orthonormal) -- given
    Eigen::MatrixXd coeff_U_given = ITHACAutilities::get_coeffs_ortho(U_1,
                                    Umodes_1);
    Eigen::MatrixXd coeff_P_given = ITHACAutilities::get_coeffs_ortho(P_1,
                                    Pmodes_1);
    // Get coefficients (Orthonormal) -- ref (unknown)
    Eigen::MatrixXd coeff_U_ref = ITHACAutilities::get_coeffs_ortho(U_2, Umodes_1);
    Eigen::MatrixXd coeff_P_ref = ITHACAutilities::get_coeffs_ortho(P_2, Pmodes_1);
    int Nsnapshots_given = (Nsnapshots + 1) / 2;
    int Nsnapshots_unknown = Nsnapshots - Nsnapshots_given;
    int mu_given = givenIndices.size() / Nsnapshots_given;
    int mu_unknown = unknownIndices.size() / Nsnapshots_unknown;
    // Vector for the given time snapshots
    Eigen::VectorXd time(Nsnapshots_given);

    for (int i = 0; i < Nsnapshots_given; i++)
    {
        time(i) = i + example.startTime;
    }

    // Vector for the time snapshots to interpolate
    Eigen::VectorXd time_interp(Nsnapshots_unknown);

    for (int i = 0; i < Nsnapshots_unknown; i++)
    {
        time_interp(i) = i + example.startTime + 0.5;
    }

    // Vector for the given mu
    Eigen::VectorXd par(mu_given);

    for (int i = 0; i < mu_given; i++)
    {
        par(i) = example.mu(2 * i, 0);
    }

    // Vector for the mu to interpolate
    Eigen::VectorXd par_interp(mu_unknown);

    for (int i = 0; i < mu_unknown; i++)
    {
        par_interp(i) = example.mu(2 * i + 1, 0);
    }

    Eigen::MatrixXd coeff_U_spl = ITHACAutilities::bsplineInterp_bivar(time, par,
                                  coeff_U_given, time_interp, par_interp, 3);
    Eigen::MatrixXd coeff_P_spl = ITHACAutilities::bsplineInterp_bivar(time, par,
                                  coeff_P_given, time_interp, par_interp, 3);
    // Export all coeffs "flattened" matrices of U & P snapshots, for the modes 1 to 15
    ITHACAstream::exportMatrix(coeff_U_given, "coeff_U_given", "eigen",
                               "ITHACAoutput/coefficients");
    ITHACAstream::exportMatrix(coeff_P_given, "coeff_P_given", "eigen",
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
    file << "\ncoeff U & P (unknown mu, unknown T) \n";
    file << "# Relative error of the Velocity coefficients in L2 norm = " <<
         interp_error_U << "\n";
    file << "# Relative error of the Pressure coefficients in L2 norm = " <<
         interp_error_P << "\n";
    file.close();
    //////////////////////////////////////////////
    // Reconstruct unknown fields (in mu and T) //
    //////////////////////////////////////////////
    // Compute interpolated Velocity Snapshots using the first 15 POD modes
    PtrList < volVectorField > U_recon;
    example.reconstruct_from_matrix(U_recon, Umodes_1, 15, coeff_U_spl);
    // Compute interpolated Pressure Snapshots using the first 15 POD modes
    PtrList < volScalarField > P_recon;
    example.reconstruct_from_matrix(P_recon, Pmodes_1, 15, coeff_P_spl);
    // Export reference fields
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
    file2 << "Error norm -- unknown mu unknown T\n";
    file2 << error_norm_U << "\n";
    file2.close();
    // \Write out the L2 relative error for the Pressure Fields
    std::ofstream file3("ITHACAoutput/Pfield_L2_error.txt");
    file3 << "Error norm -- unknown mu unknown T\n";
    file3 << error_norm_P << "\n";
    file3.close();
    exit(0);
}



/// \dir 06PODI_biv Folder of the turorial 6
/// \file
/// \brief Implementation of tutorial 6 for POD with bivariate interpolation for an unsteady Navier-Stokes problem

/// \example 06PODI_biv.C
/// \section intro_PODI_biv Introduction to tutorial 6
/// In this tutorial we implement a cubic BSpline *bivariate* interpolation of the POD 2D coefficients (time and viscosity) for an unsteady Navier-Stokes 2D problem.
/// Similar to tutorial 4, the physical problem represents an incompressible flow passing around a very long cylinder. The simulation domain is rectangular
/// with spatial bounds of [-4, 30], and [-5, 5] in the X and Y directions, respectively. The cylinder has a radius of
/// 0.5 unit length and is located at the origin. The system has a prescribed uniform inlet velocity of 5 m/s which is constant through the whole simulation.
///
/// The simulations are parametrized in both time and kinematic viscosity. Snapshots are sampled at time instances with a range (10s, 20s) and step of 1s. At the
/// same time the kinematic viscosity is also parametrized to correspond to a Reynolds number range (155, 495) with step of 20. Such snapshots are utilized for the
/// POD decomposition, from which a bivariate BSpline interpolation is performed to estimate the coefficients at new instances. The following image illustrates
/// the snapshots sampled from the parametrized simulation so as to be regarded as the "given" snapshots, whereas the "unknown" snapshots are the ones to be determined
/// using the bivariate interpolation.
/// \image html grid.png
///
/// We note that in this tutorial we performed the full order simulation on all the grid points (i.e. including both the given and unknown snapshots) so that
/// we have a reference data to be used as reference, when compared with the reconstructed snapshots from the bivariate interpolation.
///
/// \section code06 A detailed look into the code
///
/// In this section we explain the main steps necessary to construct the tutorial N°6
///
/// \subsection header ITHACA-FV header files
///
/// First of all let's have a look at the header files that need to be included and what they are responsible for.
///
/// The header files of ITHACA-FV necessary for this tutorial are: <unsteadyNS.H> for the full order unsteady NS problem,
/// <ITHACAPOD.H> for the POD decomposition, <reducedUnsteadyNS.H> for the construction of the reduced order problem,
/// and finally <ITHACAstream.H> for some ITHACA input-output operations.
///
/// \dontinclude 06PODI_biv.C
/// \skip unsteadyNS
/// \until ITHACAstream
///
/// \subsection classtutorial06 Definition of the tutorial06 class
///
/// We define the tutorial06 class as a child of the unsteadyNS class.
/// The constructor is defined with members that are the fields need to be manipulated
/// during the resolution of the full order problem using pimpleFoam. Such fields are
/// also initialized with the same initial conditions in the solver.
/// \skipline tutorial06
/// \until {}
///
/// Inside the tutorial06 class we define the offlineSolve method according to the
/// specific parametrized problem that needs to be solved. The inlet boundary is assigned with the vector (5, 0, 0) for all the simulations.
/// The method also checks if the offline solve has been previously performed so it just reads the existing snapshots from the Offline directory.
/// Otherwise it performs the offline solve by initializing the internal field, iterating and changing the parameter (i.e. kinematic viscosity of the
/// system), then solving the time dependant problem for each parameter.
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
/// First we construct the object "example" of type tutorial06:
///
/// \skipline example
///
/// Then we parse the ITHACAdict file to determine the number of modes
/// to be written out, and also the ones to be used for the projection of
/// the velocity and pressure fields onto those modes:
/// \skipline ITHACAparameters
/// \until NmodesPproj
///
/// we note that a default value can be assigned in case the parser did
/// not find the corresponding string in the ITHACAdict file.
///
/// After that we set the parameters for the time integration, so as to simulate the system from 10 seconds into 20 seconds of physical time,
/// where the step size = 0.0025 seconds, and the data are dumped every 0.5 seconds, i.e.
///
/// \skipline example.startTime
/// \until example.writeEvery
///
/// and we set the parameter (i.e. kinematic viscosity) with values that are read from the file "par", i.e.
///
/// \skipline ./par
/// \until example.mu
///
/// Now we are ready to perform the offline stage:
///
/// \skipline Solve()
///
/// We highlight that such offline procedure solves for all the possible grid points (i.e. given and required snapshots). Now we need to extract the relevant
/// data and divide into two sets: one to be used as the given set for the POD decomposition and interpolation,
/// while the other set to be used as a reference for the validations. This can be achieved by separating the offline solutions
/// (21 time snapshots for each parametrized simulation for 35 different Re, hence a total number of 735 snapshots) into one set
/// contains 198 snapshots to be used as given data, whereas the other set contains 170 snapshots to be used as a reference to validate the reconstructed
/// snapshots at same conditions via interpolation.
///
/// The set of indices that corresponds to the given fields are determined through the following:
///
/// \skipline std::vector< int > givenIndices
/// \until i += Nsnapshots;
/// \skipline }
///
/// whereas the set of indices that corresponds to the reference fields are determined through the following:
///
/// \skipline std::vector< int > unknownIndices
/// \until i += 1;
/// \skipline }
/// \until }
///
/// After that we create two PtrList corresponding to the given fields, namely "U_1" and "P_1", in addition to
/// two PtrList corresponding to the unknown fields, namely "U_2" and "P_2". Those sets of PtrList correspond to
/// the fields that have the same indices as in the previously constructed sets, i.e.
///
/// \skipline PtrList <volVectorField> U_1;
/// \until P_2.append
/// \skipline }
/// \skipline }
///
/// Now we perform POD decomposition on the given fields,
/// and we obtain the first 15 energetic spatial modes of the Velocity and the Pressure:
///
/// \skipline Umodes_1
/// \until getModes(P_1
///
/// Having the POD modes of the given dataset, we can compute the temporal coefficients (amplitudes)
/// for the given snapshots. Also we compute the temporal coefficients of the required snapshots so as
/// to use as reference to compare against the interpolated coefficients.
///
/// \skipline coeff_U_given
/// \until coeff_P_given
///
/// \skipline coeff_U_ref
/// \until coeff_P_ref
///
/// We note that due to the orthonormality of the spatial modes, the temporal coefficients
/// are computed just by performing an inner product in the Hilbert space between the POD modes and the field snapshots.
///
/// Now we proceed to the coefficients interpolation. First we create a vector for t={1, 2, 3, ..., 20} to be used as the original first input,
/// and a vector for t={1.5, 2.5, 3.5, ..., 19.5} to be used as the interpolation first input. Also we create a vector for mu corresponding to
/// Re={155, 175, 195, ..., 495} to be used as the original second input, and a vector for mu corresponding to Re={165, 185, 205, ..., 485}
/// to be used as the interpolation second input, i.e.
/// \skipline time(Nsnapshots_given)
/// \until par_interp(i)
/// \skipline }
///
/// then we perform the BSpline bivariate interpolation with degree=3:
///
/// \skipline coeff_U_spl
/// \until coeff_P_spl
///
/// The coefficients matrices that correspond to the given data, reference data, and the interpolated data are then exported:
///
/// \skipline exportMatrix
/// \until coeff_P_spl
///
/// Then we compute the interpolation relative error in the L2-norm and we write out the data into file, i.e.
///
/// \skipline interp_error_U
/// \until interp_error_P
/// \skipline ofstream
/// \until close()
///
/// Now we reconstruct the Velocity and the Pressure fields by projecting the interpolated coefficients
/// onto the first 15 modes, then we export such fields into a new directory <interp_fields>, i.e.
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
