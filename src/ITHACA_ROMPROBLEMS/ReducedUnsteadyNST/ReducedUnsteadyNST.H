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
Class
    reducedUnsteadyNST
Description
    A reduced problem for the unsteady NS plus the energy equations
SourceFiles
    reducedUnsteadyNST.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the reducedUnsteadyNST class.
/// \dir
/// Directory containing the header and source files for the reducedUnsteadyNST class

#ifndef ReducedUnsteadyNST_H
#define ReducedUnsteadyNST_H
#include "fvCFD.H"
#include "IOmanip.H"
#include "ReducedUnsteadyNS.H"
#include "unsteadyNST.H"
#include <Eigen/Dense>
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>

/// Newton object for the resolution of the reduced problem using a supremizer approach
struct newton_unsteadyNST_sup: public newton_argument<double>
{
    public:
        newton_unsteadyNST_sup() {}
        newton_unsteadyNST_sup(label Nx, label Ny,
                               unsteadyNST& problem): newton_argument<double>(Nx, Ny),
            problem(&problem),
            Nphi_u(problem.NUmodes + problem.liftfield.size() + problem.NSUPmodes),
            Nphi_p(problem.NPmodes),
            N_BC(problem.inletIndex.rows())
        {}

        label operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const;
        label df(const Eigen::VectorXd& x,  Eigen::MatrixXd& fjac) const;
        unsteadyNST* problem;
        label Nphi_u;
        label Nphi_p;
        label N_BC;
        scalar nu;
        scalar dt;
        Eigen::VectorXd y_old;
        Eigen::VectorXd BC;

};

struct newton_unsteadyNST_sup_t: public newton_argument<double>
{
    public:
        newton_unsteadyNST_sup_t() {}
        newton_unsteadyNST_sup_t(label Nx, label Ny,
                                 unsteadyNST& problem): newton_argument<double>(Nx, Ny),
            problem(&problem),
            Nphi_t(problem.NTmodes + problem.liftfieldT.size()),
            N_BC_t(problem.inletIndexT.rows())
        {}

        label operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const;
        label df(const Eigen::VectorXd& x,  Eigen::MatrixXd& fjac) const;
        unsteadyNST* problem;
        label Nphi_t;
        label N_BC_t;
        label Nphi_u;
        scalar nu;
        scalar DT;
        scalar dt;
        Eigen::VectorXd a_tmp;
        Eigen::VectorXd z_old;
        Eigen::VectorXd BC_t;

};

/*---------------------------------------------------------------------------*\
                        Class reducedProblem Declaration
\*---------------------------------------------------------------------------*/

/// Class where it is implemented a reduced problem for the <b> unsteady Navier-stokes weakly coupled with the energy equation </b> problem.
/** In this class are implemented the methods to solve a reduced version of the non-stationary Navier-stokes and energy
equations. This class must be used together with the unsteadyNST class */
class reducedUnsteadyNST: public reducedUnsteadyNS
{
    private:

    public:
        // Constructors
        /// Construct Null
        reducedUnsteadyNST();
        /// Construct Null
        ///
        /// @param[in]  problem  a full order unsteadyNS problem
        /// @param[in]  tipo     Type of pressure stabilisation method you want to use "SUP" for supremizer, "PPE" for pressure Poisson equation.
        ///
        explicit reducedUnsteadyNST(unsteadyNST& problem);

        ~reducedUnsteadyNST() {};

        /// Diffusion Term
        Eigen::MatrixXd B_matrix;

        /// Laplacian Term for temperature trasnport equation
        Eigen::MatrixXd Y_matrix;

        /// Mass Matrix Term for temperature
        Eigen::MatrixXd MT_matrix;

        /// Mass Matrix Term
        Eigen::MatrixXd M_matrix;

        /// Gradient of pressure
        Eigen::MatrixXd K_matrix;

        /// Turbulent thermal term
        List <Eigen::MatrixXd> S_matrix;

        /// Convective Term
        List <Eigen::MatrixXd> C_matrix;

        /// Convective Term for temperature transport equation
        List <Eigen::MatrixXd> Q_matrix;

        /// Divergence of velocity
        Eigen::MatrixXd P_matrix;

        /// Laplacian of pressure
        Eigen::MatrixXd D_matrix;

        /// Divergence of momentum
        List <Eigen::MatrixXd> G_matrix;

        /// List of Eigen matrices to store the online solution for temperature equation
        List < Eigen::MatrixXd> online_solutiont;

        /// List of pointers to store the modes for temperature
        PtrList<volScalarField> Tmodes;

        /// List of pointers to store the snapshots for temperature
        PtrList<volScalarField> Tsnapshots;

        /// Functor object to call the non linear solver sup. approach for velocity case
        newton_unsteadyNST_sup newton_object_sup;

        /// Functor object to call the non linear solver sup. approach for temperature case
        newton_unsteadyNST_sup_t newton_object_sup_t;

        /// Reconstructed temperature field
        PtrList<volScalarField> T_rec;

        /// Number of temperature modes
        label Nphi_t;

        /// Reconstructed temperature field
        PtrList<volScalarField> TREC;

        /// Number of parametrized boundary conditions related to temperature field
        label N_BC_t;

        /// Scalar to store the current time
        scalar time;

        /// Scalar to store the time increment
        scalar dt;

        /// Vector to store the previous temperature solution during the Newton procedure
        Eigen::VectorXd z_old;

        /// Vector to store the temperature solution during the Newton procedure
        Eigen::VectorXd z;

        /// DT
        scalar DT;

        /// Scalar to store the final time if the online simulation
        scalar finalTime;

        /// Scalar to store the final time if the online simulation
        scalar tstart;

        ///pointer to the FOM problem
        unsteadyNST* problem;

        // Functions

        /// Method to perform an online solve using a supremizer stabilisation method
        ///
        /// @param[in]  vel_now   The vector of online velocity. It is defined in
        /// with an Eigen::MatrixXd and must have one col and as many rows as the number
        /// of parametrized boundary conditions.
        /// @param[in]  temp_now   The scalar of online temperature. It is defined in
        /// with an Eigen::MatrixXd and must have one col and as many rows as the number
        /// of parametrized boundary conditions.
        /// @param[in]  startSnap The first snapshot taken from the offline snapshots
        /// and used to get the reduced initial condition.
        ///
        void solveOnline_sup(Eigen::MatrixXd& vel_now, Eigen::MatrixXd& temp_now,
                             label startSnap = 0);

        /// Method to reconstruct a solution from an online solve with a supremizer stabilisation technique.
        /// stabilisation method
        ///
        /// @param[in]  folder      The folder where you want to store the results (default is "./ITHACAOutput/online_rec").
        /// @param[in]  printevery  Variable to recover only every printevery online solutions default is 1.
        ///
        void reconstruct_sup(fileName folder = "./ITHACAOutput/online_rec",
                             label printevery = 1);

        /// Method to reconstruct a solution from an online solve with a supremizer stabilisation technique for temperature field.
        /// stabilisation method
        ///
        /// @param[in]  folder      The folder where you want to store the results (default is "./ITHACAOutput/online_rec").
        /// @param[in]  printevery  Variable to recover only every printevery online solutions default is 1.
        ///
        void reconstruct_supt(fileName folder = "./ITHACAOutput/online_rec",
                              label printevery = 1);
};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif



