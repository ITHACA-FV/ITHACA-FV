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
    reducedSteadyNS
Description
    A reduced problem for the stationary NS equations
SourceFiles
    reducedSteadyNS.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the reducedSteadyNS class.
/// \dir
/// Directory containing the header and source files for the reducedSteadyNS class

#ifndef ReducedSteadyNS_H
#define ReducedSteadyNS_H

#include "fvCFD.H"
#include "IOmanip.H"
#include "ReducedProblem.H"
#include "steadyNS.H"
#include "ITHACAutilities.H"
#include "EigenFunctions.H"
#include <Eigen/Eigen>
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>

/// Structure to implement a newton object for a stationary NS problem
struct newton_steadyNS: public newton_argument<double>
{
    public:
        newton_steadyNS() {}
        newton_steadyNS(label Nx, label Ny,
                        steadyNS& problem): newton_argument<double>(Nx,
                                    Ny),
            problem(&problem),
            Nphi_u(problem.NUmodes + problem.liftfield.size() + problem.NSUPmodes),
            Nphi_p(problem.NPmodes),
            N_BC(problem.inletIndex.rows())
        {}

        label operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const;
        label df(const Eigen::VectorXd& x,  Eigen::MatrixXd& fjac) const;

        steadyNS* problem;
        label Nphi_u;
        label Nphi_p;
        label N_BC;
        scalar nu;
        Eigen::MatrixXd tauU;
        Eigen::VectorXd BC;
};


/*---------------------------------------------------------------------------*\
                        Class reducedProblem Declaration
\*---------------------------------------------------------------------------*/

/// Class where it is implemented a reduced problem for the <b> steady Navier-stokes </b> problem.
/** In this class are implemented the methods to solve a reduced version of the stationary Navier-stokes
equations. This class must be used together with the steadyNS class */
class reducedSteadyNS: public reducedProblem
{
    private:

    public:
        // Constructors
        /// Construct Null
        reducedSteadyNS();

        /// Construct Null
        ///
        /// @param[in]  problem  a full order steadyNS problem
        /// @param[in]  tipo     Type of pressure stabilisation method you want to use "SUP" for supremizer, "PPE" for pressure Poisson equation.
        ///
        explicit reducedSteadyNS(steadyNS& problem);

        ~reducedSteadyNS() {};

        /// parameters to be read from the ITHACAdict file
        ITHACAparameters* para;

        /// Online inlet velocity vector
        Eigen::MatrixXd vel_now;

        /// Reduced matrix for tangent forces
        Eigen::MatrixXd fTau;

        /// Reduced matrix for normal forces
        Eigen::MatrixXd fN;

        /// Vector to store the previous solution during the Newton procedure
        Eigen::VectorXd y_old;

        /// Vector to store the solution during the Newton procedure
        Eigen::VectorXd y;

        /// Reduced viscosity in case of parametrized viscosity
        scalar nu;

        /// List of Eigen matrices to store the online solution
        List < Eigen::MatrixXd> online_solution;

        /// List of pointers to store the modes for velocity
        PtrList<volVectorField> Umodes;

        /// List of pointers to store the modes for pressure
        PtrList<volScalarField> Pmodes;

        /// List of pointers to store the snapshots for velocity
        PtrList<volVectorField> Usnapshots;

        /// List of pointers to store the snapshots for pressure
        PtrList<volScalarField> Psnapshots;

        /// Reconstructed pressure field
        PtrList<volScalarField> PREC;

        /// Recontructed velocity field
        PtrList<volVectorField> UREC;

        /// Reconstructed pressure fields list
        PtrList<volScalarField> pRecFields;

        /// Recontructed velocity fields list
        PtrList<volVectorField> uRecFields;

        /// Newton object used to solve the non linear problem
        newton_steadyNS newton_object;

        /// Pointer to the FOM problem
        steadyNS* problem;

        /// Number of velocity modes
        label Nphi_u;

        /// Number of pressure modes
        label Nphi_p;

        /// Number of supremizer modes
        label Nphi_sup;

        /// Number of parametrized boundary conditions
        label N_BC;

        /// Counter to count the online solutions
        label count_online_solve = 1;

        /// Penalty Factor
        Eigen::MatrixXd tauU;

        // Functions

        /// Method to perform an online solve using a PPE stabilisation method
        ///
        /// @param[in]  vel_now  The vector of online velocity. It is defined in
        /// with an Eigen::MatrixXd and must have one col and as many rows as the number
        /// of parametrized boundary conditions.
        ///
        void solveOnline_PPE(Eigen::MatrixXd vel_now);

        /// Method to perform an online solve using a supremizer stabilisation
        /// method
        ///
        /// @param[in]  vel_now  The vector of online velocity. It is defined in
        /// with an Eigen::MatrixXd and must have one col and as many rows as the number
        /// of parametrized boundary conditions.
        ///
        void solveOnline_sup(Eigen::MatrixXd vel_now);

        /// Method to reconstruct a solution from an online solve with a PPE stabilisation technique.
        /// stabilisation method
        ///
        /// @param[in]  folder      The folder where to output the solutions
        /// @param[in]  printevery  The printevery
        ///
        void reconstruct_PPE(fileName folder = "./ITHACAoutput/online_rec",
                             label printevery = 1);


        ///
        /// Method to reconstruct the solutions from an online solve.
        ///
        /// @param[in]  exportFields  A boolean variable which determines whether to export fields or not
        /// @param[in]  folder        The folder where to output the solutions in case on wants to
        /// @param[in]  printevery    The printevery variabel, when it is 1 the function outputs all stored fields
        ///
        void reconstruct(bool exportFields = false,
                         fileName folder = "./ITHACAoutput/online_rec",
                         label printevery = 1);

        /// Method to compute the reduced order forces
        ///
        /// @param      problem  a steadyNS full order problem.
        /// @param[in]  folder   The folder where to output the forces matrices
        ///
        void reconstructLiftAndDrag(steadyNS& problem, fileName folder);

        /// Method to evaluate the online inf-sup constant
        ///
        /// @return     return the reduced version of the inf-sup constant.
        double inf_sup_constant();

        ///
        /// @brief      Sets the online velocity.
        ///
        /// @param[in]  vel   The velocity
        ///
        /// @return     The rescaled velocity vector
        ///
        Eigen::MatrixXd setOnlineVelocity(Eigen::MatrixXd vel);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



#endif






