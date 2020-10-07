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
    ReducedUnsteadyNSTurbIntrusive
Description
    A reduced problem for the unsteady turbulent NS equations using fully intrusive approach
SourceFiles
    ReducedUnsteadyNSTurbIntrusive.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the ReducedUnsteadyNSTurbIntrusive class.
/// \dir
/// Directory containing the header and source files for the ReducedUnsteadyNSTurbIntrusive class

#ifndef ReducedUnsteadyNSTurbIntrusive_H
#define ReducedUnsteadyNSTurbIntrusive_H

#include "fvCFD.H"
#include "IOmanip.H"
#include "ReducedSteadyNS.H"
#include "ReducedUnsteadyNS.H"
#include "UnsteadyNSTurbIntrusive.H"
#include <Eigen/Dense>
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>


struct newtonUnsteadyNSTurbIntrusive: public newton_argument<double>
{
    public:
        newtonUnsteadyNSTurbIntrusive() {}
        newtonUnsteadyNSTurbIntrusive(label Nx, label Ny,
                                      UnsteadyNSTurbIntrusive& problem): newton_argument<double>(Nx, Ny),
            problem(&problem),
            Nphi_u(problem.nModesOnline),
            N_BC(problem.inletIndex.rows())
        {}

        label operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const;
        label df(const Eigen::VectorXd& x,  Eigen::MatrixXd& fjac) const;

        UnsteadyNSTurbIntrusive* problem;
        label Nphi_u;
        label N_BC;
        scalar nu;
        scalar dt;
        Eigen::VectorXd y_old;
        Eigen::VectorXd yOldOld;
        Eigen::VectorXd bc;
        Eigen::MatrixXd tauU;
};

struct newtonUnsteadyNSTurbIntrusivePPE: public newton_argument<double>
{
    public:
        newtonUnsteadyNSTurbIntrusivePPE() {}
        newtonUnsteadyNSTurbIntrusivePPE(label Nx, label Ny,
                                         UnsteadyNSTurbIntrusive& problem): newton_argument<double>(Nx, Ny),
            problem(&problem),
            Nphi_u(problem.NUmodes),
            Nphi_p(problem.NPmodes),
            N_BC(problem.inletIndex.rows())
        {}

        label operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const;
        label df(const Eigen::VectorXd& x,  Eigen::MatrixXd& fjac) const;

        UnsteadyNSTurbIntrusive* problem;
        label Nphi_u;
        label Nphi_p;
        label N_BC;
        scalar nu;
        scalar dt;
        Eigen::VectorXd y_old;
        Eigen::VectorXd yOldOld;
        Eigen::VectorXd bc;
        Eigen::MatrixXd tauU;
};

/*---------------------------------------------------------------------------*\
                        Class reducedProblem Declaration
\*---------------------------------------------------------------------------*/

/// Class where it is implemented a reduced problem for the <b> unsteady Navier-stokes </b> problem.
/** In this class are implemented the methods to solve a reduced version of the stationary Navier-stokes
equations. This class must be used together with the steadyNS class */
class ReducedUnsteadyNSTurbIntrusive: public reducedUnsteadyNS
{
    private:

    public:
        // Constructors
        /// Construct Null
        ReducedUnsteadyNSTurbIntrusive();
        /// Construct Null
        ///
        /// @param      problem  a full order UnsteadyNSTurbIntrusive problem
        ///
        explicit ReducedUnsteadyNSTurbIntrusive(UnsteadyNSTurbIntrusive& problem);

        ~ReducedUnsteadyNSTurbIntrusive() {};

        /// Pointer to the FOM problem
        UnsteadyNSTurbIntrusive* problem;

        /// List of pointers to store the modes for the eddy viscosity
        PtrList<volScalarField> nutModes;

        /// Function object to call the non linear solver sup approach
        newtonUnsteadyNSTurbIntrusive newtonObject;

        /// Function object to call the non linear solver sup approach
        newtonUnsteadyNSTurbIntrusivePPE newtonObjectPPE;

        /// Tha matrix containing the initial conditions for the reduced variables,
        /// in case of PPE approach it must be one column with the head being
        /// the degrees of freedom of velocity adn the tail the ones of pressure
        Eigen::MatrixXd initCond;

        // Functions

        //--------------------------------------------------------------------------
        /// @brief      Method to perform an online solve using an intrusive approach
        ///
        /// @param[in]  vel   The velocity
        ///
        void solveOnline(Eigen::MatrixXd vel);

        //--------------------------------------------------------------------------
        /// @brief      Method to perform an online solve using an intrusive approach with the usage of PPE
        ///
        /// @param[in]  velNow   The vector of online velocity. It is defined in
        /// with an Eigen::MatrixXd and must have one col and as many rows as the number
        /// of parametrized boundary conditions.
        ///
        void solveOnlinePPE(Eigen::MatrixXd vel);

        //--------------------------------------------------------------------------
        /// @brief      Method to reconstruct the solutions from an online solve with same number of modes for all variables
        ///
        /// @param[in]  folder      The folder where to output the solutions
        ///
        void reconstruct(fileName folder = "./online_rec");

        //--------------------------------------------------------------------------
        /// @brief      Method to reconstruct the solutions from an online solve with number of modes different for pressure
        ///
        /// @param[in]  folder      The folder where to output the solutions
        ///
        void reconstructPPE(fileName folder);

        //--------------------------------------------------------------------------
        /// @brief      Sets the online velocity.
        ///
        /// @param[in]  vel   The velocity
        ///
        /// @return     The rescaled velocity vector
        ///
        Eigen::MatrixXd setOnlineVelocity(Eigen::MatrixXd vel);

        //--------------------------------------------------------------------------
        /// @brief      Method to compute the reduced order forces for same number of modes of velocity and pressure
        ///
        /// @param      problem  A UnsteadyNSTurbIntrusive full order problem.
        /// @param[in]  folder   The folder where to output the forces matrices
        ///
        void reconstructLiftAndDrag(UnsteadyNSTurbIntrusive& problem,
                                    fileName folder);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



#endif






