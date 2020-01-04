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
    ReducedSteadyNSTurbIntrusive
Description
    A reduced problem for the stationary turbulent NS equations using fully intrusive approach
SourceFiles
    ReducedSteadyNSTurbIntrusive.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the ReducedSteadyNSTurbIntrusive class.
/// \dir
/// Directory containing the header and source files for the ReducedSteadyNSTurbIntrusive class

#ifndef ReducedSteadyNSTurbIntrusive_H
#define ReducedSteadyNSTurbIntrusive_H

#include "fvCFD.H"
#include "IOmanip.H"
#include "ReducedProblem.H"
#include "ReducedSteadyNS.H"
#include "SteadyNSTurbIntrusive.H"
#include "ITHACAutilities.H"
#include <Eigen/Dense>
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>

/// Structure to implement a newton object for a stationary NS problem
struct newtonSteadyNSTurbIntrusive: public newton_argument<double>
{
    public:
        newtonSteadyNSTurbIntrusive() {}
        newtonSteadyNSTurbIntrusive(int Nx, int Ny,
                                    SteadyNSTurbIntrusive& problem): newton_argument<double>(Nx, Ny),
            problem(&problem),
            Nphi_u(problem.nModesOnline),
            N_BC(problem.inletIndex.rows())
        {}

        int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const;
        int df(const Eigen::VectorXd& x,  Eigen::MatrixXd& fjac) const;

        SteadyNSTurbIntrusive* problem;
        int Nphi_u;
        int N_BC;
        scalar nu;
        Eigen::MatrixXd tauU;
        Eigen::VectorXd bc;
};


/*---------------------------------------------------------------------------*\
                        Class reducedProblem Declaration
\*---------------------------------------------------------------------------*/


/// Class where it is implemented a reduced problem for the <b> steady turbulent Navier-stokes </b> problem.
/** In this class are implemented the methods to solve a reduced version of the stationary turbulent Navier-stokes
equations. This class must be used together with the steadyNS class */
class ReducedSteadyNSTurbIntrusive: public reducedSteadyNS
{
    private:

    public:
        // Constructors
        /// Construct Null
        ReducedSteadyNSTurbIntrusive();
        /// Construct Null
        ///
        /// @param      problem  a full order SteadyNSTurb problem
        ///
        ReducedSteadyNSTurbIntrusive(SteadyNSTurbIntrusive& problem);

        ~ReducedSteadyNSTurbIntrusive() {};


        /// List of snapshots for the solution for eddy viscosity
        PtrList<volScalarField> nutFields;

        /// List of POD modes for eddy viscosity
        PtrList<volScalarField> nutModes;

        /// Reconstructed eddy viscosity field
        PtrList<volScalarField> nutRec;

        /// Number of viscosity modes
        int nphiNut;

        /// Newton Object to solve the nonlinear problem
        newtonSteadyNSTurbIntrusive newtonObject;

        /// Pointer to the FOM problem
        SteadyNSTurbIntrusive* problem;

        ///
        /// Method to perform an online solve using a supremizer stabilisation method
        ///
        /// @param[in]  velNow   The vector of online velocity. It is defined in
        /// with an Eigen::MatrixXd and must have one col and as many rows as the number
        /// of parametrized boundary conditions.
        /// @param[in]  startSnap The first snapshot taken from the offline snahpshots
        /// and used to get the reduced initial condition.
        ///
        void solveOnline(Eigen::MatrixXd velNow);

        /// Method to reconstruct the solutions from an online solve with a supremizer stabilisation technique.
        /// stabilisation method
        ///
        /// @param[in]  folder      The folder where to output the solutions
        /// @param[in]  printEvery  The printEvery
        ///
        void reconstruct(fileName folder = "./ITHACAOutput/online_rec",
                         int printEvery = 1);

        ///
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
        /// @param      problem  A SteadyNSTurbIntrusive full order problem.
        /// @param[in]  folder   The folder where to output the forces matrices
        ///
        void reconstructLiftAndDrag(SteadyNSTurbIntrusive& problem,
                                    fileName folder);


};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



#endif






