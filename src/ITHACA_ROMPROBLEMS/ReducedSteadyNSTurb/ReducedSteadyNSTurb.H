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
    ReducedSteadyNSTurb
Description
    A reduced problem for the stationary turbulent NS equations
SourceFiles
    ReducedSteadyNSTurb.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the ReducedSteadyNSTurb class.
/// \dir
/// Directory containing the header and source files for the ReducedSteadyNSTurb class

#ifndef ReducedSteadyNSTurb_H
#define ReducedSteadyNSTurb_H

#include "fvCFD.H"
#include "IOmanip.H"
#include "ReducedProblem.H"
#include "ReducedSteadyNS.H"
#include "SteadyNSTurb.H"
#include "ITHACAutilities.H"
#include <Eigen/Dense>
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>

/// Structure to implement a newton object for a stationary NS problem
struct newtonSteadyNSTurb: public newton_argument<double>
{
    public:
        newtonSteadyNSTurb() {}
        newtonSteadyNSTurb(label Nx, label Ny,
                           SteadyNSTurb& problem): newton_argument<double>(Nx, Ny),
            problem(&problem),
            Nphi_u(problem.NUmodes + problem.liftfield.size() + problem.NSUPmodes),
            nphiNut(problem.nNutModes),
            Nphi_p(problem.NPmodes),
            N_BC(problem.inletIndex.rows()),
            gNut(problem.nNutModes)
        {}

        label operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const;
        label df(const Eigen::VectorXd& x,  Eigen::MatrixXd& fjac) const;

        SteadyNSTurb* problem;
        label Nphi_u;
        label nphiNut;
        label Nphi_p;
        label N_BC;
        Eigen::VectorXd gNut;
        scalar nu;
        Eigen::MatrixXd tauU;
        Eigen::VectorXd bc;
        std::vector<SPLINTER::RBFSpline*> SPLINES;
};


/*---------------------------------------------------------------------------*\
                        Class reducedProblem Declaration
\*---------------------------------------------------------------------------*/


/// Class where it is implemented a reduced problem for the <b> steady turbulent Navier-stokes </b> problem.
/** In this class are implemented the methods to solve a reduced version of the stationary turbulent Navier-stokes
equations. This class must be used together with the steadyNS class */
class ReducedSteadyNSTurb: public reducedSteadyNS
{
    private:

    public:
        // Constructors
        /// Construct Null
        ReducedSteadyNSTurb();
        /// Construct Null
        ///
        /// @param      problem  a full order SteadyNSTurb problem
        ///
        ReducedSteadyNSTurb(SteadyNSTurb& problem);

        ~ReducedSteadyNSTurb() {};


        /// List of snapshots for the solution for eddy viscosity
        PtrList<volScalarField> nutFields;

        /// List of POD modes for eddy viscosity
        PtrList<volScalarField> nutModes;

        /// The matrix of the eddy viscosity RBF interoplated coefficients
        Eigen::MatrixXd rbfCoeffMat;

        /// Vector of eddy viscosity RBF interoplated coefficients
        Eigen::VectorXd rbfCoeff;

        /// Number of viscosity modes
        label nphiNut;

        /// Newton Object to solve the nonlinear problem
        newtonSteadyNSTurb newtonObject;

        /// Pointer to the FOM problem
        SteadyNSTurb* problem;

        /// Reconstructed eddy viscosity fields list
        PtrList<volScalarField> nutRecFields;

        /// Method to perform an online solve using a supremizer stabilisation method
        ///
        /// @param[in]  velNow   The vector of online velocity. It is defined in
        /// with an Eigen::MatrixXd and must have one col and as many rows as the number
        /// of parametrized boundary conditions.
        /// @param[in]  startSnap The first snapshot taken from the offline snahpshots
        /// and used to get the reduced initial condition.
        ///
        void solveOnlineSUP(Eigen::MatrixXd velNow);


        /// Method to reconstruct the solutions from an online solve.
        ///
        /// @param[in]  exportFields  A boolean variable which determines whether to export fields or not
        /// @param[in]  folder        The folder where to output the solutions in case on wants to
        /// @param[in]  printevery    The printevery variabel, when it is 1 the function outputs all stored fields
        ///
        void reconstruct(bool exportFields = false,
                         fileName folder = "./ITHACAoutput/online_rec",
                         label printevery = 1);

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






