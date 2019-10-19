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


struct newtonUnsteadyNSTurbSUPIntrusive: public newton_argument<double>
{
    public:
        newtonUnsteadyNSTurbSUPIntrusive() {}
        newtonUnsteadyNSTurbSUPIntrusive(int Nx, int Ny,
                                         UnsteadyNSTurbIntrusive& problem): newton_argument<double>(Nx, Ny),
            problem(&problem),
            Nphi_u(problem.NUmodes + problem.liftfield.size() + problem.NSUPmodes),
            Nphi_p(problem.NPmodes),
            N_BC(problem.inletIndex.rows())
        {}

        int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const;
        int df(const Eigen::VectorXd& x,  Eigen::MatrixXd& fjac) const;

        UnsteadyNSTurbIntrusive* problem;
        int Nphi_u;
        int Nphi_p;
        int N_BC;
        scalar nu;
        scalar dt;
        Eigen::VectorXd y_old;
        Eigen::VectorXd yOldOld;
        Eigen::VectorXd bc;
};



struct newtonUnsteadyNSTurbPPEIntrusive: public newton_argument<double>
{
    public:
        newtonUnsteadyNSTurbPPEIntrusive() {}
        newtonUnsteadyNSTurbPPEIntrusive(int Nx, int Ny,
                                         UnsteadyNSTurbIntrusive& problem): newton_argument<double>(Nx, Ny),
            problem(&problem),
            Nphi_u(problem.NUmodes + problem.liftfield.size()),
            Nphi_p(problem.NPmodes),
            N_BC(problem.inletIndex.rows())
        {}

        int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const;
        int df(const Eigen::VectorXd& x,  Eigen::MatrixXd& fjac) const;

        UnsteadyNSTurbIntrusive* problem;
        int Nphi_u;
        int Nphi_p;
        int N_BC;
        scalar nu;
        scalar dt;
        Eigen::VectorXd y_old;
        Eigen::VectorXd yOldOld;
        Eigen::VectorXd bc;
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

        /// Number of viscosity modes
        int nphiNut;

        /// List of pointers to store the modes for the eddy viscosity
        PtrList<volScalarField> nutModes;

        /// Function object to call the non linear solver sup approach
        newtonUnsteadyNSTurbSUPIntrusive newtonObjectSUP;

        /// Function object to call the non linear solver PPE approach
        newtonUnsteadyNSTurbPPEIntrusive newtonObjectPPE;

        // Functions


        /// Method to perform an online solve using a PPE stabilisation method
        ///
        /// @param[in]  velNow   The vector of online velocity. It is defined in
        /// with an Eigen::MatrixXd and must have one col and as many rows as the number
        /// of parametrized boundary conditions.
        /// @param[in]  startSnap The first snapshot taken from the offline snahpshots
        /// and used to get the reduced initial condition.
        ///
        void solveOnlinePPE(Eigen::MatrixXd velNow, label startSnap = 0);

        /// Method to perform an online solve using a supremizer stabilisation method
        ///
        /// @param[in]  velNow   The vector of online velocity. It is defined in
        /// with an Eigen::MatrixXd and must have one col and as many rows as the number
        /// of parametrized boundary conditions.
        /// @param[in]  startSnap The first snapshot taken from the offline snahpshots
        /// and used to get the reduced initial condition.
        ///
        void solveOnlineSUP(Eigen::MatrixXd velNow, label startSnap = 0);

        /// Method to reconstruct a solution from an online solve with a PPE stabilisation technique.
        /// stabilisation method
        ///
        /// @param[in]  folder      The folder where to output the solutions
        ///
        void reconstructPPE(fileName folder = "./online_rec");

        /// Method to reconstruct the solutions from an online solve with a supremizer stabilisation technique.
        /// stabilisation method
        ///
        /// @param[in]  folder      The folder where to output the solutions
        ///
        void reconstructSUP(fileName folder = "./online_rec");

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






