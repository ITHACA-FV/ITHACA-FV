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
    reducedUnsteadyNSturb

Description
    A reduced problem for the unsteady turbulent NS equations

SourceFiles
    reducedUnsteadyNSturb.C

\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the reducedUnsteadyNSturb class.
/// \dir
/// Directory containing the header and source files for the reducedUnsteadyNSturb class

#ifndef reducedUnsteadyNSturb_H
#define reducedUnsteadyNSturb_H

#include "fvCFD.H"
#include "IOmanip.H"
#include "reducedSteadyNS.H"
#include "reducedUnsteadyNS.H"
#include "unsteadyNSturb.H"
#include <Eigen/Dense>
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>


struct newton_unsteadyNSturb_sup: public newton_argument<double>
{
    public:
        newton_unsteadyNSturb_sup() {}
        newton_unsteadyNSturb_sup(int Nx, int Ny,
                                  unsteadyNSturb& problem): newton_argument<double>(Nx, Ny),
            problem(&problem),
            Nphi_u(problem.NUmodes + problem.liftfield.size() + problem.NSUPmodes),
            Nphi_nut(problem.Nnutmodes),
            Nphi_p(problem.NPmodes),
            N_BC(problem.inletIndex.rows()),
            nu_c(problem.Nnutmodes)
        {}

        int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const;
        int df(const Eigen::VectorXd& x,  Eigen::MatrixXd& fjac) const;

        int Nphi_u;
        int Nphi_nut;
        int Nphi_p;
        int N_BC;
        unsteadyNSturb* problem;
        scalar nu;
        scalar dt;
        Eigen::VectorXd y_old;
        Eigen::VectorXd BC;
        Eigen::VectorXd nu_c;
        std::vector<SPLINTER::RBFSpline*> SPLINES;
};



struct newton_unsteadyNSturb_PPE: public newton_argument<double>
{
    public:
        newton_unsteadyNSturb_PPE() {}
        newton_unsteadyNSturb_PPE(int Nx, int Ny,
                                  unsteadyNSturb& problem): newton_argument<double>(Nx, Ny),
            problem(&problem),
            Nphi_u(problem.NUmodes + problem.liftfield.size()),
            Nphi_nut(problem.Nnutmodes),
            Nphi_p(problem.NPmodes),
            N_BC(problem.inletIndex.rows()),
            nu_c(problem.Nnutmodes)
        {}

        int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const;
        int df(const Eigen::VectorXd& x,  Eigen::MatrixXd& fjac) const;

        int Nphi_u;
        int Nphi_nut;
        int Nphi_p;
        int N_BC;
        unsteadyNSturb* problem;
        scalar nu;
        scalar dt;
        Eigen::VectorXd y_old;
        Eigen::VectorXd BC;
        Eigen::VectorXd nu_c;
        std::vector<SPLINTER::RBFSpline*> SPLINES;
};



/*---------------------------------------------------------------------------*\
                        Class reducedProblem Declaration
\*---------------------------------------------------------------------------*/

/// Class where it is implemented a reduced problem for the <b> unsteady Navier-stokes </b> problem.
/** In this class are implemented the methods to solve a reduced version of the stationary Navier-stokes
equations. This class must be used together with the steadyNS class */
class reducedUnsteadyNSturb: public reducedUnsteadyNS
{
    private:

    public:
        // Constructors
        /// Construct Null
        reducedUnsteadyNSturb();
        /// Construct Null
        ///
        /// @param      problem  a full order unsteadyNSturb problem
        ///
        explicit reducedUnsteadyNSturb(unsteadyNSturb& problem);

        ~reducedUnsteadyNSturb() {};

        /// Pointer to the FOM problem
        unsteadyNSturb* problem;

        /// Reconstructed eddy viscosity field
        PtrList<volScalarField> nutREC;

        /// Number of viscosity modes
        int Nphi_nut;

        /// Function object to call the non linear solver sup approach
        newton_unsteadyNSturb_sup newton_object_sup;

        /// Function object to call the non linear solver PPE approach
        newton_unsteadyNSturb_PPE newton_object_PPE;

        /// Scalar to store the current time
        double time;

        /// Scalar to store the time increment
        double dt;

        /// Scalar to store the final time if the online simulation
        double finalTime;

        /// Scalar to store the final time if the online simulation
        double tstart;



        // Functions


        /// Method to perform an online solve using a PPE stabilisation method
        ///
        /// @param[in]  vel_now   The vector of online velocity. It is defined in
        /// with an Eigen::MatrixXd and must have one col and as many rows as the number
        /// of parametrized boundary conditions.
        /// @param[in]  startSnap The first snapshot taken from the offline snahpshots
        /// and used to get the reduced initial condition.
        ///
        void solveOnline_PPE(Eigen::MatrixXd& vel_now, label startSnap = 0);

        /// Method to perform an online solve using a supremizer stabilisation method
        ///
        /// @param[in]  vel_now   The vector of online velocity. It is defined in
        /// with an Eigen::MatrixXd and must have one col and as many rows as the number
        /// of parametrized boundary conditions.
        /// @param[in]  startSnap The first snapshot taken from the offline snahpshots
        /// and used to get the reduced initial condition.
        ///
        void solveOnline_sup(Eigen::MatrixXd& vel_now, label startSnap = 0);

        /// Method to reconstruct a solution from an online solve with a PPE stabilisation technique.
        /// stabilisation method
        ///
        /// @param[in]  folder      The folder where to output the solutions
        /// @param[in]  printevery  The printevery
        ///
        void reconstruct_PPE(fileName folder = "./online_rec", int printevery = 1);

        /// Method to reconstruct the solutions from an online solve with a supremizer stabilisation technique.
        /// stabilisation method
        ///
        /// @param[in]  folder      The folder where to output the solutions
        /// @param[in]  printevery  The printevery
        ///
        void reconstruct_sup(fileName folder = "./online_rec", int printevery = 1);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



#endif






