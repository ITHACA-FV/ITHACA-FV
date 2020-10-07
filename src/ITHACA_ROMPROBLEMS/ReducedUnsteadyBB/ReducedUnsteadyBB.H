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
    ReducedUnsteadyBB
Description
    A reduced problem for the unsteady Boussinesq equations
SourceFiles
    ReducedUnsteadyBB.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the ReducedUnsteadyBB class.
/// \dir
/// Directory containing the header and source files for the ReducedUnsteadyBB class

#ifndef ReducedUnsteadyBB_H
#define ReducedUnsteadyBB_H

#include "fvCFD.H"
#include "IOmanip.H"
#include "ReducedUnsteadyNS.H"
#include "UnsteadyBB.H"
#include <Eigen/Dense>
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>

/// Newton object for the resolution of the reduced problem using a supremizer approach
struct newton_unsteadyBB_sup: public newton_argument<double>
{
    public:
        newton_unsteadyBB_sup() {}
        newton_unsteadyBB_sup(int Nx, int Ny,
                              UnsteadyBB& problem): newton_argument<double>(Nx, Ny),
            problem(&problem),
            Nphi_u(problem.NUmodes + problem.liftfield.size() + problem.NSUPmodes),
            Nphi_t(problem.NTmodes + problem.liftfieldT.size()),
            N_BC_t(problem.inletIndexT.rows()),
            N_BC(problem.inletIndex.rows()),
            Nphi_prgh(problem.NPrghmodes)
        {}

        int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const;
        int df(const Eigen::VectorXd& x,  Eigen::MatrixXd& fjac) const;
        UnsteadyBB* problem;
        int Nphi_u;
        int Nphi_t;
        int N_BC_t;
        int N_BC;
        int Nphi_prgh;
        scalar nu;
        scalar Pr;
        scalar dt;

        Eigen::VectorXd y_old;
        Eigen::VectorXd BC_t;
        Eigen::VectorXd BC;

};

struct newton_unsteadyBB_PPE: public newton_argument<double>
{
    public:
        newton_unsteadyBB_PPE() {}
        newton_unsteadyBB_PPE(int Nx, int Ny,
                              UnsteadyBB& problem): newton_argument<double>(Nx, Ny),
            problem(&problem),
            Nphi_u(problem.NUmodes + problem.liftfield.size()),
            Nphi_t(problem.NTmodes + problem.liftfieldT.size()),
            N_BC_t(problem.inletIndexT.rows()),
            N_BC(problem.inletIndex.rows()),
            Nphi_p(problem.NPmodes),
            Nphi_prgh(problem.NPrghmodes)
        {}

        int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const;
        int df(const Eigen::VectorXd& x,  Eigen::MatrixXd& fjac) const;

        UnsteadyBB* problem;
        int Nphi_u;
        int Nphi_p;
        int Nphi_t;
        int Nphi_prgh;
        int N_BC_t;
        int N_BC;
        scalar nu;
        scalar Pr;
        scalar dt;
        int writeEvery;
        word BCmethod;

        Eigen::VectorXd y_old;
        Eigen::VectorXd BC_t;
        Eigen::VectorXd BC;
};



/*---------------------------------------------------------------------------*\
                        Class reducedProblem Declaration
\*---------------------------------------------------------------------------*/

/// Class where it is implemented a reduced problem for the <b> unsteady Navier-stokes </b> problem.
/** In this class are implemented the methods to solve a reduced version of the stationary Navier-stokes
equations. This class must be used together with the steadyNS class */
class ReducedUnsteadyBB: public reducedUnsteadyNS
{
    private:

    public:
        // Constructors
        /// Construct Null
        ReducedUnsteadyBB();
        /// Construct Null
        ///
        /// @param[in]  problem  a full order UnsteadyBB problem
        ///
        explicit ReducedUnsteadyBB(UnsteadyBB& problem);
        ~ReducedUnsteadyBB() {};

        /// List of Eigen matrices to store current online solution for temperature equation
        Eigen::MatrixXd online_solutiont;

        /// List of pointers to store the modes for temperature including lift field modes
        PtrList<volScalarField> LTmodes;

        /// List of pointers to store the modes for velocity including lift field modes
        PtrList<volVectorField> LUmodes;

        /// Function object to call the non linear solver sup approach
        newton_unsteadyBB_sup newton_object_sup;

        /// Function object to call the non linear solver PPE approach
        newton_unsteadyBB_PPE newton_object_PPE;

        /// Reconstructed temperature field
        PtrList<volScalarField> TREC;

        PtrList<volScalarField> PREC;

        /// Number of temperature modes
        int Nphi_t;

        /// Number of shifted pressure modes
        int Nphi_prgh;

        /// Number of parametrized boundary conditions related to temperature field
        int N_BC_t;

        /// DimensionedScalar Pr;
        scalar Pr;

        /// Pointer to the FOM problem
        UnsteadyBB* problem;

        // Functions

        /// Method to perform an online solve using a supremizer stabilisation method
        ///
        /// @param[in]  temp_now_BC   The scalar of online temperature. It is defined in
        /// with an Eigen::MatrixXd and must have one col and as many rows as the number
        /// of parametrized boundary conditions.
        /// @param[in] vel_now_BC   The scalar of online velocity. It is defined in
        /// with an Eigen::MatrixXd and must have one col and as many rows as the number
        /// of parametrized boundary conditions.
        /// @param[in]  NParaSet The parameter sets number for which the online solution is solved.
        /// @param[in]  startSnap The first snapshot taken from the offline snapshots
        /// and used to get the reduced initial condition default = 0.
        ///
        Eigen::MatrixXd solveOnline_sup(Eigen::MatrixXd& temp_now_BC,
                                        Eigen::MatrixXd& vel_now_BC, int NParaSet = 0, int startSnap = 0);


        /// Method to perform an online solve using a PPE stabilisation method
        ///
        /// @param[in]  temp_now_BC   The scalar of online temperature. It is defined in
        /// with an Eigen::MatrixXd and must have one col and as many rows as the number
        /// of parametrized boundary conditions.
        /// @param[in] vel_now_BC   The scalar of online velocity. It is defined in
        /// with an Eigen::MatrixXd and must have one col and as many rows as the number
        /// of parametrized boundary conditions.
        /// @param[in]  NParaSet The parameter sets number for which the online solution is solved.
        /// @param[in]  startSnap The first snapshot taken from the offline snapshots
        /// and used to get the reduced initial condition default = 0.
        ///
        Eigen::MatrixXd solveOnline_PPE(Eigen::MatrixXd& temp_now_BC,
                                        Eigen::MatrixXd& vel_now_BC, int NParaSet = 0,
                                        int startSnap = 0);

        /// Method to reconstruct a solution from an online solve with a supremizer stabilisation technique.
        /// stabilisation method
        ///
        /// @param[in]  folder      The folder where you want to store the results (default is "./ITHACAOutput/online_rec").
        /// @param[in]  printevery  Variable to recover only every printevery online solutions default is 1.
        ///
        void reconstruct_sup(fileName folder = "./ITHACAOutput/online_rec",
                             int printevery = 1);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif



