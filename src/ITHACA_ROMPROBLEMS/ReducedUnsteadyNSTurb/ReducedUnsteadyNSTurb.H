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
    ReducedUnsteadyNSTurb
Description
    A reduced problem for the unsteady turbulent NS equations
SourceFiles
    ReducedUnsteadyNSTurb.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the ReducedUnsteadyNSTurb class.
/// \dir
/// Directory containing the header and source files for the ReducedUnsteadyNSTurb class

#ifndef ReducedUnsteadyNSTurb_H
#define ReducedUnsteadyNSTurb_H

#include "fvCFD.H"
#include "IOmanip.H"
#include "ReducedSteadyNS.H"
#include "ReducedUnsteadyNS.H"
#include "UnsteadyNSTurb.H"
#include <Eigen/Dense>
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>


struct newtonUnsteadyNSTurbSUP: public newton_argument<double>
{
    public:
        newtonUnsteadyNSTurbSUP() {}
        newtonUnsteadyNSTurbSUP(int Nx, int Ny,
                                UnsteadyNSTurb& problem): newton_argument<double>(Nx, Ny),
            problem(&problem),
            Nphi_u(problem.NUmodes + problem.liftfield.size() + problem.NSUPmodes),
            nphiNut(problem.nNutModes),
            Nphi_p(problem.NPmodes),
            N_BC(problem.inletIndex.rows()),
            gNut(problem.nNutModes)
        {}

        int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const;
        int df(const Eigen::VectorXd& x,  Eigen::MatrixXd& fjac) const;

        UnsteadyNSTurb* problem;
        int Nphi_u;
        int nphiNut;
        int Nphi_p;
        int N_BC;
        scalar nu;
        scalar dt;
        Eigen::VectorXd y_old;
        Eigen::VectorXd yOldOld;
        Eigen::VectorXd bc;
        Eigen::MatrixXd tauU;
        Eigen::VectorXd gNut;
        std::vector<SPLINTER::RBFSpline*> SPLINES;
};



struct newtonUnsteadyNSTurbPPE: public newton_argument<double>
{
    public:
        newtonUnsteadyNSTurbPPE() {}
        newtonUnsteadyNSTurbPPE(int Nx, int Ny,
                                UnsteadyNSTurb& problem): newton_argument<double>(Nx, Ny),
            problem(&problem),
            Nphi_u(problem.NUmodes + problem.liftfield.size()),
            nphiNut(problem.nNutModes),
            Nphi_p(problem.NPmodes),
            N_BC(problem.inletIndex.rows()),
            gNut(problem.nNutModes)
        {}

        int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const;
        int df(const Eigen::VectorXd& x,  Eigen::MatrixXd& fjac) const;

        UnsteadyNSTurb* problem;
        int Nphi_u;
        int nphiNut;
        int Nphi_p;
        int N_BC;
        scalar nu;
        scalar dt;
        Eigen::VectorXd y_old;
        Eigen::VectorXd yOldOld;
        Eigen::VectorXd bc;
        Eigen::MatrixXd tauU;
        Eigen::VectorXd gNut;
        std::vector<SPLINTER::RBFSpline*> SPLINES;
};



/*---------------------------------------------------------------------------*\
                        Class reducedProblem Declaration
\*---------------------------------------------------------------------------*/

/// Class where it is implemented a reduced problem for the <b> unsteady Navier-stokes </b> problem.
/** In this class are implemented the methods to solve a reduced version of the stationary Navier-stokes
equations. This class must be used together with the steadyNS class */
class ReducedUnsteadyNSTurb: public reducedUnsteadyNS
{
    private:

    public:
        // Constructors
        /// Construct Null
        ReducedUnsteadyNSTurb();
        /// Construct Null
        ///
        /// @param      problem  a full order UnsteadyNSTurb problem
        ///
        explicit ReducedUnsteadyNSTurb(UnsteadyNSTurb& problem);

        ~ReducedUnsteadyNSTurb() {};

        /// Pointer to the FOM problem
        UnsteadyNSTurb* problem;

        /// Number of viscosity modes
        int nphiNut;

        /// List of pointers to store the modes for the eddy viscosity
        PtrList<volScalarField> nutModes;

        /// The matrix of the eddy viscosity RBF interoplated coefficients
        Eigen::MatrixXd rbfCoeffMat;

        /// Function object to call the non linear solver sup approach
        newtonUnsteadyNSTurbSUP newtonObjectSUP;

        /// Function object to call the non linear solver PPE approach
        newtonUnsteadyNSTurbPPE newtonObjectPPE;

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






