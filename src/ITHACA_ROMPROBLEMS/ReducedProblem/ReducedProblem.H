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
    reducedProblem
Description
    A general reduced problem class is implemented
SourceFiles
    reducedProblem.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the reducedProblem class.
/// \dir
/// Directory containing the header and source files for the reducedProblem class.
/// \dir ../reducedProblems
/// Directory containing the implementation of all the reduced problems and of the methods to solve online reduced problems.

#ifndef reducedProblem_H
#define reducedProblem_H

#include "fvCFD.H"
#include "IOmanip.H"
#include "reductionProblem.H"
#include <Eigen/Eigen>
#include "newton_argument.H"
#include "Foam2Eigen.H"


/*---------------------------------------------------------------------------*\
                        Class reducedProblem Declaration
\*---------------------------------------------------------------------------*/
/// Base class for the implementation of a reduced problem
class reducedProblem
{
    public:
        // Constructors
        /// Construct Null
        reducedProblem();
        virtual ~reducedProblem() {};

        /// Construct with reduced Problem
        ///
        /// @param  problem a reductionProblem or a child of it
        ///
        explicit reducedProblem(reductionProblem& problem);

        // Variables
        /// Viscosity
        scalar nu;

        /// Pointer to FOAM problem
        reductionProblem* problem;

        // Functions
        /// Virtual Method to perform and online Solve
        virtual void solveOnline();

        ///
        /// @brief      Linear system solver for the online problem. It can be used for any kind of variable.
        /// Boundary conditions are set to 0 as default so that the system is not constrained if no conditions
        /// are given (pressure case). The solver has to be choosen between Eigen linear solvers
        /// (fullPivLu set as default).
        ///
        /// @param[in]  LinSys      The linear system as a list of matrices.
        /// @param[in]  x           Solution of the problem at previous step.
        /// @param      residual    Residual given to let the function update it.
        /// @param[in]  bc          Boundary conditions to constrain the system.
        /// @param[in]  solverType  Solver type.
        ///
        /// @return     Updated solution to the system.
        ///
        static Eigen::MatrixXd solveLinearSys(List<Eigen::MatrixXd> LinSys, Eigen::MatrixXd x,
                                       Eigen::VectorXd& residual, const Eigen::MatrixXd& bc = Eigen::MatrixXd::Zero(0,
                                               0), const std::string solverType = "fullPivLu");
        ///
        /// @brief      Linear system solver for the online problem. It can be used for any kind of variable.
        /// Boundary conditions are set to 0 as default so that the system is not constrained if no conditions
        /// are given (pressure case). The solver has to be choosen between Eigen linear solvers
        /// (fullPivLu set as default).
        ///
        /// @param[in]  LinSys      The linear system as a list of matrices.
        /// @param[in]  x           Solution of the problem at previous step.
        /// @param      residual    Residual given to let the function update it.
        /// @param[in]  solverType  Solver type.
        ///
        /// @return     Updated solution to the system.
        ///
        static Eigen::MatrixXd solveLinearSys(List<Eigen::MatrixXd> LinSys, Eigen::MatrixXd x,
                                       Eigen::VectorXd& residual, const std::string solverType);
};


/// Class for the online evaluation of the coefficient manifold for PODI
class onlineInterp
{
    public:
        // Constructors
        /// Construct Null
        onlineInterp();
        virtual ~onlineInterp() {};

        /// Number of parameters in the mu_samples matrix in the FOM
        int Nmu_samples;

        //--------------------------------------------------------------------------
        /// @brief      Get interpolated coefficients evaluated at points from matrix "mu_interp" using the constructed parameters-coefficients manifold
        /// @param[in]  model       Type of the interpolator model, either "rbf" or spl", which must be the same as in the constructed model in the offline phase (see "buildCoeffManifold" method)
        /// @param[in]  mu_interp   The matrix of points required to be evaluated using the constructed interpolator model. The matrix has a shape (n x m) where "n" is the number of points need to be evaluated (i.e. Nsnapshots) for parametric simulations, while "m" is the number of parameters
        ///
        /// @return     Interpolated coefficient matrix
        ///
        Eigen::MatrixXd getInterpCoeffRBF(std::vector<SPLINTER::RBFSpline> rbfVec,
                                          Eigen::MatrixXd mu_interp);

        //--------------------------------------------------------------------------
        /// @brief      Get interpolated coefficients evaluated at points from matrix "mu_interp" using the constructed parameters-coefficients manifold
        /// @param[in]  model       Type of the interpolator model, either "rbf" or spl", which must be the same as in the constructed model in the offline phase (see "buildCoeffManifold" method)
        /// @param[in]  mu_interp   The matrix of points required to be evaluated using the constructed interpolator model. The matrix has a shape (n x m) where "n" is the number of points need to be evaluated (i.e. Nsnapshots) for parametric simulations, while "m" is the number of parameters
        ///
        /// @return     Interpolated coefficient matrix
        ///
        Eigen::MatrixXd getInterpCoeffSPL(std::vector<SPLINTER::BSpline> splVec,
                                          Eigen::MatrixXd mu_interp);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

