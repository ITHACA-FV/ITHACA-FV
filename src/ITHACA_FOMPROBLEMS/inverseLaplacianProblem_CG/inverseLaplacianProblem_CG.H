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
    inverseLaplacianProblem_CG
Description
    Implementation of the Alifanov's regularization for the inverse problem of
    estimating the boundary heat flux, given pointwise temperature measurements,
    in a Laplacian problem
SourceFiles
    inverseLaplacianProblem_CG.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the inverseLaplacianProblem_CG class.
/// \dir
/// Directory containing the header and source files for the inverseLaplacianProblem_CG class.


#ifndef inverseLaplacianProblem_CG_H
#define inverseLaplacianProblem_CG_H
#include "thermocouplesPlane.H"
#include "inverseLaplacianProblem.H"
#define _USE_MATH_DEFINES

using namespace SPLINTER;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                 Class inverseLaplacianProblem_CG Declaration
\*---------------------------------------------------------------------------*/

/// Implementation of a <b> Alifanov's regularization </b>to solve inverse Laplacian problems.
/** In this class, we implement the Alifanov's regularization for solving
the inverse problem of estimating the boundary heat flux, given pointwise
temperature measurements, in a Laplacian problem*/
class inverseLaplacianProblem_CG: public inverseLaplacianProblem
{

    public:
        // Constructors
        /// Null constructor
        inverseLaplacianProblem_CG();
        /// Construct with argc and argv. Reads the thermocouples dictionary
        inverseLaplacianProblem_CG(int argc, char* argv[]);
        /// Destructor
        virtual ~inverseLaplacianProblem_CG() {};


        /// Adjoint field
        autoPtr<volScalarField> _lambda;

        /// Sensibility temperature field
        autoPtr<volScalarField> _deltaT;

        /// Value fraction for the adjoint Robin BC
        List<scalar> valueFractionAdj;

        /// List of temperature solutions
        PtrList<volScalarField> Tfield;

        /// List of adjoint solutions
        PtrList<volScalarField> lambdaField;

        /// List of sensitivity solutions
        PtrList<volScalarField> deltaTfield;

        /// Flag for interpolation of the temperature measurements
        bool interpolation;

        /// Flag for definition of interpolation plane
        bool interpolationPlaneDefined;

        /// Conjugate Gradient (CG) interations counter
        int cgIter;

        /// Maximum CG iterations
        int cgIterMax;

        /// Flag to save solutions in lists
        bool saveSolInLists = 0;

        /// Vector to store the const function J
        Eigen::MatrixXd Jlist;

        /// L2 norm of the gradient of J
        double gradJ_L2norm;

        /// Absolute stopping criterion for the CG
        double Jtol;

        /// Relative stopping criterion for the CG
        double JtolRel;

        /// Conjugate coefficient
        double gamma;

        /// Denoinator of the conjugate coefficient
        double gamma_den;

        /// CG search step size
        double beta;


        /// Search direction
        List<scalar> P;

        /// Gradient of the cost function
        List<scalar> gradJ;

        /// Vector of solutions of the sensitivity problem at the thermocouples points
        Eigen::VectorXd Tsens;

        /// Temperature reduced matrix
        List<Eigen::MatrixXd> ArbT;

        /// Adjoint reduced matrix
        List<Eigen::MatrixXd> ArbLambda;

        /// Sensitivity reduced matrix
        List<Eigen::MatrixXd> ArbDeltaT;

        /// Interpolation plane
        thermocouplesPlane interpolationPlane;

        /// IDs of the cells in the interpolation plane
        Eigen::VectorXd cellsInPlane;

        // Functions

        //--------------------------------------------------------------------------
        /// Set valueFraction list values for Robin condition
        ///
        void set_valueFraction();

        //--------------------------------------------------------------------------
        /// Set BC of the adjoint problem
        ///
        void assignAdjointBC();

        //--------------------------------------------------------------------------
        /// Assign the BC for the adjoint problem and returs the source field
        ///
        /// @return  Source field
        ///
        volScalarField assignAdjointBCandSource();

        //--------------------------------------------------------------------------
        /// Set BC and IF of the sensitivity problem
        ///
        void assignSensitivityBC();

        //--------------------------------------------------------------------------
        /// Solve adjoint problem
        ///
        void solveAdjoint();

        //--------------------------------------------------------------------------
        /// Solve sensibility problem
        ///
        void solveSensitivity();

        //--------------------------------------------------------------------------
        /// Solve for direct and sensitivity
        ///
        void solve(const char* problemID);

        //--------------------------------------------------------------------------
        /// Identifies the plane defined by the thermocouples
        ///
        void defineThermocouplesPlane();

        //--------------------------------------------------------------------------
        /// Fill the Foam::vector containing the values of the sensitivity solution at
        /// the thermocouples locations
        ///
        void sensibilitySolAtThermocouplesLocations();

        //--------------------------------------------------------------------------
        /// Conjugate gradient method
        ///
        /// @return  1 if converged within cgIterMax iterations, 0 if not
        ///
        int conjugateGradient();

        //--------------------------------------------------------------------------
        /// Computes the gradient of cost function J and its L2 norm
        ///
        void computeGradJ();

        //--------------------------------------------------------------------------
        /// Computes the search direction P
        ///
        void searchDirection();

        //--------------------------------------------------------------------------
        /// Compute the search step beta
        ///
        void computeSearchStep();

        //--------------------------------------------------------------------------
        /// Updates the heat flux in the conjugate gradient iterations
        ///
        void updateHeatFlux();

        //--------------------------------------------------------------------------
        /// Convergence cher for the conjugate gradient method
        ///
        /// @return  1 if convergence criteria met, 0 if not
        ///
        int conjugateGradientConvergenceCheck();

        //--------------------------------------------------------------------------
        /// Checks if a cell crosses the interpolation plane
        ///
        /// @param[in]  cx  X coordinate of the cell center
        /// @param[in]  cy  Y coordinate of the cell center
        /// @param[in]  cz  Z coordinate of the cell center
        /// @param[in]  thermocoupleCellDim  Dimension of the cell
        ///
        /// @return  1 if the cell crosses the plane, 0 elsewise
        ///
        int isInPlane(double cx, double cy, double cz,
                      Foam::vector thermocoupleCellDim);

        //--------------------------------------------------------------------------
        /// Writes fields to file
        ///
        /// @param[in]  folder          Folder in which the fields are written
        /// @param[in]  folderNumber    Folder number
        ///
        void writeFields(label folderNumber, const char* folder);

        //--------------------------------------------------------------------------
        /// Interpolates the thermocouples measurements in the plane defined
        /// in readThermocouples() using radial basis functions.
        /// Assumes thermocouples are all at the same y coordinate.
        /// Assumes all cells have the same dimensions
        ///
        void thermocouplesInterpolation();

        //--------------------------------------------------------------------------
        /// Interpolates the thermocouples measurements in the plane defined
        /// in readThermocouples() using Radial Basis Functions (RBF)
        /// Assumes thermocouples are all at the same y coordinate.
        /// Assumes all cells have the same dimensions
        ///
        /// @param[in]  RBFweights  Weights of the RBFs
        /// @param[in]  RBFbasis    Basis for the RBF interpolation
        ///
        void thermocouplesInterpolation(DenseMatrix& RBFweights, DenseMatrix& RBFbasis);

        //--------------------------------------------------------------------------
        /// Computes the difference between direct problem solution and measure
        ///
        void differenceBetweenDirectAndMeasure();

        //--------------------------------------------------------------------------
        /// Restart fields
        ///
        /// @param[in]  fieldName   Name of the field to restart.
        ///                         Options are (T, lambda and deltaT)
        ///
        void restart(word fieldName = "all");


        //--------------------------------------------------------------------------
        /// Performs offline solution
        ///
        //void offlineSolve();

};

#endif
