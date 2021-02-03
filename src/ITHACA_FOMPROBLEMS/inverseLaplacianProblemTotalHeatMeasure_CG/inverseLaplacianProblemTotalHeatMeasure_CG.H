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
    inverseLaplacianProblemTotalHeatMeasure_CG
Description
    A general full order implementation of an inverse problem
SourceFiles
    inverseLaplacianProblemTotalHeatMeasure_CG.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the inverseLaplacianProblemTotalHeatMeasure_CG class.
/// \dir
/// Directory containing the header and source files for the inverseLaplacianProblemTotalHeatMeasure_CG class.


#ifndef inverseLaplacianProblemTotalHeatMeasure_CG_H
#define inverseLaplacianProblemTotalHeatMeasure_CG_H
#include "inverseLaplacianProblem_CG.H"
#define _USE_MATH_DEFINES

/// Class to implement a full order inverse laplacian problem
class inverseLaplacianProblemTotalHeatMeasure_CG: public
    inverseLaplacianProblem_CG
{

    public:
        // Constructors
        inverseLaplacianProblemTotalHeatMeasure_CG();
        /// Construct with argc and argv
        inverseLaplacianProblemTotalHeatMeasure_CG(int argc, char* argv[]);
        virtual ~inverseLaplacianProblemTotalHeatMeasure_CG() {};

        double gIntegral;
        double gIntegralWeight = 0;
        double gIntegral_meas = 0;


        //--------------------------------------------------------------------------

        /// Conjugate gradient method
        ///
        /// @param[out]  1 if converged within cgIterMax iterations, 0 if not
        ///
        int conjugateGradient();

        //--------------------------------------------------------------------------

        /// Computes the gradient of cost function J and its L2 norm
        ///
        void computeGradJ();


        //--------------------------------------------------------------------------

        /// Compute the search step beta
        ///
        void computeSearchStep();


        //--------------------------------------------------------------------------

        /// Convergence cher for the conjugate gradient method
        ///
        /// @param[out]  1 if convergence criteria met, 0 if not
        ///
        int conjugateGradientConvergenceCheck();

};

#endif
