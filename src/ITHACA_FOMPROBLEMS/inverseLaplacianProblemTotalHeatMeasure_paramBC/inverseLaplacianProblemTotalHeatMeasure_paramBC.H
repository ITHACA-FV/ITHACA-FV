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
    inverseLaplacianProblemTotalHeatMeasure_paramBC
Description
    A general full order implementation of an inverse problem
SourceFiles
    inverseLaplacianProblemTotalHeatMeasure_paramBC.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the inverseLaplacianProblemTotalHeatMeasure_paramBC class.
/// \dir
/// Directory containing the header and source files for the inverseLaplacianProblemTotalHeatMeasure_paramBC class.


#ifndef inverseLaplacianProblemTotalHeatMeasure_paramBC_H
#define inverseLaplacianProblemTotalHeatMeasure_paramBC_H
#include "inverseLaplacianProblem_paramBC.H"
#define _USE_MATH_DEFINES

/// Class to implement a full order inverse laplacian problem
class inverseLaplacianProblemTotalHeatMeasure_paramBC: public
    inverseLaplacianProblem_paramBC
{

    public:
        // Constructors
        inverseLaplacianProblemTotalHeatMeasure_paramBC();
        /// Construct with argc and argv
        inverseLaplacianProblemTotalHeatMeasure_paramBC(int argc, char* argv[]);
        virtual ~inverseLaplacianProblemTotalHeatMeasure_paramBC() {};

        double gIntegral;
        double gIntegralWeight = 0;
        double gIntegral_meas = 0;

        Eigen::MatrixXd Phi;
        Eigen::VectorXd phi;


        //--------------------------------------------------------------------------

        /// Performs offline computation for the parameterized BC method, if
        /// the offline directory ""./ITHACAoutputs/offlineParamBC" exists,
        /// it reads the solution from there
        ///
        /// @param[in]  force   If 1, force the offline phase to be computed
        ///
        void parameterizedBCoffline(bool force = 0);

        //--------------------------------------------------------------------------

        ///
        ///
        Eigen::VectorXd parameterizedBC(word linSys_solver = "fullPivLU",
                                        double regPar = 0);

};

#endif
