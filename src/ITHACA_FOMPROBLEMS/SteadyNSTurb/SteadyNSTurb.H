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
    SteadyNSTurb
Description
    Reduction class for a Stationary turbulent Navier-Stokes problem.
SourceFiles
    SteadyNSTurb.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the SteadyNSTurb class.
/// \dir
/// Directory containing the header and source files for the SteadyNSTurb class.

#ifndef SteadyNSTurb_H
#define SteadyNSTurb_H
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "pisoControl.H"
#include "fvOptions.H"
#include "reductionProblem.H"
#include "ITHACAstream.H"
#include <iostream>
#include <datatable.h>
#include <bspline.h>
#include <bsplinebuilder.h>
#include <rbfspline.h>
#include <spline.h>
#include <Eigen/Dense>
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                        Class SteadyNSturb Declaration
\*---------------------------------------------------------------------------*/

/// Implementation of a parametrized full order <b> steady turbulent Navier Stokes problem </b> and preparation of the the reduced matrices for the online solve.
/** In this class are implemented the methods for the offline solve of a turbulent steady NS problem and the for the generation of the reduced matrices for subsequent online solve, this class is a son
of the reduction problem class */
class SteadyNSTurb: public steadyNS
{


    public:
        // Constructors
        SteadyNSTurb();
        /// Construct with argc and argv
        SteadyNSTurb(int argc, char* argv[]);

        // Member Functions
        /// List of snapshots for the solution for eddy viscosity
        PtrList<volScalarField> nutFields;

        /// List of POD modes for eddy viscosity
        volScalarModes nutModes;

        /// Create a samples for interpolation
        std::vector<SPLINTER::DataTable*> samples;

        /// Create a RBF splines for interpolation
        std::vector<SPLINTER::RBFSpline*> rbfSplines;

        /// Turbulent viscosity matrix
        Eigen::MatrixXd btMatrix;

        /// Turbulent viscosity tensor
        List <Eigen::MatrixXd> ct2Matrix;
        Eigen::Tensor<double, 3 > ct2Tensor;

        /// Turbulent viscosity tensor
        List <Eigen::MatrixXd> ct1Matrix;
        Eigen::Tensor<double, 3 > ct1Tensor;

        /// Total Turbulent tensor
        List <Eigen::MatrixXd> cTotalMatrix;
        Eigen::Tensor<double, 3 > cTotalTensor;

        /// Total B Matrix
        Eigen::MatrixXd bTotalMatrix;

        /// The matrix of L2 projection coefficients for the eddy viscosity
        Eigen::MatrixXd coeffL2;

        /// The vector of L2 projection coefficients for the eddy viscosity snapshot
        Eigen::VectorXd nutCoeff;

        /// The way to compute the eddy viscosity snapshots
        word viscCoeff;

        /// Number of viscoisty modes used for the projection
        label nNutModes;

        /// Eddy viscosity field
        autoPtr<volScalarField> _nut;
        autoPtr<volScalarField> _nuTilda;

        //--------------------------------------------------------------------------
        /// @brief      Perform a truthsolve
        /// @param[in]  mu_now  The actual value of the parameter for this truthSolve. Used only
        /// to construct mu_interp matrix which is written out in the Offline folder, also for par file in
        /// the Parameters folder.
        ///
        void truthSolve(List<scalar> mu_now);

        //--------------------------------------------------------------------------
        /// Project using a supremizer approach
        ///
        /// @param[in]  folder     The folder name
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NPmodes    The number of pressure modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes.
        /// @param[in]  nNutModes  The number of eddy viscosity modes.
        ///
        void projectSUP(fileName folder, label NUmodes, label NPmodes, label NSUPmodes,
                        label nNutModes);

        //--------------------------------------------------------------------------
        ///
        /// @brief      bt added matrix for the turbulence treatement
        ///
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes.
        ///
        /// @return     bt matrix for turbulence treatment
        ///
        Eigen::MatrixXd btTurbulence(label NUmodes, label NSUPmodes);


        //--------------------------------------------------------------------------
        /// @brief      ct1 added matrix for the turbulence treatement
        ///
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes.
        /// @param[in]  nNutModes  The number of eddy viscosity modes.
        ///
        /// @return     ct1 matrix for turbulence treatment
        ///
        List < Eigen::MatrixXd > turbulenceTerm1(label NUmodes, label NSUPmodes,
                label nNutModes);

        //--------------------------------------------------------------------------
        /// @brief      ct1 tensor for the turbulence treatement
        ///
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes.
        /// @param[in]  nNutModes  The number of eddy viscosity modes.
        ///
        /// @return     ct1 tensor for the turbulence treatement
        ///
        Eigen::Tensor<double, 3 > turbulenceTensor1(label NUmodes,
                label NSUPmodes, label nNutModes);

        //--------------------------------------------------------------------------
        /// @brief      Method to compute one of the turbulence eddy viscosity tensors
        ///
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes.
        /// @param[in]  nNutModes  The number of eddy viscosity modes.
        ///
        /// @return     ct2 matrix for turbulence treatment
        ///
        List < Eigen::MatrixXd > turbulenceTerm2(label NUmodes, label NSUPmodes,
                label nNutModes);

        //--------------------------------------------------------------------------
        /// @brief      Method to compute one of the turbulence eddy viscosity tensors
        ///
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes.
        /// @param[in]  nNutModes  The number of eddy viscosity modes.
        ///
        /// @return     ct2 tensor for the turbulence treatement
        ///
        Eigen::Tensor<double, 3 > turbulenceTensor2(label NUmodes,
                label NSUPmodes, label nNutModes);

};

#endif









