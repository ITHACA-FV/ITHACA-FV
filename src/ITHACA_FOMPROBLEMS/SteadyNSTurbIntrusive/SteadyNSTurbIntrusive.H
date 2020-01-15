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
    SteadyNSTurbIntrusive
Description
    Reduction class for a Stationary turbulent Navier-Stokes problem using fully intrusive approach.
SourceFiles
    SteadyNSTurbIntrusive.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the SteadyNSTurbIntrusive class.
/// \dir
/// Directory containing the header and source files for the SteadyNSTurbIntrusive class.

#ifndef SteadyNSTurbIntrusive_H
#define SteadyNSTurbIntrusive_H
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
class SteadyNSTurbIntrusive: public steadyNS
{


    public:
        // Constructors
        SteadyNSTurbIntrusive();
        /// Construct with argc and argv
        SteadyNSTurbIntrusive(int argc, char* argv[]);

        // Member Functions
        /// List of snapshots for the solution for eddy viscosity
        PtrList<volScalarField> nutFields;

        /// List of POD modes for eddy viscosity
        PtrList<volScalarField> nutModes;

        /// Diffusive matrix
        Eigen::MatrixXd bMatrix;

        /// Pressure Gradient matrix
        Eigen::MatrixXd kMatrix;

        /// Turbulent viscosity matrix
        Eigen::MatrixXd btMatrix;

        /// Convective tensor
        Eigen::Tensor<double, 3 > convTensor;

        /// Turbulent viscosity tensor
        Eigen::Tensor<double, 3 > ct2Tensor;

        /// Turbulent viscosity tensor
        Eigen::Tensor<double, 3 > ct1Tensor;

        /// Total Turbulent tensor
        Eigen::Tensor<double, 3 > cTotalTensor;

        /// Total B Matrix
        Eigen::MatrixXd bTotalMatrix;

        /// Number of modes used in the online stage
        label nModesOnline;

        /// Eddy viscosity field
        autoPtr<volScalarField> _nut;
        autoPtr<volScalarField> _nuTilda;

        //--------------------------------------------------------------------------
        /// @brief      Perform a truthsolve
        ///
        /// @param[in]  mu_now  The actual value of the parameter for this truthSolve. Used only
        /// to construct mu_interp matrix which is written out in the Offline folder, also for par file in
        /// the Parameters folder.
        ///
        void truthSolve(List<scalar> mu_now);

        //--------------------------------------------------------------------------
        /// @brief      The projection method for computing the reduced order matrices
        ///
        /// @param[in]  folder  The folder used to save the reduced matrices.
        /// @param[in]  nModes  The number of modes used in the online stage for all reduced variables
        ///
        void project(fileName folder, label nModes);


        //--------------------------------------------------------------------------
        /// @brief      bt added matrix for the turbulence treatement
        ///
        /// @param[in]  nModes  The number of modes used in the online stage for all reduced variables
        ///
        /// @return     bt matrix for turbulence treatment
        ///
        Eigen::MatrixXd btTurbulence(label nModes);


        //--------------------------------------------------------------------------
        /// @brief      Method to compute one of the turbulence eddy viscosity tensors
        ///
        /// @param[in]  nModes  The number of modes used in the online stage for all reduced variables
        ///
        /// @return     ct1 tensor for the turbulence treatement
        ///
        Eigen::Tensor<double, 3 > turbulenceTensor1(label nModes);


        //--------------------------------------------------------------------------
        /// @brief      Method to compute one of the turbulence eddy viscosity tensors
        ///
        /// @param[in]  nModes  The number of modes used in the online stage for all reduced variables
        ///
        /// @return     ct2 tensor for the turbulence treatement
        ///
        Eigen::Tensor<double, 3 > turbulenceTensor2(label nModes);

        //--------------------------------------------------------------------------
        /// @brief      Diffusive Term
        ///
        /// @param[in]  nModes  The number of modes used in the online stage for all reduced variables
        ///
        /// @return     reduced matrix in Eigen::MatrixXd format for the diffusion term.
        ///
        Eigen::MatrixXd diffusiveTerm(label nModes);

        //--------------------------------------------------------------------------
        /// @brief      The method for computing the convective term tensor
        ///
        /// @param[in]  nModes  The number of modes used in the online stage for all reduced variables
        ///
        /// @return     c Tensor for the treatement of the convective term
        ///
        Eigen::Tensor<double, 3> convectiveTerm(label nModes);

        //--------------------------------------------------------------------------
        /// @brief      The method for computing the pressure gradient term with number of modes of pressure
        /// being equal to the number of modes of velocity and with modes of pressure computed using the SVD of
        /// velocity snapshots
        ///
        /// @param[in]  nModes  The number of modes used in the online stage for all reduced variables
        ///
        /// @return     A reduced order matrix which is square and represent the pressure gradient term
        ///
        Eigen::MatrixXd pressureGradientTerm(label nModes);

        //--------------------------------------------------------------------------
        /// @brief      Boundary integral modes on boundary used by the penaly method
        ///
        /// @param[in]  nModes  The number of modes used in the online stage for all reduced variables
        ///
        /// @return     A list of redueced order vectors for the penalty method which contain the boundary
        /// integral of the velocity modes.
        ///
        List< Eigen::MatrixXd > bcVelocityVec(label nModes);

        //--------------------------------------------------------------------------
        /// @brief      Boundary integral modes on boundary used by the penaly method
        ///
        /// @param[in]  nModes  The number of modes used in the online stage for all reduced variables
        ///
        /// @return     A list of redueced order matrices for the penalty method which contain the boundary
        /// integral of the scalar products of the velocity modes.
        List< Eigen::MatrixXd > bcVelocityMat(label nModes);


};

#endif









