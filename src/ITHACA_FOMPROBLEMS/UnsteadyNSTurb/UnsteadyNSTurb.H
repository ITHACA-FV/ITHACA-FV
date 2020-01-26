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
    UnsteadyNSTurb
Description
    Reduction class for a NON-Stationary turbulent NS problem
SourceFiles
    UnsteadyNSTurb.C

\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the UnsteadyNSTurb class.
/// \dir
/// Directory containing the header and source files for the UnsteadyNSTurb class.

#ifndef UnsteadyNSTurb_H
#define UnsteadyNSTurb_H
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "IOporosityModelList.H"
#include "IOMRFZoneList.H"
#include "fixedFluxPressureFvPatchScalarField.H"
#include "steadyNS.H"
#include "unsteadyNS.H"
#include <iostream>
#include <datatable.h>
#include <bspline.h>
#include <bsplinebuilder.h>
#include <rbfspline.h>
#include <spline.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                        Class SteadyNS Declaration
\*---------------------------------------------------------------------------*/

/// Implementation of a parametrized full order <b> unsteady NS problem </b> and preparation of the the reduced matrices for the online solve.
/** In this class are implemented the methods for the offline solve of a unsteady NS problem
and the for the generation of the reduced matrices for subsequent online solve, this class is a son
of the steadyNS class */
class UnsteadyNSTurb: public unsteadyNS
{
    public:
        // Constructors
        /// Construct Null
        UnsteadyNSTurb();

        /// Construct with argc and argv
        UnsteadyNSTurb(int argc, char* argv[]);

        // Specific variable for the unstationary case
        // Member Functions
        /// List of snapshots for the solution for eddy viscosity
        PtrList<volScalarField> nutFields;

        /// List of POD modes for eddy viscosity
        PtrList<volScalarField> nutModes;

        /// List of for eddy viscosity time-averaged fields
        PtrList<volScalarField> nutAve;

        /// Create a Rbf splines for interpolation
        std::vector<SPLINTER::DataTable*> samples;

        /// Create a samples for interpolation
        std::vector<SPLINTER::RBFSpline*> rbfSplines;

        /// Turbulent viscosity term
        Eigen::MatrixXd btMatrix;

        /// Time-parameter combined matrix
        Eigen::MatrixXd z;

        /// Turbulent viscosity tensor
        Eigen::Tensor<double, 3 > ct1Tensor;

        /// Turbulent average viscosity tensor for the splitting approach
        Eigen::Tensor<double, 3 > ct1AveTensor;

        /// Turbulent viscosity tensor in the PPE equation
        Eigen::Tensor<double, 3 > ct1PPETensor;

        /// Turbulent average viscosity tensor for the splitting approach in the PPE equation
        Eigen::Tensor<double, 3 > ct1PPEAveTensor;

        /// Turbulent viscosity tensor
        Eigen::Tensor<double, 3 > ct2Tensor;

        /// Turbulent average viscosity tensor for the splitting approach
        Eigen::Tensor<double, 3 > ct2AveTensor;

        /// Turbulent viscosity tensor in the PPE equation
        Eigen::Tensor<double, 3 > ct2PPETensor;

        /// Turbulent average viscosity tensor for the splitting approach in the PPE equation
        Eigen::Tensor<double, 3 > ct2PPEAveTensor;

        /// Turbulent total viscosity tensor
        Eigen::Tensor<double, 3 > cTotalTensor;

        /// Turbulent total average viscosity tensor for the splitting approach
        Eigen::Tensor<double, 3 > cTotalAveTensor;

        /// Turbulent total viscosity tensor in the PPE equation
        Eigen::Tensor<double, 3 > cTotalPPETensor;

        /// Turbulent total average viscosity tensor for the splitting approach in the PPE equation
        Eigen::Tensor<double, 3 > cTotalPPEAveTensor;

        /// Total B Matrix
        Eigen::MatrixXd bTotalMatrix;

        /// The matrix of L2 projection coefficients for the eddy viscosity
        Eigen::MatrixXd coeffL2;

        /// Velocity coefficients for RBF interpolation
        Eigen::MatrixXd velRBF;

        /// RBF functions radius
        double e = 1;

        /// Number of viscoisty modes used for the projection
        label nNutModes;

        /// Interpolation independent variable choice
        int interChoice = 1;

        /// Eddy viscosity field
        autoPtr<volScalarField> _nut;

        //--------------------------------------------------------------------------
        /// @brief      bt added matrix for the turbulence treatement
        ///
        /// @param[in]  NUmodes    The N of velocity modes.s
        /// @param[in]  NSUPmodes  The N of supremizer modes.
        ///
        /// @return     bt matrix for turbulence treatment
        ///
        Eigen::MatrixXd btTurbulence(label NUmodes, label NSUPmodes);

        //--------------------------------------------------------------------------
        /// @brief      ct1 added tensor for the turbulence treatement
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
        /// @brief      ct2 added tensor for the turbulence treatement
        ///
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes.
        /// @param[in]  nNutModes  The number of eddy viscosity modes.
        ///
        /// @return     ct2 tensor for the turbulence treatement
        ///
        Eigen::Tensor<double, 3 > turbulenceTensor2(label NUmodes,
                label NSUPmodes, label nNutModes);

        //--------------------------------------------------------------------------
        /// @brief      ct1PPE added tensor for the turbulence treatement in the PPE method
        ///
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes.
        /// @param[in]  NPmodes    The number of pressure modes.
        /// @param[in]  nNutModes  The number of eddy viscosity modes.
        ///
        /// @return     ct1PPE tensor for the turbulence treatement in the PPE method
        ///
        Eigen::Tensor<double, 3 > turbulencePPETensor1(label NUmodes,
                label NSUPmodes, label NPmodes, label nNutModes);

        //--------------------------------------------------------------------------
        /// @brief      ct2PPE added tensor for the turbulence treatement in the PPE method
        ///
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes.
        /// @param[in]  NPmodes    The number of pressure modes.
        /// @param[in]  nNutModes  The number of eddy viscosity modes.
        ///
        /// @return     ct2PPE tensor for the turbulence treatement in the PPE method
        ///
        Eigen::Tensor<double, 3 > turbulencePPETensor2(label NUmodes,
                label NSUPmodes, label NPmodes, label nNutModes);

        //--------------------------------------------------------------------------
        /// @brief      ct1Ave added tensor for approximation of the averaged part of the eddy viscosity
        ///
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes.
        ///
        /// @return     ct1Ave tensor corresponding to the approximation of the averaged part
        ///  of the eddy viscosity for the turbulence treatement
        ///
        Eigen::Tensor<double, 3> turbulenceAveTensor1(label NUmodes,
                label NSUPmodes);

        //--------------------------------------------------------------------------
        /// @brief      ct2Ave added tensor for approximation of the averaged part of the eddy viscosity
        ///
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes.
        ///
        /// @return     ct2Ave tensor corresponding to the approximation of the averaged part
        ///  of the eddy viscosity for the turbulence treatement
        ///
        Eigen::Tensor<double, 3> turbulenceAveTensor2(label NUmodes,
                label NSUPmodes);
        //--------------------------------------------------------------------------
        /// @brief      ct1PPEAve added tensor for approximation of the averaged part of the eddy viscosity
        /// with the usage of the PPE approach
        ///
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes.
        /// @param[in]  NPmodes    The number of pressure modes.
        ///
        /// @return     ct1PPEAve tensor corresponding to the approximation of the averaged part
        ///  of the eddy viscosity for the turbulence treatement with the usage of the PPE approach
        ///
        Eigen::Tensor<double, 3> turbulencePPEAveTensor1(label NUmodes,
                label NSUPmodes, label NPmodes);
        //--------------------------------------------------------------------------
        /// @brief      ct2PPEAve added tensor for approximation of the averaged part of the eddy viscosity
        /// with the usage of the PPE approach
        ///
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes.
        /// @param[in]  NPmodes    The number of pressure modes.
        ///
        /// @return     ct2PPEAve tensor corresponding to the approximation of the averaged part
        ///  of the eddy viscosity for the turbulence treatement with the usage of the PPE approach
        ///
        Eigen::Tensor<double, 3> turbulencePPEAveTensor2(label NUmodes,
                label NSUPmodes, label NPmodes);
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
        /// Project using the Poisson Equation for pressure
        ///
        /// @param[in]  folder     The folder name
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NPmodes    The number of pressure modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes which will be set to zero.
        /// @param[in]  nNutModes  The number of eddy viscosity modes.
        ///
        void projectPPE(fileName folder, label NUmodes, label NPmodes, label NSUPmodes,
                        label nNutModes);

        //--------------------------------------------------------------------------
        /// @brief      A method to compute the two matrices needed for the RBF interpolation by combining
        /// the velocity L2 projection coefficients and their time derivatives
        ///
        /// @param[in]  A            The velocity L2 projection coefficients
        /// @param[in]  G            The eddy viscoisty L2 projection coefficients
        /// @param[in]  initSnapInd  The initial snapshots indices for the different parameter samples
        /// @param[in]  timeSnap     The time rate at which snapshots were taken for the different parameter samples
        ///
        /// @return     The matrix of the combined velocity L2 projection coefficients starting from the second
        /// snapshots for each parameter sample together with their time derivatives and also the matrix of the
        /// eddy viscosity L2 projection coefficients associated with the first matrix
        ///
        List < Eigen::MatrixXd > velDerivativeCoeff(Eigen::MatrixXd A,
                Eigen::MatrixXd G,
                Eigen::VectorXd initSnapInd, Eigen::VectorXd timeSnap);

        //--------------------------------------------------------------------------
        /// @brief      A method to compute the two matrices needed for the RBF interpolation by combining
        /// the parameter samples values with the velocity L2 projection coefficients
        ///
        /// @param[in]  A            The velocity L2 projection coefficients
        /// @param[in]  G            The eddy viscoisty L2 projection coefficients
        ///
        /// @return     The matrix of the combined parameter and velocity L2 projection coefficients and also
        /// the corresponding matrix of eddy viscosity L2 projection coefficients
        ///
        List < Eigen::MatrixXd > velParCoeff(Eigen::MatrixXd A, Eigen::MatrixXd G);


        //--------------------------------------------------------------------------
        /// @brief      A method to compute the two matrices needed for the RBF interpolation by combining
        /// the parameter value and the velocity L2 projection coefficients and their time derivatives
        ///
        /// @param[in]  A            The velocity L2 projection coefficients
        /// @param[in]  G            The eddy viscoisty L2 projection coefficients
        /// @param[in]  initSnapInd  The initial snapshots indices for the different parameter samples
        /// @param[in]  timeSnap     The time rate at which snapshots were taken for the different parameter samples
        ///
        /// @return     The matrix of the combined parameter samples values with their corresponding
        /// velocity L2 projection coefficients starting from the second snapshots for each
        /// parameter sample together with their time derivatives and also the matrix of the
        /// eddy viscosity L2 projection coefficients associated with the previous matrix
        ///
        List < Eigen::MatrixXd > velParDerivativeCoeff(Eigen::MatrixXd A,
                Eigen::MatrixXd G,
                Eigen::VectorXd initSnapInd, Eigen::VectorXd timeSnap);

        //--------------------------------------------------------------------------
        /// @brief      A method to compute the matrix of the combination of the parameter sample value and
        /// the matrix of velocity projection coefficients together with their time derivatives
        /// based on the backward scheme
        ///
        /// @param[in]  A         The velocity L2 projection coefficients
        /// @param[in]  par       The parameter vector corresponding to the matrix A
        /// @param[in]  timeSnap  The time rate at which snapshots are taken in matrix A
        ///
        /// @return     the matrix of the combination of the parameter sample value and
        /// the matrix of velocity projection coefficients together with their time derivatives
        /// based on the backward scheme
        ///
        Eigen::MatrixXd velParDerivativeCoeff(Eigen::MatrixXd A,
                                              Eigen::VectorXd par, double timeSnap);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif










