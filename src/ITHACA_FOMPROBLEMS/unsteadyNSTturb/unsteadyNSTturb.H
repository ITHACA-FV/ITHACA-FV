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
    unsteadyNS

Description
    Reduction class for a NON-Stationary NS plus Energy Equation turbulent problem

SourceFiles
    unsteadyNSTturb.C

\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the unsteadyNSTturb class.
/// \dir
/// Directory containing the header and source files for the unsteadyNSTturb class.

#ifndef unsteadyNSTturb_H
#define unsteadyNSTturb_H
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "IOporosityModelList.H"
#include "IOMRFZoneList.H"
#include "fixedFluxPressureFvPatchScalarField.H"
#include "steadyNS.H"
#include "unsteadyNST.H"
#include <iostream>
#include <datatable.h>
#include <bspline.h>
#include <bsplinebuilder.h>
#include <rbfspline.h>
#include <spline.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                        Class SteadyNSTturb Declaration
\*---------------------------------------------------------------------------*/

/// Implementation of a parametrized full order <b> unsteady NST problem weakly coupled with the energy equation </b> and preparation of the the reduced matrices for the online solve.
/** In this class are implemented the methods for the offline solve of a unsteady NSTturb problem
and the for the generation of the reduced matrices for subsequent online solve, this class is a son
of the unsteadyNST class */
class unsteadyNSTturb: public unsteadyNST
{
    public:
        // Constructors
        /// Construct Null
        unsteadyNSTturb();

        /// Construct with argc and argv
        unsteadyNSTturb(int argc, char* argv[]);
        ~unsteadyNSTturb() {};

        // Specific variable for the unstationary case

        // Member Functions

        /// List of snapshots for the solution for eddy viscosity
        PtrList<volScalarField> nutFields;

        /// List of snapshots for the solution for eddy viscosity
        PtrList<volScalarField> liftfieldnut;

        /// List of snapshots to form the homogeneous eddy viscosity fields
        PtrList<volScalarField> nutomFields;

        /// List of POD modes for eddy viscosity
        PtrList<volScalarField> nuTmodes;

        /// Eddy diffusivity field
        PtrList<volScalarField> alphatfield;

        /// Scalar to store the turbulent Prandtl number
        scalar Prt;

        /// Scalar to store the Prandtl number
        scalar Pr;

        /// Number of nut field modes
        int Nphi_nut;

        /// Create a Rbf splines for interpolation
        std::vector<SPLINTER::DataTable*> SAMPLES;

        /// Create a SAMPLES for interpolation
        std::vector<SPLINTER::RBFSpline*> rbfsplines;

        /// Number of viscoisty modes used for the projection
        label Nnutmodes;

        /** @name Reduced Matrices
        *
        */
        ///@{
        //
        /// Turbulent viscosity term
        Eigen::MatrixXd BT_matrix;

        /// Turbulent viscosity term
        List <Eigen::MatrixXd> CT2_matrix;

        /// Turbulent viscosity term
        List <Eigen::MatrixXd> CT1_matrix;

        /// Turbulent diffusivity term
        List <Eigen::MatrixXd> S_matrix;

        /// Total C Matrix
        List <Eigen::MatrixXd> C_total_matrix;

        /// Total B Matrix
        Eigen::MatrixXd B_total_matrix;
        ///@}

        /// Eddy viscosity field
        autoPtr<volScalarField> _nut;

        ///Eddy viscosity for the Spalart-Allmaras turbulence model
        autoPtr<volScalarField> _nuTilda;

        /// Turbulent Prandtl number
        autoPtr<dimensionedScalar> _Prt;

        ///Prandtl number
        autoPtr<dimensionedScalar> _Pr;

        /// Turbulent thermal diffusivity
        autoPtr<volScalarField> _alphat;

        /// @brief      BT added matrix for the turbulence treatement
        ///
        /// @param[in]  NUmodes    The N of velocity modes.s
        /// @param[in]  NSUPmodes  The N of supremizer modes.
        ///
        /// @return     BT matrix for turbulence treatment
        ///
        Eigen::MatrixXd BT_turbulence(label NUmodes, label NSUPmodes);

        ///
        /// @brief      CT1 added matrix for the turbulence treatement
        ///
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes.
        /// @param[in]  Nnutmodes  The number of eddy viscosity modes.
        ///
        /// @return     CT1 matrix for turbulence treatment
        ///
        List < Eigen::MatrixXd > turbulence_term1(label NUmodes, label NSUPmodes,
                label Nnutmodes);

        ///
        /// @brief      S added matrix for the thermal turbulence treatement
        ///
        /// @param[in]  NTmodes    The number of temperature modes.
        /// @param[in]  Nnutmodes  The number of eddy viscosity modes.
        ///
        /// @return     S matrix for turbulence treatment
        ///
        List < Eigen::MatrixXd > temperature_turbulence_term(label NTmodes,
                label Nnutmodes);

        ///
        /// @brief      CT2 added matrix for the turbulence treatement
        ///
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes.
        /// @param[in]  Nnutmodes  The number of eddy viscosity modes.
        ///
        /// @return     CT2 matrix for turbulence treatment
        ///
        List < Eigen::MatrixXd > turbulence_term2(label NUmodes, label NSUPmodes,
                label Nnutmodes);

        //--------------------------------------------------------------------------
        /// @brief      Perform a truthsolve
        /// @param[in]  mu_now  The actual value of the parameter for this truthSolve. Used only
        /// to construct mu_interp matrix which is written out in the Offline folder, also for par file in
        /// the Parameters folder.
        ///
        void truthSolve(List<scalar> mu_now);

        /// Function to check if the solution must be exported.
        ///
        /// @param      timeObject  The time object of OpenFOAM.
        ///
        /// @return     1 if we must write 0 elsewhere.
        ///
        bool checkWrite(Time& timeObject);

        /// Perform a lift solve for temperature field
        void liftSolveT();

        ///
        /// Project using a supremizer approach
        ///
        /// @param[in]  folder     The folder
        /// @param[in]  NUmodes    The n umodes
        /// @param[in]  NPmodes    The n pmodes
        /// @param[in]  NSUPmodes  The nsu pmodes
        /// @param[in]  Nnutmodes  The nnutmodes
        ///
        void projectSUP(fileName folder, label NU, label NP, label NSUP, label Nnut,
                        label NT);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif










