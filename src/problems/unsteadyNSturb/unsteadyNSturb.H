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
    unsteadyNSturb

Description
    Reduction class for a NON-Stationary turbulent NS problem

SourceFiles
    unsteadyNSturb.C   

\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the unsteadyNSturb class.
/// \dir
/// Directory containing the header and source files for the unsteadyNSturb class.

#ifndef unsteadyNSturb_H
#define unsteadyNSturb_H
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
class unsteadyNSturb: public unsteadyNS
{
public:
    // Constructors
        /// Construct Null
    unsteadyNSturb();

        /// Construct with argc and argv
    unsteadyNSturb(int argc, char *argv[]);

    // Specific variable for the unstationary case
        // Member Functions
    /// List of snapshots for the solution for eddy viscosity 
    PtrList<volScalarField> nutFields;
    
    /// List of POD modes for eddy viscosity 
    PtrList<volScalarField> nuTmodes;

    // Create a Rbf splines for interpolation
    std::vector<SPLINTER::DataTable*> SAMPLES;
    // Create a SAMPLES for interpolation
    std::vector<SPLINTER::RBFSpline*> rbfsplines;
 
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

    /// Total C Matrix
    List <Eigen::MatrixXd> C_total_matrix;

    /// Total B Matrix
    Eigen::MatrixXd B_total_matrix;
    ///@}

    /// Number of viscoisty modes used for the projection
    label Nnutmodes;

    
    /// Eddy viscosity field
    autoPtr<volScalarField> _nut;
    autoPtr<volScalarField> _nuTilda;



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
    List < Eigen::MatrixXd > turbulence_term1(label NUmodes, label NSUPmodes, label Nnutmodes);

    ///
    /// @brief      CT2 added matrix for the turbulence treatement
    ///
    /// @param[in]  NUmodes    The number of velocity modes.
    /// @param[in]  NSUPmodes  The number of supremizer modes.
    /// @param[in]  Nnutmodes  The number of eddy viscosity modes.
    ///
    /// @return     CT2 matrix for turbulence treatment
    ///
    List < Eigen::MatrixXd > turbulence_term2(label NUmodes, label NSUPmodes, label Nnutmodes);

    void truthSolve(List<scalar> mu_now);


    /**
     * @brief      { function_description }
     *
     * @param[in]  folder     The folder
     * @param[in]  NUmodes    The number of velocity modes.
     * @param[in]  NPmodes    The number of pressure modes.
     * @param[in]  NSUPmodes  The number of supremizer modes.
     * @param[in]  Nnutmodes  The number of eddy viscosity modes.
     */
    void projectSUP(fileName folder, label NUmodes, label NPmodes, label NSUPmodes, label Nnutmodes);
    
    void projectPPE(fileName folder, label NUmodes, label NPmodes, label NSUPmodes, label Nnutmodes);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif










