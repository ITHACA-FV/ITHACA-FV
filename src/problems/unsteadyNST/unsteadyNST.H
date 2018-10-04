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
    Reduction class for a NON-Stationary NS problem

SourceFiles
    unsteadyNS.C

\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the unsteadyNST class.
/// \dir
/// Directory containing the header and source files for the unsteadyNS class.

#ifndef unsteadyNST_H
#define unsteadyNST_H
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pisoControl.H"
#include "fvOptions.H"
#include "IOporosityModelList.H"
#include "IOMRFZoneList.H"
#include "fixedFluxPressureFvPatchScalarField.H"
#include "unsteadyNS.H"
#include <iostream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                        Class SteadyNS Declaration
\*---------------------------------------------------------------------------*/

/// Implementation of a parametrized full order <b> unsteady NS problem </b> and preparation of the the reduced matrices for the online solve.
/** In this class are implemented the methods for the offline solve of a unsteady NS problem
and the for the generation of the reduced matrices for subsequent online solve, this class is a son
of the steadyNS class */
class unsteadyNST: public unsteadyNS
{

    public:
        /// Constructors

        /// Construct Null
        unsteadyNST();

        /// Construct with argc and argv
        unsteadyNST(int argc, char* argv[]);
        ~unsteadyNST() {};

        // Specific variable for the unstationary case
        /// Start Time (initial time to start storing the snapshots)
        scalar startTime;

        /// Final time (final time of the simulation and consequently of the acquisition of the snapshots)
        scalar finalTime;

        /// Time step of the simulation
        scalar timeStep;

        /// Time step of the writing procedure
        scalar writeEvery = timeStep;

        /// Auxiliary variable to store the next writing instant
        scalar nextWrite;

        /// adjustTimeStep
        bool adjustTimeStep;

        /// maxCourant
        scalar maxCo;

        /// maxDeltaT
        scalar maxDeltaT;

        /// List of pointers used to form the temperature snapshots matrix
        PtrList<volScalarField> Tfield;

        /// Temperature field
        autoPtr<volScalarField> _T;

        /// List of pointers used to form the temperature modes
        PtrList<volScalarField> LTmodes;

        ///dimensionedScalar DT;
        autoPtr<dimensionedScalar> _DT;

        ///dimensionedScalar nu;
        autoPtr<dimensionedScalar> _nu;

        /// List of pointers used to form the list of the temperature lifting functions
        PtrList<volScalarField> liftfieldT;

        /// List of pointers used to form the homogeneous temperature snapshots
        PtrList<volScalarField> Tomfield;

        /// pisoControl
        autoPtr<pisoControl> _piso;

        /// Reference pressure cell
        label pRefCell;

        /// Reference pressure value
        scalar pRefValue;

        /// Mesh
        mutable autoPtr<fvMesh> _mesh;

        /// Number of temperature modes used for the projection
        label NTmodes;

        /// Gradient of pressure matrix
        Eigen::MatrixXd Y_matrix;

        /// Non linear term
        List <Eigen::MatrixXd> Q_matrix;

        /// Mass Matrix T
        Eigen::MatrixXd MT_matrix;

        /// Time
        autoPtr<Time> _runTime;

        /// Perform a lift solve for velocity field
        void liftSolve();

        /// Perform a lift solve for temperature field
        void liftSolveT();

        /// Perform a truthsolve
        void truthSolve(List<scalar> mu_now);

        bool checkWrite(Time& timeObject);

        /// Specific variable for the unstationary case

        //--------------------------------------------------------------------------
        /// Project using the Poisson Equation for pressure
        ///
        /// @param[in]  folder     The folder used to save the reduced matrices.
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NPmodes    The number of pressure modes.
        /// @param[in]  NPmodes    The number of pressure modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes, 0 by default and it must remain always zero.
        ///
        void projectSUP(fileName folder, label NUmodes, label NPmodes, label NTmodes,
                        label NSUPmodes);

        //--------------------------------------------------------------------------
        /// Convective Term for Temperature
        ///
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NTmodes    The number of temperature modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes.
        ///
        /// @return     reduced third order tensor in List <Eigen::MatrixXd> format for the convective term.
        ///
        List< Eigen::MatrixXd > convective_term_temperature(label NUmodes,
                label NTmodes, label NSUPmodes);

        //--------------------------------------------------------------------------
        /// Diffusive Term for Temperature
        ///
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NTmodes    The number of temperature modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes.
        ///
        /// @return     reduced matrix in Eigen::MatrixXd format for the diffusion term for the temperature equation.
        ///
        Eigen::MatrixXd diffusive_term_temperature(label NUmodes,  label NTmodes,
                label NSUPmodes);

        //--------------------------------------------------------------------------
        /// Mass Term for Temperature
        ///
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NTmodes    The number of temperature modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes.
        ///
        /// @return     reduced matrix in Eigen::MatrixXd format for the mass matrix for the temperature equation.
        ///
        Eigen::MatrixXd mass_term_temperature(label NUmodes,  label NTmodes,
                                              label NSUPmodes);


};

#endif



