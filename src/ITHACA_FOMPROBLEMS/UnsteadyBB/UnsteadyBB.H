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
    UnsteadyBB
Description
    Reduction class for a Stationary Navier-Stokes problem.
SourceFiles
    UnsteadyBB.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the UnsteadyBB class.
/// \dir
/// Directory containing the header and source files for the UnsteadyBB class.

#ifndef UnsteadyBB_H
#define UnsteadyBB_H
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "pisoControl.H"
#include "fvOptions.H"
#include "unsteadyNS.H"
#include <iostream>


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                        Class SteadyBB Declaration
\*---------------------------------------------------------------------------*/

/// Implementation of a parametrized full order <b> unsteady Boussinesq problem </b> and preparation of the the reduced matrices for the online solve.
/** In this class are implemented the methods for the offline solve of a unsteady Boussinesq problem (two way coupling between momentum and energy equation)
and the for the generation of the reduced matrices for subsequent online solve, this class is a son
of the reduction problem class */
class UnsteadyBB: public unsteadyNS
{
    public:
        // Constructors
        /// Null constructor
        UnsteadyBB();

        /// Construct with argc and argv
        UnsteadyBB(int argc, char* argv[]);
        ~UnsteadyBB() {};

        // Member Functions

        /// List of pointers used to form the shifted pressure snapshots matrix
        PtrList<volScalarField> Prghfield;

        /// List of pointers used to form the temperature modes
        PtrList<volScalarField> Tmodes;

        /// List of pointers used to form the shifted pressure modes
        PtrList<volScalarField> Prghmodes;

        /// List of pointers used to form the temperature snapshots matrix
        PtrList<volScalarField> Tfield;

        /// List of pointers used to form the temperature snapshots matrix
        PtrList<volScalarField> Tfield_on;

        /// List of pointers used to form the temperature snapshots matrix
        PtrList<volScalarField> Pfield_on;

        /// List of pointers used to form the temperature snapshots matrix
        PtrList<volVectorField> Ufield_on;

        /// List of pointers used to form the list of lifting functions
        PtrList<volScalarField> liftfieldT;

        /// List of pointers used to form the homogeneous velocity snapshots
        PtrList<volScalarField> Tomfield;

        /// List of pointers containing the lift for temperature and the temperature field
        PtrList<volScalarField> L_T_modes;

        // Dummy variables to transform simplefoam into a class
        autoPtr<volScalarField> _UliftBC;

        /// Number of temperature modes used for the projection
        label NTmodes;

        /// Number of pressure modes used for the projection
        label NPrghmodes;

        /// Buoyancy term - momentum equation
        Eigen::MatrixXd H_matrix;

        /// Buoyancy term - PPE equation
        Eigen::MatrixXd HP_matrix;

        /// Diffusive term - energy equation
        Eigen::MatrixXd Y_matrix;

        /// Non linear convective term - energy equation
        List <Eigen::MatrixXd> Q_matrix;
        Eigen::Tensor<double, 3 > Q_tensor;

        /// Mass Matrix - energy equation
        Eigen::MatrixXd W_matrix;

        /// Mesh
        mutable autoPtr<fvMesh> _mesh;

        /// Shifted Pressure field
        autoPtr<volScalarField> _p_rgh;

        /// Temperature field
        autoPtr<volScalarField> _T;

        ///dimensionedScalar beta;
        autoPtr<dimensionedScalar> _beta;

        ///dimensionedScalar Tref;
        autoPtr<dimensionedScalar> _TRef;

        ///dimensionedScalar Pr;
        autoPtr<dimensionedScalar> _Pr;

        ///dimensionedScalar Prt;
        autoPtr<dimensionedScalar> _Prt;

        ///dimensionedScalar alphat;
        autoPtr<volScalarField> _alphat;

        ///dimensionedScalar nut;
        autoPtr<volScalarField> _nut;

        ///dimensionedScalar rhok;
        autoPtr<volScalarField> _rhok;

        ///dimensionedScalar nu;
        autoPtr<dimensionedScalar> _nu;

        /// List of pointers used to form the gravitational acceleration
        autoPtr<volScalarField> _gh;

        /// List of pointers used to form the gravitational acceleration
        autoPtr<surfaceScalarField> _ghf;

        autoPtr<dimensionedVector> _g;

        ///dimensionedScalar hRef;
        autoPtr<dimensionedScalar> _hRef;

        ///dimensionedScalar ghRef;
        autoPtr<dimensionedScalar> _ghRef;

        /// Source Term Heat
        autoPtr<volScalarField> _S;

        /// Pointers for perfoming Nested-POD Method
        /// List of pointers used to form the Weighted nested-POD modes for             velocity field
        PtrList<volVectorField> UModesWeighted;

        /// List of pointers used to form the Weighted nested-POD modes for             pressure field
        PtrList<volScalarField> PModesWeighted;

        /// List of pointers used to form the Weighted nested-POD modes for             temperature field
        PtrList<volScalarField> TModesWeighted;

        /// List of pointers used to form the Weighted nested-POD modes for             nut-field
        PtrList<volScalarField> NUTModesWeighted;

        // Functions
        /// Perform a truthsolve for parameters mu_now
        void truthSolve(List<scalar> mu_now);

        /// Perform a truthsolve for full order solution
        void truthSolve(fileName folder = "./ITHACAOutput/Offline");

        ///
        /// solve the supremizer either with the use of the pressure snaphots or the pressure modes
        ///
        /// @param[in]  type  The type of the supremizer approach, either done on the pressure snapshots or on pressure modes.
        ///
        void solvesupremizer(word type = "snapshots");

        /// Perform a lift solve for temperature
        void liftSolveT();

        /// Perform a lift solve for velocity field
        void liftSolve();

        // Wrapped Proj. Methods;

        //--------------------------------------------------------------------------
        /// Project using a supremizer approach
        ///
        /// @param[in]  folder     The folder used to save the reduced matrices.
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NPmodes    The number of pressure modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes.
        ///
        void projectSUP(fileName folder, label NUmodes, label NPmodes, label NTmodes,
                        label NSUPmodes);

        //--------------------------------------------------------------------------
        /// Project using a PPE approach
        ///
        /// @param[in]  folder     The folder used to save the reduced matrices.
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NPmodes    The number of pressure modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes.
        ///
        void projectPPE(fileName folder, label NUmodes, label NPrghmodes, label NTmodes,
                        label NSUPmodes);
        //--------------------------------------------------------------------------
        /// Gradient of pressure
        ///
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NPmodes    The number of pressure modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes.
        ///
        /// @return     reduced matrix in Eigen::MatrixXd format for the Gradient of pressure term.
        ///
        Eigen::MatrixXd pressure_gradient_term(label NUmodes, label NPrghmodes,
                                               label NSUPmodes);

        //--------------------------------------------------------------------------
        /// Diffusive Term Energy Equation
        ///
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NTmodes    The number of temperature modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes.
        ///
        /// @return     reduced matrix in Eigen::MatrixXd format for the diffusion term.
        ///
        Eigen::MatrixXd diffusive_term_temperature(label NUmodes, label NTmodes,
                label NSUPmodes);
        //--------------------------------------------------------------------------

        /// Divergence Term (supremizer approach)
        ///
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NPmodes    The number of pressure modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes.
        ///
        /// @return     reduced matrix in Eigen::MatrixXd format for the divergence term.
        ///
        Eigen::MatrixXd divergence_term(label NUmodes, label NPrghmodes,
                                        label NSUPmodes);

        //--------------------------------------------------------------------------
        /// Convective Term Energy Equation
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
        /// Mass Term Energy Equation
        ///
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NTmodes    The number of temperature modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes.
        ///
        /// @return     reduced matrix in Eigen::MatrixXd format for the mass matrix.
        ///
        Eigen::MatrixXd mass_term_temperature(label NUmodes, label NTmodes,
                                              label NSUPmodes);

        //--------------------------------------------------------------------------
        /// Buoyant Term Momentum Equation
        ///
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NTmodes    The number of temperature modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes.
        ///
        /// @return     reduced third order tensor in List <Eigen::MatrixXd> format for the convective term.
        ///
        Eigen::MatrixXd buoyant_term(label NUmodes, label NTmodes, label NSUPmodes);

        //--------------------------------------------------------------------------
        /// Buoyant Term PPE Equation
        ///
        /// @param[in]  NPrghmodes    The number of Prgh Pressure modes.
        /// @param[in]  NTmodes    The number of temperature modes.
        ///
        /// @return     reduced matrix in List <Eigen::MatrixXd> format for the buoyancy-related PPE term.
        ///
        Eigen::MatrixXd buoyant_term_poisson(label NPrghmodes, label NTmodes);

        //--------------------------------------------------------------------------
        /// Function to check if the solution must be exported.
        ///
        /// @param      timeObject  The time object of OpenFOAM.
        ///
        /// @return     1 if we must write 0 elsewhere.
        ///
        bool checkWrite(Time& timeObject);

        //--------------------------------------------------------------------------
        /// Function to change the viscosity
        ///
        /// @param[in]  mu    viscosity (scalar)
        ///
        void change_viscosity(double mu);

        //--------------------------------------------------------------------------
        /// @brief      Export convective term energy equation as a tensor
        ///
        /// @param[in]  NUmodes    The N of velocity modes
        /// @param[in]  NTmodes    The N of temperature modes
        /// @param[in]  NSUPmodes  The N of supremizer modes
        ///
        /// @return     tensor_Q
        ///
        Eigen::Tensor<double, 3 > convective_term_tens_temperature(label NUmodes,
                label NTmodes, label NSUPmodes);
};

#endif










