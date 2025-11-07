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

#include "Smagorinsky.H"
#include "kOmegaSSTBase.H"
#include "DESModelBase.H"
#include "DESModel.H"
#include "kOmegaSSTDES.H"
#include "PODParameters.H"
#include <optional>

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

        // Small construct with Parameters to access member functions linked to Smagorinsky diffusion 
        UnsteadyNSTurb(const ITHACAPOD::Parameters* myParameters);

        // Specific variable for the unstationary case
        // Member Functions
        /// List of snapshots for the solution for eddy viscosity
        PtrList<volScalarField> nutFields;

        /// List of POD modes for eddy viscosity
        volScalarModes nutModes;

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

        /// RBF shape parameters vector
        Eigen::VectorXd radii;

        /// RBF functions radius
        double e = 1;

        /// Number of viscoisty modes used for the projection
        label nNutModes;

        /// Interpolation independent variable choice
        label interChoice = 1;

        /// Eddy viscosity field
        autoPtr<volScalarField> _nut;

        /// Face velocity field
        autoPtr<surfaceVectorField> _Uf;

        /// Pressure reference cell
        label _pRefCell;

        /// Pressure reference value
        scalar _pRefValue;

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
        /// @param[in]  offlinepath Path where solution should be stored
        ///
        void truthSolve(List<scalar> mu_now, std::string& offlinepath);

        //----------------------------------------------------------------------
        /// Project using a supremizer approach
        ///
        /// @param[in]  folder     The folder name
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NPmodes    The number of pressure modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes.
        /// @param[in]  nNutModes  The number of eddy viscosity modes.
        /// @param[in]  rbfInterp  A boolean variable that determines whether to carry out the RBF interpolation or not
        ///
        void projectSUP(fileName folder, label NUmodes, label NPmodes, label NSUPmodes,
                        label nNutModes, bool rbfInterp = true);

        //----------------------------------------------------------------------
        /// Project using the Poisson Equation for pressure
        ///
        /// @param[in]  folder     The folder name
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NPmodes    The number of pressure modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes which will be
        ///                        set to zero.
        /// @param[in]  nNutModes  The number of eddy viscosity modes.
        /// @param[in]  rbfInterp  A boolean variable that determines whether to carry out the RBF interpolation or not
        ///
        void projectPPE(fileName folder, label NUmodes, label NPmodes, label NSUPmodes,
                        label nNutModes, bool rbfInterp = true);

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

private :

        /// Parameter object containing main information
        const ITHACAPOD::PODParameters* ithacaFVParameters;

        // /// Turbulence model
        // autoPtr<incompressible::turbulenceModel> turbModel;
    
 public :

    // =========================================================================
    // ============================= Smagorinsky ===============================
    // =========================================================================

        //--------------------------------------------------------------------------
        /// @brief   Compute non linear Snapshot at a given time
        /// @details Overloaded function depending on the type of the second argument
        ///          If volVectorField then it computes the Smagorinsky term
        ///
        /// @param[in]     index    label defining concerned time step
        /// @param[in,out] Smag     volVectorField in wich the Smagorinsky term will be computed
        /// @param[in]     modesU   PtrList<volVectorField>, (optional) if given then compute 
        ///                         Smag from reduced velocity
        ///
        /// @return  Smag
        ///
        void computeNonLinearSnapshot_at_time(const word& snap_time, volVectorField& Smag, 
                        std::optional<PtrList<volVectorField>> modesU = std::nullopt);

        //--------------------------------------------------------------------------
        /// @brief   Initialize Smagorinsky term with values at time 0
        /// @details Smagorinsky at time 0 is not provided by OpenFOAM, it is required to compute this quantity
        ///
        /// @return  Smag
        ///
        volVectorField initSmagFunction();

        //--------------------------------------------------------------------------
        /// @brief   Initialize Smagorinsky term of a tracer phi with values at time 0
        /// @details Smagorinsky of a tracer phi at time 0 is not provided by OpenFOAM, it is required to compute this quantity
        ///
        /// @return  Smag
        ///
        volScalarField initSmagPhiFunction(const volScalarField template_field_phi);

        //--------------------------------------------------------------------------
        /// @brief   Compute Smagorinsky Snapshot at a given time
        ///
        /// @param[in] index    label defining concerned time step term will be computed
        /// @param[in] modesU   PtrList<volVectorField> (optional) storing the velocity modes, 
        ///                     if given then compute Smag from reduced velocity
        ///
        /// @return  Smag
        ///
        volVectorField computeSmagTerm_at_time(const word& snap_time, std::optional<PtrList<volVectorField>> modesU = std::nullopt);

        //--------------------------------------------------------------------------
        /// @brief   Compute Smagorinsky Snapshot of a tracer phi at a given time
        ///
        /// @param[in] index               label defining concerned time step term will be computed
        /// @param[in] template_field_phi  volScalarField, template field of the tracer
        ///
        /// @return  Smagorinsky volScalarField
        ///
        volScalarField computeSmagTermPhi_at_time(const word& snap_time, const volScalarField template_field_phi);

        //--------------------------------------------------------------------------
        /// @brief   Compute Smagorinsky Snapshot for a given velocity field U
        ///
        /// @param[in] snapshotj volVectorField of snapshot U
        ///
        /// @return  Smagorinsky volVectorField
        ///
        volVectorField computeSmagTerm_fromU(const volVectorField& snapshotj);

        //--------------------------------------------------------------------------
        /// @brief  Compute Smagorinsky Snapshot for given velocity field U and tracer phi
        ///
        /// @param[in] snapshotj volVectorField of snapshot U
        /// @param[in] phij      volScalarField of snapshot phi
        ///
        /// @return  Smagorinsky volScalarField
        ///
        volScalarField computeSmagTermPhi_fromUPhi(const volVectorField& snapshotj, const volScalarField& phij);

        //--------------------------------------------------------------------------
        /// @brief  Compute Smagorinsky Snapshot for a given velocity field U (=computeSmagTerm_fromU)
        ///
        /// @param[in] snapshotj volVectorField of snapshot U
        /// @param[in] phij      volScalarField of snapshot phi (unsused, assumed to be = U)
        ///
        /// @return  Smagorinsky volScalarField
        ///
        volVectorField computeSmagTermPhi_fromUPhi(const volVectorField& snapshotj, const volVectorField& phij);

        //--------------------------------------------------------------------------
        /// @brief  Compute diffusion term for a given scalar field phi
        ///
        /// @param[in] coefDiff  volScalarField or volTensorField of diffusion coefficient or matrix
        /// @param[in] phi       volScalarField of snapshot phi
        ///
        /// @return  diffusion term volScalarField
        ///
        template<typename T>
        static volScalarField diffusion(const T& coefDiff, const volScalarField& phi);

        //--------------------------------------------------------------------------
        /// @brief  Compute diffusion term for a given vector field u
        ///
        /// @param[in] coefDiff  volScalarField or volTensorField of diffusion coefficient or matrix
        /// @param[in] u         volVectorField of snapshot u
        ///
        /// @return  diffusion term volVectorField
        ///
        template<typename T>
        static volVectorField diffusion(const T& coefDiff, const volVectorField& u);


    // =========================================================================
    // ======================= projFullStressFunction ==========================
    // =========================================================================

        //--------------------------------------------------------------------------
        /// @brief  Compute non linear Snapshot at a given time
        /// @details computes the projection of the Smagorinsky term
        ///          on the velocity mode given in the third argument
        ///          Other definitions to enable calling the function from a templatized class
        ///
        /// @param[in]     index      label defining concerned time step
        /// @param[in,out] phi        volScalarField in wich the projected Smagorinsky term will be computed
        /// @param[in]     modeU volVectorField containing a velocity mode to project the Smagorinsky term
        /// @param[in]     modeU_grad volVectorField containing a velocity mode to compute the symmetric gradient for Smag
        ///
        /// @return  projSmagFromNut volScalarField
        ///
        void computeNonLinearSnapshot_at_time(const word& snap_time, volScalarField& phi, volVectorField& modeU);

        void computeNonLinearSnapshot_at_time(const word& snap_time, volVectorField& phi, volVectorField& modeU)
        {
            Info << "Output field with this call to computeNonLinearSnapshot_at_time is scalar" << endl;
            abort();
        };

        //--------------------------------------------------------------------------
        /// @brief  Initialize projFullStressFunction with values at time 0
        ///
        /// @return  template field for projFullStressFunction
        ///
        volScalarField initProjSmagFunction();

        //--------------------------------------------------------------------------
        /// @brief  Compute projFullStressFunction snapshot at a given time projected on the given mode
        ///
        /// @param[in] index label defining concerned time step
        /// @param[in] mode  volVectorField containing the field on which the Smagorinsky term is projected
        ///
        /// @return  projFullStressFunction volScalarField
        ///
        volScalarField computeProjSmagTerm_at_time_fromMode(const word& snap_time, const volVectorField& mode);

        //--------------------------------------------------------------------------
        /// @brief  Compute projFullStressFunction given the Smagorinsky term and the velocity mode
        /// @details Other definitions to enable calling the function from a templatized class
        ///
        /// @param[in,out]  phi volScalarField
        /// @param[in]      Smag volVectorField containing the Smagorinsky term
        /// @param[in]      mode volVectorField containing the field on which the Smagorinsky term is projected
        ///
        /// @return  projFullStressFunction volScalarField
        ///
        void computeProjSmagTerm_fromSmag_fromMode(volScalarField& phi, const volVectorField& Smag, const volVectorField& mode);

        void computeProjSmagTerm_fromSmag_fromMode(volScalarField phi,const volScalarField& Smag, const volVectorField& mode)
        {
            Info << "Input fields for this call to computeProjSmagTerm_fromSmag_fromMode are vectors" << endl;
            abort();
        };

        void computeProjSmagTerm_fromSmag_fromMode(volVectorField phi, const volVectorField& Smag, const volVectorField& mode)
        {
            Info << "Output field for this call to computeProjSmagTerm_fromSmag_fromMode is scalar" << endl;
            abort();
        };


    // =========================================================================
    // =========================== projSmagFromNut =============================
    // =========================================================================

        //--------------------------------------------------------------------------
        /// @brief  Compute non linear Snapshot at a given time
        /// @details Computes the projection on a velocity mode of the Smagorinsky term computed from 
        ///          a nut snapshot and a velocity mode for the symmetric gradient.
        ///          Other definitions to enable calling the function from a templatized class
        ///
        /// @param[in]     index      label defining concerned time step
        /// @param[in,out] phi        volScalarField in wich the projected Smagorinsky term will be computed
        /// @param[in]     modeU_proj volVectorField containing the field on which the Smagorinsky term is projected
        /// @param[in]     modeU_grad volVectorField containing a velocity mode to compute the symmetric gradient for Smag
        ///
        /// @return  projSmagFromNut volScalarField
        ///
        void computeNonLinearSnapshot_at_time(const word& snap_time, volScalarField& phi, 
                                volVectorField& modeU_proj , volVectorField& modeU_grad);

        void computeNonLinearSnapshot_at_time(const word& snap_time, volVectorField& phi, 
                                volVectorField& modeU_proj , volVectorField& modeU_grad)
        {
            Info << "Output field with this call to computeNonLinearSnapshot_at_time is scalar" << endl;
            abort();
        };

        //--------------------------------------------------------------------------
        /// @brief  Initialize projSmagFromNut with values at time 0
        ///
        /// @return  template field for projSmagFromNut
        ///
        volScalarField initProjSmagFromNutFunction();

        //--------------------------------------------------------------------------
        /// @brief  Compute projSmagFromNut snapshot at a given time from 2 velocity modes
        /// @details Other definitions to enable calling the function from a templatized class
        ///
        /// @param[in]  index      label defining concerned time step
        /// @param[in]  modeU_proj volVectorField containing the field on which the Smagorinsky term is projected
        /// @param[in]  modeU_grad volVectorField containing a velocity mode to compute the symmetric gradient for Smag
        ///
        /// @return  projSmagFromNut volScalarField
        ///
        volScalarField computeProjSmagFromNut_at_time_fromModes(const word& snap_time, const volVectorField& modeU_proj, 
                                                           const volVectorField& modeU_grad);

        //--------------------------------------------------------------------------
        /// @brief  Compute projSmagFromNut given Nut and 2 velocity modes
        /// @details Other definitions to enable calling the function from a templatized class
        ///
        /// @param[in,out]  phi        volScalarField storing the result
        /// @param[in]      Smag       volVectorField containing the Smagorinsky term
        /// @param[in]      modeU_proj volVectorField containing the field on which the Smagorinsky term is projected
        /// @param[in]      modeU_grad volVectorField containing a velocity mode to compute the symmetric gradient for Smag
        ///
        /// @return  projSmagFromNut volScalarField
        ///
        void computeProjSmagFromNut_fromNut_fromModes(volScalarField& phi, const volScalarField& Nut, 
                                const volVectorField& modeU_proj, const volVectorField& modeU_grad);

        void computeProjSmagFromNut_fromNut_fromModes(volVectorField& phi, const volVectorField& Nut, 
                                const volVectorField& modeU_proj, const volVectorField& modeU_grad)
        {
            Info << "Output field for this call to computeProjSmagFromNut_fromNut_fromModes is a scalar field" << endl;
            abort();
        };

        void computeProjSmagFromNut_fromNut_fromModes(volScalarField& phi, const volVectorField& Nut, 
                                const volVectorField& modeU_proj, const volVectorField& modeU_grad)
        {
            Info << "Input field for this call to computeProjSmagFromNut_fromNut_fromModes is a scalar field" << endl;
            abort();
        };

        //--------------------------------------------------------------------------
        /// @brief  Compute projected diffusion term with the integration by parts (IBP) equation for a given vector field u
        ///
        /// @param[in] coefDiff volScalarField of diffusion coefficient
        /// @param[in] u volVectorField of velocity for the diffusion
        /// @param[in] v volVectorField of velocity mode for the projection
        ///
        /// @return  volScalarField of non integrated projected diffusion  
        ///
        static volScalarField projDiffusionIBP(const volScalarField& coefDiff, 
                             const volVectorField& u, const volVectorField& v);



    // =========================================================================
    // ============================= NUT =======================================
    // =========================================================================

        //--------------------------------------------------------------------------
        /// @brief  Compute non linear Snapshot at a given time
        /// @details Overloaded function depending on the type of the second argument
        ///          If volScalarField then it computes the turbulent viscosity term
        ///
        /// @param[in]     index  label defining concerned time step
        /// @param[in,out] nut    volScalarField in wich the turbulent viscosity term will be computed
        /// @param[in]     modesU PtrList<volVectorField>, optional, if given then compute nut from reduced velocity
        ///
        /// @return  nut
        ///
        void computeNonLinearSnapshot_at_time(const word& snap_time, volScalarField& nut, 
              std::optional<PtrList<volVectorField>> modesU = std::nullopt);

        //--------------------------------------------------------------------------
        /// @brief  Initialize Nut term with values at time 0
        /// @details Nut is provided by OpenFOAM, but for homogeneity purposes 
        ///          it is initialized the same way as Smagorinsky term
        ///
        /// @return  nonLinearSnapshot
        ///
        volScalarField initNutFunction();

        //--------------------------------------------------------------------------
        /// @brief  Compute non linear Snapshot at a given time
        /// @details Overloaded function depending on the type of the second argument
        ///          If volScalarField then it computes the turbulent viscosity term
        ///
        /// @param[in] index label defining concerned time step
        /// @param[in] modesU PtrList<volVectorField>, optional, if given then compute nut from reduced velocity
        ///
        /// @return  nut volScalarField
        ///
        volScalarField computeNut_at_time(const word& snap_time, std::optional<PtrList<volVectorField>> modesU = std::nullopt);

        //--------------------------------------------------------------------------
        /// @brief  Compute turbulent viscosity Snapshot for a given velocity field U
        ///
        /// @param[in] snapshotj volVectorField of snapshot U
        ///
        /// @return  nut volScalarField
        ///
        volScalarField computeNut_fromU(const volVectorField& snapshotj);

        //--------------------------------------------------------------------------
        /// @brief   Compute turbulent viscosity Snapshot for a given tensor Field S
        ///
        /// @param[in] S volTensorField  S = 0.5*(gradV+gradV^t)
        ///
        /// @return  nut volScalarField
        ///
        volScalarField computeNut_fromS(const volTensorField& S);

        // void initTurbModel();

        //--------------------------------------------------------------------------
        /// @brief   Compute S deformation tensor for a given vector field U
        ///
        /// @param[in] snapshotj volVectorField U
        ///
        /// @return  S deformation Tensor
        ///
        static volTensorField computeS_fromU(const volVectorField& snapshotj);

        //--------------------------------------------------------------------------
        /// @brief   Compute gradient for a given scalar field Phi
        ///
        /// @param[in] phij volScalarField  (omega or k)
        ///
        /// @return  S grad(phij)
        ///
        static volVectorField computeS_fromU(const volScalarField& phij);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif










