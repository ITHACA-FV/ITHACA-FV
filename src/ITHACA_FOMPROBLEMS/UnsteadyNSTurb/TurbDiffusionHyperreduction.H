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
    TurbDiffusionHyperreduction
Description
    Class to compute and hyperreduce the non-polynomial diffusion for a NON-Stationary turbulent NS problem
SourceFiles
    TurbDiffusionHyperreduction.C

\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the TurbDiffusionHyperreduction class.
/// \dir
/// Directory containing the header and source files for the UnsteadyNSTurb class.

#ifndef TurbDiffusionHyperreduction_H
#define TurbDiffusionHyperreduction_H

#include "UnsteadyNSTurb.H"
#include "PODTemplate.H"
#include "PODTemplateH1.H"
#include "hyperReduction.templates.H"


template<typename volF, typename T, typename S>
class TurbDiffusionHyperreduction
{
  protected :
    /// Parameter object containing main information
    StoredParameters* m_parameters;

    /// Parameter object to evaluate interpolated functions
    autoPtr<UnsteadyNSTurb> ithacaUnsteadyNSTurb;

    /// Hyperreduction objects to choose magic points
    autoPtr<HyperReduction<PtrList<S>>> ithacaHyperreduction;

    /// Number of snapshot, also number of time step
    label l_nSnapshot;

    /// start time of the POD decomposition
    label l_startTime;

    /// Number of Simulation snapshot, also number of time steps
    label l_nSnapshotSimulation;

    /// Path to snapshots, number, start and end time indices for each HR POD. If -1, set automatically by the POD class.
    word snapshotsFolder;
    label POD_nSnapshot = -1;
    label POD_startTime = -1;
    label POD_endTime = -1;

    /// Interpolated field template
    T template_HRInterpField;

    /// Snapshots field template
    S template_HRSnapshotsField;

    /// Pointer list of POD Modes of U
    PtrList<volF> f_spatialModesU; // TO DO : need to be duplicated (local => global) if //

    /// Pointer of mean of U
    volF* f_meanU; // TO DO : need to be duplicated (local => global) if //

    // Spatial POD-HR modes
    PtrList<S> f_spatialModesHR;

    // Temporal modes of HR field
    Eigen::MatrixXd m_temporalModesHR;

    // Test temporal modes of HR field
    Eigen::MatrixXd m_temporalModesHRSimulation;

    // Covariance matrix the resolved modes of HR field
    Eigen::MatrixXd covMatrixHR;

    /// Export folder for HR operators
    word folder_HR;

    /// Dummy Foam::Time object for mesh and submesh management
    Foam::Time runTime2;

    /// Boolean to include or not magic points' neighbourhood in the submesh
    ///  0 if only magicPoints and 1 if include magicNeighborhoods
    bool mgPoints_0_Neighborhoods_1;

    // Name of the HR method: DEIM or ECP
    word HRMethod;

    // Name of the interpolated field: fullStressFunction or nut or their reduced version
    word HRInterpField;

    // Name of the snapshots' field if ECP
    word HRSnapshotsField;

    // Name of the algorithm choice if ECP: Global or EachMode
    word ECPAlgo;

    // Centered field to perform the hyperreduction
    bool interpFieldCenteredOrNot;



  public :

        //--------------------------------------------------------------------------
        /// @brief      Class TurbDiffusionHyperreduction Constructor
        /// @details    This class is used to perform the offline part of hyperreduction on a turbulent non-stationnary diffusion
        ///
        /// @param[in]  myParameters StoredParameters object
        /// @param[in]  mgPoints_0_Neighborhoods_1 boolean  Include neighboroods in submesh or not
        /// @param[in]  template_HRInterpField volScalarField or volVectorField  Template of the interpolated field
        /// @param[in]  template_HRInterpField volScalarField or volVectorField  Template of the snapshots field
        /// @param[in]  meanU volVectorField  Mean of the velocity field
        /// @param[in]  spatialModesU PtrList<volVectorField>  POD modes of the velocity field
        TurbDiffusionHyperreduction(Parameters* myParameters, bool mgPoints_0_Neighborhoods_1,
                                    T& template_HRInterpField, S& template_HRSnapshotsField,
                                    volF*& meanU, PtrList<volF>& spatialModesU);

        //--------------------------------------------------------------------------
        /// @brief      Class TurbDiffusionHyperreduction Destructor
        ~TurbDiffusionHyperreduction();


        //--------------------------------------------------------------------------
        /// @brief   Precompute and save non polynomial field
        /// @details Precompute either the interpolated fields or the snapshot fields depending on the argument
        ///
        /// @param[in]     fieldToCompute word Name of the field to compute
        ///
        void precomputeTurbDiffusionFunctions(word& fieldToCompute);

        //--------------------------------------------------------------------------
        /// @brief   Compute the symmetric tensor of the velocity modes on the magic points submesh
        /// @details Saves the result in the StoredParameters object
        ///
        template<typename volGradF>
        void computeDefTensor(volGradF tensorMeanU);

        //--------------------------------------------------------------------------
        /// @brief   Compute the interpolated field's mean on the magic points submesh
        /// @details Saves the result in the StoredParameters object
        ///
        void common_MeanHRInterpField();
       
        //--------------------------------------------------------------------------
        /// @brief   Compute the hyperreduction by callng other methods
        /// @details Generic function working to hyperreduce any of the possible fields
        ///
        void computeTurbDiffusionHyperreduction();

};


#endif
