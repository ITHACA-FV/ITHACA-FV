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
    ITHACADMD
Description
    Implementation of a DMD of a general field
SourceFiles
    ITHACADMD.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the ITHACADMD class.
/// \dir
/// Directory containing the header and source files for the ITHACADMD class.

#ifndef ITHACADMD_H
#define ITHACADMD_H

#include "fvCFD.H"
#include "ITHACAutilities.H"
#include "ITHACAstream.H"
#include "ITHACAparameters.H"
#include "Foam2Eigen.H"
#include "EigenFunctions.H"
#include "ITHACAPOD.H"
#include <functional>
#include "Modes.H"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wnon-virtual-dtor"
#include <Spectra/GenEigsSolver.h>
#include <Spectra/SymEigsSolver.h>
#include <Eigen/Eigen>
#include <unsupported/Eigen/SparseExtra>
#include <unsupported/Eigen/MatrixFunctions>
#include <redsvd>
#pragma GCC diagnostic pop

/*---------------------------------------------------------------------------*\
                        Class reductionProblem Declaration
\*---------------------------------------------------------------------------*/

/// Class of the computation of the DMD, it exploits the SVD methods.
///
/// @tparam     Field_type  It can be scalar or vector
///
template<class Field_type>
class ITHACADMD
{
    public:

        ///
        /// @brief      Constructs the object.
        ///
        /// @param[in]  snapshots  The snapshots on which you want to perform DMD
        /// @param[in]  dt         The Time Step used to acquire the snapshots
        ///
        ITHACADMD(PtrList<GeometricField<Field_type, fvPatchField, volMesh>>& snapshots,
                  double dt);

        /// PtrList of OpenFOAM GeoometricFields where the snapshots are stored
        PtrList<GeometricField<Field_type, fvPatchField, volMesh>> snapshotsDMD;

        /// Modes object used to store the Real part of the DMD modes
        Modes<Field_type> DMDmodesReal;

        /// Modes object used to store the Imaginary part of the DMD modes
        Modes<Field_type> DMDmodesImag;

        /// Eigenvalues of the dynamics mode decomposition
        Eigen::VectorXcd eigenValues;

        /// DMD modes stored in a complex Eigen::Matrix, it is the object used for computations
        Eigen::MatrixXcd DMDEigenModes;

        /// DMD modes on the boundary stored in a List of complex Eigen::Matrix, it is the object used for computations
        List<Eigen::MatrixXcd> DMDEigenModesBC;

        /// Complex matrix used to store the POD modes, used only for compution in the projected approach
        Eigen::MatrixXcd PODm;

        /// List of complex matrices used to store the POD modes on the boundaries, used only for compution in the projected approach
        List<Eigen::MatrixXcd> PODmBC;

        /// Amplitudes of DMD
        Eigen::VectorXd Amplitudes;

        /// Complex Eigen::Matrix used to store the Dynamics of the DMD modes
        Eigen::MatrixXcd dynamics;

        /// Number of snapshots
        int NSnaps;

        /// Rank of the DMD
        int SVD_rank_public;

        /// Original time step used to acquire the snapshots
        double originalDT;

        /// If true, it uses the Randomized SVD
        bool redSVD;

        //--------------------------------------------------------------------------
        /// Get the DMD modes
        ///
        /// @param[in]  SVD_rank        The svd rank
        /// @param[in]  exact           True (default) if you want the exact DMD modes computations. False if you want the projected modes
        /// @param[in]  exportDMDmodes  A flag to export the DMD modes, both real and imaginary parts are exported
        ///
        void getModes(int SVD_rank = -1, bool exact = true,
                      bool exportDMDmodes = false);

        //--------------------------------------------------------------------------
        /// Convert the EigenModes in Matrix form into OpenFOAM GeometricFields
        ///
        void convert2Foam();

        //--------------------------------------------------------------------------
        /// Export the eigenvalues in numpy format.
        ///
        /// @param[in]  exportFolder  Export folder path
        ///
        void exportEigs(word exportFolder);

        //--------------------------------------------------------------------------
        /// Export the dynamics of DMD given an initial time step, a final one and a time step.
        ///
        /// @param[in]  tStart  Intial time
        /// @param[in]  tFinal  Final time
        /// @param[in]  dt      Time step
        ///
        void getDynamics(double tStart, double tFinal, double dt);

        //--------------------------------------------------------------------------
        /// Reconstruct and export the solution using the computed dynamics
        ///
        void reconstruct(word exportFolder, word fieldName);
};
#endif
