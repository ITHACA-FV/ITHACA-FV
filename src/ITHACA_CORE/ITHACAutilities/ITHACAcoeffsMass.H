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
Namespace
    ITHACAutilites
Description
    Utilities to compute projection coefficients and mass matrices
SourceFiles
    ITHACAcoeffsMass.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the ITHACAcoeffsMass file.
/// \dir
/// Directory containing the header and source files for the ITHACAcoeffsMass files.

#ifndef ITHACAcoeffsMass_H
#define ITHACAcoeffsMass_H

#include "fvCFD.H"
#include "IOmanip.H"
#include "freestreamFvPatchField.H"
#include <sys/stat.h>
#include <unistd.h>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#include <Eigen/Eigen>
#pragma GCC diagnostic pop
#include <functional>
#include "./colormod.H"
#include "polyMeshTools.H"
#include <chrono>
#include "mixedFvPatchFields.H"
#include "fvMeshSubset.H"
using namespace std::placeholders;
#include "Foam2Eigen.H"

/// Namespace to implement projection coefficients calculation and mass matrices
namespace ITHACAutilities
{

//--------------------------------------------------------------------------
/// Exact reconstruction using a certain number of modes for a list of
/// fields and the projection coefficients
///
/// @param[in]  modes         The modes used for reconstruction as PtrList.
/// @param[in]  Nmodes        The number of modes you want to use.
/// @param[in]  coeff_matrix  The matrix of coefficients.
///
/// @tparam     Type        vector or scalar.
/// @tparam     PatchField  fvPatchField or fvsPatchField.
/// @tparam     GeoMesh     volMesh or surfaceMesh.
///
/// @return     rec_field     The reconstructed field as PtrList.
///
template<class Type, template<class> class PatchField, class GeoMesh>
PtrList<GeometricField<Type, PatchField, GeoMesh>> reconstructFromCoeff(
            PtrList<GeometricField<Type, PatchField, GeoMesh>>& modes,
            Eigen::MatrixXd& coeff_matrix, label Nmodes);

//-----------------------------------------------------------------------------
/// Gets the mass matrix using the eigen routine
///
/// @param[in]  fields            The openFoam fields
/// @param[in]  consider_volumes  True if you want to consider the volumes
///
/// @tparam     Type        vector or scalar.
/// @tparam     PatchField  fvPatchField or fvsPatchField.
/// @tparam     GeoMesh     volMesh or surfaceMesh.
///
/// @return     The mass matrix eigen.
///
template<class Type, template<class> class PatchField, class GeoMesh>
Eigen::MatrixXd getMassMatrix(PtrList<GeometricField<Type, PatchField, GeoMesh>>&
                              fields, int nModes = 0,
                              bool consider_volumes = true);

//--------------------------------------------------------------------------
/// Gets a vector containing the volumes of each cell of the mesh.
///
/// @param[in]  snapshot   A snapshot, scalar or vector field.
///
/// @tparam     Type        vector or scalar.
/// @tparam     PatchField  fvPatchField or fvsPatchField.
/// @tparam     GeoMesh     volMesh or surfaceMesh.
///
/// @return     The FV mass matrix in Eigen::MatrixXd form.
///
template<class Type, template<class> class PatchField, class GeoMesh>
Eigen::VectorXd getMassMatrixFV(
    GeometricField<Type, PatchField, GeoMesh>& snapshot);

//------------------------------------------------------------------------------
/// Projects a snapshot on a basis function and gets
/// the coefficients of the projection
///
/// @param[in]  snapshot          The snapshot.
/// @param[in]  modes             The modes.
/// @param[in]  Nmodes            The number of modes
/// @param[in]  consider_volumes  Consider volumes or not, in not equals Frobenius Projection
///
/// @tparam     Type        vector or scalar.
/// @tparam     PatchField  fvPatchField or fvsPatchField.
/// @tparam     GeoMesh     volMesh or surfaceMesh.
///
/// @return     The coefficients of the projection.
///
template<class Type, template<class> class PatchField, class GeoMesh>
Eigen::VectorXd getCoeffs(GeometricField<Type, PatchField, GeoMesh>&
                          snapshot,
                          PtrList<GeometricField<Type, PatchField, GeoMesh>>& modes, int Nmodes = 0,
                          bool consider_volumes = true);

//------------------------------------------------------------------------------
/// Projects snapshots on a basis function and gets
/// the coefficients of the projection
///
/// @param[in]  snapshot          The snapshots.
/// @param[in]  modes             The modes.
/// @param[in]  Nmodes            The number of modes you want to use
/// @param[in]  consider_volumes  The consider volumes (if not equals Frobenius Projection)
///
/// @tparam     Type        vector or scalar.
/// @tparam     PatchField  fvPatchField or fvsPatchField.
/// @tparam     GeoMesh     volMesh or surfaceMesh.
///
/// @return     The coefficients of the projection.
///
template<class Type, template<class> class PatchField, class GeoMesh>
Eigen::MatrixXd getCoeffs(PtrList<GeometricField<Type, PatchField, GeoMesh>>&
                          snapshot,
                          PtrList<GeometricField<Type, PatchField, GeoMesh>>& modes, int Nmodes = 0,
                          bool consider_volumes = true);

//--------------------------------------------------------------------------
/// @brief      A method to compute the time-parameter combined matrix whose any single element
/// corresponds to a unique snapshot in the snapshots acquired for the offline stage
///
/// @param[in]  acquiredSnapshotsTimes  The time instants at which snapshots are acquired
/// @param[in]  parameters   The parameters for which the simulations are run ordered with respect to
/// vectors acquiredSnapshotsTimes which contain the corresponding time instants at which we acquired
/// snapshots
///
/// @return     A matrix which has in each row the parameters value and the time instants which
/// correspond to a unique snapshot in the offline folder.
///
Eigen::MatrixXd parTimeCombMat(List<Eigen::VectorXd>
                               acquiredSnapshotsTimes, Eigen::MatrixXd parameters);

}

#endif
