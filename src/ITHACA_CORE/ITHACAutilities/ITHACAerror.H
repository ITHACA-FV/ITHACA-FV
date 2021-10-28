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
\*---------------------------------------------------------------------------*/

#ifndef ITHACAerror_H
#define ITHACAerror_H

/// \file
/// Header file of the ITHACAerror file.


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

namespace ITHACAutilities
{
//--------------------------------------------------------------------------
/// Computes the relative error between two Fields in the Frobenius norm
///
/// @param      field1     The reference field to which the norm is computed
/// @param      field2     The field for which the error is computed
///
/// @tparam     TypeField  type of field can be scalar or vector
///
/// @return     The Frobenius norm of the relative error.
///

template<class Type, template<class> class PatchField, class GeoMesh>
double errorFrobRel(GeometricField<Type, PatchField, GeoMesh>& field1,
                    GeometricField<Type, PatchField, GeoMesh>& field2);

//--------------------------------------------------------------------------
/// Computes the relative error between two geometric Fields in L2 norm
///
/// @param[in]  field1  The field 1
/// @param[in]  field2  The field 2
///
/// @tparam     TypeField   type of field
///
/// @return     L2 norm of the relative error.
///
template<typename T>
double errorL2Rel(GeometricField<T, fvPatchField, volMesh>& field1,
                  GeometricField<T, fvPatchField, volMesh>& field2);

//--------------------------------------------------------------------------
/// Computes the relative error between two geometric Fields in Linf norm
///
/// @param[in]  field1  The field 1
/// @param[in]  field2  The field 2
///
/// @tparam     TypeField   type of field
///
/// @return     Linf norm of the relative error.
///
template<typename T>
double errorLinfRel(GeometricField<T, fvPatchField, volMesh>& field1,
                    GeometricField<T, fvPatchField, volMesh>& field2);

//--------------------------------------------------------------------------
/// Computes the relative error between two geometric Fields in
/// L2 norm, given the mesh volumes field
///
/// @param[in]  field1  The field 1
/// @param[in]  field2  The field 2
/// @param[in]  Volumes  The field containing the volumes of the mesh cells
///
/// @tparam     T   type of field
///
/// @return     L2 norm of the absolute error.
///
template<typename T>
double errorL2Abs(GeometricField<T, fvPatchField, volMesh>& field1,
                  GeometricField<T, fvPatchField, volMesh>& field2, volScalarField& Volumes);

//--------------------------------------------------------------------------
/// Computes the absolute error between two Fields in L2 norm
///
/// @param[in]  field1  The field 1
/// @param[in]  field2  The field 2
///
/// @tparam     TypeField   type of field
///
/// @return     L2 norm of the absolute error.
///
template<typename T>
double errorL2Abs(GeometricField<T, fvPatchField, volMesh>& field1,
                  GeometricField<T, fvPatchField, volMesh>& field2);

//--------------------------------------------------------------------------
/// Computes the relative error in L2 norm between two lists of fields
///
/// @param      fields1  The reference fields to which the norm is computed
/// @param      fields2  The fields for which the error is computed
///
/// @tparam     TypeField   type of field
///
/// @return     Column vector, in each row the L2 norm of the relative error.
///
template<typename T>
Eigen::MatrixXd errorL2Rel(PtrList<GeometricField<T, fvPatchField, volMesh>>&
                           fields1,
                           PtrList<GeometricField<T, fvPatchField, volMesh>>& fields2);

//--------------------------------------------------------------------------
/// Computes the relative error in the Frobenius norm between two lists of fields
///
/// @param      fields1  The reference fields to which the norm is computed
/// @param      fields2  The fields for which the error is computed
///
/// @tparam     TypeField   type of field
///
/// @return     Column vector, in each row the Frobenius norm of the relative error.
///
template<typename T>
Eigen::MatrixXd errorFrobRel(PtrList<GeometricField<T, fvPatchField, volMesh>>&
                             fields1,
                             PtrList<GeometricField<T, fvPatchField, volMesh>>& fields2);

//--------------------------------------------------------------------------
/// Computes the relative error in L2 norm between two lists of Fields
/// given the list of mesh volumes fields
///
/// @param[in]  fields1  The fields 1, scalar or vector
/// @param[in]  fields2  The fields 2, scalar or vector
/// @param[in]  Volumes  The list of the mesh volumes fields
///
/// @tparam     TypeField   type of field
///
/// @return     Column vector, in each row the L2 norm of the error.
///
template<typename T>
Eigen::MatrixXd errorL2Abs(
    PtrList<GeometricField<T, fvPatchField, volMesh>>& fields1,
    PtrList<GeometricField<T, fvPatchField, volMesh>>& fields2,
    PtrList<volScalarField>& Volumes);

//--------------------------------------------------------------------------
/// Computes the absolute error in L2 norm between two lists of Fields
///
/// @param[in]  fields1  The fields 1, scalar or vector
/// @param[in]  fields2  The fields 2, scalar or vector
///
/// @tparam     TypeField   type of field
///
/// @return     olumn vector, in each row the L2 norm of the absolute error.
///
template<typename T>
Eigen::MatrixXd errorL2Abs(PtrList<GeometricField<T, fvPatchField, volMesh>>&
                           fields1,
                           PtrList<GeometricField<T, fvPatchField, volMesh>>& fields2);

//------------------------------------------------------------------------------
/// Evaluate the L2 norm of a geometric field
///
/// @param[in]  field  The field.
///
/// @tparam     T      Type of field, scalar or vector
///
/// @return     L2 norm of the volVectorField.
///
template<class T>
double L2Norm(GeometricField<T, fvPatchField, volMesh>& field);

//------------------------------------------------------------------------------
/// Evaluate the Linf norm of a geometric field
///
/// @param[in]  field  The field.
///
/// @tparam     T      Type of field, scalar or vector
///
/// @return     Linf norm of the volVectorField.
///
template<class T>
double LinfNorm(GeometricField<T, fvPatchField, volMesh>& field);

//------------------------------------------------------------------------------
/// Evaluate the H1 seminorm of a geoemtri
///
/// @param[in]  field  The field.
///
/// @tparam     T      Type of field, scalar or vector
///
/// @return     H1 seminorm of the volVectorField.
///
template<class T>
double H1Seminorm(GeometricField<T, fvPatchField, volMesh>& field);

//------------------------------------------------------------------------------
/// Evaluate the Frobenius norm of a field
///
/// @param[in]  field  The field.
///
/// @tparam     T      Type of field, scalar or vector
///
/// @return     Frobenius norm of the field.
///
template<class Type, template<class> class PatchField, class GeoMesh>
double frobNorm(GeometricField<Type, PatchField, GeoMesh>& field);

//--------------------------------------------------------------------------
/// Evaluate the L2 norm of a field on a boundary patch
///
/// @param[in]  mesh   The mesh.
/// @param[in]  field  The field.
/// @param[in]  patch  Name of the patch.
///
/// @return     L2 norm of the field on the patch.
///
double L2normOnPatch(fvMesh& mesh, volScalarField& field, word patch);

//--------------------------------------------------------------------------
/// Evaluate the L2 inner product between two scalarLists
///
/// @param[in]  mesh   The mesh.
/// @param[in]  field1 The first field.
/// @param[in]  field2 The second field.
/// @param[in]  patch  Name of the patch.
///
/// @return     L2 norm of the field on the patch.
///
double L2productOnPatch(fvMesh& mesh, List<scalar>& field1,
                        List<scalar>& field2, word patch);

//--------------------------------------------------------------------------
/// Evaluate the Linf norm of a field on a boundary patch
///
/// @param[in]  mesh   The mesh.
/// @param[in]  field  The field.
/// @param[in]  patch  Name of the patch.
///
/// @return     Linf norm of the field on the patch.
///
double LinfNormOnPatch(fvMesh& mesh, volScalarField& field, word patch);

//--------------------------------------------------------------------------
/// Evaluate the integral on a patch
///
/// @param[in]  mesh   The mesh.
/// @param[in]  field  The field.
/// @param[in]  patch  Name of the patch.
///
/// @return     Integral of the field on the patch.
///
double integralOnPatch(fvMesh& mesh, volScalarField& field, word patch);

//--------------------------------------------------------------------------
/// Evaluate the integral on a patch
///
/// @param[in]  mesh   The mesh.
/// @param[in]  field  A list of values of the field on the patch.
/// @param[in]  patch  Name of the patch.
///
/// @return     Integral of the field on the patch.
///
double integralOnPatch(fvMesh& mesh, List<scalar> field, word patch);

}


#endif
