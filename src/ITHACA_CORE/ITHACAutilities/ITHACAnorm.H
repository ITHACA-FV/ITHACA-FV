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

#ifndef ITHACAnorm_H
#define ITHACAnorm_H

/// \file
/// Header file of the ITHACAnorm file.

#include "ITHACAparameters.H"
#include "fvCFD.H"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#include <Eigen/Eigen>
#pragma GCC diagnostic pop
#include "ITHACAstream.H"

namespace ITHACAutilities
{

//------------------------------------------------------------------------------
/// @brief      Perform the dot product of v and w for a given hilbertSpace

/// @param      v a volVectorField
/// @param      w a volVectorField
/// @param      hilbertSpacePOD hibert space use to compute the dot product, can be L2, L2wBC, H1, wH1
/// @param      weightBC in the case of L2wBC, set the weight of the dot product at the boundary
/// @param      patchBC  in the case of L2wBC, select which boundary patch to apply the dot product on pathBC
/// @param      weightH1 coefficient in front of the dot product of the gradient 
double dot_product_POD(const volVectorField& v, const volVectorField& w,
                        const word& hilbertSpacePOD,
                        const double& weightBC=0., const word& patchBC="inlet",
                        const double& weightH1 = 1.0);

//------------------------------------------------------------------------------
/// @brief      Perform the dot product of v and w for a given hilbertSpace

/// @param      v a volScalarField
/// @param      w a volScalarField
/// @param      hilbertSpacePOD hibert space use to compute the dot product, can be L2, L2wBC, H1, wH1
/// @param      weightBC in the case of L2wBC, set the weight of the dot product at the boundary
/// @param      patchBC  in the case of L2wBC, select which boundary patch to apply the dot product on pathBC
/// @param      weightH1 coefficient in front of the dot product of the gradient 
double dot_product_POD(const volScalarField& v, const volScalarField& w,
                        const word& hilbertSpacePOD,
                        const double& weightBC=0., const word& patchBC="inlet",
                        const double& weightH1 = 1.0);


//------------------------------------------------------------------------------
/// @brief      Perform the dot product of v and w for a given hilbertSpace

/// @param      v a volTensorField
/// @param      w a volTensorField
/// @param      hilbertSpacePOD hibert space use to compute the dot product, can be L2, L2wBC, H1, wH1
/// @param      weightBC in the case of L2wBC, set the weight of the dot product at the boundary
/// @param      patchBC  in the case of L2wBC, select which boundary patch to apply the dot product on pathBC
/// @param      weightH1 coefficient in front of the dot product of the gradient 
double dot_product_POD(const volTensorField& v, const volTensorField& w,
                        const word& hilbertSpacePOD,
                        const double& weightBC=0., const word& patchBC="inlet",
                        const double& weightH1 = 1.0);


//------------------------------------------------------------------------------
/// @brief      Perform the dot product of two PtrList of type T

/// @tparam     T type of v and w, vol[Scalar|Vector|Tensor]Field
/// @param      v field of type T
/// @param      w field of type T
/// @param      hilbertSpacePOD hibert space use to compute the dot product, can be L2, L2wBC, H1, wH1
/// @param      weightBC in the case of L2wBC, set the weight of the dot product at the boundary
/// @param      patchBC  in the case of L2wBC, select which boundary patch to apply the dot product on pathBC
/// @param      weightH1 coefficient in front of the dot product of the gradient 
template<typename T>
Eigen::MatrixXd dot_product_POD(PtrList<T>& v, PtrList<T>& w,
                                const word& hilbertSpacePOD,
                                const double& weightBC,
                                const word& patchBC,
                                const double& weightH1 = 1.0);

//------------------------------------------------------------------------------
/// @brief      Perform the L2 dot product of v and w

/// @param      v a volVectorField
/// @param      w a volVectorField  
/// @return     Result of the dot product
double dot_product_L2(const volVectorField& v, const volVectorField& w);


//------------------------------------------------------------------------------
/// @brief      Perform the L2 dot product of v and w

/// @param      v a volScalarField
/// @param      w a volScalarField  
/// @return     Result of the dot product
double dot_product_L2(const volScalarField& v, const volScalarField& w);

//------------------------------------------------------------------------------
/// @brief      Perform the L2 dot product of v and w

/// @param      v a volScalarField
/// @param      w a volScalarField  
/// @return     Result of the dot product
double dot_product_L2(const volTensorField& v, const volTensorField& w);

//------------------------------------------------------------------------------
/// @brief      Perform the H1 dot product of v and w

/// @tparam     T type of v and w, can be vol[Scalar|Vector]Field
/// @param      v field of type T
/// @param      w field of type T  
/// @return     Result of the dot product
template<typename T>
double dot_product_H1(const T& v, const T& w, const double& weightH1 = 1.0);

//------------------------------------------------------------------------------
/// @brief      Perform the L2 dot product of v and w at a given boundary patch

/// @tparam     T type of v & w, can be vol[Scalar|Vector|Tensor]Field
/// @param      f1BC_i field of type T
/// @param      f2BC_i field of type T
/// @param      AreaFace set the weight of the dot product at the boundary
/// @param      patchBC  select which boundary patch to apply the dot product on pathBC
/// @param      d dimension of the field f1BC_i

/// @return     result of the dot product
double dot_product_patch(const Eigen::VectorXd& f1BC_i,const Eigen::VectorXd& f2BC_i, const scalarField& AreaFace, const int& d);
template<typename T>

//------------------------------------------------------------------------------
/// @brief      Perform the L2 dot product of v and w at a given boundary 

/// @tparam     T type of v & w, can be vol[Scalar|Vector|Tensor]Field
/// @param      v field of type T
/// @param      w field of type T
/// @param      patchBC  select which boundary patch to apply the dot product

/// @return     result of the dot product
double dot_product_boundary(const T& v,const T& w, const word& patchBC);

//------------------------------------------------------------------------------
/// @brief      Perform the L2 dot product of v and w at the boundaries

/// @tparam     T type of v & w, can be vol[Scalar|Vector|Tensor]Field
/// @param      v field of type T
/// @param      w field of type T
/// @param      weightBC set the weight of the dot product at the boundary
/// @param      patchBC  select which boundary patch to apply the dot product on pathBC

/// @return     result of the dot product
template<typename T>
double dot_product_L2wBC(const T& v, const T& w, const double& weightBC, const word& patchBC="inlet");

//------------------------------------------------------------------------------
/// @brief      Compute the norm using the L2 dot product of v at the boundaries

/// @tparam     T type of v & w, can be vol[Scalar|Vector|Tensor]Field
/// @param      v field of type T
/// @param      w field of type T
/// @param      weightBC set the weight of the dot product at the boundary
/// @param      patchBC  select which boundary patch to apply the dot product on pathBC

/// @return     result of the dot product
template<typename T>
double norm_L2wBC(const T& v, const double& weightBC, const word& patchBC="inlet");

//------------------------------------------------------------------------------
/// @brief      Compute the norm of v

/// @tparam     T type of v, vol[Scalar|Vector]Field
/// @param      v field of type T
/// @param      hilbertSpacePOD hibert space use to compute the dot product, can be L2, L2wBC, H1, wH1
/// @param      weightBC in the case of L2wBC, set the weight of the dot product at the boundary
/// @param      patchBC  in the case of L2wBC, select which boundary patch to apply the dot product on pathBC
/// @param      weightH1 coefficient in front of the dot product of the gradient 

/// @return     Norm of v
template<typename T>
double norm_POD(const T& v,
                        const word& hilbertSpacePOD,
                        const double& weightBC=0., const word& patchBC="inlet",
                        const double& weightH1 = 1.0);

//------------------------------------------------------------------------------
/// @brief      Compute the L2 norm of v

/// @tparam     T type of v, vol[Scalar|Vector|Tensor]Field
/// @param      v field of type T

/// @return     Norm of v
template<typename T>
double L2Norm(const T v);

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
