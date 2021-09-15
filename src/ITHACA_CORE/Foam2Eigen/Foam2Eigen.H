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
     foam2eigen
 Description
     Converts foam objects into Eigen and viceversa
 SourceFiles
     foam2eigen.C
 \*---------------------------------------------------------------------------*/

/// \file
/// Header file of the Foam2Eigen class.
/// \dir
/// Directory containing the header and source files for the Foam2Eigen class.

#ifndef Foam2Eigen_H
#define Foam2Eigen_H

#include "fvCFD.H"
#include "IOmanip.H"
#include "ITHACAassert.H"
#include <tuple>
#include <sys/stat.h>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#include <Eigen/Eigen>
#pragma GCC diagnostic pop
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


/*---------------------------------------------------------------------------*\
  Class foam2eigen Declaration
\*---------------------------------------------------------------------------*/

/// Class to converts OpenFOAM objects into Eigen and viceversa
class Foam2Eigen
{
    private:

    public:

        //----------------------------------------------------------------------
        /// @brief      Convert a FvMatrix OpenFOAM matrix (Linear System) into
        ///             a Eigen Matrix A and a source vector b
        ///
        /// @param[in]  foam_matrix       The foam matrix can be fvScalarMatrix
        ///                               or fvVectorMatrix
        /// @param[out] A                 The matrix or the list of matrices,
        ///                               can be dense or sparse
        /// @param[out] b                 The source term vector or the list of
        ///                               source term vectors, always dense
        ///
        /// @tparam     type_foam_matrix  The type of foam matrix can be scalar
        ///                               or vector
        /// @tparam     type_A            The type of matrix or the list of
        ///                               matrices, can be dense
        ///                               (Eigen::MatrixXd) or sparse
        ///                               (Eigen::SparseMatrix<double>) or a
        ///                               List of them
        /// @tparam     type_B            The type source term vector or the
        ///                               list can be Eigen::VectorXd or
        ///                               List<Eigen::VectorXd>
        ///
        template <class type_foam_matrix, class type_A, class type_B>
        static void fvMatrix2Eigen(fvMatrix<type_foam_matrix> foam_matrix, type_A& A,
                                   type_B& b);

        //----------------------------------------------------------------------
        /// @brief      Convert a ldu OpenFOAM matrix into a Eigen Matrix A
        ///
        /// @param[in]  foam_matrix       The foam matrix can be fvScalarMatrix
        ///                               or fvVectorMatrix
        /// @param[out] A                 The matrix or the list of matrices,
        ///                               can be dense or sparse
        ///
        /// @tparam     type_foam_matrix  The type of foam matrix can be
        ///                               fvScalarMatrix or fvVectorMatrix
        /// @tparam     type_A            The type of matrix or the list of
        ///                               matrices, can be dense
        ///                               (Eigen::MatrixXd) or sparse
        ///                               (Eigen::SparseMatrix<double>) or a
        ///                               List of them
        ///
        template <class type_foam_matrix, class type_A>
        static void fvMatrix2EigenM(fvMatrix<type_foam_matrix>& foam_matrix, type_A& A);

        //----------------------------------------------------------------------
        /// @brief      Convert a ldu OpenFOAM matrix into a source vector b
        ///
        /// @param[in]  foam_matrix       The foam matrix can be fvScalarMatrix
        ///                               or fvVectorMatrix
        /// @param[out] b                 The source term vector or the list of
        ///                               source term vectors, always dense
        ///
        /// @tparam     type_foam_matrix  The type of foam matrix can be
        ///                               fvScalarMatrix or fvVectorMatrix
        /// @tparam     type_B            The type source term vector or the
        ///                               list can be Eigen::VectorXd or
        ///                               List<Eigen::VectorXd>
        ///
        template <class type_foam_matrix, class type_B>
        static void fvMatrix2EigenV(fvMatrix<type_foam_matrix>& foam_matrix, type_B& b);


        //----------------------------------------------------------------------
        /// @brief      Convert a PtrList of snapshots to Eigen matrix (only
        ///             internal field)
        ///
        /// @param[in]  fields      The fields can be a PtrList<volScalarField>
        ///                         or PtrList<volVectorField>
        /// @param[in]  Nfields     The nfields
        ///
        /// @tparam     Type        vector or scalar.
        /// @tparam     PatchField  fvPatchField or fvsPatchField.
        /// @tparam     GeoMesh     volMesh or surfaceMesh.
        ///
        /// @return     An Eigen matrix containing as columns the snapshots
        ///
        template<class Type, template<class> class PatchField, class GeoMesh>
        static Eigen::MatrixXd PtrList2Eigen(
            PtrList<GeometricField<Type, PatchField, GeoMesh>>& fields,
            label Nfields = -1);

        //----------------------------------------------------------------------
        /// @brief      Convert a vector OpenFOAM field into an Eigen Vector
        ///
        /// @param[in]  field       The field
        ///
        /// @tparam     PatchField  fvPatchField or fvsPatchField.
        /// @tparam     GeoMesh     volMesh or surfaceMesh.
        ///
        /// @return     Dense Eigen vector
        ///
        template<template<class> class PatchField, class GeoMesh>
        static Eigen::VectorXd field2Eigen(GeometricField<vector, PatchField, GeoMesh>&
                                           field);

        //----------------------------------------------------------------------
        /// @brief      Convert a scalar OpenFOAM field into an Eigen Vector
        ///
        /// @param[in]  field       The field
        ///
        /// @tparam     PatchField  fvPatchField or fvsPatchField.
        /// @tparam     GeoMesh     volMesh or surfaceMesh.
        ///
        /// @return     Dense Eigen vector
        ///
        template<template<class> class PatchField, class GeoMesh>
        static Eigen::VectorXd field2Eigen(GeometricField<scalar, PatchField, GeoMesh>&
                                           field);

        //----------------------------------------------------------------------
        /// @brief      Convert an OpenFOAM field into an Eigen Vector
        ///
        /// @param[in]  field    The field
        ///
        /// @tparam     Type     vector or scalar.
        /// @tparam     GeoMesh  volMesh or surfaceMesh.
        ///
        /// @return     Dense Eigen vector
        ///
        template<class Type, class GeoMesh>
        static Eigen::VectorXd field2Eigen(const
                                           DimensionedField<Type, GeoMesh>& field);

        //----------------------------------------------------------------------
        /// @brief      Convert an OpenFOAM field into an Eigen Vector
        ///
        /// @param[in]  field  The field
        ///
        /// @tparam     Type   Type of the field can be volScalarField or
        ///                    volVectorField
        ///
        /// @return     Dense   Eigen vector
        ///
        template<class Type>
        static Eigen::VectorXd field2Eigen(const Field<Type>& field);

        //----------------------------------------------------------------------
        /// @brief      Convert an OpenFOAM scalar field to a List of Eigen
        ///             Vectors, one for each boundary
        ///
        /// @param[in]  field       The field
        ///
        /// @tparam     PatchField  fvPatchField or fvsPatchField.
        /// @tparam     GeoMesh     volMesh or surfaceMesh.
        ///
        /// @return     Dense Eigen vector
        ///
        template<template<class> class PatchField, class GeoMesh>
        static List<Eigen::VectorXd> field2EigenBC(
            GeometricField<scalar, PatchField, GeoMesh>& field);

        //----------------------------------------------------------------------
        /// @brief      Convert an OpenFOAM vector field to a List of Eigen
        ///             Vectors, one for each boundary
        ///
        /// @param[in]  field       The field
        ///
        /// @tparam     PatchField  fvPatchField or fvsPatchField.
        /// @tparam     GeoMesh     volMesh or surfaceMesh.
        ///
        /// @return     Dense Eigen vector
        ///
        template<template<class> class PatchField, class GeoMesh>
        static List<Eigen::VectorXd> field2EigenBC(
            GeometricField<vector, PatchField, GeoMesh>& field);


        //----------------------------------------------------------------------
        /// @brief      Convert an OpenFOAM scalar field to a List of Eigen
        ///             Vectors, one for each boundary
        ///
        /// @param[in]  fields      The field
        /// @param[in]  Nfields     The nfields
        ///
        /// @tparam     PatchField  fvPatchField or fvsPatchField.
        /// @tparam     GeoMesh     volMesh or surfaceMesh.
        ///
        /// @return     Dense Eigen matrix
        ///
        template<template<class> class PatchField, class GeoMesh>
        static List<Eigen::MatrixXd> PtrList2EigenBC(
            PtrList<GeometricField<scalar, PatchField, GeoMesh>>& fields,
            label Nfields = -1);

        //----------------------------------------------------------------------
        /// @brief      Convert an OpenFOAM vector field to a List of Eigen
        ///             Vectors, one for each boundary
        ///
        /// @param[in]  fields      The field
        /// @param[in]  Nfields     The nfields
        ///
        /// @tparam     PatchField  fvPatchField or fvsPatchField.
        /// @tparam     GeoMesh     volMesh or surfaceMesh.
        ///
        /// @return     Dense Eigen matrix
        ///
        template<template<class> class PatchField, class GeoMesh>
        static List<Eigen::MatrixXd> PtrList2EigenBC(
            PtrList<GeometricField<vector, PatchField, GeoMesh>>& fields,
            label Nfields = -1);

        //----------------------------------------------------------------------
        /// @brief         Convert a vector in Eigen format into an OpenFOAM
        ///                scalar GeometricField
        ///
        /// @param[in/out] field         OpenFOAM GeometricField
        /// @param[in]     eigen_vector  Vector in Eigen format
        ///
        /// @tparam        PatchField    fvPatchField or fvsPatchField.
        /// @tparam        GeoMesh       volMesh or surfaceMesh.
        ///
        /// @return        OpenFOAM GeometricField
        ///
        template<template<class> class PatchField, class GeoMesh>
        static GeometricField<scalar, PatchField, GeoMesh> Eigen2field(
            GeometricField<scalar, PatchField, GeoMesh>& field,
            Eigen::VectorXd& eigen_vector, bool correctBC=true);

        //----------------------------------------------------------------------
        /// @brief         Convert a vector in Eigen format into an OpenFOAM
        ///                vector GeometricField
        ///
        /// @param[in/out] field         OpenFOAM GeometricField
        /// @param[in]     eigen_vector  Vector in Eigen format
        ///
        /// @tparam        PatchField    fvPatchField or fvsPatchField.
        /// @tparam        GeoMesh       volMesh or surfaceMesh.
        ///
        /// @return        OpenFOAM GeometricField
        ///
        template<template<class> class PatchField, class GeoMesh>
        static GeometricField<vector, PatchField, GeoMesh> Eigen2field(
            GeometricField<vector, PatchField, GeoMesh>& field,
            Eigen::VectorXd& eigen_vector, bool correctBC=true);

        //----------------------------------------------------------------------
        /// @brief         Converts a matrix in Eigen format into an OpenFOAM
        ///                Field
        ///
        /// @param[in/out] field   OpenFOAM Field
        /// @param[in]     matrix  Matrix in Eigen format
        ///
        /// @tparam        Type    vector or scalar.
        ///
        /// @return        OpenFOAM Field
        ///
        template <class Type>
        static Field<Type> Eigen2field(
            Field<Type>& field, Eigen::MatrixXd& matrix, bool correctBC=true);

        //----------------------------------------------------------------------
        /// @brief      Obtain the Mass Matrix from a list of reduced basis
        ///
        /// @param      matrix      The matrix
        /// @param[in]  modes       The OpenFOAM modes contained in a PtrList
        /// @param[in]  Nmodes      Number of requested modes
        ///
        /// @tparam     Type        vector or scalar.
        /// @tparam     PatchField  fvPatchField or fvsPatchField.
        /// @tparam     GeoMesh     volMesh or surfaceMesh.
        ///
        /// @return     An Eigen matrix containing the mass matrix
        ///
        template<class Type, template<class> class PatchField, class GeoMesh>
        static std::tuple<Eigen::MatrixXd, Eigen::VectorXd> projectFvMatrix(
            fvMatrix<Type>& matrix,
            PtrList<GeometricField<Type, PatchField, GeoMesh>>& modes, label Nmodes);


        //----------------------------------------------------------------------
        /// @brief      Perform the projection of an OpenFOAM field onto a set
        ///             of modes using the Eigen matrix multiplication
        ///
        ///             It is the standard Galerkin projection of a field onto a
        ///             reduced basis:
        ///
        ///             @f[  \mathbf{f} \in R^{Nx1} , \mathbf{f_r} \in R^{nx1} ,
        ///             \mathbf{M} \in R^{Nxn}
        /// @f]
        ///
        ///             @f[ \mathbf{f_r} = \mathbf{M^T f}
        /// @f]
        ///
        /// @param[in]  field       The f OpenFOAM field
        /// @param[in]  modes       The M OpenFOAM modes contained in a PtrList
        /// @param[in]  Nmodes      The n number of modes you want to use (The
        ///                         size of input modes can be bigger than the
        ///                         requested number of modes)
        ///
        /// @tparam     Type        vector or scalar.
        /// @tparam     PatchField  fvPatchField or fvsPatchField.
        /// @tparam     GeoMesh     volMesh or surfaceMesh.
        ///
        /// @return     Eigen vector with the projected field
        ///
        template<class Type, template<class> class PatchField, class GeoMesh>
        static Eigen::VectorXd projectField(
            GeometricField<Type, PatchField, GeoMesh>& field,
            PtrList<GeometricField<Type, PatchField, GeoMesh>>& modes,
            label Nmodes);

        //----------------------------------------------------------------------
        /// @brief      Obtain the Mass Matrix from a list of reduced basis
        ///
        /// @param      modes       The OpenFOAM modes contained in a PtrList
        /// @param[in]  Nmodes      The number of requested modes
        ///
        /// @tparam     Type        vector or scalar.
        /// @tparam     PatchField  fvPatchField or fvsPatchField.
        /// @tparam     GeoMesh     volMesh or surfaceMesh.
        ///
        /// @return     Eigen matrix containing the mass matrix
        ///
        template<class Type, template<class> class PatchField, class GeoMesh>
        static Eigen::MatrixXd MassMatrix(
            PtrList<GeometricField<Type, PatchField, GeoMesh>>& modes, label Nmodes);

        //----------------------------------------------------------------------
        /// @brief      Convert a PtrList of OpenFOAM fvMatrix into a tuple of
        ///             lists of Eigen Sparse Matrices and source vectors
        ///
        /// @param[in]  MatrixList  The matrix list
        ///
        /// @tparam     Type        { description }
        /// @tparam     [in]  Type  type of fvMatrix can be either scalar or vectorial
        ///
        /// @return     a tuple containing a list of sparse matrices and a list
        ///             of source vectors
        ///
        template<class Type>
        static std::tuple<List<Eigen::SparseMatrix<double>>, List<Eigen::VectorXd>>
                LFvMatrix2LSM(PtrList<fvMatrix<Type>>& MatrixList);

        //--------------------------------------------------------------------------
        /// @brief      Convert a Foam List into an Eigen matrix with one column
        ///
        /// @param[in]  list        The Foam List
        ///
        /// @tparam     type_matrix The type of data in the List (int or double)
        ///
        /// @return     The Eigen matrix
        ///
        template <class type_matrix>
        static Eigen::Matrix<type_matrix, Eigen::Dynamic, Eigen::Dynamic>
        List2EigenMatrix ( List<type_matrix> list );

        //--------------------------------------------------------------------------
        /// @brief      Convert an Eigen matrix with one column into a Foam List
        ///
        /// @param[in]  matrix       The Eigen matrix
        ///
        /// @tparam     type_matrix The type of data in the List (int or double)
        ///
        /// @return     The Foam List
        ///
        template <class type_matrix>
        static List<type_matrix> EigenMatrix2List (
            Eigen::Matrix<type_matrix, Eigen::Dynamic, Eigen::Dynamic> matrix );
};
#endif
