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
    ITHACAutilities
Description
    Utilities to manage files and check existence of folders
SourceFiles
    ITHACAutilities.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the ITHACAutilities namespace.
/// \dir
/// Directory containing the header and source files for the ITHACAutilities class.

#ifndef ITHACAutilities_H
#define ITHACAutilities_H

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
#include "ITHACAgeometry.H"
#include "ITHACAsystem.H"
#include "ITHACAerror.H"
#include "ITHACAassign.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


/*---------------------------------------------------------------------------*\
                    Namespace ITHACAutilities Declaration
\*---------------------------------------------------------------------------*/

/// Namespace to implement some useful functions inside ITHACA-FV
namespace ITHACAutilities
{

//--------------------------------------------------------------------------
/// @brief      Generates random matrix with random values in an interval
///
/// @param[in]  rows  Number of rows
/// @param[in]  cols  Number of cols
/// @param[in]  min   Minimum of the interval
/// @param[in]  max   Maximum of the interval
///
/// @return     the matrix
///
Eigen::MatrixXd rand(int rows, int cols, double min, double max);

//--------------------------------------------------------------------------
/// @brief      Generates random matrix with random values in an interval
///
/// @param[in]  rows    Number of rows
/// @param[in]  minMax  The min and max Matrix, it has as many rows as columns of the output matrix.
///
/// @return     Random matrix
///
Eigen::MatrixXd rand(int rows, Eigen::MatrixXd minMax);


//--------------------------------------------------------------------------
/// This function checks if ratio is an integer
///
/// @param[in]  ratio  The ratio
///
/// @return     1 if ratio is an integer
///
bool isInteger(double ratio);

//--------------------------------------------------------------------------
/// This function checks if the case is turbulent
///
/// @return     1 if the case is turbulent
///
bool isTurbulent();


//--------------------------------------------------------------------------
/// Exact reconstruction using a certain number of modes for a list of
/// fields and the projection coefficients
///
/// @param[in]  modes         The modes used for reconstruction as PtrList.
/// @param[in]  Nmodes        The number of modes you want to use.
/// @param[in]  coeff_matrix  The matrix of coefficients.
///
/// @tparam     TypeField   type of field
///
/// @return     rec_field     The reconstructed field as PtrList.
///
template<class TypeField>
PtrList<TypeField> reconstruct_from_coeff(PtrList<TypeField>& modes,
        Eigen::MatrixXd& coeff_matrix, label Nmodes);


//--------------------------------------------------------------------------
/// Computes a Mass Matrix from a list of Basis Functions (vectorial) for L2 projection
///
/// @param[in]  modes  The modes in PtrList form.
///
/// @return     The mass matrix in Eigen::MatrixXd form.
///
Eigen::MatrixXd get_mass_matrix(PtrList<volVectorField> modes,
                                int Nmodes = 0);

//--------------------------------------------------------------------------
/// Computes a Mass Matrix from a list of Basis Functions (scalar) for L2 projection
///
/// @param[in]  modes  The modes in PtrList form.
///
/// @return     The mass matrix in Eigen::MatrixXd form.
///
Eigen::MatrixXd get_mass_matrix(PtrList<volScalarField> modes,
                                int Nmodes = 0);

//-----------------------------------------------------------------------------
/// Gets the mass matrix using the eigen routine
///
/// @param[in]  fields            The openFoam fields
/// @param[in]  consider_volumes  True if you want to consider the volumes
///
/// @tparam     T                 type of field
///
/// @return     The mass matrix eigen.
///
template<typename T>
Eigen::MatrixXd get_mass_matrix_Eigen(PtrList<T>& fields,
                                      bool consider_volumes = true);

//--------------------------------------------------------------------------
/// Gets a vector containing the volumes of each cell of the mesh.
///
/// @param[in]  snapshot   A snapshot, scalar or vector field.
///
/// @tparam     TypeField   type of field
///
/// @return     The FV mass matrix in Eigen::MatrixXd form.
///
template<class TypeField>
Eigen::VectorXd get_mass_matrix_FV(
    GeometricField<TypeField, fvPatchField, volMesh>& snapshot);

//--------------------------------------------------------------------------
/// Projects a snapshot vector field on a non-orthogonal basis function and
/// gets the coefficients of the projection
///
/// @param[in]  snapshot  The snapshots.
/// @param[in]  modes     The modes.
///
/// @return     The coefficients of the projection.
///
Eigen::VectorXd get_coeffs(volVectorField snapshot,
                           PtrList<volVectorField>& modes, int Nmodes = 0);

//--------------------------------------------------------------------------
/// Project a snapshot vector field on a non-orthogonal basis function
/// and get the coefficients of the projection
///
/// @param[in]  snapshot  The snapshots.
/// @param[in]  modes     The modes.
/// @param[in]  Nmodes    Number of modes.
///
/// @return     The coefficients of the projection.
///
Eigen::VectorXd get_coeffs(volScalarField snapshot,
                           PtrList<volScalarField>& modes, int Nmodes = 0);

//--------------------------------------------------------------------------
/// Project a snapshot scalar field on a non-orthogonal basis function and get the coefficients of the projection
///
/// @param[in]  snapshots  The snapshots matrix.
/// @param[in]  modes     The modes.
/// @param[in]  Nmodes    Number of modes.
///
/// @return     The coefficients of the projection.
///
template<class TypeField>
Eigen::MatrixXd get_coeffs(PtrList<TypeField> snapshots,
                           PtrList<TypeField> modes,  int Nmodes = 0);

//--------------------------------------------------------------------------
/// Project a snapshot scalar field on an orthogonal basis function
/// and get the coefficients of the projection
///
/// @param[in]  snapshot  The snapshot.
/// @param[in]  modes     The modes.
/// @param[in]  Nmodes    Number of modes, 0 for all.
///
/// @tparam     TypeField   type of field
///
/// @return     The coefficients of the projection.
///
template<class TypeField>
Eigen::VectorXd get_coeffs_ortho(
    GeometricField<TypeField, fvPatchField, volMesh> snapshot,
    PtrList<GeometricField<TypeField, fvPatchField, volMesh>>& modes,
    int Nmodes = 0);

//--------------------------------------------------------------------------
/// @brief      Gets the coeffs ortho.
///
/// @param[in]  snapshots  The snapshots
/// @param      modes      The modes
/// @param[in]  Nmodes     Number of modes, 0 for all.
///
/// @tparam     TypeField   type of field
///
/// @return     The coeffs ortho.
///
template<class TypeField>
Eigen::MatrixXd get_coeffs_ortho(PtrList<TypeField> snapshots,
                                 PtrList<TypeField>& modes, int Nmodes = 0);

//--------------------------------------------------------------------------
/// @brief      Gets the coeffs of Frobenius projection.
///
/// @param[in]  snapshots  The snapshots
/// @param      modes      The modes
/// @param[in]  Nmodes     Number of modes, 0 for all.
///
/// @tparam     TypeField   type of field
///
/// @return     The coeffs of Frobenius projection.
///
template<class TypeField>
Eigen::MatrixXd getCoeffsFrobenius(PtrList<TypeField> snapshots,
                                   PtrList<TypeField>& modes, int nModes = 0);

//--------------------------------------------------------------------------
/// @brief      A method to comput the time-parameter combined matrix whose any single element
/// corresponds to a unique snapshot in the snapshots acquired for the offline stage
///
/// @param[in]  acquiredSnapshotsTimes  The time instants at which snapshots are acquired
/// @param[in]  parameters   The parameters for which the simulations are run oredered with respect to
/// vectors acquiredSnapshotsTimes which contain the corresponding time instants at which we acquired
/// snapshots
///
/// @return     A matrix which has in each row the parameters value and the time instants which
/// correspond to a unique snapshot in the offline folder.
///
Eigen::MatrixXd parTimeCombMat(List<Eigen::VectorXd>
                               acquiredSnapshotsTimes, Eigen::MatrixXd parameters);


//--------------------------------------------------------------------------
/// @brief      Returns a scalarField that containes the non-orthogonality value of a given mesh
///
/// @param      mesh      The mesh for which you want to compute the non-orthogonality value
/// @param      NonOrtho  The scalarField containing the non-orthogonality
///
/// @return     The scalarField containing the non-orthogonality
///
volScalarField meshNonOrtho(fvMesh& mesh, volScalarField& NonOrtho);

//--------------------------------------------------------------------------
/// @brief      Rotate a list of points in clockwise direction given an axis of rotation and an angle in degrees.
///
/// @param[in]  originalPoints    The original points
/// @param[in]  CentreOfRotation  The Axis Of Rotation
/// @param[in]  AngleOfRotation   The angle of rotation in degrees
///
/// @return     The Rotated Points
///
List<vector> rotatePoints(const List<vector>& originalPoints,
                          vector AxisOfRotation, double AngleOfRotation);



};
#endif







