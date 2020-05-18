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
File
    ITHACAgeometry
Description
    set of methods for geometrical parametrization
SourceFiles
    ITHACAgeometry.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the geometry namespace. It contains functions to
/// \dir
/// Directory containing the header and source files for the geometry namespace.

#ifndef ITHACAgeometry_H
#define ITHACAgeometry_H

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
///
/// @brief      Gives the indices conteined into a defined box
///
/// @param[in]  mesh  The mesh
///
/// @param[in]  indices  The list where indices have to be stored
///
/// @param[in]  Box  The box defined by two opposite points
///
/// @param[in]  points2Move  The list containing the coordines of the points into the box
///
/// @return     the labels list of the points
///
labelList getIndicesFromBox(fvMesh& mesh, List<label> indices,
                            Eigen::MatrixXd Box, List<vector>& points2Move);

//--------------------------------------------------------------------------
/// @brief      Gets the indices of the cells around a certain cell.
///
/// @param      mesh    The mesh
/// @param[in]  index   The index of the considered cell
/// @param[in]  layers  The number of layers to be considered
///
/// @return     The indices.
///
List<int> getIndices(fvMesh& mesh, int index, int layers);

//--------------------------------------------------------------------------
/// @brief      Gets the indices of the cells around a certain cell.
///
/// @param      mesh    The mesh
/// @param[in]  index_row   The index of the considered cell
/// @param[in]  index_col   The index of the considered cell
/// @param[in]  layers  The number of layers to be considered
///
/// @return     The indices.
///
List<int> getIndices(fvMesh& mesh, int index_row, int index_col,
                     int layers);

//--------------------------------------------------------------------------
/// @brief      Get the point coordinates and indices from patch.
///
/// @param[in]  mesh     The mesh
/// @param[in]  ind      The index of the patch from which you want to extract the points
/// @param[out] points   The coordinates of the points on the patch
/// @param[out] indices  The indices of the points on the patch
///
void getPointsFromPatch(fvMesh& mesh, label ind, List<vector>& points,
                        labelList& indices);

//--------------------------------------------------------------------------
/// @brief      Get position of displaced segment of points given the displacements of the end points
///
/// @param[in]  x0    The points of the segment before displacement
/// @param[in]  mux1  The x displacement of the 1st end of the segment
/// @param[in]  mux2  The x displacement of the 2nd end of the segment
/// @param[in]  muy1  The y displacement of the 1st end of the segment
/// @param[in]  muy2  The y displacement of the 2nd end of the segment
///
/// @return     The points coordinates of the displaced segment
///
List<vector> displacedSegment(List<vector> x0, double mux1, double mux2,
                              double muy1, double muy2);
//----------------------------------------------------------------------
/// @brief      Displace a Point belonging to a given segment
///
/// @param[in]  x0       The initial coordinate of the point
/// @param[in]  x_low    The position of the lower end of the segment
/// @param[in]  x_up     The posizion of the upper end of the segment
/// @param[in]  mux_low  The x displacement of the lower end of the segment
/// @param[in]  mux_up   The x displacement of the upper end of the segment
/// @param[in]  muy_low  The y displacement of the lower end of the segment
/// @param[in]  muy_up   The y displacement of the upper end of the segment
/// @param[in]  muz_low  The z displacement of the lower end of the segment
/// @param[in]  muz_up   The z displacement of the upper end of the segment
///
/// @return     The displaced point
///
vector displacePoint(vector x0, vector x_low, vector x_up,
                     double mux_low, double mux_up,
                     double muy_low, double muy_up, double muz_low = 0.0, double muz_up = 0.0);

//--------------------------------------------------------------------------
/// @brief      A function that rotates the mesh rigidly, by a fixed specified angle for the points
/// lying within a disc of radius r1, and by a an angle which diminishes from the former angle to zero
/// either linearly or following a Sinusoidal relation or Sigmoid function in a disc of radius
/// r2 larger than r1.
///
/// @param      mesh                  The mesh which you want to rotate rigidly
/// @param[in]  r1                    The inner disc radius
/// @param[in]  r2                    The outer disc radius
/// @param[in]  axis                  The axis of rotation
/// @param[in]  alpha                 The angle of rotation in degrees
/// @param[in]  movingPointsIDs       The list of id's of the points to be rotated
/// @param[in]  radii                 The distances of the points to be rotated to the origin of the rotation
/// @param[in]  angleVariationMethod  The angle variation method either Linear or Sinusoidal or Sigmoid functions
/// @param[in]  v                     The tuning parameter in case of using the Sigmoid function
///
/// @return     The vector field of the newly rotated points
///
Field<vector> rotateMesh(fvMesh& mesh, double r1, double r2, vector axis,
                         double alpha, labelList movingPointsIDs,
                         List<double> radii, word angleVariationMethod = "Linear", double v = 3);

//--------------------------------------------------------------------------
/// @brief      Functions that return a Rotation Matrix given an axis of rotation and an angle in degrees
///
/// @param[in]  AxisOfRotation   The axis of rotation
/// @param[in]  AngleOfRotation  The angle of rotation in degrees
///
/// @return     The rotation matrix in Eigen format
///
Eigen::MatrixXd rotationMatrix(vector AxisOfRotation,
                               double AngleOfRotation);
//--------------------------------------------------------------------------
/// @brief      Compute the distance between the boundary face center and the boundary cell center
///
/// @param[in]  BC_ind   Index of the boundary patch
/// @param[in]  mesh     Mesh
///
/// @return     Vector of the face to cell center distance in Eigen format
///
Eigen::VectorXd boudaryFaceToCellDistance(fvMesh& mesh, label BC_ind);

//--------------------------------------------------------------------------
/// @brief      Gets the indices of the points which are inside a disc of a specified radius before
/// carrying out a rotation around an axis.
///
/// @param      mesh    The mesh for which you want to find the indices of the points inside a disc
/// @param[in]  radius  The radius of the disc
/// @param[in]  origin  The origin of the disc
/// @param[in]  axis    The axis of the rotation
/// @param      radii   The distances of the points lying inside the disc to the origin
///
/// @return     The indices of the points lying inside the specified the disc.
///
List<int> getIndicesFromDisc(fvMesh& mesh, double radius,
                             vector origin, vector axis, List<double>& radii);

//--------------------------------------------------------------------------
/// @brief      Gets the indices from coordinates of a given box.
///
/// @param      field   A geometric Field from OpenFoam field
/// @param[in]  Box     The box
///
/// @tparam     type_f  type of the Field can be Scalar or Vector
///
/// @return     The indices from box.
///
template<typename type_f>
List<int> getIndicesFromBox(
    GeometricField<type_f, fvPatchField, volMesh>& field, Eigen::MatrixXd Box);

//--------------------------------------------------------------------------
/// @brief      Gets a subMesh from a box of coordinates and a given field (used only for the mesh).
///
/// @param      field   The field passed
/// @param[in]  Box     The box of coordinates 3*2 matrix
///
/// @tparam     type_f  type of the Field can be Scalar or Vector
///
/// @return     The subMesh from the box. It can later used to interpolate any Field
///
template<typename type_f>
fvMeshSubset* getSubMeshFromBox(
    GeometricField<type_f, fvPatchField, volMesh>& field, Eigen::MatrixXd Box);

}

#endif
