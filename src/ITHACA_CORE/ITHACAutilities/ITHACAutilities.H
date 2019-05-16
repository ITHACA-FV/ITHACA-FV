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
/// Header file of the ITHACAutilities class.
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
using namespace std::placeholders;

#include "Foam2Eigen.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


/*---------------------------------------------------------------------------*\
                        Class ITHACAutilities Declaration
\*---------------------------------------------------------------------------*/

/// Class to implement some useful functions inside ITHACA-FV
class ITHACAutilities
{

    private:

    public:

        //--------------------------------------------------------------------------
        /// @brief      Gets the indices of the cells around a certain cell.
        ///
        /// @param      mesh    The mesh
        /// @param[in]  index   The index of the considered cell
        /// @param[in]  layers  The number of layers to be considered
        ///
        /// @return     The indices.
        ///
        static List<int> getIndices(fvMesh& mesh, int index, int layers);

        //--------------------------------------------------------------------------
        /// @brief      Gets the indices of the cells around a certain cell.
        ///
        /// @param      mesh    The mesh
        /// @param[in]  index   The index of the considered cell
        /// @param[in]  layers  The number of layers to be considered
        ///
        /// @return     The indices.
        ///
        static List<int> getIndices(fvMesh& mesh, int index_row, int index_col,
                                    int layers);

        //--------------------------------------------------------------------------
        /// @brief      Creates a symbolic links to 0, system and constant
        ///
        /// @param[in]  folder  The folder where to create the symbolic links
        ///
        static void createSymLink(word folder);

        //--------------------------------------------------------------------------
        /// @brief      generate random matrix with random values in an interval
        ///
        /// @param[in]  rows  The rows
        /// @param[in]  cols  The cols
        /// @param[in]  min   The minimum
        /// @param[in]  max   The maximum
        ///
        /// @return     the matrix
        ///
        static Eigen::MatrixXd rand(int rows, int cols, double min, double max);

        //--------------------------------------------------------------------------
        /// @brief      generate random matrix with random values in an interval
        ///
        /// @param[in]  rows    The rows
        /// @param[in]  minMax  The min and max Matrix, as mant rows as columns of the output matrix.
        ///
        /// @return     the matrix
        ///
        static Eigen::MatrixXd rand(int rows, Eigen::MatrixXd minMax);

        //--------------------------------------------------------------------------
        /// Check if a folder exists
        ///
        /// @param[in]  folder  The folder you want to check.
        ///
        /// @return     1 if the folder exists 0 elsewhere.
        ///
        static bool check_folder(word folder);

        //--------------------------------------------------------------------------
        /// Check if the offline data folder exists
        ///
        /// @return     1 if the folder exists 0 elsewhere.
        ///
        static bool check_off();

        //--------------------------------------------------------------------------
        /// @brief      function that returns true if a file exists
        ///
        /// @param[in]  fileName  The name of the file
        ///
        /// @return     1 if the file exists 0 elsewhere.
        ///
        static bool check_file(std::string fileName);

        //--------------------------------------------------------------------------
        /// Check if the POD data folder exists
        ///
        /// @return     1 if the folder exists 0 elsewhere.
        ///
        static bool check_pod();

        //--------------------------------------------------------------------------
        /// Check if the supremizer folder exists
        ///
        /// @return     1 if the folder exists 0 elsewhere.
        ///
        static bool check_sup();

        //--------------------------------------------------------------------------
        /// Exact reconstruction using a certain number of modes for a list of
        /// fields and the projection coefficients
        ///
        /// @param[in]  modes         The modes used for reconstruction as PtrList.
        /// @param[in]  Nmodes        The number of modes you want to use.
        /// @param[in]  coeff_matrix  The matrix of coefficients.
        /// @return     rec_field     The reconstructed field as PtrList.
        ///
        template<typename T>
        static PtrList<T> reconstruct_from_coeff(PtrList<T>& modes,
                Eigen::MatrixXd& coeff_matrix, label Nmodes);

        //--------------------------------------------------------------------------
        /// Function to compute the relative error between two volVectorFields in L2 norm
        ///
        /// @param[in]  field1  The field 1
        /// @param[in]  field2  The field 2
        ///
        /// @return     L2 norm of the relative error.
        ///
        static double error_fields(volVectorField& field1, volVectorField& field2);

        //--------------------------------------------------------------------------
        /// Function to compute the absolute error between two volVectorFields in L2 norm
        ///
        /// @param[in]  field1  The field 1
        /// @param[in]  field2  The field 2
        ///
        /// @return     L2 norm of the absolute error.
        ///
        static double error_fields_abs(volVectorField& field1, volVectorField& field2);

        //--------------------------------------------------------------------------
        /// Function to compute the relative error between two volScalarFields in L2 norm
        ///
        /// @param[in]  field1  The field 1
        /// @param[in]  field2  The field 2
        ///
        /// @return     L2 norm of the absolute error.
        ///
        static double error_fields(volScalarField& field1, volScalarField& field2);

        //--------------------------------------------------------------------------
        /// Function to compute the absolute error between two volScalarFields in L2 norm
        ///
        /// @param[in]  field1  The field 1
        /// @param[in]  field2  The field 2
        ///
        /// @return     L2 norm of the absolute error.
        ///
        static double error_fields_abs(volScalarField& field1, volScalarField& field2);

        //--------------------------------------------------------------------------
        /// Function to compute the relative error in L2 norm between two lists of volVectorFields
        ///
        /// @param[in]  fields1  The fields 1
        /// @param[in]  fields2  The fields 2
        ///
        /// @return     A vector containing the error between the two fields in Eigen::MatrixXd form.
        ///
        static Eigen::MatrixXd error_listfields(PtrList<volVectorField>& fields1,
                                                PtrList<volVectorField>& fields2);

        //--------------------------------------------------------------------------
        /// Function to compute the relative error in L2 norm between two lists of volScalarFields
        ///
        /// @param[in]  fields1  The fields 1
        /// @param[in]  fields2  The fields 2
        ///
        /// @return     A vector containing the error between the two fields in Eigen::MatrixXd form.
        ///
        static Eigen::MatrixXd error_listfields(PtrList<volScalarField>& fields1,
                                                PtrList<volScalarField>& fields2);

        //--------------------------------------------------------------------------
        /// Function to compute the absolute error in L2 norm between two lists of volScalarFields
        ///
        /// @param[in]  fields1  The fields 1
        /// @param[in]  fields2  The fields 2
        ///
        /// @return     A vector containing the error between the two fields in Eigen::MatrixXd form.
        ///
        static Eigen::MatrixXd error_listfields_abs(PtrList<volScalarField>& fields1,
                PtrList<volScalarField>& fields2);

        //--------------------------------------------------------------------------
        /// Function to compute the absolute error in L2 norm between two lists of volVectorFields
        ///
        /// @param[in]  fields1  The fields 1
        /// @param[in]  fields2  The fields 2
        ///
        /// @return     A vector containing the error between the two fields in Eigen::MatrixXd form.
        ///
        static Eigen::MatrixXd error_listfields_abs(PtrList<volVectorField>& fields1,
                PtrList<volVectorField>& fields2);


        //--------------------------------------------------------------------------
        /// Function to compute a Mass Matrix from a list of Basis Functions (vectorial) for L2 projection
        ///
        /// @param[in]  modes  The modes in PtrList form.
        ///
        /// @return     The mass matrix in Eigen::MatrixXd form.
        ///
        static Eigen::MatrixXd get_mass_matrix(PtrList<volVectorField> modes,
                                               int Nmodes = 0);

        //--------------------------------------------------------------------------
        /// Function to compute a Mass Matrix from a list of Basis Functions (scalar) for L2 projection
        ///
        /// @param[in]  modes  The modes in PtrList form.
        ///
        /// @return     The mass matrix in Eigen::MatrixXd form.
        ///
        static Eigen::MatrixXd get_mass_matrix(PtrList<volScalarField> modes,
                                               int Nmodes = 0);

        //--------------------------------------------------------------------------
        /// Function to compute the FV Mass Matrix from a vector or scalar field: diagonal matrix where non zero elements are the volumes of the cells.
        ///
        /// @param[in]  snapshot  A snapshot in volScalarField or volVectorField form.
        ///
        /// @return     The FV mass matrix in Eigen::MatrixXd form.
        ///
        template<class type_f>
        static Eigen::MatrixXd get_mass_matrix_FV(
            GeometricField<type_f, fvPatchField, volMesh>& snapshot);

        //--------------------------------------------------------------------------
        /// Project a snapshot scalar field on a non-orthogonal basis function and get the coefficients of the projection
        ///
        /// @param[in]  snapshots  The snapshots matrix.
        /// @param[in]  modes     The modes.
        ///
        /// @return     The coefficients of the projection.
        ///
        static Eigen::MatrixXd get_coeffs(PtrList<volScalarField> snapshots,
                                          PtrList<volScalarField> modes,  int Nmodes = 0);

        //--------------------------------------------------------------------------
        /// Project a snapshot vector field on a non-orthogonal basis function and get the coefficients of the projection
        ///
        /// @param[in]  snapshots  The snapshots matrix.
        /// @param[in]  modes     The modes.
        ///
        /// @return     The coefficients of the projection.
        ///
        static Eigen::MatrixXd get_coeffs(PtrList<volVectorField> snapshots,
                                          PtrList<volVectorField> modes, int Nmodes = 0);

        //--------------------------------------------------------------------------
        /// Project a snapshot vector field on a non-orthogonal basis function and get the coefficients of the projection
        ///
        /// @param[in]  snapshot  The snapshots.
        /// @param[in]  modes     The modes.
        ///
        /// @return     The coefficients of the projection.
        ///
        static Eigen::VectorXd get_coeffs(volVectorField snapshot,
                                          PtrList<volVectorField>& modes, int Nmodes = 0);

        //--------------------------------------------------------------------------
        /// Project a snapshot vector field on a non-orthogonal basis function
        /// and get the coefficients of the projection
        ///
        /// @param[in]  snapshot  The snapshots.
        /// @param[in]  modes     The modes.
        ///
        /// @return     The coefficients of the projection.
        ///
        static Eigen::VectorXd get_coeffs(volScalarField snapshot,
                                          PtrList<volScalarField>& modes, int Nmodes = 0);

        //--------------------------------------------------------------------------
        /// @brief      Gets the coeffs ortho.
        ///
        /// @param[in]  snapshots  The snapshots
        /// @param      modes      The modes
        ///
        /// @return     The coeffs ortho.
        ///
        static Eigen::MatrixXd get_coeffs_ortho(PtrList<volScalarField> snapshots,
                                                PtrList<volScalarField>& modes, int Nmodes = 0);

        //--------------------------------------------------------------------------
        /// @brief      Gets the coeffs ortho.
        ///
        /// @param[in]  snapshots  The snapshots
        /// @param      modes      The modes
        ///
        /// @return     The coeffs ortho.
        ///
        static Eigen::MatrixXd get_coeffs_ortho(PtrList<volVectorField> snapshots,
                                                PtrList<volVectorField>& modes, int Nmodes = 0);

        //--------------------------------------------------------------------------
        /// Project a snapshot scalar field on an orthogonal basis function
        /// and get the coefficients of the projection
        ///
        /// @param[in]  snapshot  The snapshot.
        /// @param[in]  modes     The modes.
        ///
        /// @return     The coefficients of the projection.
        ///
        static Eigen::VectorXd get_coeffs_ortho(volScalarField snapshot,
                                                PtrList<volScalarField>& modes, int Nmodes = 0);

        //--------------------------------------------------------------------------
        /// Project a snapshot scalar field on an orthogonal basis function
        /// and get the coefficients of the projection
        ///
        /// @param[in]  snapshot  The snapshot.
        /// @param[in]  modes     The modes.
        ///
        /// @return     The coefficients of the projection.
        ///
        static Eigen::VectorXd get_coeffs_ortho(volVectorField snapshot,
                                                PtrList<volVectorField>& modes, int Nmodes = 0);

        //--------------------------------------------------------------------------
        /// Assign internal field condition
        ///
        /// @param[in,out]  s      field where you want to assign the internal field condition
        /// @param          value  The value you want to assign
        ///
        /// @tparam         T      type of field (volVectorField or volScalarField)
        /// @tparam         G      type of value you want to assign (OpenFOAM vector or scalar)
        ///
        template<typename T, typename G>
        static void assignIF(T& s, G& value);

        //--------------------------------------------------------------------------
        /// Assign one to volScalarField
        ///
        /// @param[in,out] s     field where you want to assign the internal field
        ///                      condition
        /// @param         L     List of indices
        ///
        static void assignONE(volScalarField& s, List<int>& L);

        //--------------------------------------------------------------------------
        /// Assign Boundary Condition to a volVectorField
        ///
        /// @param[in]  s       field where you want to assign the BC.
        /// @param[in]  BC_ind  The BC index.
        /// @param[in]  value   The value you want to assign (it must be an OpenFOAM vector).
        ///
        static void assignBC(volVectorField& s, label BC_ind, Vector<double>& value);

        //--------------------------------------------------------------------------
        /// Assign Boundary Condition to a volScalarField
        ///
        /// @param[in]  s       field where you want to assign the BC.
        /// @param[in]  BC_ind  The BC index.
        /// @param[in]  value   The value you want to assign (it must be a double).
        ///
        static void assignBC(volScalarField& s, label BC_ind, double& value);

        //--------------------------------------------------------------------------
        /// Assign Boundary Condition to a volVectorField
        ///
        /// @param[in]  s       field where you want to assign the BC.
        /// @param[in]  BC_ind  The BC index.
        /// @param[in]  valueVec   The value you want to assign (it must be an OpenFOAM vector).
        ///
        static void assignBC(volVectorField& s, label BC_ind, Eigen::MatrixXd valueVec);

        //--------------------------------------------------------------------------
        /// Assign Boundary Condition to a volScalarField
        ///
        /// @param[in]  s       field where you want to assign the BC.
        /// @param[in]  BC_ind  The BC index.
        /// @param[in]  valueVec   The value you want to assign (it must be a double).
        ///
        static void assignBC(volScalarField& s, label BC_ind, Eigen::MatrixXd valueVec);

        //--------------------------------------------------------------------------
        /// Evaluate the L2 norm of a volScalarField
        ///
        /// @param[in]  field  The field.
        ///
        /// @return     L2 norm of the volScalarField.
        ///
        static double L2norm(volScalarField field);

        //--------------------------------------------------------------------------
        /// Evaluate the L2 norm of a volVectorField
        ///
        /// @param[in]  field  The field.
        ///
        /// @return     L2 norm of the volVectorField.
        ///
        static double L2norm(volVectorField field);

        //--------------------------------------------------------------------------
        /// Evaluate the H1 seminorm of a volScalarField
        ///
        /// @param[in]  field  The field.
        ///
        /// @return     H1 seminorm of the volScalarField.
        ///
        static double H1seminorm(volScalarField field);

        //--------------------------------------------------------------------------
        /// Evaluate the H1 seminorm of a volVectorField
        ///
        /// @param[in]  field  The field.
        ///
        /// @return     H1 seminorm of the volVectorField.
        ///
        static double H1seminorm(volVectorField field);

        //--------------------------------------------------------------------------
        /// Convert to List of volVectorField snapshots to eigen matrix (only internalfield)
        ///
        /// @param      fields  The fields in PtrList format.
        ///
        /// @return     { description_of_the_return_value }
        ///
        static Eigen::MatrixXd foam2eigen(PtrList<volVectorField>& fields);

        //--------------------------------------------------------------------------
        /// Convert to List of volScalarField snapshots to eigen matrix (only internalfield)
        ///
        /// @param      fields  The fields in PtrList format.
        ///
        /// @return     { description_of_the_return_value }
        ///
        static Eigen::MatrixXd foam2eigen(PtrList<volScalarField>& fields);

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcomment"
        //--------------------------------------------------------------------------
        /// Set value of a volScalarField to a constant inside a given box
        ///
        /// @details the Box must be defined with a 2*3 Eigen::MatrixXd in the following way
        ///
        /// \\[
        /// Box=\\begin{bmatrix}
        /// x_1 & y_1 & z_1 \\\\ \
        /// x_2 & y_2 & z_2 \\\\
        /// \\end{bmatrix}
        /// \\]
        ///
        /// where \f$ x_1, y_1, z_1 \f$ and \f$ x_2, y_2, z_2 \f$ are the coordinates of the two corners defining the box.
        ///
        /// @param[in]  field  The field.
        /// @param[in]  Box    The box.
        /// @param[in]  value  The value you want to give to the volScalarField
        ///
        static void setBoxToValue(volScalarField& field, Eigen::MatrixXd Box,
                                  double value);
#pragma GCC diagnostic pop

        //-----------------------------------------------------------------------------
        /// @brief      Gets the mass matrix using the eigen routine
        ///
        /// @param[in]  fields            The openFoam fields
        /// @param[in]  consider_volumes  True if you want to consider the volumes
        ///
        /// @tparam     T                 type of field
        ///
        /// @return     The mass matrix eigen.
        ///
        template<typename T>
        static Eigen::MatrixXd get_mass_matrix_Eigen(PtrList<T>& fields,
                bool consider_volumes = true);

        //--------------------------------------------------------------------------
        /// @brief      Change the boundary condition type for a GeometricField
        ///
        /// @param      field   The field
        /// @param[in]  BCtype  The The BC type you want to set
        /// @param[in]  BC_ind  The BC index
        ///
        /// @tparam     type_f  The type of field it can be scalar or vector
        ///
        template<typename type_f>
        static void changeBCtype(GeometricField<type_f, fvPatchField, volMesh>& field,
                                 word BCtype, label BC_ind);

        //--------------------------------------------------------------------------
        /// @brief      Get the point coordinates and indices from patch.
        ///
        /// @param[in]  mesh     The mesh
        /// @param[in]  ind      The index of the patch from which you want to extract the points
        /// @param[out] points   The coordinates of the points on the patch
        /// @param[out] indices  The indices of the points on the patch
        ///
        static void getPointsFromPatch(fvMesh& mesh, label ind, List<vector>& points,
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
        static List<vector> displacedSegment(List<vector> x0, double mux1, double mux2,
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
        static vector displacePoint(vector x0, vector x_low, vector x_up,
                                    double mux_low, double mux_up,
                                    double muy_low, double muy_up, double muz_low = 0.0, double muz_up = 0.0);

        //--------------------------------------------------------------------------
        /// @brief      Sets some given Indices of a list of objects to given values
        ///
        /// @param      ind2set       The indices of the list that you want to assign
        /// @param      value2set     The values you want to set
        /// @param      movingIDS     The index of the list in the original list that you want to assign
        /// @param      originalList  The original list
        ///
        /// @tparam     T             Type of list
        ///
        template<typename T>
        static void setIndices2Value(labelList& ind2set, List<T>& value2set,
                                     labelList& movingIDS, List<T>& originalList);

        //--------------------------------------------------------------------------
        /// @brief      Returns a scalarField that containes the non-orthogonality value of a given mesh
        ///
        /// @param      mesh      The mesh for which you want to compute the non-orthogonality value
        /// @param      NonOrtho  The scalarField containing the non-orthogonality
        ///
        /// @return     NonOrtho  The scalarField containing the non-orthogonality
        ///
        static volScalarField meshNonOrtho(fvMesh& mesh, volScalarField& NonOrtho);

        //--------------------------------------------------------------------------
        /// @brief      Rotate a list of points in clockwise direction given an axis of rotation and an angle in degrees.
        ///
        /// @param[in]  originalPoints    The original points
        /// @param[in]  CentreOfRotation  The Axis Of Rotation
        /// @param[in]  AngleOfRotation   The angle of rotation in degrees
        ///
        /// @return     The Rotated Points
        ///
        static List<vector> rotatePoints(const List<vector>& originalPoints,
                                         vector AxisOfRotation, double AngleOfRotation);

        //--------------------------------------------------------------------------
        /// @brief      Functions that return a Rotation Matrix given an axis of rotation and an angle in degrees
        ///
        /// @param[in]  AxisOfRotation   The axis of rotation
        /// @param[in]  AngleOfRotation  The angle of rotation in degrees
        ///
        /// @return     The rotation matrix in Eigen format
        ///
        static Eigen::MatrixXd rotationMatrix(vector AxisOfRotation,
                                              double AngleOfRotation);
        
	//--------------------------------------------------------------------------
        /// @brief      Compute the distance between the boundary face center and the boundary cell center
        ///
        /// @param[in]  BC_ind   Index of the boundary patch
        /// @param[in]  mesh     Mesh
        ///
	/// @return     Vector of the face to cell center distance in Eigen format
        ///
        static Eigen::VectorXd boudaryFaceToCellDistance(fvMesh& mesh, label BC_ind);


        //--------------------------------------------------------------------------
        /// @brief      Assign value of a boundary condition of type "mixed"
        ///
        /// @param[in,out]   field      Field
        /// @param[in]       BC_ind     Index of the boundary patch
        /// @param[in]       value      Reference fixed value
        /// @param[in]       grad       Reference patch normal gradient
        /// @param[in]       valueFrac  Weight field
        ///
        /// @tparam          type_f     Type of the field, it can be scalar or vector.
        ///
        template<typename type_f>
        static void assignMixedBC(GeometricField<type_f, fvPatchField, volMesh>& field,
                                  label BC_ind, Eigen::MatrixXd& value,
                                  Eigen::MatrixXd& grad, Eigen::MatrixXd& valueFrac);
};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//

template<>
void ITHACAutilities::assignMixedBC(
    GeometricField<scalar, fvPatchField, volMesh>& field, label BC_ind,
    Eigen::MatrixXd& value,
    Eigen::MatrixXd& grad, Eigen::MatrixXd& valueFrac);

template<>
void ITHACAutilities::assignMixedBC(
    GeometricField<vector, fvPatchField, volMesh>& field, label BC_ind,
    Eigen::MatrixXd& value,
    Eigen::MatrixXd& grad, Eigen::MatrixXd& valueFrac);

template<typename T, typename G>
void ITHACAutilities::assignIF(T& s, G& value)
{
    for (label i = 0; i < s.internalField().size(); i++)
    {
        s.ref()[i] = value;
    }
}

template<typename T>
Eigen::MatrixXd ITHACAutilities::get_mass_matrix_Eigen(PtrList<T>& fields,
        bool consider_volumes)
{
    Eigen::MatrixXd F = Foam2Eigen::PtrList2Eigen(fields);
    Eigen::VectorXd V = Foam2Eigen::field2Eigen(fields[0].mesh().V());
    Eigen::MatrixXd VM = V.asDiagonal();
    Eigen::MatrixXd M;

    if (consider_volumes)
    {
        M = F.transpose() * VM * F;
    }
    else
    {
        M = F.transpose() * F;
    }

    return M;
}

template<>
Eigen::MatrixXd ITHACAutilities::get_mass_matrix_FV(
    GeometricField<scalar, fvPatchField, volMesh>& snapshot);


template<>
PtrList<volVectorField> ITHACAutilities::reconstruct_from_coeff(
    PtrList<volVectorField>& modes, Eigen::MatrixXd& coeff_matrix, label Nmodes);

template<>
PtrList<volScalarField> ITHACAutilities::reconstruct_from_coeff(
    PtrList<volScalarField>& modes, Eigen::MatrixXd& coeff_matrix, label Nmodes);

template<>
void ITHACAutilities::changeBCtype(
    GeometricField<scalar, fvPatchField, volMesh>& field, word BCtype,
    label BC_ind);

template<>
void ITHACAutilities::changeBCtype(
    GeometricField<vector, fvPatchField, volMesh>& field, word BCtype,
    label BC_ind);

template<>
void ITHACAutilities::setIndices2Value(labelList& ind2set,
                                       List<vector>& value2set, labelList& movingIDS, List<vector>& originalList);



#endif







