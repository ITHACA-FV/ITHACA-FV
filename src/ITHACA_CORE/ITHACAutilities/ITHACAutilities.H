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
#include "ITHACAstream.H"
#include "ITHACAparameters.H"
#include "ITHACAgeometry.H"
#include "ITHACAsystem.H"
#include "ITHACAerror.H"
#include "ITHACAassign.H"
#include "ITHACAcoeffsMass.H"


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
Eigen::MatrixXd rand(label rows, label cols, double min, double max);

//--------------------------------------------------------------------------
/// @brief      Generates random matrix with random values in an interval
///
/// @param[in]  rows    Number of rows
/// @param[in]  minMax  The min and max Matrix, it has as many rows as columns of the output matrix.
///
/// @return     Random matrix
///
Eigen::MatrixXd rand(label rows, Eigen::MatrixXd minMax);


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

//------------------------------------------------------------------------------
/// @brief      Combine a list of list into a single list with unique elements
///
/// @param      doubleList  The list of list
///
/// @tparam     T           type of element label, double, etc...
///
/// @return     list in unique order
///
template<typename T>
List<T> combineList(List<List<T >> & doubleList);

//------------------------------------------------------------------------------
/// @brief      Using the Eigen library, using the SVD decomposition method to solve the
///             matrix pseudo-inverse, the default error er is 0
///
/// @param      origin  Matrix to invert
/// @param      er      Error
///
/// @return     Inversed matrix
///
Eigen::MatrixXd pinv_eigen_based(Eigen::MatrixXd& origin, const float er = 0);

//------------------------------------------------------------------------------
/// @brief      Invert a matrix given the method name in the ITHACAdict
///
/// @param[in]      matrixToInvert  Matrix to invert
/// @param[in]      methodName      Method
///
/// @return         Inversed matrix
///
Eigen::MatrixXd invertMatrix(Eigen::MatrixXd& matrixToInvert,
                             const word inversionMethod);



//------------------------------------------------------------------------------
/// @brief      Multiplication between a field of type vol[Scalar|Vector|Tensor]Field
///             and a double
///
/// @param      f1  field 
/// @param      alpha double value
///
template<typename T>
void multField(T& f1, double alpha);

//------------------------------------------------------------------------------
/// @brief      Multiplication between field of type vol[Scalar|Vector|Tensor]Field
///             and a vector
///
/// @param      f1  field 
/// @param      alphaVec vector filled with double
///
template<typename T>
void multField(T &f1, const Eigen::VectorXd alphaVec);

//------------------------------------------------------------------------------
/// @brief      Multiplication between a list of field of type
///             vol[Scalar|Vector|Tensor]Field
///
///
/// @param      f1 PtrList of field 
/// @param      f2 vector filled with double
///
template<typename T>
void multField(PtrList<T> &f1, const Eigen::VectorXd alphaVec);


//------------------------------------------------------------------------------
/// @brief      Perform the following operation 
///             f1 + f2 * alpha with f1 and f2 being of type
///             vol[Scalar|Vector|Tensor]Field and alpha a double
///
///
/// @param      f1 field 
/// @param      f2 field
/// @param      alpha double
///
template<typename T>
void addFields(T& f1, const T& f2c, double alpha = 1.0);

//------------------------------------------------------------------------------
/// @brief      Perform the substraction (f1 - f2) between two fields 
///              of type vol[Scalar|Vector|Tensor]Field and alpha a double
///
/// @param      f1 field 
/// @param      f2 field
///
template<typename T>
void subtractFields(T& f1, const T& f2);

//------------------------------------------------------------------------------
/// @brief      Set a field of type vol[Scalar|Vector|Tensor]Field to 0
///
/// @param      f1 field to set to 0 
///
template<typename T>
void setToZero(T& f1);

//------------------------------------------------------------------------------
/// @brief      Tensor field product between a volScalarField and a volTensorField
///
/// @param      coef field
/// @param      f1 field 
///
volTensorField tensorFieldProduct(const volScalarField& coef, const volTensorField& S);

//------------------------------------------------------------------------------
/// @brief      Tensor field product between two volTensorFields
///
/// @param      coef field 
/// @param      S field 
///
volTensorField tensorFieldProduct(const volTensorField& coef, const volTensorField& S);

//------------------------------------------------------------------------------
/// @brief      Return the dimension of a volTensorField
///
/// @param      v 
int dimensionField(const volTensorField& v);

//------------------------------------------------------------------------------
/// @brief      Return the dimension of a volVectorField
///
/// @param      v 
int dimensionField(const volVectorField& v);

//------------------------------------------------------------------------------
/// @brief      Return the dimension of a volScalarField
///
/// @param      v 
int dimensionField(const volScalarField& v);



//------------------------------------------------------------------------------
/// @brief      Trim a string cf. 
///             https://stackoverflow.com/questions/25829143/trim-whitespace-from-a-string
///
/// @param      String s
std::string str_trim(std::string const& s);

//------------------------------------------------------------------------------
/// @brief      Format a string
///
/// @param      s string
/// @param      Nmax max length of the string
void str_format_io( std::string const& s, unsigned int Nmax );

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
double norm_L2(const T v);

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

};

#endif
