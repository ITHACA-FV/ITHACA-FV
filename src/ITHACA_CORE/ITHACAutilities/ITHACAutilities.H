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
List<T> combineList(List<List<T>>& doubleList);

//------------------------------------------------------------------------------
/// @brief      Using the Eigen library, using the SVD decomposition method to solve the
///             matrix pseudo-inverse, the default error er is 0
///
/// @param      origin  Matrix to invert
/// @param      er      Error
///
/// @return     Inversed matrix
///
Eigen::MatrixXd pinv_eigen_based(Eigen::MatrixXd &origin, const float er = 0);

//------------------------------------------------------------------------------
/// @brief      Invert a matrix given the method name in the ITHACAdict
///
/// @param[in]      matrixToInvert  Matrix to invert
/// @param[in]      methodName      Method 
///
/// @return         Inversed matrix
///
Eigen::MatrixXd invertMatrix(Eigen::MatrixXd& matrixToInvert, const word inversionMethod);



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




};
#endif
