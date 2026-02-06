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

#ifndef ITHACAfieldsOperations_H
#define ITHACAfieldsOperations_H

/// \file
/// Header file of the ITHACAfieldsOperations file.

#include "fvCFD.H"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#include <Eigen/Eigen>
#pragma GCC diagnostic pop
using namespace std::placeholders;
#include "ITHACAstream.H"
#include "ITHACAparameters.H"
#include "ITHACAgeometry.H"
#include "ITHACAsystem.H"
#include "ITHACAassign.H"


namespace ITHACAutilities
{

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


//-----------------------------------------------------------------------------
/// @brief      A function to compute time-averaged fields for a set of different parameter samples
/// and also the fields with the corresponding averaged subtracted
///
/// @param[in]  fields     The fields from which the time-averaged fields have to be computed
/// @param[in]  ind        The indices at which fields for different samples of the parameter start
/// @param      ave        The computed time-averaged fields
///
/// @tparam     TypeField  Type of field
///
/// @return     A list of fields which correspond to the original fields subtracted by
/// the time-averaged part
///
template<class TypeField>
PtrList<TypeField> averageSubtract(PtrList<TypeField>
                                   fields, Eigen::MatrixXd ind, PtrList<TypeField>& ave);

 
//-----------------------------------------------------------------------------
/// @brief      Calculates the average of a list of fields
///
/// @param      fields     The fields
///
/// @tparam     TypeField  { description }
///
/// @return     The average.
///
template<class TypeField>
TypeField computeAverage(PtrList<TypeField>& fields);


//--------------------------------------------------------------------------
/// @brief      Normalize list of Geometric fields
///
/// @param      fields  The fields to be normalized
///
/// @tparam     Type  type of the Field can be Scalar or Vector
///
template<typename Type>
void normalizeFields(
    PtrList<GeometricField<Type, fvPatchField, volMesh >>& fields);


//------------------------------------------------------------------------------
/// @brief      /// @brief      Gets the values from a list of indices
///
/// @param      field    The field from which you want to extract the values
/// @param      indices  The indices at the point
/// @param      xyz      Definition of the x, y, z coordinate of the value you want to extract 0 for x, 1 for y, 2 for z
///
/// @tparam     Type     type of the Field can be Scalar or Vector
///
/// @return     The values.
///
template<typename Type>
Eigen::MatrixXd getValues(GeometricField<Type, fvPatchField,
                          volMesh>& field, labelList& indices, labelList* xyz = NULL);


//------------------------------------------------------------------------------
/// @brief      Gets the values from a list of indices
///
/// @param      field    The field
/// @param      indices  The indices
/// @param      xyz      The xyz
///
/// @tparam     Type     { description }
///
/// @return     The values.
///
template<typename Type>
Eigen::MatrixXd getValues(PtrList<GeometricField<Type, fvPatchField,
                          volMesh >>& field, labelList& indices, labelList* xyz = NULL);

};

#endif
