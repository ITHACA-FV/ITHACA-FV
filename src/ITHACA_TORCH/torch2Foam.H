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

/// \file
/// Header file of the torch2Foam namespace. It contains functions to
/// \dir
/// Directory containing the header and source files for the torch2Foam namespace.

#ifndef torch2Foam_H
#define torch2Foam_H

#include <iostream>
#include <torch/script.h>
#include <torch/torch.h>
#include <Eigen/Eigen>
#include "ITHACAassert.H"
#include "fvCFD.H"

namespace ITHACAtorch
{
namespace torch2Foam
{
//--------------------------------------------------------------------------
/// @brief      Convert an OpenFOAM field to a torch tensor, in case of
///             field<vector> the tensor is arranged as t = [fx;fy;fz]
///
/// @param      field   The field
///
/// @tparam     type_f  type can be scalar or vector
///
/// @return     the torch tensor
///
template <class type_f>
torch::Tensor field2Torch(Field<type_f>& field);

//--------------------------------------------------------------------------
/// @brief      Convert an Torch TensorOpenFOAM to an OpenFoam Field
///
/// @param      tTensor  The Torch tensor
///
/// @tparam     type_f  type can be scalar or vector
///
/// @return     the OF field
///
template <class type_f>
Field<type_f> torch2Field(torch::Tensor& tTensor);

template <class type_f>
torch::Tensor ptrList2Torch(PtrList<type_f>& ptrList);

template <class type_f>
PtrList<Field<type_f>> torch2PtrList(torch::Tensor& tTensor);
}

}

#endif
