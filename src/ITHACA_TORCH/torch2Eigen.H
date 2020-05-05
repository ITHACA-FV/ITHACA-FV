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
    torch2Eigen
Description
    set of functions to convert Eigen object to Torch objects and viceversa
SourceFiles
    torch2Eigen.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the torch2Eigen namespace. It contains functions to
/// \dir
/// Directory containing the header and source files for the torch2Eigen namespace.

#ifndef torch2Eigen_H
#define torch2Eigen_H

#include <iostream>
#include <torch/script.h>
#include <torch/torch.h>
#include <Eigen/Eigen>
#include "ITHACAassert.H"

namespace ITHACAtorch
{
namespace torch2Eigen
{
///
/// @brief      Convert an eigen Matrix to a torch tensor
///
/// @param[in]  eigenMatrix  The eigen matrix
///
/// @tparam     type         Can be double, float, int
///
/// @return     a matrix in Eigen float format
///
template<class type>
torch::Tensor eigenMatrix2torchTensor(
    Eigen::Matrix<type, Eigen::Dynamic, Eigen::Dynamic> eigenMatrix);


///
/// @brief      Convert a torch tensor to an eigen Matrix
///
/// @param      torchTensor  The torch tensor
///
/// @tparam     type         Can be double, float, int
///
/// @return     a torch tensor
///
template<class type>
Eigen::Matrix<type, Eigen::Dynamic, Eigen::Dynamic> torchTensor2eigenMatrix(
    torch::Tensor& torchTensor);
}

}

#endif