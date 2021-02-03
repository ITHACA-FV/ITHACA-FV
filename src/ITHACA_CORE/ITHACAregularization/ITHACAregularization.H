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
    ITHACAregularization
Description
    Methods for regularizing ill-conditione linear systems
SourceFiles
    ITHACAregularization.C
    ITHACAregularizationTemplates.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the ITHACAregularization class, it contains the implementation of
/// several methods for input output operations.
/// \dir
/// Directory containing the header, source and template files for the ITHACAregularization class, ITHACAparameters and Assert Functions.

#ifndef ITHACAregularization_H
#define ITHACAregularization_H

#include "fvCFD.H"
#include "IOmanip.H"
#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>
#include <algorithm>
#include <fstream>
#include <string>
#include <stdexcept>
#include <sstream>
#include <vector>
#include <cstdio>
#include <typeinfo>
#include <iostream>
#include <cassert>
#include <zlib.h>
#include <map>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wignored-attributes"
#include <Eigen/Eigen>
#include "EigenFunctions.H"
#include "ITHACAutilities.H"


#include <unsupported/Eigen/CXX11/Tensor>
#pragma GCC diagnostic pop

#define MAXBUFSIZE (static_cast<int> (1e6))
#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


/*---------------------------------------------------------------------------*\
                        Class reductionProblem Declaration
\*---------------------------------------------------------------------------*/
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic pop

/// Namespace for regularization of ill-conditione linear system
namespace ITHACAregularization
{

//--------------------------------------------------------------------------
/// @brief      Truncated Singular Value regularization
///
/// @param[in]  A      Matrix of coefficient
/// @param[in]  b      Column vector of solutions
/// @param[in]  filter Number of singular values to keep
///
/// @return     Column vector of variables
///
Eigen::VectorXd  TSVD(Eigen::MatrixXd A, Eigen::MatrixXd b,
                      int filter);

//--------------------------------------------------------------------------
/// @brief      Truncated Singular Value regularization
///
/// @param[in]  A                Matrix of coefficient
/// @param[in]  b                Column vector of solutions
/// @param[in]  parameterMethod  Regularization parameter selection method
///
/// @return     Column vector of variables
///
Eigen::VectorXd  TSVD(Eigen::MatrixXd A, Eigen::MatrixXd b,
                      double noiseVariance, word parameterMethod);

//--------------------------------------------------------------------------
/// @brief      Truncated Singular Value regularization
///
/// @param[in]  A      Matrix of coefficient
/// @param[in]  b      Column vector of solutions
///
/// @return     Column vector of variables
///
Eigen::VectorXd  Tikhonov(Eigen::MatrixXd A, Eigen::MatrixXd b,
                          double regularizationParameter);

};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



#endif
