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
};
#endif
