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
Namespace
    ITHACAutilites
Description
    Utilities to manage files and check existence of folders (system operations)
SourceFiles
    ITHACAsystem.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the ITHACAsystem file.
/// \dir
/// Directory containing the header and source files for the ITHACAsystem files.

#ifndef ITHACAsystem_H
#define ITHACAsystem_H

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

/// Namespace to implement some useful system calls inside ITHACA-FV
namespace ITHACAutilities
{
//--------------------------------------------------------------------------
/// @brief      Creates a symbolic links to 0, system and constant
///
/// @param[in]  folder  The folder where to create the symbolic links
///
void createSymLink(word folder);

//--------------------------------------------------------------------------
/// Checks if a folder exists
///
/// @param[in]  folder  The folder you want to check.
///
/// @return     1 if the folder exists 0 elsewhere.
///
bool check_folder(word folder);

//--------------------------------------------------------------------------
/// Check if the offline data folder "./ITHACAoutput/Offline" exists
///
/// @return     1 if the folder exists 0 elsewhere.
///
bool check_off();

//--------------------------------------------------------------------------
/// Function that returns true if a file exists
///
/// @param[in]  fileName  The name of the file
///
/// @return     1 if the file exists 0 elsewhere.
///
bool check_file(std::string fileName);

//--------------------------------------------------------------------------
/// Check if the POD data folder "./ITHACAoutput/POD" exists
///
/// @return     1 if the folder exists 0 elsewhere.
///
bool check_pod();

//--------------------------------------------------------------------------
/// Check if the supremizer folder exists
///
/// @return     1 if the folder exists 0 elsewhere.
///
bool check_sup();

}

#endif
