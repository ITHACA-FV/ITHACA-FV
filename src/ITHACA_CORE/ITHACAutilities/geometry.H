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
    geometry
Description
    Reduction class for a Stationary Navier-Stokes problem.
SourceFiles
    geometry.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the geometry namespace. It contains functions to
/// \dir
/// Directory containing the header and source files for the geometry namespace.

#ifndef geometry_H
#define geometry_H

#include "fvCFD.H"
#include <iostream>
#include <Eigen/Eigen>

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

}

#endif
