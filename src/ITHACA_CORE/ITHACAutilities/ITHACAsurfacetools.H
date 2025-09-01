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

#ifndef functionObjects_ITHACAsurfacetools_H
#define functionObjects_ITHACAsurfacetools_H

#include <Eigen/Eigen>
#include "volPointInterpolation.H"
#include "Foam2Eigen.H"
#include "fvMesh.H"
#include "fvMeshSubset.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace ITHACAutilities
{
namespace ITHACAsurfacetools
{

template<typename T>
List<label> surfaceIndexInt(T& field, const label patchInt,
                            const label patchExt);

template<typename T>
List<label> surfaceIndexExt(T& field, const label patchInt,
                            const label patchExt);

template<typename T, typename V>
void surfaceValuesInt(T& field, const label patchInt, const label patchExt,
                      List<V> & result);

template<typename T, typename V>
void surfaceValuesExt(T& field, const label patchInt, const label patchExt,
                      List<V> & result);

template<typename T>
Foam::Vector<scalar> surfaceFindMirrorPoint(T& field, const label patchInt,
        const label patchExt, const label cellID);

template<typename T>
label surfaceFindClosest(T& field, const label patchInt, const label patchExt,
                         Foam::Vector<scalar> point);

template<typename T, typename V>
void surfaceAverage(T& field, const label patchInt, const label patchExt,
                    List<V> & result);

template<typename T, typename V>
void surfaceJump(T& field, const label patchInt, const label patchExt,
                 List<V> & result);

}; // End namespace ITHACAsurfacetools
} // End namespace ITHACAutilities

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
