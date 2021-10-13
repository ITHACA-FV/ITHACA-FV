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
    ConvLayer
Description
    set of functions to convert Eigen object to Torch objects and viceversa
SourceFiles
    ConvLayer.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the ConvLayer class. It contains functions to

#ifndef ConvLayer_H
#define ConvLayer_H

#include <iostream>
#include <torch/script.h>
#include <torch/torch.h>
#include <Eigen/Eigen>
#include "ITHACAassert.H"
#include "fvCFD.H"
#include "boxToCell.H"
#include "cellSet.H"
#include "Filter.H"


namespace ITHACAtorch
{
template<class Type, template<class> class PatchField, class GeoMesh>
class ConvLayer
{
    public:
        //--------------------------------------------------------------------------
        /// @brief      Construct using Time as functionObject
        ///
        /// @param[in]  db    The database
        ///
        ConvLayer(PtrList<GeometricField<Type, PatchField, GeoMesh>>& snapshots);

        PtrList<GeometricField<Type, PatchField, GeoMesh>>& _snapshots;
        const fvMesh& mesh;

        Vector<label> domainDivision;
        Vector<scalar> filterSize;
        Vector<scalar> domainSize;
        Vector<scalar> ds;

        List<point> convPoints;

        void setFilterSize(double dx, double dy, double dz);
        void setDomainDivision(label Nx, label Ny, label Nz);

        torch::Tensor filter();
        treeBoundBox box;

        List<labelList> cellsInBoxes;
        List<scalarList> weights;

        bool isDomainDivisionSet = false;
        bool isFilterSizeSet = false;

        autoPtr<Filter> flt;
        autoPtr<IOdictionary> convDict;
};

typedef ConvLayer<scalar, fvPatchField, volMesh> volScalarLayer;
typedef ConvLayer<vector, fvPatchField, volMesh> volVectorLayer;
typedef ConvLayer<scalar, fvsPatchField, surfaceMesh> surfaceScalarLayer;

}
#endif
