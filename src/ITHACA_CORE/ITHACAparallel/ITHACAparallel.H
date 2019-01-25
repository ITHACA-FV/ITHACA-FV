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
    ITHACAparallel
Description
    Methods for to deal with parallel computations
SourceFiles
    ITHACAparallel.C
\*---------------------------------------------------------------------------*/
#ifndef ITHACAparallel_H
#define ITHACAparallel_H

#include "fvCFD.H"

class ITHACAparallel
{
    public:
        /// Delete empty constructor
        ITHACAparallel() = delete;

        ///
        /// @brief      Constructs the object.
        ///
        /// @param      field  A field can be either a volvectorfield or a volscalarfield.
        ///
        ITHACAparallel(fvMesh& mesh);


        /// Function to sumpend MPI
        static void suspendMPI();

        /// Function to resume MPI
        static void resumeMPI();

        /// Totoal number of internal field cells
        int N_IF_glob;

        int N_BF;

        /// Cell proc addressing
        autoPtr<labelIOList> indices;

        /// Face proc addressing
        autoPtr<labelIOList> indicesF;

        /// Boundary Field variables representing the total number of faces for each processor
        autoPtr<labelList> Gsize_BF;

        /// Indexing
        autoPtr< List<labelList>> IndFaceLocal;

        /// Local size at the boundary
        autoPtr<labelList> Lsize_BF;

        /// StartFace on the gloabl Mesh
        autoPtr<labelList> Start;

        //-
        static List<label> oldProcIDs_;

        //-
        static List<label> newProcIDs_;

};

#endif


