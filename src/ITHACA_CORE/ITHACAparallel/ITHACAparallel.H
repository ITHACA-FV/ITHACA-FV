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

template<class type_of_field>
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
        ITHACAparallel(GeometricField<type_of_field, fvPatchField, volMesh>& field);

        /// Function to sumpend MPI
        void suspendMPI();

        /// Function to resume MPI
        static void resumeMPI();

        /// Totoal number of internal field cells
        int N_IF_glob;

        /// Cell proc addressing
        autoPtr<labelIOList> indices;

        // Boundary Field variables representing the total number of faces for each processor
        autoPtr< List<labelList>> Gsize_BF;

        /// It gives for each boundary the starting face number
        autoPtr<List< labelList >> startFaces;

        /// ID of the old process
        static List<label> oldProcIDs_;

        /// ID of the new process
        static List<label> newProcIDs_;
};

template<class type_of_field>
List<label> ITHACAparallel<type_of_field>::oldProcIDs_(0);

template<class type_of_field>
List<label> ITHACAparallel<type_of_field>::newProcIDs_(0);

template<class type_of_field>
ITHACAparallel<type_of_field>::ITHACAparallel(
    GeometricField<type_of_field, fvPatchField, volMesh>& field)
{
    // IF construction
    indices = new labelIOList(
        IOobject
        (
            "cellProcAddressing",
            field.mesh().facesInstance(),
            field.mesh().meshSubDir,
            field.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    N_IF_glob = (field.size());
    reduce(N_IF_glob, sumOp<int>());
    // BF construction
    startFaces = new List<labelList> (field.boundaryFieldRef().size(),
                                      labelList(Pstream::nProcs(), 0));
    Gsize_BF = new List<labelList> (field.boundaryFieldRef().size(),
                                    labelList(Pstream::nProcs(), 0));

    for (int i = 0; i < field.boundaryFieldRef().size(); i++)
    {
        Gsize_BF()[i][Pstream::myProcNo()] = field.boundaryFieldRef()[i].size();
        reduce(Gsize_BF()[i], sumOp<labelList>());
    }

    for (int i = 0; i < field.boundaryFieldRef().size(); i++)
    {
        for (label k = 1; k < Pstream::nProcs(); k++)
        {
            startFaces()[i][k] = startFaces()[i][k - 1] + Gsize_BF()[i][k - 1];
        }
    }
}

template<class type_of_field>
void ITHACAparallel<type_of_field>::suspendMPI()
{
    Pstream::parRun() = false;
    label comm        = Pstream::worldComm;
    oldProcIDs_       = Pstream::procID(comm);
    newProcIDs_       = List<label> (1, 0);
    Pstream::procID(comm) = newProcIDs_;
}

template<class type_of_field>
void ITHACAparallel<type_of_field>::resumeMPI()
{
    label comm  = Pstream::worldComm;
    Pstream::procID(comm) = oldProcIDs_;
    Pstream::parRun() = true;
}
#endif


