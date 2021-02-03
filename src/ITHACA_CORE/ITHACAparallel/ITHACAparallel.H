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
#include "ITHACAassert.H"

/// Class for parallel handling, it has several functions to deal
/// with parallel problems, the class is implemented as a singletone
class ITHACAparallel
{
        ///
        /// @brief      Constructs the object.
        ///
        /// @param      mesh  A mesh of the problem
        ///
        ITHACAparallel(fvMesh& mesh, Time& localTime);

        /// Pointer to the ITHACAparallel class
        static ITHACAparallel* instance;

    public:
        ///
        /// @brief      Gets an instance of ithacaparallel, to be used if the instance is not existing.
        ///
        /// @param      mesh       The mesh
        /// @param      localTime  The runtime object
        ///
        /// @return     The instance.
        ///
        static ITHACAparallel* getInstance(fvMesh& mesh, Time& localTime);

        ///
        /// @brief      Gets an instance of ithacaparallel, to be used if the instance is already existing.
        ///
        /// @return     The instance.
        ///
        static ITHACAparallel* getInstance();

        /// Delete empty constructor
        ITHACAparallel() = delete;

        /// Function to sumpend MPI
        static void suspendMPI();

        /// Function to resume MPI
        static void resumeMPI();

        /// Function to get a global field from a parallel one
        template<class Type>
        List<Field <Type>> combineFields(GeometricField<Type, fvPatchField, volMesh>&
                                         field);

        template<class type>
        GeometricField<type, fvPatchField, volMesh> constructGlobalField(
            GeometricField<type, fvPatchField, volMesh> field);

        /// Totoal number of internal field cells
        label N_IF_glob;

        /// Number of Boundary Patches without considering processors
        label N_BF;

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

        /// ID of the OLD process
        static List<int> oldProcIDs_;

        /// ID of the NEW process
        static List<int> newProcIDs_;

        /// runTime defined locally
        Time& runTime;

        /// Mesh object defined locally
        fvMesh& mesh;

};

template<>
List<Field<scalar>> ITHACAparallel::combineFields(
                     GeometricField<scalar, fvPatchField, volMesh>& field);

template<class type>
GeometricField<type, fvPatchField, volMesh>
ITHACAparallel::constructGlobalField(GeometricField<type, fvPatchField, volMesh>
                                     field)
{
    autoPtr<Time> globalTimePtr_;
    autoPtr<fvMesh> globalMeshPtr_;
    Time& localTime = runTime;
    fvMesh& localMesh = mesh;
    const word globalConstant = localTime.rootPath() + "/" +
                                localTime.globalCaseName() + "/constant";
    const fileName gRootPath = localTime.rootPath();
    const fileName gCaseName = localTime.globalCaseName();
    globalTimePtr_.reset
    (
        new Time
        (
            gRootPath,
            gCaseName
        )
    );
    Time& globalTime = globalTimePtr_();
    globalMeshPtr_.reset
    (
        new fvMesh
        (
            IOobject
            (
                fvMesh::defaultRegion,
                globalTime.timeName(),
                globalTime,
                IOobject::MUST_READ
            )
        )
    );
    GeometricField<type, fvPatchField, volMesh> F_glob
    (
        IOobject
        (
            field.name(),
            globalTimePtr_().timeName(),
            globalMeshPtr_(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        globalMeshPtr_()
    );
    return F_glob;
}
#endif


