/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::dynamicFvMesh

Description
    Abstract base class for geometry and/or topology changing fvMesh.

    Supports optional update controls that may be used by custom solvers:
    \table
        Property | Description                              | Required | Default
        updateControl   | See time controls below           | no  | timeStep
        updateInterval  | Steps/time between update phases  | no  | 1
    \endtable

See also
    Foam::timeControl

SourceFiles
    dynamicFvMesh.C
    dynamicFvMeshNew.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicFvMesh_H
#define dynamicFvMesh_H

#include "fvMesh.H"
#include "timeControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class dynamicFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicFvMesh
:
    public fvMesh
{
    // Private Data

        //- Optional update control
        timeControl timeControl_;


    // Private Member Functions

        //- Read the updateControl/updateInterval from dynamicMeshDict
        void readDict();

        //- No copy construct
        dynamicFvMesh(const dynamicFvMesh&) = delete;

        //- No copy assignment
        void operator=(const dynamicFvMesh&) = delete;


public:

    //- Runtime type information
    TypeName("dynamicFvMesh");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            dynamicFvMesh,
            IOobject,
            (const IOobject& io),
            (io)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            dynamicFvMesh,
            doInit,
            (
                const IOobject& io,
                const bool doInit
            ),
            (io, doInit)
        );


    // Constructors

        //- Construct from an IOobject
        explicit dynamicFvMesh(const IOobject& io, const bool doInit=true);

        //- Construct from components without boundary.
        //  Boundary is added using addFvPatches() member function
        dynamicFvMesh(const IOobject& io, const zero, const bool syncPar=true);

        //- Construct from components without boundary.
        //  Boundary is added using addFvPatches() member function
        dynamicFvMesh
        (
            const IOobject& io,
            pointField&& points,
            faceList&& faces,
            labelList&& allOwner,
            labelList&& allNeighbour,
            const bool syncPar = true
        );

        //- Construct without boundary from cells rather than owner/neighbour.
        //  Boundary is added using addPatches() member function
        dynamicFvMesh
        (
            const IOobject& io,
            pointField&& points,
            faceList&& faces,
            cellList&& cells,
            const bool syncPar = true
        );


    // Selectors

        //- Select, construct and return the dynamicFvMesh
        //  If the constant/dynamicMeshDict does not exist
        //  a staticFvMesh is returned
        static autoPtr<dynamicFvMesh> New(const IOobject& io);


        //- Select, construct and return the dynamicFvMesh
        //  If the constant/dynamicMeshDict does not exist
        //  a staticFvMesh is returned
        static autoPtr<dynamicFvMesh> New
        (
            const argList& args,
            const Time& runTime
        );


    //- Destructor
    virtual ~dynamicFvMesh() = default;


    // Member Functions

        //- Initialise all non-demand-driven data
        virtual bool init(const bool doInit);

        //- Is mesh dynamic
        virtual bool dynamic() const
        {
            return true;
        }

        //- Update the mesh if controller permits
        virtual bool controlledUpdate();

        //- Update the mesh for both mesh motion and topology change
        virtual bool update() = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
