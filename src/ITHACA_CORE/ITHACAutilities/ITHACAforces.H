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

#ifndef functionObjects_ITHACAforces_H
#define functionObjects_ITHACAforces_H

#include "fvMeshFunctionObject.H"
#include "logFiles.H"
#include "coordinateSystem.H"
#include "volFieldsFwd.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class forces Declaration
\*---------------------------------------------------------------------------*/

class ITHACAforces
    :
    public fvMeshFunctionObject,
    public logFiles
{

    protected:

        // Protected data

        //- Enumeration for ensuring the right file is accessed
        enum fileID
        {
            MAIN_FILE = 0,
            BINS_FILE = 1
        };

        //- Pressure, viscous and porous force per bin
        List<Field<vector>> force_;

        //- Pressure, viscous and porous moment per bin
        List<Field<vector>> moment_;


        // Read from dictionary

        //- Patches to integrate forces over
        labelHashSet patchSet_;

        //- Name of pressure field
        word pName_;

        //- Name of velocity field
        word UName_;

        //- Name of density field (optional)
        word rhoName_;

        //- Is the force density being supplied directly?
        Switch directForceDensity_;

        //- The name of the force density (fD) field
        word fDName_;

        //- Reference density needed for incompressible calculations
        scalar rhoRef_;

        //- Reference pressure
        scalar pRef_;

        //- Coordinate system used when evaluting forces/moments
        coordinateSystem coordSys_;

        //- Flag to indicate whether we are using a local co-ordinate sys
        bool localSystem_;

        //- Flag to include porosity effects
        bool porosity_;


        // Bin information

        //- Number of bins
        label nBin_;

        //- Direction used to determine bin orientation
        vector binDir_;

        //- Distance between bin divisions
        scalar binDx_;

        //- Minimum bin bounds
        scalar binMin_;

        //- Bin positions along binDir
        List<point> binPoints_;

        //- Should bin data be cumulative?
        bool binCumulative_;


        //- Initialised flag
        bool initialised_;


        // Protected Member Functions

        //- Create file names for forces and bins
        wordList createFileNames(const dictionary& dict) const;

        //- Output file header information
        virtual void writeFileHeader(const label i);

        //- Initialise the fields
        void initialise();

        //- Return the effective viscous stress (laminar + turbulent).
        tmp<volSymmTensorField> devRhoReff() const;

        //- Dynamic viscosity field
        tmp<volScalarField> mu() const;

        //- Return rho if specified otherwise rhoRef
        tmp<volScalarField> rho() const;

        //- Return rhoRef if the pressure field is dynamic, i.e. p/rho
        //  otherwise return 1
        scalar rho(const volScalarField& p) const;

        //- Accumulate bin data
        void applyBins
        (
            const vectorField& Md,
            const vectorField& fN,
            const vectorField& fT,
            const vectorField& fP,
            const vectorField& d
        );

        //- Helper function to write force data
        void writeForces();

        //- Helper function to write bin data
        void writeBins();

        //- Disallow default bitwise copy construct
        ITHACAforces(const ITHACAforces&);

        //- Disallow default bitwise assignment
        void operator=(const ITHACAforces&);


    public:

        //- Runtime type information
        TypeName("ITHACAforces");


        // Constructors

        //- Construct from Time and dictionary
        ITHACAforces
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Construct from objectRegistry and dictionary
        ITHACAforces
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary&
        );


        //- Destructor
        virtual ~ITHACAforces();


        // Member Functions

        //- Read the forces data
        virtual bool read(const dictionary&);

        //- Calculate the forces and moments
        virtual void calcForcesMoment();

        //- Return the total force
        virtual vector forceEff() const;

        //- Return the total moment
        virtual vector momentEff() const;

        //- Return the viscous forces
        virtual vector forceTau() const;

        //- Return the pressure forces
        virtual vector forcePressure() const;

        //- Return the porous forces
        virtual vector forcePorous() const;

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the forces
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
