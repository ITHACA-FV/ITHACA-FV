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
#include "writeFile.H"
#include "cartesianCS.H"
#include "volFieldsFwd.H"
#include "HashSet.H"
#include "Tuple2.H"
#include "OFstream.H"
#include "Switch.H"

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
    public writeFile
{

    protected:

        // Protected data

        //- Pressure, viscous and porous force per bin
        List<Field<vector>> force_;

        //- Pressure, viscous and porous moment per bin
        List<Field<vector>> moment_;

        // File streams

        //- Forces
        autoPtr<OFstream> forceFilePtr_;

        //- Moments
        autoPtr<OFstream> momentFilePtr_;

        //- Force bins
        autoPtr<OFstream> forceBinFilePtr_;

        //- Moment bins
        autoPtr<OFstream> momentBinFilePtr_;

        //- Local force
        autoPtr<OFstream> localForceFilePtr_;

        //- Local moment
        autoPtr<OFstream> localMomentFilePtr_;

        //- Local force bins
        autoPtr<OFstream> localForceBinFilePtr_;

        //- Local moment bins
        autoPtr<OFstream> localMomentBinFilePtr_;


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
        coordSystem::cartesian coordSys_;

        //- Flag to indicate whether we are using a local coordinates
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


        //- Write fields flag
        bool writeFields_;

        //- Initialised flag
        bool initialised_;


        // Protected Member Functions

        //- Create a field name
        word fieldName(const word& name) const;

        //- Create the output files
        void createFiles();

        //- Write header for integrated data
        void writeIntegratedHeader(const word& header, Ostream& os) const;

        //- Write header for binned data
        void writeBinHeader(const word& header, Ostream& os) const;

        //- Initialise the fields
        void initialise();

        //- Initialise the collection bins
        void initialiseBins();

        //- Reset the fields prior to accumulation of force/moments
        void resetFields();

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

        //- Add patch contributions to force and moment fields
        void addToFields
        (
            const label patchi,
            const vectorField& Md,
            const vectorField& fN,
            const vectorField& fT,
            const vectorField& fP
        );

        //- Add cell contributions to force and moment fields
        void addToFields
        (
            const labelList& cellIDs,
            const vectorField& Md,
            const vectorField& fN,
            const vectorField& fT,
            const vectorField& fP
        );

        //- Helper function to write integrated forces and moments
        void writeIntegratedForceMoment
        (
            const string& descriptor,
            const vectorField& fm0,
            const vectorField& fm1,
            const vectorField& fm2,
            autoPtr<OFstream>& osPtr
        ) const;

        //- Write force data
        void writeForces();

        //- Helper function to write binned forces and moments
        void writeBinnedForceMoment
        (
            const List<Field<vector>>& fm,
            autoPtr<OFstream>& osPtr
        ) const;

        //- Write binned data
        void writeBins();

        //- No copy construct
        ITHACAforces(const ITHACAforces&) = delete;

        //- No copy assignment
        void operator=(const ITHACAforces&) = delete;


    public:

        //- Runtime type information
        TypeName("ITHACAforces");


        // Constructors

        //- Construct from Time and dictionary
        ITHACAforces
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict,
            const bool readFields = true
        );

        //- Construct from objectRegistry and dictionary
        ITHACAforces
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict,
            const bool readFields = true
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
