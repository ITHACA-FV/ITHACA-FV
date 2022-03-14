/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
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
    Foam::weaklyCoupledFsi


\*---------------------------------------------------------------------------*/

#ifndef weaklyCoupledFsi_H
#define weaklyCoupledFsi_H


#include "forces.H"
#include "volFieldsFwd.H"
#include "Tuple2.H"
#include "OFstream.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

        // Forward declaration of classes
        class objectRegistry;
        class dictionary;
        class polyMesh;
        class mapPolyMesh;

        namespace functionObjects
        {

                /*---------------------------------------------------------------------------*\
                                           Class weaklyCoupledFsi Declaration
                \*---------------------------------------------------------------------------*/

                class weaklyCoupledFsi: public forces
                {

                public:

                        //- cylinder mass
                        scalar M_;
                        
                        //- damping coefficient
                        scalar C_;
                        
                        //- rigidity coefficient
                        scalar K_;
                        
                        //- ratio of real cylinder length to depth of domain
                        scalar R_;

                        //- maximum amplitude of displacement
                        scalar Ymax_;

                        //- true if after restart data should be appended to previous log 
                        //  false if log should be overwritten
                        bool append_;

                        //- current state of system (y,Vy)
                        Pair<scalar> Y_; 
                        
                        //- old state of system (y, Vy)
                        Pair<scalar> Yold_;
                        
                        //- pointer to stream with results
                        autoPtr<OFstream> fsiResPtr_;

                public:

                        //- Disallow default bitwise copy construct
                        weaklyCoupledFsi(const weaklyCoupledFsi&);

                        //- Disallow default bitwise assignment
                        void operator=(const weaklyCoupledFsi&);
                        
                        //-
                        void createFsiOutFile(const dictionary& dict);
                        
                        //- Creates out new file, appending previous results if needed
                        void createFsiOutFile(const dictionary& dict, bool app);


                public:

                    //- Runtime type information
                    TypeName("weaklyCoupledFsi");


                    // Constructors

                        //- Construct from given Time and dictionary.
                        weaklyCoupledFsi
                        (
                            const word& name,
                            const Time&,
                            const dictionary&
                        );
                        
                        //- Construct from given objectRegistry and dictionary.
                        weaklyCoupledFsi
                        (
                            const word& name,
                            const objectRegistry&,
                            const dictionary&
                        );

                    //- Destructor
                    virtual ~weaklyCoupledFsi();


                    // Member Functions

                        //- Read the weaklyCoupledFsi data
                        virtual bool read(const dictionary&);
                        
                        //-
                        virtual bool execute();
                        
                        //- Write the weaklyCoupledFsi
                        virtual bool write();
                        
                        //- Distributes displacements between slave processes
                        // and sets cellDisplacement field Y component on patch
                        void setDisplacements(volVectorField& yDispl);
                };


        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        } // End namespace functionObjects

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
