/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

Application
    icoFoam

Group
    grpIncompressibleSolvers

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids.

    \heading Solver details
    The solver uses the PISO algorithm to solve the continuity equation:

        \f[
            \div \vec{U} = 0
        \f]

    and momentum equation:

        \f[
            \ddt{\vec{U}}
          + \div \left( \vec{U} \vec{U} \right)
          - \div \left(\nu \grad \vec{U} \right)
          = - \grad p
        \f]

    Where:
    \vartable
        \vec{U} | Velocity
        p       | Pressure
    \endvartable

    \heading Required fields
    \plaintable
        U       | Velocity [m/s]
        p       | Kinematic pressure, p/rho [m2/s2]
    \endplaintable

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pisoControl.H"
#include "Foam2Eigen.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for incompressible, laminar flow"
        " of Newtonian fluids."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    Foam::Info<< "Create time\n" << Foam::endl;
    Foam::Time runTime(Foam::Time::controlDictName, args);
    Info << "Create mesh for time = "
     << runTime.timeName() << nl << endl;
     Foam::autoPtr<Foam::fvMesh> meshPtr(nullptr);

meshPtr = autoPtr<fvMesh>
        (
            new fvMesh
            (
                IOobject
                (
                    fvMesh::defaultRegion,
                    runTime.timeName(),
                    runTime,
                    IOobject::MUST_READ
                )
            )
        );

Foam::fvMesh& mesh = meshPtr();

    pisoControl piso(mesh);

    #include "createFields.H"
    
    Info << p << endl;
    
    Eigen::MatrixXd p_eigen;
    p_eigen = Foam2Eigen::field2Eigen(p);
    std::cerr << p_eigen << std::endl;
    exit(0);
    
    Info << U << endl;
    Info << R << endl;


        // Info << *(&U.ref()[0][0]) << endl;
        // Info << *(&U.ref()[0][0]+1) << endl;
        // Info << *(&U.ref()[0][0]+2) << endl;
        // exit(0);
       

   
        // Info << R[0] << endl;

        

        // Info << U << endl;
        // Info << *(&U.ref()[0][0]) << endl;
        // Info << *(&U.ref()[0][0]+1) << endl;
        // Info << *(&U.ref()[0][0]+2) << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
