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
    
    ////////// SCALAR FIELD TEST //////////

    Eigen::VectorXd p_eigen;
    p_eigen = Foam2Eigen::field2Eigen(p);

    Info << "\nOriginal p:\n" << endl;
    Info << p << endl;

    Info << "\nFirst two elements of p, by reference:\n" << endl;
    std::cerr << *(&p.ref()[0]) << std::endl;
    std::cerr << *(&p.ref()[0]+1) << std::endl;

    Info << "\nFirst two elements of p_eigen:\n" << endl;
    std::cerr << p_eigen(0) << std::endl;
    std::cerr << p_eigen(1) << std::endl;

    // Change p_eigen
    p_eigen(0) = 0.0;
    p = Foam2Eigen::Eigen2field(p, p_eigen, false);

    Info << "\nConverting back to field, after setting R_eigen(9) = 0.0 (First row, first column, second element):\n" << endl;
    Info << *(&p.ref()) << endl;

    // uncomment lines below to test error handling for wrong size
    // Info << "\nConverting back to field with a resized p (3×3). Eigen2Field should throw an error:\n" << endl;
    // Eigen::MatrixXd p_eigen_r = Eigen::Map<Eigen::MatrixXd>((p_eigen.data()),3,3);
    // Foam2Eigen::Eigen2field(pp, p_eigen_r, false);

    // Testing for BCs
    List<Eigen::VectorXd> p_bc = Foam2Eigen::field2EigenBC(p);
    Info << "\np_bc (Eigen Map):\n" << endl;
    Info << p_bc << endl;
    // exit(0);

    ////////// VECTOR FIELD TEST //////////

    Info << "\nOriginal U:\n" << endl;
    Info << U << endl;
    Eigen::VectorXd U_eigen;
    U_eigen = Foam2Eigen::field2Eigen(U);

    std::cerr << "\nFirst 4 elements of U, by reference:\n" << std::endl;
    Info << *(&U.ref()[0][0]) << endl;
    Info << *(&U.ref()[0][0]+1) << endl;
    Info << *(&U.ref()[0][0]+2) << endl;
    Info << *(&U.ref()[0][0]+3) << endl;

    std::cerr << "\nFirst 4 elements of U_eigen:\n" << std::endl;
    std::cerr << U_eigen(0) << std::endl;
    std::cerr << U_eigen(1) << std::endl;
    std::cerr << U_eigen(2) << std::endl;
    std::cerr << U_eigen(3) << std::endl;
    
    // Change U_eigen
    U_eigen(3) = 0.0;
    U = Foam2Eigen::Eigen2field(U, U_eigen, false);
    std::cerr << "\nConverting back to field:, after setting U_eigen(3) = 0.0 (second row, first column)\n" << std::endl;
    Info << U << endl;

    // uncomment lines below to test error handling for wrong size
    // Info << "\nConverting back to field with a resized U (9×3). Eigen2Field should throw an error:\n" << endl;
    // Eigen::MatrixXd U_eigen_r = Eigen::Map<Eigen::MatrixXd>((U_eigen.data()),3 ,9);
    // U_eigen_r.transposeInPlace();
    // std::cerr << U_eigen_r << std::endl;
    // volVectorField UU(U);
    // Foam2Eigen::Eigen2field(UU, U_eigen_r, false);

    // Testing for BCs
    List<Eigen::VectorXd> U_bc;
    U_bc = Foam2Eigen::field2EigenBC(U);
    Info << "\nU_bc (Eigen Map):\n" << endl;
    Info << U_bc << endl;
    // exit(0);

    ////////// TENSOR FIELD TEST //////////

    Eigen::VectorXd R_eigen;
    R_eigen = Foam2Eigen::field2Eigen(R);

    Info << "\nOriginal R\n" << endl;
    Info << *(&R.ref()) << endl;
    Info << *(R.ref().data()) << endl;

    Info << "\nFirst row of first element:\n" << endl;
    Info << *(&R.ref()[0][0]) << endl;
    Info << *(&R.ref()[0][0]+1) << endl;
    Info << *(&R.ref()[0][0]+2) << endl;

    Info << "\nSecond row of first element:\n" << endl;
    Info << *(&R.ref()[0][0]+3) << endl;
    Info << "\nFirst row of second element:\n" << endl;
    Info << *(&R.ref()[0][0]+9) << endl;

    Info << "\nR_eigen\n" << endl;
    Info << "\nFirst three elements:\n" << endl;
    Info << R_eigen(0) << endl;
    Info << R_eigen(1) << endl;
    Info << R_eigen(2) << endl;
    Info << "\n4th element:\n" << endl;
    Info << R_eigen(3) << endl;
    Info << "\n9th element:\n" << endl;
    Info << R_eigen(9) << endl;

    // Change R_eigen
    R_eigen(9) = 0.0;
    std::cerr << "\nConverting back to field, after setting R_eigen(9) = 0.0 (First row, first column, second element):\n" << std::endl;
    R = Foam2Eigen::Eigen2field(R, R_eigen, false);
    Info << R << endl;
    // exit(0);
    
    // uncomment lines below to test error handling for wrong size
    // Eigen::MatrixXd R_eigen_r = Eigen::Map<Eigen::MatrixXd>((R_eigen.data()),9 ,9);
    // R_eigen_r.transposeInPlace();
    // Info << "\nResized R_eigen:\n" << endl;
    // std::cerr << R_eigen_r << std::endl;
    // volTensorField RR(R);
    // Info << "\nConverting back to field with a resized U (9×9). Eigen2Field should throw an error:\n" << endl;
    // Foam2Eigen::Eigen2field(RR, R_eigen_r, false);

    // Testing for BCs
    List<Eigen::VectorXd> R_bc;
    R_bc = Foam2Eigen::field2EigenBC(R);
    Info << "\nR_bc (Eigen Map):\n" << endl;
    Info << R_bc << endl;

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
