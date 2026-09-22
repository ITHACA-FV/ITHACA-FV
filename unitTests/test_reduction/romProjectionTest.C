#include "fvCFD.H"
#include "PstreamReduceOps.H"
#include "Modes.H"

#include <iomanip>
#include <string>

using namespace Foam;


// -----------------------------------------------------------------------------
// Helpers
// -----------------------------------------------------------------------------

scalar relativeL2Error
(
    const volScalarField& reference,
    const volScalarField& approximation
)
{
    scalar err2 = 0.0;
    scalar ref2 = 0.0;

    const scalarField& volumes = reference.mesh().V();

    forAll(reference, celli)
    {
        const scalar diff =
            reference[celli]
          - approximation[celli];

        err2 +=
            volumes[celli]
           *sqr(diff);

        ref2 +=
            volumes[celli]
           *sqr(reference[celli]);
    }

    reduce(err2, sumOp<scalar>());
    reduce(ref2, sumOp<scalar>());

    return Foam::sqrt(err2/ref2);
}


scalar relativeL2Error
(
    const volVectorField& reference,
    const volVectorField& approximation
)
{
    scalar err2 = 0.0;
    scalar ref2 = 0.0;

    const scalarField& volumes = reference.mesh().V();

    forAll(reference, celli)
    {
        err2 +=
            volumes[celli]
           *magSqr
            (
                reference[celli]
              - approximation[celli]
            );

        ref2 +=
            volumes[celli]
           *magSqr(reference[celli]);
    }

    reduce(err2, sumOp<scalar>());
    reduce(ref2, sumOp<scalar>());

    return Foam::sqrt(err2/ref2);
}


void printTestResult
(
    const word& fieldType,
    const word& projectionType,
    const Eigen::MatrixXd& reducedMatrix,
    const Eigen::MatrixXd& reducedRhs,
    const Eigen::MatrixXd& coefficient,
    const scalar relativeError,
    const scalar tolerance,
    label& failures
)
{
    const scalar coeffError =
        Foam::mag(coefficient(0,0) - 1.0);

    const bool passed =
        coeffError < tolerance
     && relativeError < tolerance;

    if (!passed)
    {
        ++failures;
    }

    if (Pstream::master())
    {
        std::cout
            << "\n----------------------------------------\n"
            << " " << fieldType
            << " / " << projectionType
            << "\n----------------------------------------\n"
            << "Ar =\n"
            << reducedMatrix
            << "\n\nbr =\n"
            << reducedRhs
            << "\n\na =\n"
            << coefficient
            << "\n\n|a - 1| = "
            << coeffError
            << "\nRelative L2 error = "
            << relativeError
            << "\nResult: "
            << (passed ? "PASS" : "FAIL")
            << "\n"
            << std::endl;
    }
}


// -----------------------------------------------------------------------------
// Main
// -----------------------------------------------------------------------------

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    Foam::Time runTime
    (
        Foam::Time::controlDictName,
        args
    );

    Foam::fvMesh mesh
    (
        Foam::IOobject
        (
            Foam::fvMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        )
    );

    ITHACAparameters* para =
        ITHACAparameters::getInstance(mesh, runTime);

    (void)para;


    if (Pstream::master())
    {
        std::cout
            << std::setprecision(16)
            << "\n========================================\n"
            << " ITHACA-FV fvMatrix projection test\n"
            << " MPI ranks = "
            << Pstream::nProcs()
            << "\n========================================\n"
            << std::endl;
    }


    // -------------------------------------------------------------------------
    // Global geometry information.
    //
    // Coordinates are used to define the source so that the same physical
    // problem is obtained independently of the domain decomposition.
    // -------------------------------------------------------------------------

    scalar xmin = GREAT;
    scalar xmax = -GREAT;

    scalar ymin = GREAT;
    scalar ymax = -GREAT;

    scalar zmin = GREAT;
    scalar zmax = -GREAT;

    forAll(mesh.C(), celli)
    {
        const vector& C =
            mesh.C()[celli];

        xmin = min(xmin, C.x());
        xmax = max(xmax, C.x());

        ymin = min(ymin, C.y());
        ymax = max(ymax, C.y());

        zmin = min(zmin, C.z());
        zmax = max(zmax, C.z());
    }

    reduce(xmin, minOp<scalar>());
    reduce(xmax, maxOp<scalar>());

    reduce(ymin, minOp<scalar>());
    reduce(ymax, maxOp<scalar>());

    reduce(zmin, minOp<scalar>());
    reduce(zmax, maxOp<scalar>());

    const scalar Lx =
        max(xmax - xmin, SMALL);

    const scalar Ly =
        max(ymax - ymin, SMALL);

    const scalar Lz =
        max(zmax - zmin, SMALL);


    // A deliberately strict tolerance. The exact-space tests should normally
    // be close to machine precision in both serial and parallel.
    const scalar tolerance = 1e-10;

    label failures = 0;


    // =========================================================================
    //
    //                         SCALAR EXACT-SPACE TEST
    //
    // =========================================================================

    if (Pstream::master())
    {
        std::cout
            << "\n========================================\n"
            << " SCALAR EXACT-SPACE TEST\n"
            << "========================================\n"
            << std::endl;
    }


    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );


    // -------------------------------------------------------------------------
    // Solve scalar FOM:
    //
    //                    A T = b
    // -------------------------------------------------------------------------

    fvScalarMatrix scalarAFOM
    (
        -fvm::laplacian(T)
    );

    forAll(scalarAFOM.source(), celli)
    {
        const vector& C =
            mesh.C()[celli];

        const scalar x =
            (C.x() - xmin)/Lx;

        const scalar y =
            (C.y() - ymin)/Ly;

        const scalar z =
            (C.z() - zmin)/Lz;

        scalarAFOM.source()[celli] =
              1.0
            + 0.5*x
            + 0.25*y
            + 0.125*z;
    }

    scalarAFOM.solve();


    volScalarField TFOM
    (
        IOobject
        (
            "TFOM",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        T
    );


    // -------------------------------------------------------------------------
    // Exact reduced space:
    //
    //                    phi_0 = TFOM
    //
    // Hence both Galerkin and Petrov-Galerkin must return a = 1.
    // -------------------------------------------------------------------------

    PtrList<volScalarField> scalarModes(1);

    scalarModes.set
    (
        0,
        new volScalarField
        (
            IOobject
            (
                "scalarMode_0",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            TFOM
        )
    );

    volScalarModes scalarITHACAModes;
    scalarITHACAModes = scalarModes;


    // Reassemble exactly the same discrete system.
    fvScalarMatrix scalarA
    (
        -fvm::laplacian(T)
    );

    forAll(scalarA.source(), celli)
    {
        const vector& C =
            mesh.C()[celli];

        const scalar x =
            (C.x() - xmin)/Lx;

        const scalar y =
            (C.y() - ymin)/Ly;

        const scalar z =
            (C.z() - zmin)/Lz;

        scalarA.source()[celli] =
              1.0
            + 0.5*x
            + 0.25*y
            + 0.125*z;
    }


    // -------------------------------------------------------------------------
    // Scalar Galerkin
    // -------------------------------------------------------------------------

    List<Eigen::MatrixXd> scalarG =
        scalarITHACAModes.project
        (
            scalarA,
            1,
            "G"
        );

    Eigen::MatrixXd aScalarG =
        scalarG[0]
       .fullPivLu()
       .solve(scalarG[1]);


    volScalarField TScalarG
    (
        IOobject
        (
            "TScalarG",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        TFOM
    );

    scalarITHACAModes.reconstruct
    (
        TScalarG,
        aScalarG,
        "TScalarG"
    );

    const scalar scalarGError =
        relativeL2Error
        (
            TFOM,
            TScalarG
        );

    printTestResult
    (
        "SCALAR",
        "G",
        scalarG[0],
        scalarG[1],
        aScalarG,
        scalarGError,
        tolerance,
        failures
    );


    // -------------------------------------------------------------------------
    // Scalar Petrov-Galerkin
    // -------------------------------------------------------------------------

    List<Eigen::MatrixXd> scalarPG =
        scalarITHACAModes.project
        (
            scalarA,
            1,
            "PG"
        );

    Eigen::MatrixXd aScalarPG =
        scalarPG[0]
       .fullPivLu()
       .solve(scalarPG[1]);


    volScalarField TScalarPG
    (
        IOobject
        (
            "TScalarPG",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        TFOM
    );

    scalarITHACAModes.reconstruct
    (
        TScalarPG,
        aScalarPG,
        "TScalarPG"
    );

    const scalar scalarPGError =
        relativeL2Error
        (
            TFOM,
            TScalarPG
        );

    printTestResult
    (
        "SCALAR",
        "PG",
        scalarPG[0],
        scalarPG[1],
        aScalarPG,
        scalarPGError,
        tolerance,
        failures
    );


    // =========================================================================
    //
    //                         VECTOR EXACT-SPACE TEST
    //
    // =========================================================================

    if (Pstream::master())
    {
        std::cout
            << "\n========================================\n"
            << " VECTOR EXACT-SPACE TEST\n"
            << "========================================\n"
            << std::endl;
    }


    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );


    // -------------------------------------------------------------------------
    // Solve vector FOM:
    //
    //                    A U = b
    // -------------------------------------------------------------------------

    fvVectorMatrix vectorAFOM
    (
        -fvm::laplacian(U)
    );

    forAll(vectorAFOM.source(), celli)
    {
        const vector& C =
            mesh.C()[celli];

        const scalar x =
            (C.x() - xmin)/Lx;

        const scalar y =
            (C.y() - ymin)/Ly;

        const scalar z =
            (C.z() - zmin)/Lz;

        vectorAFOM.source()[celli] =
            vector
            (
                1.0 + x + 0.10*y,
                2.0 + y + 0.15*z,
                3.0 + z + 0.20*x
            );
    }

    vectorAFOM.solve();


    volVectorField UFOM
    (
        IOobject
        (
            "UFOM",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U
    );


    // -------------------------------------------------------------------------
    // Exact reduced space:
    //
    //                    phi_0 = UFOM
    //
    // Hence both Galerkin and Petrov-Galerkin must return a = 1.
    // -------------------------------------------------------------------------

    PtrList<volVectorField> vectorModes(1);

    vectorModes.set
    (
        0,
        new volVectorField
        (
            IOobject
            (
                "vectorMode_0",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            UFOM
        )
    );

    volVectorModes vectorITHACAModes;
    vectorITHACAModes = vectorModes;


    // Reassemble exactly the same vector system.
    fvVectorMatrix vectorA
    (
        -fvm::laplacian(U)
    );

    forAll(vectorA.source(), celli)
    {
        const vector& C =
            mesh.C()[celli];

        const scalar x =
            (C.x() - xmin)/Lx;

        const scalar y =
            (C.y() - ymin)/Ly;

        const scalar z =
            (C.z() - zmin)/Lz;

        vectorA.source()[celli] =
            vector
            (
                1.0 + x + 0.10*y,
                2.0 + y + 0.15*z,
                3.0 + z + 0.20*x
            );
    }


    // -------------------------------------------------------------------------
    // Vector Galerkin
    // -------------------------------------------------------------------------

    List<Eigen::MatrixXd> vectorG =
        vectorITHACAModes.project
        (
            vectorA,
            1,
            "G"
        );

    Eigen::MatrixXd aVectorG =
        vectorG[0]
       .fullPivLu()
       .solve(vectorG[1]);


    volVectorField UVectorG
    (
        IOobject
        (
            "UVectorG",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        UFOM
    );

    vectorITHACAModes.reconstruct
    (
        UVectorG,
        aVectorG,
        "UVectorG"
    );

    const scalar vectorGError =
        relativeL2Error
        (
            UFOM,
            UVectorG
        );

    printTestResult
    (
        "VECTOR",
        "G",
        vectorG[0],
        vectorG[1],
        aVectorG,
        vectorGError,
        tolerance,
        failures
    );


    // -------------------------------------------------------------------------
    // Vector Petrov-Galerkin
    // -------------------------------------------------------------------------

    List<Eigen::MatrixXd> vectorPG =
        vectorITHACAModes.project
        (
            vectorA,
            1,
            "PG"
        );

    Eigen::MatrixXd aVectorPG =
        vectorPG[0]
       .fullPivLu()
       .solve(vectorPG[1]);


    volVectorField UVectorPG
    (
        IOobject
        (
            "UVectorPG",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        UFOM
    );

    vectorITHACAModes.reconstruct
    (
        UVectorPG,
        aVectorPG,
        "UVectorPG"
    );

    const scalar vectorPGError =
        relativeL2Error
        (
            UFOM,
            UVectorPG
        );

    printTestResult
    (
        "VECTOR",
        "PG",
        vectorPG[0],
        vectorPG[1],
        aVectorPG,
        vectorPGError,
        tolerance,
        failures
    );


    // =========================================================================
    //
    //                              SUMMARY
    //
    // =========================================================================

    reduce(failures, maxOp<label>());

    if (Pstream::master())
    {
        std::cout
            << "\n========================================\n"
            << " TEST SUMMARY\n"
            << "========================================\n"
            << "Execution mode : "
            << (Pstream::parRun() ? "PARALLEL" : "SERIAL")
            << "\nMPI ranks      : "
            << Pstream::nProcs()
            << "\nTolerance      : "
            << tolerance
            << "\nFailures       : "
            << failures
            << "\nOverall result : "
            << (failures == 0 ? "PASS" : "FAIL")
            << "\n========================================\n"
            << std::endl;
    }


    return failures == 0 ? 0 : 1;
}
