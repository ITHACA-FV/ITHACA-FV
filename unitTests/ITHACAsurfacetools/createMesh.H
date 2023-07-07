/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Description
    Create a fvMesh (specified region or defaultRegion) with
    additional handling of -dry-run and -dry-run-write options.

Required Variables
    - args [argList]
    - runTime [Time]

Provided Variables
    - regionName [word]
    - mesh [fvMesh]
    - meshPtr [autoPtr<fvMesh>]

\*---------------------------------------------------------------------------*/

Foam::autoPtr<Foam::fvMesh> meshPtr(nullptr);
Foam::word regionName(Foam::polyMesh::defaultRegion);

if (args.found("dry-run") || args.found("dry-run-write"))
{
    Foam::Info
        << "Operating in 'dry-run' mode: case will run for 1 time step.  "
        << "All checks assumed OK on a clean exit" << Foam::endl;

    // Allow region in combination with dry-run
    args.readIfPresent("region", regionName);

    Foam::FieldBase::allowConstructFromLargerSize = true;

    // Create a simplified 1D mesh and attempt to re-create boundary conditions
    meshPtr.reset
    (
        new Foam::simplifiedMeshes::columnFvMesh(runTime, regionName)
    );

    // Stop after 1 iteration of the simplified mesh

    if (args.found("dry-run-write"))
    {
        // Using saWriteNow triggers function objects execute(), write()
        runTime.stopAt(Foam::Time::saWriteNow);

        // Make sure mesh gets output to the current time (since instance
        // no longer constant)
        meshPtr().setInstance(runTime.timeName());
    }
    else
    {
        // Using saNoWriteNow triggers function objects execute(),
        // but not write()
        runTime.stopAt(Foam::Time::saNoWriteNow);
    }

    Foam::functionObject::outputPrefix = "postProcessing-dry-run";
}
else
{
    if (args.readIfPresent("region", regionName))
    {
        Foam::Info << "Create mesh " << regionName;
    }
    else
    {
        Foam::Info << "Create mesh";
    }
    Foam::Info << " for time = " << runTime.timeName() << Foam::nl;

    meshPtr.reset
    (
        new Foam::fvMesh
        (
            Foam::IOobject
            (
                regionName,
                runTime.timeName(),
                runTime,
                Foam::IOobject::MUST_READ
            ),
            false
        )
    );
    meshPtr().init(true);   // initialise all (lower levels and current)

    Foam::Info << Foam::endl;
}

Foam::fvMesh& mesh = meshPtr();


// ************************************************************************* //
