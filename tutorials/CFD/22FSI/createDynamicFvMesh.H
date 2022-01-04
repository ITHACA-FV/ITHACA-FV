Foam::word regionName = Foam::polyMesh::defaultRegion;

if (args.readIfPresent("region", regionName))
{
Foam::Info
    << "Create dynamic mesh " << regionName << " for time = "
    << runTime.timeName() << Foam::nl << Foam::endl;
}
else
{
Foam::Info
    << "Create mesh for time = "
    << runTime.timeName() << Foam::nl << Foam::endl;
}

autoPtr<dynamicFvMesh> meshPtr
(
    dynamicFvMesh::New
    (
        IOobject
        (
            regionName,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    )
);

dynamicFvMesh& mesh = meshPtr();

//autoPtr<dynamicFvMesh> meshPtr(dynamicFvMesh::New(args, runTime));

//dynamicFvMesh& mesh = meshPtr();
