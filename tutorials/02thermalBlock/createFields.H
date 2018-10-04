volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
volScalarField yPos = mesh.C().component(vector::Y);
volScalarField xPos = mesh.C().component(vector::X);

volScalarField S
(
    IOobject
    (
        "S",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    T.mesh(),
    dimensionedScalar("zero", dimensionSet(0, 0, -1, 1, 0, 0, 0), 0)
);



Info << "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);


dimensionedScalar DT
(
    "DT",
    dimensionSet( 0, 2, -1, 0, 0, 0, 0),
    scalar(1)
);
