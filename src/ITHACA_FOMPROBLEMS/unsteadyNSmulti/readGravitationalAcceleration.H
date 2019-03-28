Info<< "\nReading g" << endl;

_g = autoPtr<uniformDimensionedVectorField>
(
    new uniformDimensionedVectorField 
    (
        IOobject
    (
        "g",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
    )
);
uniformDimensionedVectorField& g = _g();
