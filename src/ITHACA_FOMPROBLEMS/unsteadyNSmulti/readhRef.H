Info<< "\nReading hRef" << endl;
   
_hRef = autoPtr<uniformDimensionedScalarField>
(
    new uniformDimensionedScalarField
    (
        IOobject
        (
            "hRef",
            runTime.constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        dimensionedScalar("hRef", dimLength, 0)
    )
);

uniformDimensionedScalarField& hRef = _hRef();
