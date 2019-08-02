Info << "Reading decay heat properties" << endl;
IOdictionary decayheatProperties
(
    IOobject
    (
        "decayheatProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

_decLam1 = autoPtr<dimensionedScalar>
           (
               new dimensionedScalar
               (
                   decayheatProperties.lookup("decLam1")
               )
           );
dimensionedScalar& decLam1 = _decLam1();

_decLam2 = autoPtr<dimensionedScalar>
           (
               new dimensionedScalar
               (
                   decayheatProperties.lookup("decLam2")
               )
           );
dimensionedScalar& decLam2 = _decLam2();

_decLam3 = autoPtr<dimensionedScalar>
           (
               new dimensionedScalar
               (
                   decayheatProperties.lookup("decLam3")
               )
           );
dimensionedScalar& decLam3 = _decLam3();

_decBeta1 = autoPtr<dimensionedScalar>
            (
                new dimensionedScalar
                (
                    decayheatProperties.lookup("decBeta1")
                )
            );
dimensionedScalar& decBeta1 = _decBeta1();

_decBeta2 = autoPtr<dimensionedScalar>
            (
                new dimensionedScalar
                (
                    decayheatProperties.lookup("decBeta2")
                )
            );
dimensionedScalar& decBeta2 = _decBeta2();

_decBeta3 = autoPtr<dimensionedScalar>
            (
                new dimensionedScalar
                (
                    decayheatProperties.lookup("decBeta3")
                )
            );
dimensionedScalar& decBeta3 = _decBeta3();

