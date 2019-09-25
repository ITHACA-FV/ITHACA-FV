Info << "Reading transport properties" << endl;
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

_betaTE = autoPtr<dimensionedScalar>
          (
              new dimensionedScalar
              (
                  transportProperties.lookup("beta")
              )
          );
dimensionedScalar& betaTE = _betaTE();

_nu = autoPtr<dimensionedScalar>
      (
          new dimensionedScalar
          (
              transportProperties.lookup("nu")
          )
      );
dimensionedScalar& nu = _nu();

_Pr = autoPtr<dimensionedScalar>
      (
          new dimensionedScalar
          (
              transportProperties.lookup("Pr")
          )
      );
dimensionedScalar& Pr = _Pr();

_Prt = autoPtr<dimensionedScalar>
       (
           new dimensionedScalar
           (
               transportProperties.lookup("Prt")
           )
       );
dimensionedScalar& Prt = _Prt();


_rhoRef = autoPtr<dimensionedScalar>
          (
              new dimensionedScalar
              (
                  transportProperties.lookup("rhoRef")
              )
          );
dimensionedScalar& rhoRef = _rhoRef();

_CpRef = autoPtr<dimensionedScalar>
         (
             new dimensionedScalar
             (
                 transportProperties.lookup("CpRef")
             )
         );
dimensionedScalar& CpRef = _CpRef();

_Tref = autoPtr<dimensionedScalar>
        (
            new dimensionedScalar
            (
                transportProperties.lookup("Tref")
            )
        );
dimensionedScalar& Tref = _Tref();

IOdictionary precursorTransportProperties
(
    IOobject
    (
        "precursorTransportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

_Sc = autoPtr<dimensionedScalar>
      (
          new dimensionedScalar
          (
              precursorTransportProperties.lookup("Sc")
          )
      );
dimensionedScalar& Sc = _Sc();

_Sct = autoPtr<dimensionedScalar>
       (
           new dimensionedScalar
           (
               precursorTransportProperties.lookup("Sct")
           )
       );
dimensionedScalar& Sct = _Sct();
