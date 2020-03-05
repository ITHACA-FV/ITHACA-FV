Info << "Reading transportProperties\n" << endl;

_transportProperties = autoPtr<IOdictionary>
                       (
                           new  IOdictionary
                           (
                               IOobject
                               (
                                   "transportProperties",
                                   runTime.constant(),
                                   mesh,
                                   IOobject::MUST_READ_IF_MODIFIED,
                                   IOobject::NO_WRITE
                               )
                           )
                       );

IOdictionary& transportProperties = _transportProperties();

_nu = autoPtr<dimensionedScalar>
      (
          new dimensionedScalar
          (
              "nu",
              dimViscosity,
              transportProperties
          )
      );

Info << "Reading field p\n" << endl;
_p = autoPtr<volScalarField>
     (
         new volScalarField
         (
             IOobject
             (
                 "p",
                 runTime.timeName(),
                 mesh,
                 IOobject::MUST_READ,
                 IOobject::AUTO_WRITE
             ),
             mesh
         )
     );
volScalarField& p = _p();

volScalarField p0(p);

_p0 = autoPtr<volScalarField>
      (
          new volScalarField(p0)
      );


Info << "Reading field U\n" << endl;
_U = autoPtr<volVectorField>
     (
         new volVectorField
         (
             IOobject
             (
                 "U",
                 runTime.timeName(),
                 mesh,
                 IOobject::MUST_READ,
                 IOobject::AUTO_WRITE
             ),
             mesh
         )
     );
volVectorField& U = _U();

volVectorField U0(U);

_U0 = autoPtr<volVectorField>
      (
          new volVectorField(U0)
      );


Info << "Reading/calculating face flux field phi\n" << endl;

_phi = autoPtr<surfaceScalarField>
       (
           new surfaceScalarField
           (
               IOobject
               (
                   "phi",
                   runTime.timeName(),
                   mesh,
                   IOobject::READ_IF_PRESENT,
                   IOobject::AUTO_WRITE
               ),
               linearInterpolate(U) & mesh.Sf()
           )
       );
surfaceScalarField& phi = _phi();
surfaceScalarField phi0(phi);

_phi0 = autoPtr<surfaceScalarField>
        (
            new surfaceScalarField(phi0)
        );


pRefCell = 0;
pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());