/*---------------------------------------------------------------------------*\

Global
    createUf

Description
    Creates and initialises the velocity velocity field Uf.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//from /scratch/mathlab/packages/source/OpenFOAM/OpenFOAM-v2106/src/finiteVolume/cfdTools/incompressible
/*#ifndef createUf_H
#define createUf_H

Info<< "Reading/calculating face velocity Uf\n" << endl;

_Uf = autoPtr<surfaceVectorField>
     (
         new surfaceVectorField
         (
             IOobject
             (
                 "Uf",
                 runTime.timeName(),
                 mesh,
                 IOobject::MUST_READ,
                 IOobject::AUTO_WRITE
             ),
             fvc::interpolate(U)
         )
     );
surfaceVectorField& Uf = _Uf();

surfaceVectorField U0(Uf);

_U0 = autoPtr<surfaceVectorField>
      (
          new surfaceVectorField(U0)
      );

#endif */




/*Info << "Reading field p\n" << endl;
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



*/
surfaceVectorField Uf
(
    IOobject
    (
        "Uf",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::interpolate(U)
);




// ************************************************************************* //
