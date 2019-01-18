/*---------------------------------------------------------------------------*\
     ██╗████████╗██╗  ██╗ █████╗  ██████╗ █████╗       ███████╗██╗   ██╗
     ██║╚══██╔══╝██║  ██║██╔══██╗██╔════╝██╔══██╗      ██╔════╝██║   ██║
     ██║   ██║   ███████║███████║██║     ███████║█████╗█████╗  ██║   ██║
     ██║   ██║   ██╔══██║██╔══██║██║     ██╔══██║╚════╝██╔══╝  ╚██╗ ██╔╝
     ██║   ██║   ██║  ██║██║  ██║╚██████╗██║  ██║      ██║      ╚████╔╝
     ╚═╝   ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝      ╚═╝       ╚═══╝

 * In real Time Highly Advanced Computational Applications for Finite Volumes
 * Copyright (C) 2017 by the ITHACA-FV authors
-------------------------------------------------------------------------------

  License
  This file is part of ITHACA-FV

  ITHACA-FV is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ITHACA-FV is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with ITHACA-FV. If not, see <http://www.gnu.org/licenses/>.

 \*---------------------------------------------------------------------------*/
Info << "Reading thermophysical properties\n" << endl;

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

// Laminar viscocity [m2/s]
_nu = autoPtr<dimensionedScalar>
      (
          new dimensionedScalar
          (
              "nu",
              dimViscosity,
              transportProperties.lookup("nu")
          )
      );
dimensionedScalar& nu = _nu();

// Turbulent viscocity [m2/s]
_nut = autoPtr<volScalarField>
       (
           new volScalarField
           (
               IOobject
               (
                   "nut",
                   runTime.timeName(),
                   mesh,
                   IOobject::MUST_READ,
                   IOobject::AUTO_WRITE
               ),
               mesh
           )
       );
volScalarField& nut = _nut();

_UliftBC = autoPtr<volScalarField>
           (
               new volScalarField
               (
                   IOobject
                   (
                       "UliftBC",
                       runTime.timeName(),
                       mesh,
                       IOobject::MUST_READ,
                       IOobject::AUTO_WRITE
                   ),
                   mesh
               )
           );
volScalarField& UliftBC = _UliftBC();

// Laminar Prandtl number
_Pr = autoPtr<dimensionedScalar>
      (
          new dimensionedScalar
          (
              "Pr",
              dimless,
              transportProperties.lookup("Pr")
          )
      );
dimensionedScalar& Pr = _Pr();

// Turbulent Prandtl number
_Prt = autoPtr<dimensionedScalar>
       (
           new dimensionedScalar
           (
               "Prt",
               dimless,
               transportProperties.lookup("Prt")
           )
       );
dimensionedScalar& Prt = _Prt();

// Thermal expansion coefficient [1/K]
_beta = autoPtr<dimensionedScalar>
        (
            new dimensionedScalar
            (
                "beta",
                dimless / dimTemperature,
                transportProperties.lookup("beta")
            )
        );
dimensionedScalar& beta = _beta();

// Reference temperature [K]
_TRef = autoPtr<dimensionedScalar>
        (
            new dimensionedScalar
            (
                "TRef",
                dimTemperature,
                transportProperties.lookup("TRef")
            )
        );
dimensionedScalar& TRef = _TRef();

// Gravitational constant
Info << "\nReading g" << endl;
_g = autoPtr<dimensionedVector>
     (
         new dimensionedVector
         (
             "g",
             dimensionSet(0, 1, -2, 0, 0),
             vector(0, -9.81, 0)
         )
     );
dimensionedVector& g = _g();

Info << "\nReading hRef" << endl;
_hRef = autoPtr<dimensionedScalar>
        (
            new dimensionedScalar
            (
                "hRef",
                dimLength,
                0
            )
        );
dimensionedScalar& hRef = _hRef();

Info << "Calculating field ghRef\n" << endl;
_ghRef = autoPtr<dimensionedScalar>
         (
             new dimensionedScalar
             (
                 "ghRef",
                 dimAcceleration* dimLength,
                 0
             )
         );
dimensionedScalar& ghRef = _ghRef();

Info << "Reading field gh\n" << endl;
_gh = autoPtr<volScalarField>
      (
          new volScalarField
          (
              IOobject
              (
                  "gh",
                  runTime.timeName(),
                  mesh
              ),
              (g& mesh.C()) - ghRef
          )
      );
volScalarField& gh = _gh();

Info << "Reading field ghf\n" << endl;
_ghf = autoPtr<surfaceScalarField>
       (
           new surfaceScalarField
           (
               IOobject
               (
                   "ghf",
                   runTime.timeName(),
                   mesh
               ),
               (g& mesh.Cf()) - ghRef)
       );
surfaceScalarField& ghf = _ghf();

// Reading fields
Info << "Reading field T\n" << endl;
_T = autoPtr<volScalarField>
     (
         new volScalarField
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
         )
     );
volScalarField& T = _T();

Info << "Reading field p_rgh\n" << endl;
_p_rgh = autoPtr<volScalarField>
         (
             new volScalarField
             (
                 IOobject
                 (
                     "p_rgh",
                     runTime.timeName(),
                     mesh,
                     IOobject::MUST_READ,
                     IOobject::AUTO_WRITE
                 ),
                 mesh
             )
         );
volScalarField& p_rgh = _p_rgh();

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


#include "createPhi.H"

_laminarTransport = autoPtr<singlePhaseTransportModel>
                    (
                        new singlePhaseTransportModel( U, phi )
                    );
singlePhaseTransportModel& laminarTransport = _laminarTransport();

Info << "Creating turbulence model\n" << endl;
turbulence = autoPtr<incompressible::turbulenceModel>
             (
                 incompressible::turbulenceModel::New(U, phi, laminarTransport)
             );

// Kinematic density for buoyancy force
_rhok = autoPtr<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "rhok",
                    runTime.timeName(),
                    mesh
                ),
                1.0 - beta * (T - TRef)
            )
        );
volScalarField& rhok = _rhok();

// kinematic turbulent thermal conductivity m2/s
Info << "Reading field alphat\n" << endl;
_alphat = autoPtr<volScalarField>
          (
              new volScalarField
              (
                  IOobject
                  (
                      "alphat",
                      runTime.timeName(),
                      mesh,
                      IOobject::MUST_READ,
                      IOobject::AUTO_WRITE
                  ),
                  mesh
              )
          );
volScalarField& alphat = _alphat();

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
                 IOobject::NO_READ,
                 IOobject::AUTO_WRITE
             ),
             p_rgh + rhok* gh
         )
     );
volScalarField& p = _p();


pRefCell = 0;
pRefValue = 0.0;

pimpleControl& pimple = _pimple();

setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
         (
             "p",
             p.dimensions(),
             pRefValue - getRefCellValue(p, pRefCell)
         );
}

mesh.setFluxRequired(p_rgh.name());

_MRF = autoPtr<IOMRFZoneList>
       (
           new IOMRFZoneList(mesh)
       );
