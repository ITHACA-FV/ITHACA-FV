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

_dt = autoPtr<dimensionedScalar>
      (
          new dimensionedScalar
          (
              "dt",
              dimensionSet(0, 0, 1, 0, 0),
              1.0
          )
      );

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
               mesh
           )
       );
surfaceScalarField phi = _phi();

surfaceScalarField phi0(phi);

_phi0 = autoPtr<surfaceScalarField>
        (
            new surfaceScalarField(phi0)
        );

pRefCell = 0;
pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("SIMPLE"), pRefCell, pRefValue);

_laminarTransport = autoPtr<singlePhaseTransportModel>
                    (
                        new singlePhaseTransportModel( U, phi )
                    );
singlePhaseTransportModel& laminarTransport = _laminarTransport();

turbulence = autoPtr<incompressible::turbulenceModel>
             (
                 incompressible::turbulenceModel::New(U, phi, laminarTransport)
             );

_MRF = autoPtr<IOMRFZoneList>
       (
           new IOMRFZoneList(mesh)
       );
