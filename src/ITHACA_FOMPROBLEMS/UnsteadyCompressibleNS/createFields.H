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

Info<< "Reading thermophysical properties\n" << endl;

_pThermo = autoPtr<psiThermo>
(
    psiThermo::New(mesh)
);
psiThermo& thermo = _pThermo();

volScalarField& e = thermo.he();

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

_rho = autoPtr<volScalarField>
     (
         new volScalarField
         (
             IOobject
             (
                 "rho",
                 runTime.timeName(),
                 mesh,
                 IOobject::NO_READ,
                 IOobject::AUTO_WRITE
             ),
             thermo.rho()
         )
     );

volScalarField& rho = _rho();

_rhoU = autoPtr<volVectorField>
      (
          new volVectorField
          (
              IOobject
              (
                  "rhoU",
                  runTime.timeName(),
                  mesh,
                  IOobject::NO_READ,
                  IOobject::NO_WRITE
              ),
              rho*U
          )
      );

volVectorField& rhoU = _rhoU();

_rhoE = autoPtr<volScalarField>
      (
          new volScalarField
          (
              IOobject
              (
                  "rhoE",
                  runTime.timeName(),
                  mesh,
                  IOobject::NO_READ,
                  IOobject::NO_WRITE
              ),
              rho*(e + 0.5*magSqr(U))
          )
      );

_pos = autoPtr<surfaceScalarField>
      (
          new surfaceScalarField
          (
              IOobject
              (
                  "pos",
                  runTime.timeName(),
                  mesh
              ),
              mesh,
              dimensionedScalar("pos", dimless, 1.0)
          )
      );

_neg = autoPtr<surfaceScalarField>
      (
          new surfaceScalarField
          (
              IOobject
              (
                  "neg",
                  runTime.timeName(),
                  mesh
              ),
              mesh,
              dimensionedScalar("neg", dimless, -1.0)
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
                   mesh
              ),
              fvc::flux(rhoU)
           )
       );

surfaceScalarField& phi = _phi();

Info << "Creating turbulence model\n" << endl;
turbulence = autoPtr<compressible::turbulenceModel>
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);
