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

pThermo = autoPtr<fluidThermo>
          (
              fluidThermo::New(mesh)
          );

fluidThermo& thermo = pThermo();

thermo.validate(args.executable(), "h", "e");

_p = autoPtr<volScalarField>
     (
         new volScalarField(thermo.p())
     );

volScalarField& p = _p();

volScalarField p0(thermo.p());

_p0 = autoPtr<volScalarField>
      (
          new volScalarField(p0)
      );

_rho = autoPtr<volScalarField>
       (
           new volScalarField
           (
               IOobject
               (
                   "rho",
                   runTime.timeName(),
                   mesh,
                   IOobject::READ_IF_PRESENT,
                   IOobject::AUTO_WRITE
               ),
               thermo.rho()
           )
       );

volScalarField& rho = _rho();

volScalarField rho0(rho);

_rho0 = autoPtr<volScalarField>
      (
          new volScalarField(rho0)
      );

_E = autoPtr<volScalarField>
     (
         new volScalarField(thermo.he())
     );

volScalarField E0(_E);

_E0 = autoPtr<volScalarField>
      (
          new volScalarField(E0)
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

#include "compressibleCreatePhi.H"

_pressureControl = autoPtr<pressureControl>
                   (
                       new pressureControl(p, rho, simple.dict())
                   );

//pressureControl& pressureControl = _pressureControl();

mesh.setFluxRequired(p.name());

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

_initialMass = autoPtr<dimensionedScalar>
               (
                   new dimensionedScalar(fvc::domainIntegrate(rho))
               );

_MRF = autoPtr<IOMRFZoneList>
       (
           new IOMRFZoneList(mesh)
       );

_psi = autoPtr<volScalarField>
       (
           new volScalarField(thermo.psi())
       );

_fvOptions = autoPtr<fv::options>(new fv::options(mesh));

//turbulence->validate();

/*_laminarTransport = autoPtr<singlePhaseTransportModel>
                    (
                        new singlePhaseTransportModel( U, phi )
                    );
singlePhaseTransportModel& laminarTransport = _laminarTransport();*/