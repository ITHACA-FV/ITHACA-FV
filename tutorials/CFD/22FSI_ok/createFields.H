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
std::cout << "*****create fields of a  dynamic case****\n"<< std::endl; 
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
#include "createPhi.H"

label pRefCell = 0;
scalar pRefValue = 0.0;
//setRefCell(p, mesh.solutionDict().subDict("SIMPLE"), pRefCell, pRefValue);
setRefCell(p, mesh.solutionDict().subDict("PIMPLE"), pRefCell, pRefValue);
//setRefCell(p, _pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

_laminarTransport = autoPtr<singlePhaseTransportModel>
                    (
                        new singlePhaseTransportModel( U, phi )
                    );
singlePhaseTransportModel& laminarTransport = _laminarTransport();

turbulence = autoPtr<incompressible::turbulenceModel>
             (
                 incompressible::turbulenceModel::New(U, phi, laminarTransport)
             );
//#include "createMRF.H"
_MRF = autoPtr<IOMRFZoneList>
       (
           new IOMRFZoneList(mesh)
       );
#include "createFvOptions.H"
