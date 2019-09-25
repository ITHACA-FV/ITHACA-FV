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


_tau = autoPtr<dimensionedScalar>
       (
           new dimensionedScalar
           (
               transportProperties.lookup("tau")
           )
       );
dimensionedScalar& tau = _tau();

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


