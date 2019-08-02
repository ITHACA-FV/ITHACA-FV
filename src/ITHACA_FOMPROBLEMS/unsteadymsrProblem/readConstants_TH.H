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

