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

volScalarField T0(T);

_T0 = autoPtr<volScalarField>
      (
          new volScalarField(T0)
      );




_dec1 = autoPtr<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "dec1",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );
volScalarField& dec1 = _dec1();

volScalarField dec10(dec1);

_dec10 = autoPtr<volScalarField>
         (
             new volScalarField(dec10)
         );



_dec2 = autoPtr<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "dec2",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );
volScalarField& dec2 = _dec2();


volScalarField dec20(dec2);

_dec20 = autoPtr<volScalarField>
         (
             new volScalarField(dec20)
         );



_dec3 = autoPtr<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "dec3",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );
volScalarField& dec3 = _dec3();


volScalarField dec30(dec3);

_dec30 = autoPtr<volScalarField>
         (
             new volScalarField(dec30)
         );


