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
_flux = autoPtr<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "flux",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );
volScalarField& flux = _flux();

volScalarField flux0(flux);

_flux0 = autoPtr<volScalarField>
         (
             new volScalarField(flux0)
         );



_prec1 = autoPtr<volScalarField>
         (
             new volScalarField
             (
                 IOobject
                 (
                     "prec1",
                     runTime.timeName(),
                     mesh,
                     IOobject::MUST_READ,
                     IOobject::AUTO_WRITE
                 ),
                 mesh
             )
         );
volScalarField& prec1 = _prec1();

volScalarField prec10(prec1);

_prec10 = autoPtr<volScalarField>
          (
              new volScalarField(prec10)
          );


_prec2 = autoPtr<volScalarField>
         (
             new volScalarField
             (
                 IOobject
                 (
                     "prec2",
                     runTime.timeName(),
                     mesh,
                     IOobject::MUST_READ,
                     IOobject::AUTO_WRITE
                 ),
                 mesh
             )
         );
volScalarField& prec2 = _prec2();

volScalarField prec20(prec2);

_prec20 = autoPtr<volScalarField>
          (
              new volScalarField(prec20)
          );


_prec3 = autoPtr<volScalarField>
         (
             new volScalarField
             (
                 IOobject
                 (
                     "prec3",
                     runTime.timeName(),
                     mesh,
                     IOobject::MUST_READ,
                     IOobject::AUTO_WRITE
                 ),
                 mesh
             )
         );
volScalarField& prec3 = _prec3();

volScalarField prec30(prec3);

_prec30 = autoPtr<volScalarField>
          (
              new volScalarField(prec30)
          );



_prec4 = autoPtr<volScalarField>
         (
             new volScalarField
             (
                 IOobject
                 (
                     "prec4",
                     runTime.timeName(),
                     mesh,
                     IOobject::MUST_READ,
                     IOobject::AUTO_WRITE
                 ),
                 mesh
             )
         );
volScalarField& prec4 = _prec4();

volScalarField prec40(prec4);

_prec40 = autoPtr<volScalarField>
          (
              new volScalarField(prec40)
          );



_prec5 = autoPtr<volScalarField>
         (
             new volScalarField
             (
                 IOobject
                 (
                     "prec5",
                     runTime.timeName(),
                     mesh,
                     IOobject::MUST_READ,
                     IOobject::AUTO_WRITE
                 ),
                 mesh
             )
         );
volScalarField& prec5 = _prec5();

volScalarField prec50(prec5);

_prec50 = autoPtr<volScalarField>
          (
              new volScalarField(prec50)
          );



_prec6 = autoPtr<volScalarField>
         (
             new volScalarField
             (
                 IOobject
                 (
                     "prec6",
                     runTime.timeName(),
                     mesh,
                     IOobject::MUST_READ,
                     IOobject::AUTO_WRITE
                 ),
                 mesh
             )
         );
volScalarField& prec6 = _prec6();

volScalarField prec60(prec6);

_prec60 = autoPtr<volScalarField>
          (
              new volScalarField(prec60)
          );



_prec7 = autoPtr<volScalarField>
         (
             new volScalarField
             (
                 IOobject
                 (
                     "prec7",
                     runTime.timeName(),
                     mesh,
                     IOobject::MUST_READ,
                     IOobject::AUTO_WRITE
                 ),
                 mesh
             )
         );
volScalarField& prec7 = _prec7();

volScalarField prec70(prec7);

_prec70 = autoPtr<volScalarField>
          (
              new volScalarField(prec70)
          );



_prec8 = autoPtr<volScalarField>
         (
             new volScalarField
             (
                 IOobject
                 (
                     "prec8",
                     runTime.timeName(),
                     mesh,
                     IOobject::MUST_READ,
                     IOobject::AUTO_WRITE
                 ),
                 mesh
             )
         );
volScalarField& prec8 = _prec8();

volScalarField prec80(prec8);

_prec80 = autoPtr<volScalarField>
          (
              new volScalarField(prec80)
          );




