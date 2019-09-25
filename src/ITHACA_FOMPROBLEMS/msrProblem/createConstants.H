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
#include "readConstants.H"
#include "readConstants_N.H"
#include "readConstants_TH.H"


_betaTot = autoPtr<dimensionedScalar>
           (
               new dimensionedScalar
               (
                   beta1 + beta2 + beta3 + beta4 + beta5 + beta6 + beta7 + beta8
               )
           );
dimensionedScalar& betaTot = _betaTot();

_decbetaTot = autoPtr<dimensionedScalar>
              (
                  new dimensionedScalar
                  (
                      decBeta1 + decBeta2 + decBeta3
                  )
              );
dimensionedScalar& decbetaTot = _decbetaTot();

_v = autoPtr<volScalarField>
     (
         new volScalarField
         (
             IOobject
             (
                 "v",
                 runTime.timeName(),
                 mesh,
                 IOobject::NO_READ,
                 IOobject::AUTO_WRITE
             ),
             mesh,
             dimensionedScalar("v", dimensionSet(-1, 3, 0, 0, 0, 0, 0), 1 / rhoRef.value())
         )
     );
volScalarField& v = _v();

volScalarField v0(v);

_v0 = autoPtr<volScalarField>
      (
          new volScalarField(v0)
      );

_D = autoPtr<volScalarField>
     (
         new volScalarField
         (
             IOobject
             (
                 "D",
                 runTime.timeName(),
                 mesh,
                 IOobject::NO_READ,
                 IOobject::AUTO_WRITE
             ),
             mesh,
             dimensionedScalar("D", dimensionSet(0, 1, 0, 0, 0, 0, 0), D1_0.value())
         )
     );
volScalarField& D = _D();

volScalarField D0(D);

_D0 = autoPtr<volScalarField>
      (
          new volScalarField(D0)
      );

_NSF = autoPtr<volScalarField>
       (
           new volScalarField
           (
               IOobject
               (
                   "NSF",
                   runTime.timeName(),
                   mesh,
                   IOobject::NO_READ,
                   IOobject::AUTO_WRITE
               ),
               mesh,
               dimensionedScalar("NSF", dimensionSet(0, -1, 0, 0, 0, 0, 0), NSF1_0.value())
           )
       );
volScalarField& NSF = _NSF();

volScalarField NSF0(NSF);

_NSF0 = autoPtr<volScalarField>
        (
            new volScalarField(NSF0)
        );

_A = autoPtr<volScalarField>
     (
         new volScalarField
         (
             IOobject
             (
                 "A",
                 runTime.timeName(),
                 mesh,
                 IOobject::NO_READ,
                 IOobject::AUTO_WRITE
             ),
             mesh,
             dimensionedScalar("A", dimensionSet(0, -1, 0, 0, 0, 0, 0), A1_0.value())
         )
     );
volScalarField& A = _A();

volScalarField A0(A);

_A0 = autoPtr<volScalarField>
      (
          new volScalarField(A0)
      );


_SP = autoPtr<volScalarField>
      (
          new volScalarField
          (
              IOobject
              (
                  "SP",
                  runTime.timeName(),
                  mesh,
                  IOobject::NO_READ,
                  IOobject::AUTO_WRITE
              ),
              mesh,
              dimensionedScalar("SP", dimensionSet(1, 1, -2, 0, 0, 0, 0), SP1_0.value())
          )
      );
volScalarField& SP = _SP();

volScalarField SP0(SP);

_SP0 = autoPtr<volScalarField>
       (
           new volScalarField(SP0)
       );

_TXS = autoPtr<volScalarField>
       (
           new volScalarField
           (
               IOobject
               (
                   "TXS",
                   runTime.timeName(),
                   mesh,
                   IOobject::NO_READ,
                   IOobject::AUTO_WRITE
               ),
               v* SP
           )
       );
volScalarField& TXS = _TXS();

volScalarField TXS0(TXS);

_TXS0 = autoPtr<volScalarField>
        (
            new volScalarField(TXS0)
        );

_logT = autoPtr<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "logT",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                log(max(0.5, log(T / TrefXS)))
            )
        );
volScalarField& logT = _logT();

_alphat = autoPtr<volScalarField>
          (
              new volScalarField
              (
                  IOobject
                  (
                      "alphat",
                      runTime.timeName(),
                      mesh,
                      IOobject::NO_READ,
                      IOobject::AUTO_WRITE
                  ),
                  turbulence->nut() / Prt
              )
          );
volScalarField& alphat = _alphat();

_difft = autoPtr<volScalarField>
         (
             new volScalarField
             (
                 IOobject
                 (
                     "difft",
                     runTime.timeName(),
                     mesh,
                     IOobject::NO_READ,
                     IOobject::AUTO_WRITE
                 ),
                 turbulence->nut() / Sct
             )
         );
volScalarField& difft = _difft();
