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
#include "createRDeltaT.H"

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

_mixture = autoPtr<immiscibleIncompressibleTwoPhaseMixture>
	   (
		new immiscibleIncompressibleTwoPhaseMixture(U, phi)
	   );
immiscibleIncompressibleTwoPhaseMixture& mixture = _mixture();

//need autoPtr?
volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();


// Need to store rho for ddt(rho, U)
_rho = autoPtr<volScalarField>
( 
new volScalarField
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1 + alpha2*rho2
)
);
volScalarField& rho = _rho();
rho.oldTime();

// Mass flux

_rhoPhi = autoPtr<surfaceScalarField>
(
    new
surfaceScalarField 
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
)
);
surfaceScalarField& rhoPhi = _rhoPhi();


turbulence = autoPtr<incompressible::turbulenceModel>
             (
                 incompressible::turbulenceModel::New(U, phi, mixture)
             );

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

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
    p_rgh + rho*gh
)
    );
volScalarField& p = _p();

label pRefCell = 0;
scalar pRefValue = 0.0;
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
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

_MRF = autoPtr<IOMRFZoneList>
       (
           new IOMRFZoneList(mesh)
       );
#include "createUfIfPresent.H"