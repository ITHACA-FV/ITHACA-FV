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


/// \file
/// Source file of the steadyNS class.

#include "steadyNS_simple.H"
#include "viscosityModel.H"

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructor
steadyNS_simple::steadyNS_simple() {}

steadyNS_simple::steadyNS_simple(int argc, char* argv[])
    :
    steadyNS(argc, argv)
{
    Info << offline << endl;
}

fvVectorMatrix steadyNS_simple::get_Umatrix(volVectorField& U,
        volScalarField& p)
{
    IOMRFZoneList& MRF = _MRF();
    surfaceScalarField& phi = _phi();
    fv::options& fvOptions = _fvOptions();
    MRF.correctBoundaryVelocity(U);
    fvVectorMatrix Ueqn
    (
        fvm::div(phi, U)
        + MRF.DDt(U)
        + turbulence->divDevReff(U)
        ==
        fvOptions(U)
    );
    Ueqn.relax();
    fvOptions.constrain(Ueqn);
    Ueqn_global = &Ueqn;
    return Ueqn;
}

fvScalarMatrix steadyNS_simple::get_Pmatrix(volVectorField& U,
        volScalarField& p, scalar& presidual)
{
    IOMRFZoneList& MRF = _MRF();
    surfaceScalarField& phi = _phi();
    simpleControl& simple = _simple();
    fv::options& fvOptions = _fvOptions();
    MRF.correctBoundaryVelocity(U);
    fvMesh& mesh = _mesh();
    volScalarField rAU(1.0 / Ueqn_global->A());
    volVectorField HbyA(constrainHbyA(rAU * Ueqn_global->H(), U, p));
    surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
    MRF.makeRelative(phiHbyA);
    adjustPhi(phiHbyA, U, p);
    tmp<volScalarField> rAtU(rAU);

    if (simple.consistent())
    {
        rAtU = 1.0 / (1.0 / rAU - Ueqn_global->H1());
        phiHbyA +=
            fvc::interpolate(rAtU() - rAU) * fvc::snGrad(p) * mesh.magSf();
        HbyA -= (rAU - rAtU()) * fvc::grad(p);
    }

    constrainPressure(p, U, phiHbyA, rAtU(), MRF);
    int i = 0;

    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA)
        );
        pEqn.setReference(pRefCell, pRefValue);

        if (i == 0)
        {
            presidual = pEqn.solve().initialResidual();
        }
        else
        {
            pEqn.solve().initialResidual();
        }

        if (simple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - pEqn.flux();
        }

        i++;
    }

    //p.storePrevIter(); // Perché ho dovuto metterlo se nel solver non c'è???
    p.relax();
    U = HbyA - rAtU() * fvc::grad(p);
    U.correctBoundaryConditions();
    fvOptions.correct(U);
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA)
    );
    return pEqn;
}

void steadyNS_simple::truthSolve2(List<scalar> mu_now)
{
    Time& runTime = _runTime();
    volScalarField& p = _p();
    volVectorField& U = _U();
    fv::options& fvOptions = _fvOptions();
    simpleControl& simple = _simple();
    singlePhaseTransportModel& laminarTransport = _laminarTransport();
#include "NLsolve.H"
    ITHACAstream::exportSolution(U, name(counter), "./ITHACAoutput/Offline/");
    ITHACAstream::exportSolution(p, name(counter), "./ITHACAoutput/Offline/");
    Ufield.append(U);
    Pfield.append(p);
    counter++;
    writeMu(mu_now);
    // --- Fill in the mu_samples with parameters (mu) to be used for the POD sample points
    mu_samples.conservativeResize(mu_samples.rows() + 1, mu_now.size());

    for (int i = 0; i < mu_now.size(); i++)
    {
        mu_samples(mu_samples.rows() - 1, i) = mu_now[i];
    }

    // Resize to Unitary if not initialized by user (i.e. non-parametric problem)
    if (mu.cols() == 0)
    {
        mu.resize(1, 1);
    }

    if (mu_samples.rows() == mu.cols())
    {
        ITHACAstream::exportMatrix(mu_samples, "mu_samples", "eigen",
                                   "./ITHACAoutput/Offline");
    }
}
