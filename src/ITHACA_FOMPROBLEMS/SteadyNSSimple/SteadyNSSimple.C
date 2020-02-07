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

#include "SteadyNSSimple.H"
#include "viscosityModel.H"

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructor
SteadyNSSimple::SteadyNSSimple() {}

SteadyNSSimple::SteadyNSSimple(int argc, char* argv[])
    :
    steadyNS(argc, argv)
{
    Info << offline << endl;
}

fvVectorMatrix SteadyNSSimple::get_Umatrix(volVectorField& U,
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

fvScalarMatrix SteadyNSSimple::get_Pmatrix(volVectorField& U,
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

void SteadyNSSimple::truthSolve2(List<scalar> mu_now, word Folder)
{
    Time& runTime = _runTime();
    volScalarField& p = _p();
    volVectorField& U = _U();
    fvMesh& mesh = _mesh();
    surfaceScalarField& phi = _phi();
    simpleControl& simple = _simple();
    singlePhaseTransportModel& laminarTransport = _laminarTransport();
    scalar residual = 1;
    scalar uresidual = 1;
    Vector<double> uresidual_v(0, 0, 0);
    scalar presidual = 1;
    scalar csolve = 0;
    // Variable that can be changed
    turbulence->read();
    std::ofstream res_os;
    res_os.open("./ITHACAoutput/Offline/residuals", std::ios_base::app);
#if OFVER == 6

    while (simple.loop(runTime) && residual > tolerance && csolve < maxIter )
#else
    while (simple.loop() && residual > tolerance && csolve < maxIter )
#endif
    {
        Info << "Time = " << runTime.timeName() << nl << endl;
        volScalarField nueff = turbulence->nuEff();
        fvVectorMatrix UEqn
        (
            fvm::div(phi, U)
            - fvm::laplacian(nueff, U)
            - fvc::div(nueff * dev2(T(fvc::grad(U))))
        );
        UEqn.relax();
        UEqn == - fvc::grad(p);

        if (simple.momentumPredictor())
        {
            uresidual_v = solve(UEqn).initialResidual();
        }

        scalar C = 0;

        for (label i = 0; i < 3; i++)
        {
            if (C < uresidual_v[i])
            {
                C = uresidual_v[i];
            }
        }

        uresidual = C;
        volVectorField HbyA(constrainHbyA(1.0 / UEqn.A() * UEqn.H(), U, p));
        surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
        int i = 0;

        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn
            (
                fvm::laplacian(1.0 / UEqn.A(), p) == fvc::div(phiHbyA)
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

#include "continuityErrs.H"
        p.relax();
        // Momentum corrector
        U = HbyA - 1.0 / UEqn.A() * fvc::grad(p);
        U.correctBoundaryConditions();
        residual = max(presidual, uresidual);
        Info << "Time = " << runTime.timeName() << nl << endl;
        laminarTransport.correct();
        turbulence->correct();
    }

    res_os << residual << std::endl;
    res_os.close();
    runTime.setTime(runTime.startTime(), 0);
    ITHACAstream::exportSolution(U, name(counter), Folder);
    ITHACAstream::exportSolution(p, name(counter), Folder);
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
                                   Folder);
    }
}