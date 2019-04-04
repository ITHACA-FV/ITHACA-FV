rho = thermo.rho();

volScalarField rAU(1.0 / Ueqn.A());
volScalarField rAtU(1.0 / (1.0 / rAU - Ueqn.H1()));
volVectorField HbyA(constrainHbyA(rAU* Ueqn.H(), U, p));
//volScalarField rAU(1.0/Ueqn_global->A()); // For the new "getUmatrix" formulation
//volScalarField rAtU(1.0/(1.0/rAU - Ueqn_global->H1())); // For the new "getUmatrix" formulation
//volVectorField HbyA(constrainHbyA(rAU*Ueqn_global->H(), U, p)); // For the new "getUmatrix" formulation
//tUEqn.clear();

bool closedVolume = false;

surfaceScalarField phiHbyA("phiHbyA", fvc::interpolate(rho)*fvc::flux(HbyA));
MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

volScalarField rhorAtU("rhorAtU", rho* rAtU);

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, rho, U, phiHbyA, rhorAtU, MRF);

if (simple.transonic())
{
    surfaceScalarField phid
    (
        "phid",
        (fvc::interpolate(psi) / fvc::interpolate(rho))*phiHbyA
    );
    phiHbyA +=
        fvc::interpolate(rho * (rAtU - rAU)) * fvc::snGrad(p) * mesh.magSf()
        - fvc::interpolate(psi * p) * phiHbyA / fvc::interpolate(rho);
    HbyA -= (rAU - rAtU) * fvc::grad(p);

    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvc::div(phiHbyA)
            + fvm::div(phid, p)
            - fvm::laplacian(rhorAtU, p)
            ==
            fvOptions(psi, p, rho.name())
        );
        // Relax the pressure equation to maintain diagonal dominance
        pEqn.relax();
        pEqn.setReference
        (
            pressureControl.refCell(),
            pressureControl.refValue()
        );
        presidual = pEqn.solve().initialResidual();

        if (simple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + pEqn.flux();
        }
    }
}
else
{
    closedVolume = adjustPhi(phiHbyA, U, p);
    phiHbyA += fvc::interpolate(rho * (rAtU - rAU)) * fvc::snGrad(p) * mesh.magSf();
    HbyA -= (rAU - rAtU) * fvc::grad(p);

    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvc::div(phiHbyA)
            - fvm::laplacian(rhorAtU, p)
            ==
            fvOptions(psi, p, rho.name())
        );
        pEqn.setReference
        (
            pressureControl.refCell(),
            pressureControl.refValue()
        );
        presidual = pEqn.solve().initialResidual();

        if (simple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + pEqn.flux();
        }
    }
}

// The incompressibe form of the continuity error check is appropriate for
// steady-state compressible also.
#include "incompressible/continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

U = HbyA - rAtU * fvc::grad(p);
U.correctBoundaryConditions();
fvOptions.correct(U);

bool pLimited = pressureControl.limit(p);

// For closed-volume cases adjust the pressure and density levels
// to obey overall mass continuity
if (closedVolume)
{
    p += (initialMass - fvc::domainIntegrate(psi * p))
         / fvc::domainIntegrate(psi);
}

if (pLimited || closedVolume)
{
    p.correctBoundaryConditions();
}

rho = thermo.rho();

if (!simple.transonic())
{
    rho.relax();
}
