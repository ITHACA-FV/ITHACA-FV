volScalarField rAU(1.0 /
                   Ueqn.A()); // Inverse of the diagonal part of the U equation matrix
//volScalarField rAU(1.0/Ueqn_global->A()); // For the new "getUmatrix" formulation
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho* rAU));
volVectorField HbyA(constrainHbyA(rAU* Ueqn.H(), U,
                                  p)); // H is the extra diagonal part summed to the r.h.s. of the U equation
//volVectorField HbyA(constrainHbyA(rAU*Ueqn_global->H(), U, p)); // For the new "getUmatrix" formulation
//tUEqn.clear();

bool closedVolume = false;

surfaceScalarField phiHbyA("phiHbyA", fvc::interpolate(rho)*fvc::flux(HbyA));
MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, rho, U, phiHbyA, rhorAUf, MRF);

if (simple.transonic())
{
    surfaceScalarField phid
    (
        "phid",
        (fvc::interpolate(psi) / fvc::interpolate(rho))*phiHbyA
    );
    phiHbyA -= fvc::interpolate(psi * p) * phiHbyA / fvc::interpolate(rho);

    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvc::div(phiHbyA)
            + fvm::div(phid, p)
            - fvm::laplacian(rhorAUf, p)
            ==
            fvOptions(psi, p, rho.name())
        );
        // Relax the pressure equation to ensure diagonal-dominance
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
    //Passa sempre da qui!!
    closedVolume = adjustPhi(phiHbyA, U, p);
    //p.storePrevIter();

    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvc::div(phiHbyA)
            - fvm::laplacian(rhorAUf, p)
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

#include "incompressible/continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

U = HbyA - rAU * fvc::grad(p);
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
