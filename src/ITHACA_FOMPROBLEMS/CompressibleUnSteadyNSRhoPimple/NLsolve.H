
//Time& runTime = _runTime();
dynamicFvMesh& mesh = meshPtr();
dimensionedScalar& initialMass = _initialMass();
surfaceScalarField& phi = _phi();
volScalarField& rho = _rho();
fluidThermo& thermo = pThermo();
pressureControl& pressureControl = _pressureControl();
volScalarField& psi = _psi();
bool closedVolume = false;
IOMRFZoneList& MRF = _MRF();
fv::options& fvOptions = _fvOptions();
const dimensionedScalar rhoMax("rhoMax", dimDensity, GREAT, pimple.dict());
const dimensionedScalar rhoMin("rhoMin", dimDensity, Zero, pimple.dict());


Info << "\nStarting time loop\n" << endl;
#include "createRhoUfIfPresent.H"
#include "createDpdt.H"
#include "createK.H"

while (runTime.run())
{
    //#include "readDyMControls.H"
    // Store divrhoU from the previous mesh so that it can be mapped
    // and used in correctPhi to ensure the corrected phi has the
    // same divergence
    autoPtr<volScalarField> divrhoU;

    if (correctPhi)
    {
        divrhoU.reset
        (
            new volScalarField
            (
                "divrhoU",
                fvc::div(fvc::absolute(phi, rho, U))
            )
        );
    }

#include "compressibleCourantNo.H"
    //#include "setDeltaT.H"
    runTime.setEndTime(finalTime);
    ++runTime;
    Info << "Time = " << runTime.timeName() << nl << endl;

    // --- Pressure-velocity PIMPLE corrector loop
    while (pimple.loop())
    {
        if (pimple.firstIter() || moveMeshOuterCorrectors)
        {
            // Store momentum to set rhoUf for introduced faces.
            autoPtr<volVectorField> rhoU;

            if (rhoUf.valid())
            {
                rhoU.reset(new volVectorField("rhoU", rho * U));
            }

            // Do any mesh changes
            mesh.controlledUpdate();

            if (mesh.changing())
            {
                MRF.update();

                if (correctPhi)
                {
                    // Calculate absolute flux
                    // from the mapped surface velocity
                    phi = mesh.Sf() & rhoUf();
#include "correctPhi.H"
                    // Make the fluxes relative to the mesh-motion
                    fvc::makeRelative(phi, rho, U);
                }

                if (checkMeshCourantNo)
                {
#include "meshCourantNo.H"
                }
            }
        }

        if (pimple.firstIter() && !pimple.SIMPLErho())
        {
#include "rhoEqn.H"
        }

#include "UEqn.H"
#include "EEqn.H"

        // --- Pressure corrector loop
        while (pimple.correct())
        {
#include "pEqn.H"
        }

        if (pimple.turbCorr())
        {
            turbulence->correct();
        }
    }

    rho = thermo.rho();

    if (checkWrite(runTime))
    {
        ITHACAstream::exportSolution(U, name(counter), folder);
        ITHACAstream::exportSolution(p, name(counter), folder);
        ITHACAstream::exportSolution(E, name(counter), folder);
        ITHACAstream::writePoints(meshPtr().points(), folder,
                                  name(counter) + "/polyMesh/");
        //ITHACAstream::exportSolution(_nut, name(counter), "./ITHACAoutput/Offline/");
        std::ofstream of(folder + name(counter) + "/" + runTime.timeName());
        Ufield.append(U.clone());
        Pfield.append(p.clone());
        Efield.append(E.clone());
        nextWrite += writeEvery;
        counter++;
        //nutFields.append(_nut);
    }
}
