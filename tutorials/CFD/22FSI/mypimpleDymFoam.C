int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
 
    pimpleControl pimple(mesh);
 
    #include "createFields.H"
    #include "createUf.H"
    #include "createFvOptions.H"
    #include "readTimeControls.H"
    #include "createPcorrTypes.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"
 
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
    Info<< "\nStarting time loop\n" << endl; // Time loop as per controlDict
 
    while (runTime.run()) // while time < endTime
    {
        #include "readControls.H"  // Read adjustTimeStep, maxCo, maxDeltaT
        #include "CourantNo.H" // Calculate Mean and Max Courant Numbers
 
        #include "setDeltaT.H" // Get the deltaT defined in controlDict
        if(adjustTimeStep){
          scalar maxDeltaTFact = maxCo/(CoNum+SMALL);
          scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);
          runTime.setDeltaT(
            min(deltaTFact*runTime.deltaTValue(), maxDeltaT)
          );
        }
 
        runTime++;
 
        Info<< "Time = " << runTime.timeName() << nl << endl;
 
        mesh.update();
 
        // Calculate absolute flux from the mapped surface velocity
        phi = mesh.Sf() & Uf;
 
        if (mesh.changing() && correctPhi)
        {
            #include "correctPhi.H" // Flux field phi of velocity U at cell faces
        }
 
        // Make the flux relative to the mesh motion
        fvc::makeRelative(phi, U);
 
        if (mesh.changing() && checkMeshCourantNo)
        {
            #include "meshCourantNo.H"
        }
 
        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            // Set-up the linear algebra for momentum equation and UEqn.H
            // corresponds to the momentum predictor, and pEqn.H to the 
            // corrector step(s).
            #include "UEqn.H"
 
            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H" // include the pressure equation
            }
 
            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }
 
        runTime.write();
        
        //Print on the screen information-computation and clock time
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }
 
    Info<< "End\n" << endl;
 
    return 0;
}
