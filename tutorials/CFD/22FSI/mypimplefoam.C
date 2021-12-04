#include "fvCFD.H"

#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

int main(int argc, char *argv[])   {
  
  // Check and set root case folder: $FOAM_SRC/OpenFOAM/include -Checks the 
  // basic folder structure, verifies there is a control dict present ...
  // also deals with parsing command line arguments and options.
#include "setRootCase.H"
#include "listOptions.H"
  Foam::argList args(argc, argv);
    if (!args.checkRootCase()) {
      Foam::FatalError.exit();  
	//exit() is a member function of class 'FatalError' in namespace 'Foam'
    }
#include "listOutput.H"

// Create time directories: $FOAM_SRC/OpenFOAM/include
// Create the time system (instance called runTime) 
#include "createTime.H"
    Foam::Info<< "Create time\n" << Foam::endl;
    Foam::Time runTime(Foam::Time::controlDictName,args);

// Read mesh: $FOAM_SRC/OpenFOAM/include
  // Create fvMesh (instances of objects (or classes) called mesh). mesh.C() 
  // and .Cf() return vector fields denoting centres of each cell and internal
  // face. Calling the mesh.C().size() method yields total size of the mesh.
  #include "createMesh.H"
    Foam::Info
      << "Create mesh for time = "
      << runTime.timeName() << Foam::nl << Foam::endl;
      Foam::fvMesh mesh 
        (Foam::IOobject 
          (Foam::fvMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ  //Must read mesh every timestep
          )
        );

#include "createFields.H"
    Info<< "Reading field p\n" << endl;
    volScalarField p  (        //define p volume scalar field
      IOobject  (
        "p",                   //name of the dictionary, field
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,   //must read 'p' field every timestep
        IOobject::AUTO_WRITE
      ),
      mesh
    );

    Info<< "Reading field U\n" << endl;
    volVectorField U  (        //define U volume vector field
      IOobject (
        "U",                   //name of the dictionary
        runTime.timeName(),    // look for dictionary "U"
        mesh,
        IOobject::MUST_READ,   //must read 'U' field every timestep
        IOobject::AUTO_WRITE
      ),
     mesh
    );

#include "createPhi.H"    // Flux field phi of velocity U at cell faces

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("PIMPLE"), pRefCell, pRefValue);

    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::turbulenceModel> turbulence  (
      incompressible::turbulenceModel::New(U, phi, laminarTransport)
    );
    
    //Declare and initialise the cumulative continuity error.
    //$FOAM_SRC/finiteVolume/cfdTools/general/include	
    #include "initContinuityErrs.H"
    #ifndef initContinuityErrs_H
    #define initContinuityErrs_H
    scalar cumulativeContErr = 0;
    #endif
pimpleControl pimple(mesh);
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;  //Time loop as per controlDict
    while (runTime.run())  {         //While time < endTime
      #include "readTimeControls.H"  // Read adjustTimeStep, maxCo, maxDeltaT
      
      #include "CourantNo.H"         // Calculate Mean and Max Courant Numbers 

      /*$FOAM_SRC/finiteVolume/cfdTools/general/include:  Reset the timestep 
        to maintain a constant maximum courant Number. Reduction of time-step 
        is immediate, but increase is damped to avoid unstable oscillations */
      #include "setDeltaT.H"         // Get the deltaT defined in controlDict
      if (adjustTimeStep) {
        scalar maxDeltaTFact = maxCo/(CoNum + SMALL);
        scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);
        runTime.setDeltaT (
         min (
          deltaTFact*runTime.deltaTValue(),
          maxDeltaT
         )
        );
      }
      runTime++;

      Info<< "Time = " << runTime.timeName() << nl << endl;

      // --- Pressure-velocity PIMPLE corrector loop:  as many outer correctors 
      // as nOuterCorrectors  and as many inner correctors as nCorrectors defined
      // in fvSolution dictionary.
      while (pimple.loop())   {    
        
		//Set-up the linear algebra for momentum equation and solve. UEqn.H 
		//corresponds to the momentum predictor, and pEqn.H corresponds to the
		//corrector step(s).
		#include "UEqn.H"

        // --- Pressure corrector loop
        while (pimple.correct()) {
          #include "pEqn.H"  //Include the pressure Equation
        }
        if (pimple.turbCorr())  {
          turbulence->correct();
        }
      }
      //Print on the screen information - computational and clock time
      runTime.write();
      Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
          << "  ClockTime = " << runTime.elapsedClockTime() << " s"
          << nl << endl;
    }
    Info<< "End\n" << endl;
    return 0;
}

