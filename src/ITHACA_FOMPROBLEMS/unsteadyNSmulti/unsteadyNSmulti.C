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

#include "unsteadyNSmulti.H"

/// \file
/// Source file of the unsteadyNSmulti class.

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Construct Null
unsteadyNSmulti::unsteadyNSmulti(){};

// Construct from zero
unsteadyNSmulti::unsteadyNSmulti(int argc, char* argv[])
{
 _args = autoPtr<argList>
 (
   new argList(argc, argv)
   );
 if (!_args->checkRootCase())
 {
   Foam::FatalError.exit();
 }

 argList& args = _args();
 #include "createTime.H"
 #include "createMesh.H"
 _pimple = autoPtr<pimpleControl>
 (
   new pimpleControl
   (
     mesh
     )
   );
 pimpleControl& pimple = _pimple();
 #include "createTimeControls.H"               
 correctPhi = pimple.dict().lookupOrDefault("correctPhi", mesh.dynamic());
 checkMeshCourantNo = pimple.dict().lookupOrDefault("checkMeshCourantNo", false);
 moveMeshOuterCorrectors  = pimple.dict().lookupOrDefault("moveMeshOuterCorrectors", false);
 #include "createFields.H"
 #include "createAlphaFluxes.H"
 #include "createFvOptions.H"
 #include "initCorrectPhi.H"
 para = new ITHACAparameters;
 offline = ITHACAutilities::check_off();
 podex = ITHACAutilities::check_pod();
 supex = ITHACAutilities::check_sup();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void unsteadyNSmulti::truthSolve(List<scalar> mu_now)
{
#include "initContinuityErrs.H"
    Time& runTime = _runTime();
    surfaceScalarField& phi = _phi();
    fvMesh& mesh = _mesh();
    fv::options& fvOptions = _fvOptions();
    pimpleControl& pimple = _pimple();
    volScalarField& p = _p();
    volScalarField& p_rgh = _p_rgh();
    volScalarField& rho = _rho();
    volVectorField& U = _U();
    surfaceScalarField& rhoPhi = _rhoPhi();
    immiscibleIncompressibleTwoPhaseMixture& mixture = _mixture();
    volScalarField& alpha1(mixture.alpha1());
    volScalarField& alpha2(mixture.alpha2());
    const dimensionedScalar& rho1(mixture.rho1());
    const dimensionedScalar& rho2(mixture.rho2());
    uniformDimensionedVectorField& g = _g();
    uniformDimensionedScalarField& hRef = _hRef();
    volScalarField& gh = _gh();
    surfaceScalarField& ghf = _ghf();
    surfaceScalarField& alphaPhi10 = _alphaPhi10();
    tmp<surfaceScalarField>& talphaPhi1Corr0 = _talphaPhi1Corr0();
    IOMRFZoneList& MRF = _MRF();

    instantList Times = runTime.times();
    runTime.setEndTime(finalTime);
    // Perform a TruthSolve
    runTime.setTime(Times[1], 1);
    runTime.setDeltaT(timeStep);
    nextWrite = startTime;
    // Initialize Nsnapshots
    int nsnapshots = 0;

    
        while (runTime.run())
    {
        #include "readDyMControls.H"
        #include "CourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"
              runTime++;
        Info << "Time = " << runTime.timeName() << nl << endl;

        while (pimple.loop())
        {
            #include "alphaControls.H"
            #include "alphaEqnSubCycle.H"
            mixture.correct();
            #include "UEqn.H"
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

                runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

//     Time& runTime = _runTime();
//     surfaceScalarField& phi = _phi();
//     fvMesh& mesh = _mesh();
//     fv::options& fvOptions = _fvOptions();
//     pimpleControl& pimple = _pimple();
//     volScalarField p = _p();
//     volVectorField U = _U();
//     IOMRFZoneList& MRF = _MRF();
//     singlePhaseTransportModel& laminarTransport = _laminarTransport();
//     instantList Times = runTime.times();
//     runTime.setEndTime(finalTime);
//     // Perform a TruthSolve
//     runTime.setTime(Times[1], 1);
//     runTime.setDeltaT(timeStep);
//     nextWrite = startTime;
//     // Initialize Nsnapshots
//     int nsnapshots = 0;

//     // Start the time loop
//     while (runTime.run())
//     {
// #include "readTimeControls.H"
// #include "CourantNo.H"
// #include "setDeltaT.H"
//         runTime.setEndTime(finalTime + timeStep);
//         Info << "Time = " << runTime.timeName() << nl << endl;

//         // --- Pressure-velocity PIMPLE corrector loop
//         while (pimple.loop())
//         {

// #include "alphaControls.H"
// #include "alphaEqnSubCycle.H"

// mixture.correct();

// #include "UEqn.H"

//             // --- Pressure corrector loop
//             while (pimple.correct())
//             {
// #include "pEqn.H"
//             }

//             if (pimple.turbCorr())
//             {
//                 //laminarTransport.correct(); //not present in interFoam (?)
//                 turbulence->correct();
//             }
//         }

//         Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
//              << "  ClockTime = " << runTime.elapsedClockTime() << " s"
//              << nl << endl;

//         if (checkWrite(runTime))
//         {
//             nsnapshots += 1;
//             ITHACAstream::exportSolution(U, name(counter), "./ITHACAoutput/Offline/");
//             ITHACAstream::exportSolution(p, name(counter), "./ITHACAoutput/Offline/");
//             std::ofstream of("./ITHACAoutput/Offline/" + name(counter) + "/" +
//                              runTime.timeName());
//             Ufield.append(U);
//             Pfield.append(p);
//             counter++;
//             nextWrite += writeEvery;
//             writeMu(mu_now);
//             // --- Fill in the mu_samples with parameters (time, mu) to be used for the PODI sample points
//             mu_samples.conservativeResize(mu_samples.rows() + 1, mu_now.size() + 1);
//             mu_samples(mu_samples.rows() - 1, 0) = atof(runTime.timeName().c_str());

//             for (int i = 0; i < mu_now.size(); i++)
//             {
//                 mu_samples(mu_samples.rows() - 1, i + 1) = mu_now[i];
//             }
//         }

//         runTime++;
//     }

//     // Resize to Unitary if not initialized by user (i.e. non-parametric problem)
//     if (mu.cols() == 0)
//     {
//         mu.resize(1, 1);
//     }

//     if (mu_samples.rows() == nsnapshots * mu.cols())
//     {
//         ITHACAstream::exportMatrix(mu_samples, "mu_samples", "eigen",
//                                    "./ITHACAoutput/Offline");
//     }
}
