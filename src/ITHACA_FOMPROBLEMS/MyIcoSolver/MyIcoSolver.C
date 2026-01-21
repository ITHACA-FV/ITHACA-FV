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

#include "MyIcoSolver.H"

/// \file
/// Source file of the MyIcoSolver class.

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Construct Null
MyIcoSolver::MyIcoSolver() {}

// Construct from zero
MyIcoSolver::MyIcoSolver(int argc, char* argv[])
    :
    UnsteadyProblem()
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
    _piso = autoPtr<pisoControl>
            (
                new pisoControl
                (
                    mesh
                )
            );
    ITHACAdict = new IOdictionary
    (
        IOobject
        (
            "ITHACAdict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
#include "createFields.H"
#include "createFvOptions.H"
    para = ITHACAparameters::getInstance(mesh, runTime);
    offline = ITHACAutilities::check_off();
    podex = ITHACAutilities::check_pod();
    //supex = ITHACAutilities::check_sup();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void MyIcoSolver::truthSolve(List<scalar> mu_now, fileName folder)
{
    Time& runTime = _runTime();
    surfaceScalarField& phi = _phi();
    fvMesh& mesh = _mesh();
#include "initContinuityErrs.H"
    fv::options& fvOptions = _fvOptions();
    pisoControl& piso = _piso();
    volScalarField& p = _p();
    volVectorField& U = _U();
    IOMRFZoneList& MRF = _MRF();
    singlePhaseTransportModel& laminarTransport = _laminarTransport();
    instantList Times = runTime.times();
    runTime.setEndTime(finalTime);
    // Perform a TruthSolve
    runTime.setTime(Times[1], 1);
    runTime.setDeltaT(timeStep);
    nextWrite = startTime;
    /*
        // Set time-dependent velocity BCs for initial condition
        if (timedepbcMethod == "yes")
        {
            for (label i = 0; i < inletPatch.rows(); i++)
            {
                Vector<double> inl(0, 0, 0);

                for (label j = 0; j < inl.size(); j++)
                {
                    inl[j] = timeBCoff(i* inl.size() + j, 0);
                }

                assignBC(U, inletPatch(i, 0), inl);
            }
        }

        // Export and store the initial conditions for velocity and pressure
        ITHACAstream::exportSolution(U, name(counter), folder);
        ITHACAstream::exportSolution(p, name(counter), folder);
        std::ofstream of(folder + name(counter) + "/" +
                         runTime.timeName());
        Ufield.append(U.clone());
        Pfield.append(p.clone());
        counter++;
        nextWrite += writeEvery;*/

    // Start the time loop
    while (runTime.run())
    {
#include "readTimeControls.H"
#include "CourantNo.H"
#include "setDeltaT.H"
        runTime.setEndTime(finalTime);
        runTime++;
        Info << "Time = " << runTime.timeName() << nl << endl;
        /*
        // Set time-dependent velocity BCs
        if (timedepbcMethod == "yes")
        {
            for (label i = 0; i < inletPatch.rows(); i++)
            {
                Vector<double> inl(0, 0, 0);

                for (label j = 0; j < inl.size(); j++)
                {
                    inl[j] = timeBCoff(i* inl.size() + j, counter2);
                }

                assignBC(U, inletPatch(i, 0), inl);
            }

            counter2 ++;
        }*/

        // --- Pressure-velocity PIMPLE corrector loop
        while (piso.loop())
        {
#include "UEqn.H"

            // --- Pressure corrector loop
            while (piso.correct())
            {
#include "pEqn.H"
            }

            if (piso.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;

        if (checkWrite(runTime))
        {
            ITHACAstream::exportSolution(U, name(counter), folder);
            ITHACAstream::exportSolution(p, name(counter), folder);
            Ufield.append(U.clone());
            Pfield.append(p.clone());
            counter++;
            nextWrite += writeEvery;
            /*
            writeMu(mu_now);
            // --- Fill in the mu_samples with parameters (time, mu) to be used for the PODI sample points
            mu_samples.conservativeResize(mu_samples.rows() + 1, mu_now.size() + 1);
            mu_samples(mu_samples.rows() - 1, 0) = atof(runTime.timeName().c_str());

            for (label i = 0; i < mu_now.size(); i++)
            {
                mu_samples(mu_samples.rows() - 1, i + 1) = mu_now[i];
            }*/
        }
    }

    /*
    // Resize to Unitary if not initialized by user (i.e. non-parametric problem)
    if (mu.cols() == 0)
    {
        mu.resize(1, 1);
    }

    if (mu_samples.rows() == counter* mu.cols())
    {
        ITHACAstream::exportMatrix(mu_samples, "mu_samples", "eigen",
                                   folder);
    }*/
}

void MyIcoSolver::restart()
{
    //_runTime().objectRegistry::clear();
    //_mesh().objectRegistry::clear();
    // _mesh.clear();
    // _runTime.clear();
    //_simple.clear();
    _p.clear();
    _U.clear();
    _phi.clear();
    turbulence.clear();
    _fvOptions.clear();
    argList& args = _args();
    Time& runTime = _runTime();
    runTime.setTime(0, 1);
    Foam::fvMesh& mesh = _mesh();
    pisoControl& piso = _piso();
    Info << "ReReading field p\n" << endl;
    _p = autoPtr<volScalarField>
         (
             new volScalarField
             (
                 IOobject
                 (
                     "p",
                     runTime.timeName(),
                     mesh,
                     IOobject::MUST_READ,
                     IOobject::AUTO_WRITE
                 ),
                 mesh
             )
         );
    volScalarField& p = _p();
    Info << "ReReading field U\n" << endl;
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
    Info << "ReReading/calculating face flux field phi\n" << endl;
    _phi = autoPtr<surfaceScalarField>
           (
               new surfaceScalarField
               (
                   IOobject
                   (
                       "phi",
                       runTime.timeName(),
                       mesh,
                       IOobject::READ_IF_PRESENT,
                       IOobject::AUTO_WRITE
                   ),
                   linearInterpolate(U) & mesh.Sf()
               )
           );
    surfaceScalarField& phi = _phi();
    pRefCell = 0;
    pRefValue = 0.0;
    setRefCell(p, piso.dict(), pRefCell, pRefValue);
    _laminarTransport = autoPtr<singlePhaseTransportModel>
                        (
                            new singlePhaseTransportModel( U, phi )
                        );
    singlePhaseTransportModel& laminarTransport = _laminarTransport();
    turbulence = autoPtr<incompressible::turbulenceModel>
                 (
                     incompressible::turbulenceModel::New(U, phi, laminarTransport)
                 );
    _MRF = autoPtr<IOMRFZoneList>
           (
               new IOMRFZoneList(mesh)
           );
    _fvOptions = autoPtr<fv::options>(new fv::options(mesh));
    turbulence->validate();
}
