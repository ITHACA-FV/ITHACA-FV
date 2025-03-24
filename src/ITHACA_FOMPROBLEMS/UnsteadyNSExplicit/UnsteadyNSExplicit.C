/*---------------------------------------------------------------------------*\
     ██╗████████╗██╗  ██╗ █████╗  ██████╗ █████╗       ███████╗██╗   ██╗
     ██║╚══██╔══╝██║  ██║██╔══██╗██╔════╝██╔══██╗      ██╔════╝██║   ██║
     ██║   ██║   ███████║███████║██║     ███████║█████╗█████╗  ██║   ██║
     ██║   ██║   ██╔══██║██╔══██║██║     ██╔══██║╚════╝██╔══╝  ╚██╗ ██╔╝
     ██║   ██║   ██║  ██║██║  ██║╚██████╗██║  ██║      ██║      ╚████╔╝
     ╚═╝   ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝      ╚═╝       ╚═══╝

 * In real Time Highly Advanced Computational Applications for Finite Volumes
 * Copyright (C) 2020 by the ITHACA-FV authors
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

#include "UnsteadyNSExplicit.H"
#include "fvCFD.H"

/// \file
/// Source file of the UnsteadyNSExplicit class.

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Construct Null
UnsteadyNSExplicit::UnsteadyNSExplicit() {}

// Construct from zero
UnsteadyNSExplicit::UnsteadyNSExplicit(int argc, char* argv[])
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
    para = ITHACAparameters::getInstance(mesh, runTime);
    bcMethod = ITHACAdict->lookupOrDefault<word>("bcMethod", "none");
    M_Assert(bcMethod == "lift" || bcMethod == "penalty" || bcMethod == "none",
             "The BC method must be set to lift or penalty or none in ITHACAdict");
    fluxMethod = ITHACAdict->lookupOrDefault<word>("fluxMethod", "inconsistent");
    M_Assert(fluxMethod == "inconsistent" || fluxMethod == "consistent",
             "The flux method must be set to inconsistent or consistent in ITHACAdict");
    offline = ITHACAutilities::check_off();
    podex = ITHACAutilities::check_pod();
    supex = ITHACAutilities::check_sup();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void UnsteadyNSExplicit::truthSolve(List<scalar> mu_now, fileName folder)
{
    Time& runTime = _runTime();
    fvMesh& mesh = _mesh();
    volScalarField& p = _p();
    volVectorField& U = _U();
    surfaceScalarField& phi = _phi();

    if (fluxMethod == "inconsistent")
    {
        phi = fvc::flux(U);
    }

#include "initContinuityErrs.H"
    dimensionedScalar nu = _nu() * mu_now[0];
    dimensionedScalar dt = timeStep * _dt();
    IOMRFZoneList& MRF = _MRF();
    singlePhaseTransportModel& laminarTransport = _laminarTransport();
    instantList Times = runTime.times();
    runTime.setEndTime(finalTime);
    // Perform a TruthSolve
    runTime.setTime(Times[1], 1);
    runTime.setDeltaT(timeStep);
    nextWrite = startTime;
    // Export and store the initial conditions for velocity, pressure and flux
    ITHACAstream::exportSolution(U, name(counter), folder);
    ITHACAstream::exportSolution(p, name(counter), folder);
    ITHACAstream::exportSolution(phi, name(counter), folder);
    std::ofstream of(folder + name(counter) + "/" +
                     runTime.timeName());
    Ufield.append(U.clone());
    Pfield.append(p.clone());
    Phifield.append(phi.clone());
    counter++;
    nextWrite += writeEvery;

    // Start the time loop
    while (runTime.run())
    {
#include "readTimeControls.H"
#include "CourantNo.H"
#include "setDeltaT.H"
        runTime.setEndTime(finalTime);
        runTime++;
        Info << "Time = " << runTime.timeName() << nl << endl;

        if (fluxMethod == "inconsistent")
        {
#include "IFM.H"
        }
        else if (fluxMethod == "consistent")
        {
#include "CFM.H"
        }

#include "initContinuityErrs.H"
        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;

        if (checkWrite(runTime))
        {
            ITHACAstream::exportSolution(U, name(counter), folder);
            ITHACAstream::exportSolution(p, name(counter), folder);
            ITHACAstream::exportSolution(phi, name(counter), folder);
            std::ofstream of(folder + name(counter) + "/" +
                             runTime.timeName());
            Ufield.append(U.clone());
            Pfield.append(p.clone());
            Phifield.append(phi.clone());
            counter++;
            nextWrite += writeEvery;
        }
    }
}






