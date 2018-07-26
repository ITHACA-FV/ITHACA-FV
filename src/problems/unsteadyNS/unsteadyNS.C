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

#include "unsteadyNS.H"

/// \file
/// Source file of the unsteadyNS class.

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Construct Null
unsteadyNS::unsteadyNS() {}

// Construct from zero
unsteadyNS::unsteadyNS(int argc, char *argv[])
{
#include "setRootCase.H"
#include "createTime.H"
#include "createMesh.H"
		_pimple = autoPtr<pimpleControl>
			  (
				  new pimpleControl
				  (
					  mesh
				  )
			  );
#include "createFields.H"
#include "createFvOptions.H"

	supex = ITHACAutilities::check_sup();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void unsteadyNS::truthSolve(List<scalar> mu_now)
{
	#include "initContinuityErrs.H"
	Time& runTime = _runTime();
	surfaceScalarField& phi = _phi();
	fvMesh& mesh = _mesh();
	fv::options& fvOptions = _fvOptions();
	pimpleControl& pimple = _pimple();
	volScalarField p = _p();
	volVectorField U = _U();
	IOMRFZoneList& MRF = _MRF();
	singlePhaseTransportModel& laminarTransport = _laminarTransport();

	instantList Times = runTime.times();
	runTime.setEndTime(finalTime);
	// Perform a TruthSolve
	runTime.setTime(Times[1], 1);
	runTime.setDeltaT(timeStep);
	nextWrite = startTime;

	// Start the time loop
	while (runTime.run())
	{
#include "readTimeControls.H"
#include "CourantNo.H"
#include "setDeltaT.H"
		runTime.setEndTime(finalTime+timeStep);
		Info << "Time = " << runTime.timeName() << nl << endl;

		// --- Pressure-velocity PIMPLE corrector loop
		while (pimple.loop())
		{
#include "UEqn.H"
			// --- Pressure corrector loop
			while (pimple.correct())
			{
#include "pEqn.H"
			}
			if (pimple.turbCorr())
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
			exportSolution(U, name(counter), "./ITHACAoutput/Offline/");
			exportSolution(p, name(counter), "./ITHACAoutput/Offline/");
			std::ofstream of("./ITHACAoutput/Offline/"+name(counter)+"/"+runTime.timeName());
			Ufield.append(U);
			Pfield.append(p);
			counter++;
			nextWrite += writeEvery;
			writeMu(mu_now);
		}
		runTime++;
	}
}

bool unsteadyNS::checkWrite(Time& timeObject)
{
	scalar diffnow = mag(nextWrite - atof(timeObject.timeName().c_str()));

	scalar diffnext = mag(nextWrite - atof(timeObject.timeName().c_str()) - timeObject.deltaTValue());
	if ( diffnow < diffnext)
	{
		return true;
	}
	else
	{
		return false;
	}
}




