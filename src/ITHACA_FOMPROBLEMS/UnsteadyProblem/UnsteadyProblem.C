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

#include "UnsteadyProblem.H"

/// \file
/// Source file of the UnsteadyProblem class.
///

bool UnsteadyProblem::checkWrite(Time& timeObject)
{
    scalar diffnow = mag(nextWrite - atof(timeObject.timeName().c_str()));
    scalar diffnext = mag(nextWrite - atof(timeObject.timeName().c_str()) -
                          timeObject.deltaTValue());

    if ( diffnow < diffnext)
    {
        return true;
    }
    else
    {
        return false;
    }
}

void UnsteadyProblem::setTimes(Time& timeObject)
{
    ITHACAparameters* para(ITHACAparameters::getInstance());
    startTime = para->ITHACAdict->lookupOrDefault<scalar>("startTime",
                timeObject.startTime().value());
    finalTime = para->ITHACAdict->lookupOrDefault<scalar>("finalTime",
                timeObject.endTime().value());
    timeStep = para->ITHACAdict->lookupOrDefault<scalar>("timeStep",
               timeObject.deltaT().value());
    writeEvery = para->ITHACAdict->lookupOrDefault<scalar>("writeEvery", timeStep);
    M_Assert(finalTime > startTime,
             "The finalTime needs to be bigger than the startTime");
    M_Assert(finalTime - startTime > timeStep,
             "The timeStep needs to be bigger than the entire simulation Time");
    M_Assert(writeEvery >= timeStep,
             "The writeEvery needs to larger or equal to the timeStep");
    timeObject.setEndTime(finalTime);
    timeObject.setDeltaT(timeStep);
}
