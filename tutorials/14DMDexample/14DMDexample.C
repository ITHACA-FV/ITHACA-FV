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
Description
    Example of an unsteady NS problem reproduced with DMD
SourceFiles
    14DMDexample.C
\*---------------------------------------------------------------------------*/

#include "unsteadyNS.H"
#include "ITHACAPOD.H"
#include "ITHACADMD.H"
#include "ReducedUnsteadyNS.H"
#include "ITHACAstream.H"
#include <chrono>
#include<math.h>
#include<iomanip>
#include "redsvd"

class tutorial14: public unsteadyNS
{
    public:
        explicit tutorial14(int argc, char* argv[])
            :
            unsteadyNS(argc, argv),
            U(_U()),
            p(_p())
        {}

        // Fields To Perform
        volVectorField& U;
        volScalarField& p;

        void offlineSolve()
        {
            List<scalar> mu_now(1);

            if (offline)
            {
                ITHACAstream::read_fields(Ufield, U, "./ITHACAoutput/Offline/");
                ITHACAstream::read_fields(Pfield, p, "./ITHACAoutput/Offline/");
            }
            else
            {
                for (label i = 0; i < mu.cols(); i++)
                {
                    mu_now[0] = mu(0, i);
                    change_viscosity( mu(0, i));
                    truthSolve(mu_now);
                }
            }
        }
};

/*---------------------------------------------------------------------------*\
                               Starting the MAIN
\*---------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
    // Construct the tutorial04 object
    tutorial14 example(argc, argv);
    // Read parameters from ITHACAdict file
    ITHACAparameters* para = ITHACAparameters::getInstance(example._mesh(),
                             example._runTime());
    // Read some input parameters for DMD
    double dtDMD = readScalar(para->ITHACAdict->lookup("dtDMD"));
    double finalTimeDMD = readScalar(para->ITHACAdict->lookup("finalTimeDMD"));
    double startTimeDMD = readScalar(para->ITHACAdict->lookup("startTimeDMD"));
    int numberOfModesDMD = readInt(para->ITHACAdict->lookup("numberOfModesDMD"));
    bool exactDMD = readBool(para->ITHACAdict->lookup("exactDMD"));
    bool exportDMDModes = readBool(para->ITHACAdict->lookup("exportDMDModes"));
    word exportFolder = string(para->ITHACAdict->lookup("exportFolder"));
    /// Set the number of parameters
    example.Pnumber = 1;
    /// Set samples (for DMD you set only one sample)
    example.Tnumber = 1;
    /// Set the parameters infos
    example.setParameters();
    // Set the parameter ranges
    example.mu_range(0, 0) = 0.01;
    example.mu_range(0, 1) = 0.01;
    // Generate equispaced samples inside the parameter range
    example.genEquiPar();
    // Time parameters
    example.startTime = 60;
    example.finalTime = 85;
    example.timeStep = 0.01;
    example.writeEvery = 0.1;
    // Read the snapshots from Offline folder or compute new snapshots
    example.offlineSolve();
    // Create The DMD class for the velocity field
    word exportFieldNameU = example.Ufield[0].name() + "_" + name(numberOfModesDMD);
    ITHACADMD<vector> DMDv(example.Ufield, example.writeEvery);
    DMDv.getModes(numberOfModesDMD, exactDMD, exportDMDModes);
    DMDv.getDynamics(startTimeDMD, finalTimeDMD, dtDMD);
    DMDv.reconstruct(exportFolder, exportFieldNameU);
    // Create The DMD class for the pressure field
    word exportFieldNameP = example.Pfield[0].name() + "_" + name(numberOfModesDMD);
    ITHACADMD<scalar> DMDp(example.Pfield, example.writeEvery);
    DMDp.getModes(numberOfModesDMD, exactDMD, exportDMDModes);
    DMDp.getDynamics(startTimeDMD, finalTimeDMD, dtDMD);
    DMDp.reconstruct(exportFolder, exportFieldNameP);
    return 0;
}

