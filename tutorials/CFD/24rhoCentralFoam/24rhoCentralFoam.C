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
    Example of a compressible flow Reduction Problem
SourceFiles
    24rhoCentralFoam.C
\*---------------------------------------------------------------------------*/

#include "UnsteadyCompressibleNS.H"
#include "ITHACAPOD.H"

/// \brief Class where the tutorial number 24 is implemented.
/// \details It is a child of the UnsteadyCompressibleNS class and some of its
/// functions are overridden to be adapted to the specific case.
class tutorial24: public UnsteadyCompressibleNS
{
    public:
        explicit tutorial24(int argc, char* argv[])
            :
            UnsteadyCompressibleNS(argc, argv),
            U(_U()),
            rho(_rho())
        {}

        // Fields to perform
        volScalarField& rho;
        volVectorField& U;

        psiThermo& thermo = _pThermo();

        volScalarField& T = thermo.T();
        volScalarField& p = thermo.p();

        void offlineSolve(word folder = "./ITHACAoutput/Offline/")
        {
            if (offline)
            {
                ITHACAstream::read_fields(Pfield, p, folder);
                ITHACAstream::read_fields(Tfield, T, folder);
                ITHACAstream::read_fields(Ufield, U, folder);
                ITHACAstream::read_fields(rhofield, rho, folder);
            }
            else
            {
                truthSolve();
            }
        }
};


int main(int argc, char* argv[])
{
    // Create the example object of the tutorial24 type
    tutorial24 example(argc, argv);
    // Read some parameters from file
    ITHACAparameters* para = ITHACAparameters::getInstance(example._mesh(),
                             example._runTime());
    int nModes = para->ITHACAdict->lookupOrDefault<int>("nModes", 10);
    example.startTime = para->ITHACAdict->lookupOrDefault<scalar>("startTime", 0);
    example.finalTime = para->ITHACAdict->lookupOrDefault<scalar>("finalTime", 7e-03);
    example.timeStep  = para->ITHACAdict->lookupOrDefault<scalar>("timeStep", 1e-06);
    example.writeEvery = para->ITHACAdict->lookupOrDefault<scalar>("writeEvery", 5e-05);
    // Perform an Offline Solve
    example.offlineSolve();
    // Get Modes
    ITHACAPOD::getModes(example.Pfield, example.Pmodes, example.p.name(),
                        example.podex, 0, 0,
                        nModes);
    ITHACAPOD::getModes(example.Tfield, example.Tmodes, example.T.name(),
                        example.podex, 0, 0,
                        nModes);
    ITHACAPOD::getModes(example.Ufield, example.Umodes, example._U().name(),
                        example.podex, 0, 0,
                        nModes);
    ITHACAPOD::getModes(example.rhofield, example.rhomodes, example._rho().name(),
                        example.podex, 0, 0,
                        nModes);
    // Project modes and reconstruct solution
    PtrList<volScalarField> projectedSnapshotsP, projectedSnapshotsT, projectedSnapshotsRho;
    PtrList<volVectorField> projectedSnapshotsU;

    example.Pmodes.projectSnapshots(example.Pfield, projectedSnapshotsP, nModes);
    example.Tmodes.projectSnapshots(example.Tfield, projectedSnapshotsT, nModes);
    example.Umodes.projectSnapshots(example.Ufield, projectedSnapshotsU, nModes);
    example.rhomodes.projectSnapshots(example.rhofield, projectedSnapshotsRho, nModes);

    ITHACAstream::exportFields(projectedSnapshotsP, "./ITHACAoutput/Reconstruction", "p");
    ITHACAstream::exportFields(projectedSnapshotsT, "./ITHACAoutput/Reconstruction", "T");
    ITHACAstream::exportFields(projectedSnapshotsU, "./ITHACAoutput/Reconstruction", "U");
    ITHACAstream::exportFields(projectedSnapshotsRho, "./ITHACAoutput/Reconstruction", "rho");

    // Compute the error on the testing set (velocity)
    Eigen::MatrixXd error = ITHACAutilities::errorL2Rel(example.Ufield,
                            projectedSnapshotsU);
}
