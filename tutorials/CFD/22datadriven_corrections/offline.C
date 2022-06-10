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
    Example of NS-Stokes Reduction Problem

\*---------------------------------------------------------------------------*/
#include "steadyNS.H"
#include "UnsteadyNSTurb.H"
#include <math.h>

#include "ITHACAutilities.H"
#include "reductionProblem.H"

#include "ReducedSteadyNS.H"
#include "ReducedUnsteadyNSTurb.H"
#include "ReducedUnsteadyNS.H"

#include <Eigen/Dense>
#include <Eigen/Core>
#include "forces.H"
#include "IOmanip.H"
#include "forces.H"
#include "forceCoeffs.H"
#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include "Foam2Eigen.H"
#include <chrono>
#include "ITHACAstream.H"
#include "ITHACAPOD.H"
#include "EigenFunctions.H"
#include <chrono>
#include <Eigen/SVD>
#include <Eigen/SparseLU>

#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>
#include <algorithm>
#include <fstream>
#include <string>
#include <stdexcept>
#include <sstream>
#include <vector>
#include <cstdio>
#include <typeinfo>
#include <iostream>
#include <cassert>
#include <zlib.h>


class tutorial22 : public UnsteadyNSTurb
{
    public:
        explicit tutorial22(int argc, char *argv[])
            :
            UnsteadyNSTurb(argc, argv),
            U(_U()),
            p(_p()),
            nut(_nut())
        {}

        // Relevant Fields
        volVectorField& U;
        volScalarField& p;
        volScalarField& nut;

        void offlineSolve(std::string offlinepath)
        {
            Vector<double> inl(0, 0, 0);
            List<scalar> mu_now(1);
            mu_now[0] =  1e-05;

            if ((offline) && (ITHACAutilities::check_folder(offlinepath) == true))
            {
                ITHACAstream::read_fields(Ufield, U, offlinepath);
                ITHACAstream::read_fields(Pfield, p, offlinepath);
                ITHACAstream::read_fields(nutFields, nut, offlinepath);
            }

            else
            {
                truthSolve(mu_now, offlinepath);
            }
        }
        Eigen::MatrixXd vectorTensorMult(Eigen::VectorXd g, Eigen::Tensor<double, 3> c ,
                                         Eigen::VectorXd a)
        {
            int prodDim = c.dimension(0);
            Eigen::MatrixXd prod;
            prod.resize(prodDim, 1);

            for (int i = 0; i < prodDim; i++)
            {
                prod(i, 0) = g.transpose() *
                             Eigen::SliceFromTensor(c, 0, i) * a;
            }

            return prod;
        }
};


void supremizer_approach(tutorial22& example);
void poisson_approach(tutorial22& example);

/*---------------------------------------------------------------------------*\
                               Starting the MAIN
\*---------------------------------------------------------------------------*/

int main(int argc, char *argv[])
{
    if (argc == 1)
    {
        std::cout << "Pass 'supremizer' or 'poisson' as first arguments."
                  << std::endl;
        exit(0);
    }

    // process arguments removing "offline" or "online" keywords
    int argc_proc = argc - 1;
    char* argv_proc[argc_proc];
    argv_proc[0] = argv[0];

    if (argc > 2)
    {
        std::copy(argv + 2, argv + argc, argv_proc + 1);
    }

    argc--;
    tutorial22 example(argc, argv);
    if (std::strcmp(argv[1], "supremizer") == 0)
    {
        // perform the offline stage, extracting the modes from the snapshots' dataset corresponding to parOffline
        supremizer_approach(example);
    }

    else if (std::strcmp(argv[1], "poisson") == 0)
    {
        poisson_approach(example);
    }
    else
    {
        std::cout << "Pass supremizer, poisson" << std::endl;
    }

    exit(0);
}

void supremizer_approach(tutorial22& example)

{
    word filename("./par");
    Eigen::VectorXd par;
    example.inletIndex.resize(1, 2);
    example.inletIndex << 0, 0;
    example.inletIndexT.resize(1, 1);
    example.inletIndexT << 1;
    ITHACAparameters* para = ITHACAparameters::getInstance(example._mesh(),
                             example._runTime());
    int NmodesU = para->ITHACAdict->lookupOrDefault<int>("NmodesU", 5);
    int NmodesP = para->ITHACAdict->lookupOrDefault<int>("NmodesP", 5);
    int NmodesSUP = para->ITHACAdict->lookupOrDefault<int>("NmodesSUP", 5);
    int NmodesNUT = para->ITHACAdict->lookupOrDefault<int>("NmodesNUT", 5);
    int NmodesProject = para->ITHACAdict->lookupOrDefault<int>("NmodesProject", 5);
    int NmodesMatrixRec = para->ITHACAdict->lookupOrDefault<int>("NmodesMatrixRec",
                          5);
    double penaltyFactor =
        para->ITHACAdict->lookupOrDefault<double>("penaltyFactor", 5);
    double U_BC = para->ITHACAdict->lookupOrDefault<double>("U_BC", 0.001);
    double romStartTime = para->ITHACAdict->lookupOrDefault<double>("romStartTime",
                          0);
    double romEndTime = para->ITHACAdict->lookupOrDefault<double>("romEndTime", 3);
    double romTimeStep = para->ITHACAdict->lookupOrDefault<double>("romTimeStep",
                         0.001);
    double e = para->ITHACAdict->lookupOrDefault<double>("RBFradius", 1);
    example.startTime = 79.992;
    example.finalTime = 99.996;
    example.timeStep = 0.0002;
    example.writeEvery = 0.004;
    example.offlineSolve("./ITHACAoutput/Offline/");
    example.solvesupremizer();

    ITHACAPOD::getModes(example.nutFields, example.nutModes, example._nut().name(),
                        example.podex, 0, 0, NmodesProject);
    ITHACAPOD::getModes(example.Ufield, example.Umodes, example._U().name(),
                        example.podex, 0, 0, NmodesProject);
    ITHACAPOD::getModes(example.Pfield, example.Pmodes, example._p().name(),
                        example.podex, 0, 0, NmodesProject);
    ITHACAPOD::getModes(example.supfield, example.supmodes, example._U().name(),
                        example.podex,
                        example.supex, 1, NmodesProject);
    example.projectSUP("./Matrices", NmodesU, NmodesP, NmodesSUP, NmodesNUT);

    Eigen::MatrixXd coeefs = ITHACAutilities::getCoeffs(example.Ufield, example.L_U_SUPmodes);
    Eigen::MatrixXd coeefsNut = ITHACAutilities::getCoeffs(example.nutFields, example.nutModes);
    Eigen::MatrixXd coeefsP = ITHACAutilities::getCoeffs(example.Pfield, example.Pmodes);

    cnpy::save(coeefs, "./ITHACAoutput/Matrices/coeefs.npy");
    cnpy::save(coeefsNut, "./ITHACAoutput/Matrices/coeefsNut.npy");
    cnpy::save(coeefsP, "./ITHACAoutput/Matrices/coeefsP.npy");
}


void poisson_approach(tutorial22& example)
{
    word filename("./par");
    Eigen::VectorXd par;
    example.inletIndex.resize(1, 2);
    example.inletIndex << 0, 0;
    example.inletIndexT.resize(1, 1);
    example.inletIndexT << 1;
    ITHACAparameters* para = ITHACAparameters::getInstance(example._mesh(),
                             example._runTime());
    int NmodesU = para->ITHACAdict->lookupOrDefault<int>("NmodesU", 5);
    int NmodesP = para->ITHACAdict->lookupOrDefault<int>("NmodesP", 5);
    int NmodesSUP = para->ITHACAdict->lookupOrDefault<int>("NmodesSUP", 5);
    int NmodesNUT = para->ITHACAdict->lookupOrDefault<int>("NmodesNUT", 5);
    int NmodesProject = para->ITHACAdict->lookupOrDefault<int>("NmodesProject", 5);
    int NmodesMatrixRec = para->ITHACAdict->lookupOrDefault<int>("NmodesMatrixRec",
                          5);
    double penaltyFactor =
        para->ITHACAdict->lookupOrDefault<double>("penaltyFactor", 5);
    double U_BC = para->ITHACAdict->lookupOrDefault<double>("U_BC", 0.001);
    double romStartTime = para->ITHACAdict->lookupOrDefault<double>("romStartTime",
                          0);
    double romEndTime = para->ITHACAdict->lookupOrDefault<double>("romEndTime", 3);
    double romTimeStep = para->ITHACAdict->lookupOrDefault<double>("romTimeStep",
                         0.001);
    double e = para->ITHACAdict->lookupOrDefault<double>("RBFradius", 1);
    example.startTime = 79.992;
    example.finalTime = 99.996;
    example.timeStep = 0.0002;
    example.writeEvery = 0.004;
    example.offlineSolve("./ITHACAoutput/Offline");
    example.solvesupremizer();

    //modNut = ITHACAPOD::getModes(example.nutFields, example.nutModes, example._nut().name(),
    //                    example.podex, 0, 0, NmodesProject);
    ITHACAPOD::getModes(example.Ufield, example.Umodes, example._U().name(),
                        example.podex, 0, 0, NmodesProject);
    ITHACAPOD::getModes(example.Pfield, example.Pmodes, example._p().name(),
                        example.podex, 0, 0, NmodesProject);
    ITHACAPOD::getModes(example.supfield, example.supmodes, example._U().name(),
                        example.podex,
                        example.supex, 1, NmodesProject);
    
    example.projectPPE("./Matrices", NmodesU, NmodesP, NmodesSUP, NmodesNUT);

    Eigen::MatrixXd coeefs = ITHACAutilities::getCoeffs(example.Ufield, example.L_U_SUPmodes);
    Eigen::MatrixXd coeefsNut = ITHACAutilities::getCoeffs(example.nutFields, example.nutModes);
    Eigen::MatrixXd coeefsP = ITHACAutilities::getCoeffs(example.Pfield, example.Pmodes);

    cnpy::save(coeefs, "./ITHACAoutput/Matrices/coeefs.npy");
    cnpy::save(coeefsNut, "./ITHACAoutput/Matrices/coeefsNut.npy");
    cnpy::save(coeefsP, "./ITHACAoutput/Matrices/coeefsP.npy");

}

