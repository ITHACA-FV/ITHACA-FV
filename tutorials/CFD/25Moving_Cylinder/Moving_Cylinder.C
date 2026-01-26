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
    Example of an unsteady NS Reduction Problem
SourceFiles
    22fsi.C
\*---------------------------------------------------------------------------*/

#include "fsiBasic.H"
#include "ITHACAPOD.H"
//#include "ReducedSimpleSteadyNS.H"
#include "ITHACAstream.H"
#include "dynamicFvMesh.H"
#include "ReducedProblem.H"
#include "ReducedFsi.H"
#include <chrono>
#include<math.h>
#include<iomanip>
#include "pointConstraints.H"
#include "mathematicalConstants.H"
//#include "zoneMotion.H"


class tutorial22: public fsiBasic
{
    public:
        tutorial22(int argc, char* argv[])
            : fsiBasic(argc, argv),
              U(_U()),
              p(_p()),
              pd(_pointDisplacement())
        {
            //point0 = meshPtr().points();
        }
        /// Velocity field
        volVectorField& U;
        /// Pressure field
        volScalarField& p;
        pointVectorField& pd;
        /// Initial coordinates of the grid points
        //vectorField point0;
        void offlineSolve(word folder = "./ITHACAoutput/Offline/")
        {
            if (offline )
            {
                ITHACAstream::readMiddleFields(Ufield, U, folder);
                ITHACAstream::readMiddleFields(Pfield, p, folder);
                ITHACAstream::readMiddleFields(Dfield, pd, folder);
                //ITHACAstream::readMiddleFields(Efield, E, folder);
                //mu_samples = ITHACAstream::readMatrix("./parsOff_mat.txt");
            }
            else if (!offline)
            {
                for (label i = 0; i < mu.rows(); i++)
                {
                    //std::clock_t startOff;
                    //startOff = std::clock();
                    //folderN = i;
                    std::cout << "Current mu = " << mu(i, 0) << std::endl;
                    //change_viscosity(mu(i, 0));
                    //change_stiffness(mu(i,0));
                    //meshPtr().movePoints(Mesh0->points());
                    word param_name = "damping";
                    updateStiffnessAndRebuildSolver(mu(i, 0), param_name);
                    //exit(0);
                    //assignIF(_U(), Uinl);
                    //folder = folder + "/" +  name(i+1);
                    startTime  = 0;
                    finalTime  = 30;
                    timeStep   = 0.003; //0.0025; //4e-03; // ok dt=0.001
                    writeEvery = 1e-01;
                    truthSolve(i, folder);
                    word localFolder = folder + "../" + "/DataFromFoam_" + name(i + 1);
                    prepareFoamData(localFolder);
                    restart();
                    /// Clear Data members
                    fomforcex.clear();
                    fomforcey.clear();
                    //centerofmassx.clear();
                    centerofmassy.clear();
                    //centerofmassz.clear();
                    //this->Ufield.last().checkIn(this->meshes[i]);
                    //this->Ufield.last().checkIn();
                    // ITHACAPOD::getModes(this->Ufield, this->Umodes, this->_U().name(),
                    //                      this->podex, 0, 0, 10);
                    //     exit(0);
                }
            }
        }
};




/*------------------------------------------------------------------------------------------*\
                               Starting the MAIN
\*------------------------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
    // Construct the tutorial object
    tutorial22 example(argc, argv);
    //std::clock_t startOff;
    //double durationOff;
    // Read some parameters from file
    ITHACAparameters* para = ITHACAparameters::getInstance(example.meshPtr(),
                             example._runTime());
    //Eigen::MatrixXd parOff;
    // Read the par file where the parameters are stored
    std::ifstream exFileOff("./parsOff_mat.txt");

    if (exFileOff)
    {
        example.mu  = ITHACAstream::readMatrix("./parsOff_mat.txt");
    }

    int NmodesUout = para->ITHACAdict->lookupOrDefault<int>("NmodesUout", 40);
    int NmodesPout = para->ITHACAdict->lookupOrDefault<int>("NmodesPout", 40);
    int NmodesDout = para->ITHACAdict->lookupOrDefault<int>("NmodesDout", 40);
    int NmodesUproj = para->ITHACAdict->lookupOrDefault<int>("NmodesUproj", 15);
    int NmodesPproj = para->ITHACAdict->lookupOrDefault<int>("NmodesPproj", 5);
    int NmodesDproj = para->ITHACAdict->lookupOrDefault<int>("NmodesDproj", 1);
    // //Perform the offline solve
    //startOff= std::clock();
    //Info << example.Mesh0->checkMesh(true) << endl;
    //exit(0);
    //std::cerr << "debug point 3" << std::endl;
    example.offlineSolve();
    //std::cerr << "debug point 4" << std::endl;
    //example.meshPtr().movePoints( example.Mesh0->points());

    // bool meshesDiffer = false;
    // Info << example.Mesh0->points()[50] << endl;
    // dynamicFvMesh& newmesh = example.meshPtr();
    // //Info << example.Ufield[0].mesh().points()[0] << endl;
    // pointField newPoints  = newmesh.points()+ example.Dfield[50];
    // newmesh.movePoints(newPoints );

    // Info << newPoints[50] << endl;

    // forAll(newPoints, i)
    // {
    //     if (newmesh.points()[i] != example.Mesh0->points()[i])  // or a tolerance like 1e-10
    //     {
    //         meshesDiffer = true;
    //         break;
    //     }
    // }
    // if (meshesDiffer)
    //     Info << "Meshes differ in point locations." << endl;
    // else
    // Info << "Meshes are geometrically identical." << endl;
    // //std::cout << "The Offline phase  duration  is  equal  to " << durationOff << std::endl;
    // exit(0);

    //example.meshPtr().movePoints(example.Dfield[0]);
    if (example.podex == 0 )
    {
        ITHACAPOD::getModes(example.Ufield, example.Umodes, example._U().name(),
                            example.podex, 0, 0, NmodesUout);
        //ITHACAPOD::getModes(UfieldEnrich, online.Umodes, online._U().name(),
        //            example.podex, 0, 0, NmodesUout);
        ITHACAPOD::getModes(example.Pfield, example.Pmodes, example._p().name(),
                            example.podex, 0, 0, NmodesPout);
        ITHACAPOD::getModes(example.Dfield, example.Dmodes,
                            example._pointDisplacement().name(),
                            example.podex, 0, 0, NmodesDout);
    }
    else
    {
        ITHACAstream::read_fields(example.Umodes, example._U(),
                                  "./ITHACAoutput/POD/");
        ITHACAstream::read_fields(example.Pmodes, example._p(),
                                  "./ITHACAoutput/POD/");
        ITHACAstream::read_fields(example.Dmodes, example._pointDisplacement(),
                                  "./ITHACAoutput/POD/");
    }

    example.loadCentreOfMassY();
    //Info << example.Dfield.size() << endl;
    /// Getting the POD coeffs modes for pointDisplacement
    example.coeffL2 = ITHACAutilities::getCoeffs(example.Dfield,
                      example.Dmodes,
                      NmodesDproj, false);
    //exit(0);
    Eigen::MatrixXd parsOn;
    std::ifstream exFileOn("./parsOn_mat.txt");

    if (exFileOn)
    {
        parsOn  = ITHACAstream::readMatrix("./parsOn_mat.txt");
    }

    /// Generate data for testing
    word test_folder = "./ITHACAoutput/TestingOff/";

    /// include Symbolinks
    if (!ITHACAutilities::check_folder(test_folder))
    {
        tutorial22 testkOff(argc, argv);
        testkOff.mu  = parsOn;
        //Perform the offline solve
        testkOff.offline = false;

        for (label k = 0; k < parsOn.rows(); k++)
        {
            //std::cout << "Current mu = " << parsOn(k, 0) << std::endl;
            //example.meshPtr().movePoints( example.Mesh0->points());
            word param_name = "damping";
            testkOff.updateStiffnessAndRebuildSolver(parsOn(k, 0), param_name);
            testkOff.startTime  = 0;
            testkOff.finalTime  = 30;
            testkOff.timeStep   = 0.003; //0.0025; //4e-03; // ok dt=0.001
            testkOff.writeEvery = 1e-01;
            testkOff.truthSolve(k, test_folder);
            word localFolder = test_folder + "/DataFromFoam_" + name(k + 1);
            testkOff.prepareFoamData(localFolder);
            testkOff.restart();
            testkOff.Ufield.clear();
            testkOff.Pfield.clear();
            testkOff.Dfield.clear();
            testkOff.fomforcex.clear();
            testkOff.fomforcey.clear();
            //testkOff.centerofmassx.clear();
            testkOff.centerofmassy.clear();
            //testkOff.centerofmassz.clear();
        }

        ITHACAutilities::createSymLink("./0", test_folder);
        ITHACAutilities::createSymLink("./system", test_folder);
        ITHACAutilities::createSymLink("./constant", test_folder);
        testkOff.offline = true;
    }

    /*-----------------------------------------------------------*\
                               Starting the online part
    \*------------------------------------------------------------*/
    /// ############### contruct the reduced the class object ###################
    ReducedFsi reduced(example);

    //reducedBasicFsi reduced(example)
    // //exit(0);
    // reduced.startTime = example.startTime;
    // reduced.finalTime = example.finalTime;
    // reduced.timeStep = example.timeStep;
    // reduced.writeEvery = example.writeEvery;
    //Perform the online solutions
    //example.restart();
    // std::clock_t startOn;
    // double durationOn;
    // startOn = std::clock();
    if (exFileOn)
    {
        parsOn  = ITHACAstream::readMatrix("./parsOn_mat.txt");
    }

    word folder = "./ITHACAoutput/Online/";
    ITHACAutilities::createSymLink("./0", folder);
    ITHACAutilities::createSymLink("./system", folder);
    ITHACAutilities::createSymLink("./constant", folder);

    for (label i = 0; i < parsOn.rows(); i++)
    {
        // std::cout << "Current mu = "
        //           << parsOn(i,0) << std::endl;
        //example.meshPtr().movePoints( example.Mesh0->points());
        word param_name = "damping";
        example.updateStiffnessAndRebuildSolver(parsOn(i, 0), param_name);
        reduced.startTime  = 0; //example.startTime;
        reduced.finalTime  = 30; //example.finalTime;;
        reduced.timeStep   = 0.003; //example.timeStep;
        reduced.writeEvery = 1e-1; //example.writeEvery;
        /// Solving the reduced problem
        reduced.solveOnline_Pimple(NmodesUproj,
                                   NmodesPproj,
                                   NmodesDproj,
                                   folder);
        word localFolder = folder +  name(i + 1);

        for (int k = 0; k < reduced.UredFields.size(); ++k)
        {
            ITHACAstream::exportSolution(reduced.UredFields[k], name(k + 1),
                                         localFolder);
            ITHACAstream::exportSolution(reduced.PredFields[k], name(k + 1),
                                         localFolder);
            ITHACAstream::exportSolution(reduced.Dfield[k],     name(k + 1),
                                         localFolder);
            ITHACAstream::writePoints(reduced.ListOfpoints[k], localFolder,
                                      name(k + 1) + "/polyMesh/");
        }

        word DataRom = folder + "../" + "/DataFromRom_" + name(i + 1);
        reduced.prepareRomData(DataRom);
        /// Restart the simulation
        example.restart();
        reduced.UredFields.clear();
        reduced.PredFields.clear();
        reduced.Dfield.clear();
        reduced.romforcey.clear();
        reduced.romforcex.clear();
        reduced.centerofmassy.clear();
    }

    //durationOn = (std::clock() - startOn);
    //std::cout << "The Online  phase  duration  is  equal  to " << durationOn << std::endl;
    exit(0);
    Eigen::MatrixXd errL2U = ITHACAutilities::errorL2Rel(example.Ufield,
                             reduced.UredFields);
    std::cout <<
    "======================= errL2U completed================================" <<
    "\n";
    Eigen::MatrixXd errL2P = ITHACAutilities::errorL2Rel(example.Pfield,
                             reduced.PredFields);
    cnpy::save(errL2U, "./ITHACAoutput/DataFromRom/errL2U_" + name(
                   NmodesUproj) + "_" + name(NmodesPproj) + ".npy");
    cnpy::save(errL2P, "./ITHACAoutput/DataFromRom/errL2P_" + name(
                   NmodesUproj) + "_" + name(NmodesPproj) + ".npy");
    exit(0);
}
