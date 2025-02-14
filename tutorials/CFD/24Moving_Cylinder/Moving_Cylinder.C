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
    Tutorial of flow around a moving cylinder
\*---------------------------------------------------------------------------*/

#include "fsiBasic.H"
#include "ITHACAPOD.H"
#include "ReducedSimpleSteadyNS.H"
#include "ITHACAstream.H"
#include "dynamicFvMesh.H"
#include "ReducedProblem.H"
#include "ReducedFsi.H"
#include <chrono>
#include<math.h>
#include<iomanip>
#include "pointConstraints.H"
#include "mathematicalConstants.H"


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
    void offlineSolve(word folder="./ITHACAoutput/Offline/")
    {
        List<scalar> mu_now(1);

        if (offline)
        {
            ITHACAstream::read_fields(Ufield, U, folder);
            ITHACAstream::read_fields(Pfield, p, folder);
            ITHACAstream::read_fields(Dfield, pd, folder);
             // mu_samples =
             //    ITHACAstream::readMatrix("./ITHACAoutput/Offline/mu_samples_mat.txt");
        }
        else
        {
            mu_now[0] = mu(0, 0); 
            truthSolve(mu_now, folder);
            //restart();

        }

    }

}; 
/*----------------------------------------------------------------------------------------------------------*\
                               Starting the MAIN
\*-----------------------------------------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
    // Construct the tutorial object
    tutorial22 example(argc, argv);
    tutorial22 online(argc,argv);
    std::clock_t startOff;
    double durationOff;
    // Read some parameters from file
    ITHACAparameters* para = ITHACAparameters::getInstance(example.meshPtr(),example._runTime());
    // int NmodesUout = para->ITHACAdict->lookupOrDefault<int>("NmodesUout", 41);
    // int NmodesPout = para->ITHACAdict->lookupOrDefault<int>("NmodesPout", 41);
    // int NmodesDout = para->ITHACAdict->lookupOrDefault<int>("NmodesDout", 41);

    int NmodesUout  =  readInt(para->ITHACAdict->lookup("NmodesUout"));
    int NmodesPout  =  readInt(para->ITHACAdict->lookup("NmodesPout"));
    int NmodesDout  =  readInt(para->ITHACAdict->lookup("NmodesDout"));

    int NmodesUproj  = readInt(para->ITHACAdict->lookup("NmodesUproj"));
    int NmodesPproj  = readInt(para->ITHACAdict->lookup("NmodesPproj"));
    int NmodesDproj  = readInt(para->ITHACAdict->lookup("NmodesDproj"));
    // int NmodesUproj = para->ITHACAdict->lookupOrDefault<int>("NmodesUproj", 8);
    // int NmodesPproj = para->ITHACAdict->lookupOrDefault<int>("NmodesPproj", 16);
    // int NmodesDproj = para->ITHACAdict->lookupOrDefault<int>("NmodesDproj", 1);

    // Read the par file where the parameters are stored
    word filename("./par");
    example.mu = ITHACAstream::readMatrix(filename);
    // Time parameters: We can use Ioodictionnary to access time parameters
    example.startTime  = 0;
    example.finalTime  = 5;
    example.timeStep   = 0.002; 
    example.writeEvery = 1e-01;

    // //Perform the offline solve
    startOff= std::clock();
    example.offlineSolve();
    
    durationOff = (std::clock() - startOff);
    std::cout << "The Offline phase  duration  is  equal  to " << durationOff << std::endl;
    //exit(0);
    if(!ITHACAutilities::check_folder("./ITHACAoutput/DataFromFoam"))
    {

        mkDir("ITHACAoutput/DataFromFoam");
        Eigen::VectorXd fomforcex = Foam2Eigen::field2Eigen(example.fomforcex);
        cnpy::save(fomforcex, "./ITHACAoutput/DataFromFoam/fomforcex.npy");
        Eigen::VectorXd fomforcey = Foam2Eigen::field2Eigen(example.fomforcey);
        cnpy::save(fomforcey, "./ITHACAoutput/DataFromFoam/fomforcey.npy");

        Eigen::MatrixXd CentreOfMassY = Foam2Eigen::field2Eigen(example.centerofmassy);
        cnpy::save(CentreOfMassY, "./ITHACAoutput/DataFromFoam/CentreOfMassY.npy");
    }

    if(std::ifstream("./ITHACAoutput/DataFromFoam/CentreOfMassY.npy"))
    {
        Info << "################ Reading pitch ##############" << endl;
        cnpy::load(online.CylDispl, "./ITHACAoutput/DataFromFoam/CentreOfMassY.npy");
    }
     //std::cout << "example.podex = " << example.podex << std::endl;
 
    if(example.podex==0 )
    {
       ITHACAPOD::getModes(example.Ufield, online.Umodes, online._U().name(),
                    example.podex, 0, 0, NmodesUout);
       ITHACAPOD::getModes(example.Pfield, online.Pmodes, example._p().name(),
                            example.podex, 0, 0,NmodesPout);
       ITHACAPOD::getModes(example.Dfield, online.Dmodes, online._pointDisplacement().name(),
                            example.podex, 0, 0, NmodesDout);
    }
    else
    {

     //std::cout << "================================" << "\n";   
    //  ITHACAPOD::getModes(example.Ufield, online.Umodes, online._U().name(),
    //                 example.podex, 0, 0, NmodesUout);
    //  ITHACAPOD::getModes(example.Pfield, online.Pmodes, example._p().name(),
    //                         example.podex, 0, 0,NmodesPout);
    // ITHACAPOD::getModes(example.Dfield, online.Dmodes, online._pointDisplacement().name(),
    //                         example.podex, 0, 0, NmodesDout); 
      ITHACAstream::read_fields(online.Umodes, online._U(), "./ITHACAoutput/POD/");
      ITHACAstream::read_fields(online.Pmodes, online._p(), "./ITHACAoutput/POD/");
      ITHACAstream::read_fields(online.Dmodes, online._pointDisplacement(), "./ITHACAoutput/POD/");

    }
   online.coeffL2 = ITHACAutilities::getCoeffs(example.Dfield, online.Dmodes, NmodesDproj, false);
   if( (!ITHACAutilities::check_folder("./ITHACAoutput/Matrices/")) )
   {
            mkDir("ITHACAoutput/Matrices");
        Eigen::MatrixXd CoeffsPd = ITHACAutilities::getCoeffs(example.Dfield, online.Dmodes, NmodesDproj, false);
        Eigen::MatrixXd CoeffsP = ITHACAutilities::getCoeffs(example.Pfield, online.Pmodes, NmodesPproj, true);
        Eigen::MatrixXd CoeffsU = ITHACAutilities::getCoeffs(example.Ufield, online.Umodes, NmodesUproj, true);
        /// Export the matrix coefficient of pointDisplacement Field
        cnpy::save(CoeffsPd, "./ITHACAoutput/Matrices/CoeffsPd.npy");
        cnpy::save(CoeffsP, "./ITHACAoutput/Matrices/CoeffsP.npy");
        cnpy::save(CoeffsU, "./ITHACAoutput/Matrices/CoeffsU.npy");
    }
    

    /// ############### contruct the reduced the class object ###################
    ReducedFsi reduced(online);
    reduced.startTime = example.startTime;
    reduced.finalTime = example.finalTime;
    reduced.timeStep = example.timeStep;
    reduced.writeEvery = example.writeEvery;
    //Perform the online solutions
 
    scalar mu_now = example.mu(0, 0);
    std::clock_t startOn;
    double durationOn;
    startOn = std::clock();
    reduced.solveOnline_Pimple(mu_now, NmodesUproj, NmodesPproj, NmodesDproj);
    durationOn = (std::clock() - startOn);
    std::cout << "The Online  phase  duration  is  equal  to " << durationOn << std::endl;

    std::cout << "======================= ONLINE PHASE COMPLETED ================================" << "\n";    
    if (!ITHACAutilities::check_folder("./ITHACAoutput/DataFromRom"))
    {
        mkDir("ITHACAoutput/DataFromRom");
        Eigen::VectorXd RedCentersOfMassy = Foam2Eigen::field2Eigen(reduced.centerofmassy);
        cnpy::save(RedCentersOfMassy, "./ITHACAoutput/DataFromRom/RedCentersOfMassy_" + name(NmodesUproj) + "_" + name(NmodesPproj) + ".npy");

        Eigen::VectorXd romforcey = Foam2Eigen::field2Eigen(reduced.romforcey);
        cnpy::save(romforcey, "./ITHACAoutput/DataFromRom/romforcey_" + name(NmodesUproj) + "_" + name(NmodesPproj) + ".npy");

        Eigen::VectorXd romforcex = Foam2Eigen::field2Eigen(reduced.romforcex);
        cnpy::save(romforcex, "./ITHACAoutput/DataFromRom/romforcex_" + name(NmodesUproj) + "_" + name(NmodesPproj) + ".npy");

        Eigen::VectorXd pdcoeffrbf = Foam2Eigen::field2Eigen(reduced.pdcoeffrbf);
        ITHACAstream::exportMatrix(pdcoeffrbf, "CoeffsRbf", "python","./ITHACAoutput/DataFromRom/RedPdCoeff/");
       /// Export reduced coefficients vectors
        ITHACAstream::exportMatrix(reduced.CoeffP,"CoeffsP", "python", "./ITHACAoutput/DataFromRom/ReducedCoeffs/");
        ITHACAstream::exportMatrix(reduced.CoeffU,"CoeffsU", "python", "./ITHACAoutput/DataFromRom/ReducedCoeffs/");

    }
    
    Eigen::MatrixXd errL2U = ITHACAutilities::errorL2Rel(example.Ufield, reduced.UredFields);
    std::cout << "======================= errL2U completed================================" << "\n";
    Eigen::MatrixXd errL2P = ITHACAutilities::errorL2Rel(example.Pfield, reduced.PredFields);
    cnpy::save(errL2U, "./ITHACAoutput/DataFromRom/errL2U_" + name(NmodesUproj) + "_" + name(NmodesPproj) + ".npy");
    cnpy::save(errL2P, "./ITHACAoutput/DataFromRom/errL2P_" + name(NmodesUproj) + "_" + name(NmodesPproj) + ".npy");

    exit(0);
}

