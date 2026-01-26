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
Tutorial of compressible and unsteady flow around a moving airfoil
\*---------------------------------------------------------------------------*/

#include "CompressibleUnSteadyRhoPimple.H"
#include "ITHACAPOD.H"
#include "ITHACAstream.H"
#include "Foam2Eigen.H"
#include "DEIM.H"
//#include "ReducedProblem.H"
#include <chrono>
#include "HyperReducedCompressibleUnSteadyNS.H"

template<typename O>
class HROM : public DEIM<O>
{
    public:
        using DEIM<O>::DEIM;
        PtrList<O> fields;
        autoPtr<O> subField;

};


class tutorial26: public CompressibleUnSteadyRhoPimple
{
    public:
        tutorial26(int argc, char* argv[])
            : CompressibleUnSteadyRhoPimple(argc, argv),
              U(_U()),
              p(_p()),
              E(_E())
              //,pd(_pointDisplacement())
        {
            //point0 = meshPtr().points();
        }
        /// Velocity field
        volVectorField& U;
        /// Pressure field
        volScalarField& p;
        /// Energy field
        volScalarField& E;
        /// Hyper-reduced object
        //autoPtr<HROM<O>> DeimObject;
        ///grid nodes field
        //pointVectorField& pd;
        /// Initial coordinates of the grid points
        void offlineSolve(word folder = "./ITHACAoutput/Offline/")
        {
            //List<scalar> mu_now(1);
            if (offline)
            {
                ITHACAstream::read_fields(Ufield, U, folder);
                ITHACAstream::read_fields(Pfield, p, folder);
                ITHACAstream::read_fields(Efield, E, folder);
            }
            else
            {
                truthSolve(folder);
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
    tutorial26 example(argc, argv);
    std::clock_t startOff;
    double durationOff;
    // Read some parameters from file
    ITHACAparameters* para = ITHACAparameters::getInstance(example.meshPtr(),
                             example._runTime());
    // Read the par file where the parameters are stored
    int NmodesUout  =  readInt(para->ITHACAdict->lookup("NmodesUout"));
    int NmodesPout  =  readInt(para->ITHACAdict->lookup("NmodesPout"));
    int NmodesEout  =  readInt(para->ITHACAdict->lookup("NmodesEout"));
    int NmodesUproj  = readInt(para->ITHACAdict->lookup("NmodesUproj"));
    int NmodesPproj  = readInt(para->ITHACAdict->lookup("NmodesPproj"));
    int NmodesEproj  = readInt(para->ITHACAdict->lookup("NmodesEproj"));
    // word filename("./par");
    // example.mu = ITHACAstream::readMatrix(filename);
    // Time parameters: We can use Ioodictionnary to access time parameters
    example.startTime  = 0;
    example.finalTime  = 0.15;
    example.timeStep   = 2e-06;
    example.writeEvery = 4e-04;
    // //Perform the offline solve
    startOff = std::clock();
    example.offlineSolve();
    //exit(0);
    durationOff = (std::clock() - startOff);
    std::cout << "The Offline phase  duration  is  equal  to " << durationOff <<
              std::endl;

    if (example.podex == 0 )
    {
        ITHACAPOD::getModes(example.Ufield, example.Umodes, example._U().name(),
                            example.podex, 0, 0, NmodesUout);
        ITHACAPOD::getModes(example.Pfield, example.Pmodes, example._p().name(),
                            example.podex, 0, 0, NmodesPout);
        ITHACAPOD::getModes(example.Efield, example.Emodes, example.E().name(),
                            example.podex, 0, 0, NmodesEout);
    }
    else
    {
        ITHACAstream::read_fields(example.Umodes, example._U(), "./ITHACAoutput/POD/");
        ITHACAstream::read_fields(example.Pmodes, example._p(), "./ITHACAoutput/POD/");
        ITHACAstream::read_fields(example.Emodes, example._E(), "./ITHACAoutput/POD/");
    }

    HyperReducedCompressibleUnSteadyNS hyperreduced(example);
    // Info << hyperreduced.Umodes.size() << endl;
    // Info << hyperreduced.Pmodes.size() << endl;
    // Info << hyperreduced.Emodes.size() << endl;
    hyperreduced.startTime = example.startTime;
    hyperreduced.finalTime = example.finalTime;
    hyperreduced.timeStep = example.timeStep;
    hyperreduced.writeEvery = example.writeEvery;
    /// Solving the hyper-reduced problem
    hyperreduced.SolveHyperReducedSys(NmodesUproj, NmodesPproj, NmodesEproj);
    exit(0);
    // online.coeffL2 = ITHACAutilities::getCoeffs(example.Dfield, online.Dmodes, NmodesDproj, false);
    // if( (!ITHACAutilities::check_folder("./ITHACAoutput/Matrices/")) )
    // {
    //          mkDir("ITHACAoutput/Matrices");
    //      Eigen::MatrixXd CoeffsPd = ITHACAutilities::getCoeffs(example.Dfield, online.Dmodes, NmodesDproj, false);
    //      Eigen::MatrixXd CoeffsP = ITHACAutilities::getCoeffs(example.Pfield, online.Pmodes, NmodesPproj, true);
    //      Eigen::MatrixXd CoeffsU = ITHACAutilities::getCoeffs(example.Ufield, online.Umodes, NmodesUproj, true);
    //      /// Export the matrix coefficient of pointDisplacement Field
    //      cnpy::save(CoeffsPd, "./ITHACAoutput/Matrices/CoeffsPd.npy");
    //      cnpy::save(CoeffsP, "./ITHACAoutput/Matrices/CoeffsP.npy");
    //      cnpy::save(CoeffsU, "./ITHACAoutput/Matrices/CoeffsU.npy");
    //  }
    //  /// ############### contruct the reduced the class object ###################
    //  ReducedFsi reduced(online);
    //  reduced.startTime = example.startTime;
    //  reduced.finalTime = example.finalTime;
    //  reduced.timeStep = example.timeStep;
    //  reduced.writeEvery = example.writeEvery;
    //  //Perform the online solutions
    //  scalar mu_now = example.mu(0, 0);
    //  std::clock_t startOn;
    //  double durationOn;
    //  startOn = std::clock();
    //  reduced.solveOnline_Pimple(mu_now, NmodesUproj, NmodesPproj, NmodesDproj);
    //  durationOn = (std::clock() - startOn);
    //  std::cout << "The Online  phase  duration  is  equal  to " << durationOn << std::endl;
    //  std::cout << "======================= ONLINE PHASE COMPLETED ================================" << "\n";
    //  if (!ITHACAutilities::check_folder("./ITHACAoutput/DataFromRom"))
    //  {
    //      mkDir("ITHACAoutput/DataFromRom");
    //      Eigen::VectorXd RedCentersOfMassy = Foam2Eigen::field2Eigen(reduced.centerofmassy);
    //      cnpy::save(RedCentersOfMassy, "./ITHACAoutput/DataFromRom/RedCentersOfMassy_" + name(NmodesUproj) + "_" + name(NmodesPproj) + ".npy");
    //      Eigen::VectorXd romforcey = Foam2Eigen::field2Eigen(reduced.romforcey);
    //      cnpy::save(romforcey, "./ITHACAoutput/DataFromRom/romforcey_" + name(NmodesUproj) + "_" + name(NmodesPproj) + ".npy");
    //      Eigen::VectorXd romforcex = Foam2Eigen::field2Eigen(reduced.romforcex);
    //      cnpy::save(romforcex, "./ITHACAoutput/DataFromRom/romforcex_" + name(NmodesUproj) + "_" + name(NmodesPproj) + ".npy");
    //      Eigen::VectorXd pdcoeffrbf = Foam2Eigen::field2Eigen(reduced.pdcoeffrbf);
    //      ITHACAstream::exportMatrix(pdcoeffrbf, "CoeffsRbf", "python","./ITHACAoutput/DataFromRom/RedPdCoeff/");
    //     /// Export reduced coefficients vectors
    //      ITHACAstream::exportMatrix(reduced.CoeffP,"CoeffsP", "python", "./ITHACAoutput/DataFromRom/ReducedCoeffs/");
    //      ITHACAstream::exportMatrix(reduced.CoeffU,"CoeffsU", "python", "./ITHACAoutput/DataFromRom/ReducedCoeffs/");
    //  }
    // Eigen::MatrixXd errL2U = ITHACAutilities::errorL2Rel(example.Ufield, reduced.UredFields);
    // std::cout << "======================= errL2U completed================================" << "\n";
    // Eigen::MatrixXd errL2P = ITHACAutilities::errorL2Rel(example.Pfield, reduced.PredFields);
    // cnpy::save(errL2U, "./ITHACAoutput/DataFromRom/errL2U_" + name(NmodesUproj) + "_" + name(NmodesPproj) + ".npy");
    // cnpy::save(errL2P, "./ITHACAoutput/DataFromRom/errL2P_" + name(NmodesUproj) + "_" + name(NmodesPproj) + ".npy");
    exit(0);
}


