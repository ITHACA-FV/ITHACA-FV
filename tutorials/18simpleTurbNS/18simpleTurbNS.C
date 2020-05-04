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
    Example of steady NS Reduction Problem
SourceFiles
    03steadyNS.C
\*---------------------------------------------------------------------------*/

#include "SteadyNSSimple.H"
#include "ITHACAstream.H"
#include "ITHACAPOD.H"
#include "ReducedSimpleSteadyNS.H"
#include "forces.H"
#include "IOmanip.H"


class tutorial18 : public SteadyNSSimple
{
    public:
        /// Constructor
        explicit tutorial18(int argc, char* argv[])
            :
            SteadyNSSimple(argc, argv),
            U(_U()),
            p(_p()),
            phi(_phi())
        {}

        /// Velocity field
        volVectorField& U;
        /// Pressure field
        volScalarField& p;
        ///
        surfaceScalarField& phi;

        /// Perform an Offline solve
        void offlineSolve()
        {
            Vector<double> inl(0, 0, 0);
            List<scalar> mu_now(1);

            // if the offline solution is already performed read the fields
            if (offline && ITHACAutilities::isTurbulent())
            {
                ITHACAstream::readMiddleFields(Ufield, U, "./ITHACAoutput/Offline/");
                ITHACAstream::readMiddleFields(Pfield, p, "./ITHACAoutput/Offline/");
                auto nut = _mesh().lookupObject<volScalarField>("nut");
                //readNut(nutFields, nut, "./ITHACAoutput/Offline/");
                ITHACAstream::readConvergedFields(nutFields, nut, "./ITHACAoutput/Offline/");
                mu_samples = ITHACAstream::readMatrix("./ITHACAoutput/Offline/mu_samples_mat.txt");
            }
            else if (offline)
            {
                ITHACAstream::read_fields(Ufield, U, "./ITHACAoutput/Offline/");
                ITHACAstream::read_fields(Pfield, p, "./ITHACAoutput/Offline/");
                mu_samples =
                    ITHACAstream::readMatrix("./ITHACAoutput/Offline/mu_samples_mat.txt");
            }
            else
            {
                Vector<double> Uinl(1, 0, 0);
                label BCind = 0;

                for (label i = 0; i < mu.rows(); i++)
                {
                    mu_now[0] = mu(i, 0);
                    change_viscosity(mu_now[0]);
                    assignIF(U, Uinl);
                    truthSolve2(mu_now);
                }
            }
        }

};

int main(int argc, char* argv[])
{
    // Construct the tutorial object
    tutorial18 example(argc, argv);
    // Read some parameters from file
    ITHACAparameters* para = ITHACAparameters::getInstance(example._mesh(),
                             example._runTime());
    // Read the par file where the parameters are stored
    std::ifstream exFileOff("./parsOff_mat.txt");
    if (exFileOff)
    {
        example.mu  = ITHACAstream::readMatrix("./parsOff_mat.txt");
    }

    else
    {
        example.mu  = Eigen::VectorXd::LinSpaced(50, 1.00e-04, 1.00e-05);
        ITHACAstream::exportMatrix(example.mu , "parsOff", "eigen", "./");
    }

    Eigen::MatrixXd parOn;
    std::ifstream exFileOn("./parsOn_mat.txt");
    if (exFileOn)
    {
        parOn = ITHACAstream::readMatrix("./parsOn_mat.txt");
    }

    else
    {
        parOn = ITHACAutilities::rand(20, 1, 1.00e-04, 1.00e-05);
        ITHACAstream::exportMatrix(parOn, "parsOn", "eigen", "./");
    }
    //word filename("./par");
    //example.mu = ITHACAstream::readMatrix(filename);
    // Set the inlet boundaries patch 0 directions x and y
    example.inletIndex.resize(1, 2);
    example.inletIndex(0, 0) = 0;
    example.inletIndex(0, 1) = 0;

    // Set the maximum iterations number for the offline phase
    example.maxIter = para->ITHACAdict->lookupOrDefault<int>("maxIter", 2000);

    // Perform the offline solve
    example.offlineSolve();
    ITHACAstream::read_fields(example.liftfield, example.U, "./lift/");
    // Homogenize the snapshots
    example.computeLift(example.Ufield, example.liftfield, example.Uomfield);
    // Perform POD on velocity and pressure and store the first 10 modes
    ITHACAPOD::getModes(example.Uomfield, example.Umodes, example._U().name(),
                        example.podex, 0, 0,
                        example.NUmodesOut);
    ITHACAPOD::getModes(example.Pfield, example.Pmodes, example._p().name(),
                        example.podex, 0, 0,
                        example.NPmodesOut);
    if(ITHACAutilities::isTurbulent())
    {
        ITHACAPOD::getModes(example.nutFields, example.nutModes, "nut",
                        example.podex, 0, 0, example.NNutModesOut);
        // Create the RBF for turbulence
        example.getTurbRBF(example.NNutModes);
    }

    // Create the reduced object
    reducedSimpleSteadyNS reduced(example);
    PtrList<volVectorField> U_rec_list;
    PtrList<volScalarField> P_rec_list;
    // Reads inlet volocities boundary conditions.
    word vel_file(para->ITHACAdict->lookup("online_velocities"));
    Eigen::MatrixXd vel = ITHACAstream::readMatrix(vel_file);

    // Set the maximum iterations number for the online phase
    reduced.maxIterOn = para->ITHACAdict->lookupOrDefault<int>("maxIterOn", 2000);

    //Perform the online solutions
    for (label k = 0; k < parOn.size(); k++)
    {
        scalar mu_now = parOn(k, 0);
        example.restart();
        example.change_viscosity(mu_now);
        reduced.setOnlineVelocity(vel);
        if(ITHACAutilities::isTurbulent())
        {
            std::cout << "sono turb" << std::endl;
            reduced.solveOnline_Simple(mu_now, example.NUmodes, example.NPmodes, example.NNutModes);
        }
        else
        {
            std::cout << "non sono turb" << std::endl;
            reduced.solveOnline_Simple(mu_now, example.NUmodes, example.NPmodes, example.NNutModes);
   
        }
    }

    exit(0);
}
