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
    Example of turbulent steady NS Reduction Problem solved by the use of the SIMPLE algorithm
SourceFiles
    019simpleTurbGeom.C
\*---------------------------------------------------------------------------*/

#include "SteadyNSSimple.H"
#include "ITHACAstream.H"
#include "ITHACAPOD.H"
#include "ReducedSimpleSteadyNS.H"
#include "forces.H"
#include "IOmanip.H"
#include "RBFMotionSolver.H"


class tutorial19 : public SteadyNSSimple
{
    public:
        /// Constructor
        explicit tutorial19(int argc, char* argv[])
            :
            SteadyNSSimple(argc, argv),
            U(_U()),
            p(_p()),
            phi(_phi()),
            mesh(_mesh()),
            curX(mesh.points()),
            point0(mesh.points())
        {}

        /// Velocity field
        volVectorField& U;
        /// Pressure field
        volScalarField& p;
        /// Flux field
        surfaceScalarField& phi;
        /// The mesh
        fvMesh& mesh;
        /// Initial coordinates of the grid points
        vectorField point0;
        /// List to store the moved coordinates
        List<vector> curX;

        /// Perform an Offline solve
        void offlineSolve(Eigen::MatrixXd Box, List<label> patches)
        {
            Vector<double> inl(0, 0, 0);
            List<scalar> mu_now(1);

            // if the offline solution is already performed read the fields
            if (offline && ITHACAutilities::isTurbulent())
            {
                ITHACAstream::readMiddleFields(Ufield, U, "./ITHACAoutput/Offline/");
                ITHACAstream::readMiddleFields(Pfield, p, "./ITHACAoutput/Offline/");
                auto nut = _mesh().lookupObject<volScalarField>("nut");
                ITHACAstream::readConvergedFields(nutFields, nut, "./ITHACAoutput/Offline/");
            }
            else if (offline)
            {
                ITHACAstream::read_fields(Ufield, U, "./ITHACAoutput/Offline/");
                ITHACAstream::read_fields(Pfield, p, "./ITHACAoutput/Offline/");
            }
            else
            {
                Vector<double> Uinl(1, 0, 0);

                for (label i = 0; i < mu.rows(); i++)
                {
                    _mesh().movePoints(point0);
                    List<vector> points2Move;
                    labelList boxIndices = ITHACAutilities::getIndicesFromBox(_mesh(), patches, Box,
                                           points2Move);
                    mu_now[0] = mu(i, 0);
                    linearMovePts(mu_now[0], points2Move);

                    for (int j = 0; j < boxIndices.size(); j++)
                    {
                        curX[boxIndices[j]] = points2Move[j];
                    }

                    Field<vector> curXV(curX);
                    _mesh().movePoints(curXV);
                    ITHACAstream::writePoints(_mesh().points(), "./ITHACAoutput/Offline/",
                                              name(i + 1) + "/polyMesh/");
                    assignIF(U, Uinl);
                    truthSolve2(mu_now);
                    int middleF = 1;

                    while (ITHACAutilities::check_folder("./ITHACAoutput/Offline/" + name(
                            i + 1) + "/" + name(middleF)))
                    {
                        word command("ln -s  $(readlink -f ./ITHACAoutput/Offline/" + name(
                                         i + 1) + "/polyMesh/) ./ITHACAoutput/Offline/" + name(i + 1) + "/" + name(
                                         middleF) + "/" + " >/dev/null 2>&1");
                        std::cout.setstate(std::ios_base::failbit);
                        system(command);
                        std::cout.clear();
                        middleF++;
                    }
                }
            }
        }

        void linearMovePts(double angle, List<vector>& points2Move)
        {
            double sMax;
            scalarList x;
            scalarList y;

            for (label i = 0; i < points2Move.size(); i++)
            {
                x.append(points2Move[i].component(0));
                y.append(points2Move[i].component(1));
            }

            double xMin = min(x);
            double xMax = max(x);
            double yMin = min(y);
            double yMax = max(y);
            sMax = (yMax - yMin) * std::tan(M_PI * angle / 180);

            for (label i = 0; i < points2Move.size(); i++)
            {
                points2Move[i].component(0) = points2Move[i].component(0) +
                                              (yMax - points2Move[i].component(1)) / (yMax - yMin) * (xMax -
                                                      points2Move[i].component(0)) / (xMax - xMin) * sMax;
            }
        }
};

int main(int argc, char* argv[])
{
    // Construct the tutorial object
    tutorial19 example(argc, argv);
    // Read some parameters from file
    ITHACAparameters* para = ITHACAparameters::getInstance(example._mesh(),
                             example._runTime());
    // Read the files where the parameters are stored
    std::ifstream exFileOff("./angOff_mat.txt");

    if (exFileOff)
    {
        example.mu  = ITHACAstream::readMatrix("./angOff_mat.txt");
    }
    else
    {
        example.mu  = Eigen::VectorXd::LinSpaced(50, 0, 60);
        ITHACAstream::exportMatrix(example.mu, "angOff", "eigen", "./");
    }

    Eigen::MatrixXd angOn;
    std::ifstream exFileOn("./angOn_mat.txt");

    if (exFileOn)
    {
        angOn = ITHACAstream::readMatrix("./angOn_mat.txt");
    }
    else
    {
        angOn = ITHACAutilities::rand(20, 1, 2, 58);
        ITHACAstream::exportMatrix(angOn, "angOn", "eigen", "./");
    }

    // Set the inlet boundaries patch 0 directions x and y
    example.inletIndex.resize(1, 2);
    example.inletIndex(0, 0) = 0;
    example.inletIndex(0, 1) = 0;
    //Set the box including the step
    Eigen::MatrixXd Box(2, 3);
    Box << 1.98, 0.01, 0.11,
        7.02, -0.666669, -0.01;
    //Select the patches to be moved
    List<label> movPat;
    movPat.append(3);
    movPat.append(5);
    // Set the maximum iterations number for the offline phase
    example.maxIter = para->ITHACAdict->lookupOrDefault<int>("maxIter", 2000);
    // Perform the offline solve
    example.offlineSolve(Box, movPat);
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

    if (ITHACAutilities::isTurbulent())
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
    for (label k = 0; k < angOn.size(); k++)
    {
        scalar mu_now = angOn(k, 0);
        example.restart();
        reduced.setOnlineVelocity(vel);

        if (ITHACAutilities::isTurbulent())
        {
            example._mesh().movePoints(example.point0);
            List<vector> points2Move;
            labelList boxIndices = ITHACAutilities::getIndicesFromBox(example._mesh(),
                                   movPat, Box, points2Move);
            example.linearMovePts(mu_now, points2Move);

            for (int j = 0; j < boxIndices.size(); j++)
            {
                example.curX[boxIndices[j]] = points2Move[j];
            }

            Field<vector> curXV(example.curX);
            example._mesh().movePoints(curXV);
            ITHACAstream::writePoints(example._mesh().points(),
                                      "./ITHACAoutput/Reconstruct/", name(k + 1) + "/polyMesh/");
            reduced.solveOnline_Simple(mu_now, example.NUmodes, example.NPmodes,
                                       example.NNutModes);
        }
        else
        {
            example._mesh().movePoints(example.point0);
            List<vector> points2Move;
            labelList boxIndices = ITHACAutilities::getIndicesFromBox(example._mesh(),
                                   movPat, Box, points2Move);
            example.linearMovePts(mu_now, points2Move);

            for (int j = 0; j < boxIndices.size(); j++)
            {
                example.curX[boxIndices[j]] = points2Move[j];
            }

            Field<vector> curXV(example.curX);
            example._mesh().movePoints(curXV);
            ITHACAstream::writePoints(example._mesh().points(),
                                      "./ITHACAoutput/Reconstruct/", name(k + 1) + "/polyMesh/");
            reduced.solveOnline_Simple(mu_now, example.NUmodes, example.NPmodes,
                                       example.NNutModes);
        }
    }

    exit(0);
}
