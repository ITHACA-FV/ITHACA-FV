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
Class
    ReducedBasicFsi
Description
    A reduced problem for the Basic Fsi with mesh motion equations
SourceFiles
    ReducedBasicFsi.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the ReducedBasicFsi class.
/// \dir
/// Directory containing the header and source files for the ReducedBasicFsi class

#ifndef ReducedBasicFsi_H
#define ReducedBasicFsi_H

#include "fvCFD.H"
#include "IOmanip.H"
#include "ReducedProblem.H"
#include "fsiBasic.H"
#include "ITHACAutilities.H"
#include "ReducedSimpleSteadyNS.H"
#include "EigenFunctions.H"
#include <Eigen/Eigen>
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>
#include "Modes.H"


class ReducedBasicFsi : public reducedSimpleSteadyNS
{
    public:
        /// Constructor
        ReducedBasicFsi() {}
        explicit ReducedBasicFsi(fsiBasic& FOMproblem);
        ~ReducedBasicFsi() {}

        /// Variables

        /// Lifted velocity modes.
        volVectorModes ULmodes;
        /// Pointer to the Full problem.
        fsiBasic* problem;
        /// Scalar to store the current time
        scalar time;
        /// Scalar to store the time increment
        double dt;

        /// Scalar to store the final time if the online simulation
        scalar finalTime;

        /// Scalar to store the initial time if the online simulation
        scalar tstart;

        /// A variable for storing the reduced coefficients
        double storeEvery;

        /// A variable for exporting the fields
        double exportEvery;

        /// Projected gradient of the pressure modes.
        Eigen::MatrixXd projGradModP;

        /// Imposed boundary conditions.
        Eigen::MatrixXd vel_now;

        /// Maximum iterations number for the online step
        int maxIterOn = 1000;

        /// Counter.
        int counter = 0;

        int UprojN;
        int PprojN;
        // Set the initial time
        // instantList Times = runTime.times();
        // runTime.setEndTime(finalTime);

        // // Perform a TruthSolve
        // runTime.setTime(Times[1], 1);
        // runTime.setDeltaT(timeStep);
        // nextWrite = tstart; // timeStep initialization

        // Function to perform the online phase
        void solveOnline_Pimple(scalar mu_now, int NmodesUproj, int NmodesPproj, int NmodesSup = 0, fileName  folder = "./ITHACAoutput/Reconstruct/");
        ///
        /// @brief      It checks if the number of imposed boundary conditions is correct
        /// and set the inlet velocity equal to the given one.
        ///
        /// @param[in]  vel   Inlet velocity
        ///
        void setOnlineVelocity(Eigen::MatrixXd vel);
};

