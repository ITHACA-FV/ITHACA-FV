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
    reducedSteadyNS
Description
    A reduced problem for the stationary NS equations
SourceFiles
    reducedSteadyNS.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the reducedSteadyNS class.
/// \dir
/// Directory containing the header and source files for the reducedSteadyNS class

#ifndef ReducedCompressibleSteadyNS_H
#define ReducedCompressibleSteadyNS_H

#include "fvCFD.H"
#include "IOmanip.H"
#include "ReducedProblem.H"
#include "SteadyNSSimple.H"
#include "ReducedSimpleSteadyNS.H"
#include "ReducedSteadyNSTurb.H"
#include "CompressibleSteadyNS.H"
#include "ITHACAutilities.H"
#include "EigenFunctions.H"
#include <Eigen/Eigen>
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>
#include "Modes.H"


/*---------------------------------------------------------------------------*\
                        Class reducedProblem Declaration
\*---------------------------------------------------------------------------*/

/// Class where it is implemented a reduced problem for the <b> steady Navier-stokes </b> problem.
/** In this class are implemented the methods to solve a reduced version of the stationary Navier-stokes
equations. This class must be used together with the steadyNS class */
class ReducedCompressibleSteadyNS: public ReducedSteadyNSTurb
{
    private:

    public:
        // Constructors
        /// Construct Null
        ReducedCompressibleSteadyNS();

        /// Construct Null
        ///
        /// @param[in]  problem  a full order steadyNS problem
        /// @param[in]  tipo     Type of pressure stabilisation method you want to use "SUP" for supremizer, "PPE" for pressure Poisson equation.
        ///
        explicit ReducedCompressibleSteadyNS(CompressibleSteadyNS& problem);

        ~ReducedCompressibleSteadyNS() {};

        // Functions

        /// Method to perform an online solve using a PPE stabilisation method
        ///
        /// @param[in]  mu_now  The online viscosity.
        /// @param[in]  NmodesUproj  The number of modes one wants to use for velocity projection.
        /// @param[in]  NmodesPproj  The number of modes one wants to use for pressure projection.
        /// @param[in]  NmodesEproj  The number of modes one wants to use for energy projection.
        ///
        void solveOnlineCompressible(scalar mu_now, int NmodesUproj, int NmodesPproj, int NmodesEproj);

        ///
        /// @brief      It checks if the number of imposed boundary conditions is correct
        /// and set the inlet velocity equal to the given one.
        ///
        /// @param[in]  vel   Inlet velocity
        ///
        void setOnlineVelocity(Eigen::MatrixXd vel);

        ///
        /// @brief      It assembles the reduced oeprators using the modes
        ///
        /// @param[in]  NmodesUproj  The number of modes one wants to use for velocity projection.
        /// @param[in]  NmodesPproj  The number of modes one wants to use for pressure projection.
        /// @param[in]  NmodesEproj  The number of modes one wants to use for energy projection.
        ///
        void projectReducedOperators(int NmodesUproj, int NmodesPproj, int NmodesEproj);

        // Variables

        /// Lifted velocity modes.
        volVectorModes ULmodes;

        // /// List of snapshots for the solution for eddy viscosity
        // PtrList<volScalarField> nutFields;

        // /// List of POD modes for eddy viscosity
        // volScalarModes nutModes;

        /// Full problem.
        CompressibleSteadyNS* problem;

        /// Imposed boundary conditions.
        Eigen::MatrixXd vel_now;

        /// Counter.
        int counter = 0;

        //Reduced operators

        /// Projected gradient of the pressure modes.
        Eigen::MatrixXd projGradModP;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



#endif






