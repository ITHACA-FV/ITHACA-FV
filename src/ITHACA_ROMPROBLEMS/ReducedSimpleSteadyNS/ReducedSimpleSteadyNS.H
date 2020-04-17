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

#ifndef ReducedSimpleSteadyNS_H
#define ReducedSimpleSteadyNS_H

#include "fvCFD.H"
#include "IOmanip.H"
#include "ReducedProblem.H"
#include "SteadyNSSimple.H"
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
class reducedSimpleSteadyNS: public reducedProblem
{
    private:

    public:
        // Constructors
        /// Construct Null
        reducedSimpleSteadyNS();

        /// Construct Null
        ///
        /// @param[in]  problem  a full order steadyNS problem
        /// @param[in]  tipo     Type of pressure stabilisation method you want to use "SUP" for supremizer, "PPE" for pressure Poisson equation.
        ///
        explicit reducedSimpleSteadyNS(SteadyNSSimple& problem);

        ~reducedSimpleSteadyNS() {};

        // Functions

        /// Method to perform an online solve using a PPE stabilisation method
        ///
        /// @param[in]  vel_now  The vector of online velocity. It is defined in
        /// with an Eigen::MatrixXd and must have one col and as many rows as the number
        /// of parametrized boundary conditions.
        ///
        /// @param[in]  NmodesUproj The number of modes required for velocity projection
        ///
        /// @param[in]  NmodesPproj The number of modes required for pressure projection
        ///
        /// @param[in]  Folder The folder where we want to esport solutions (default is "./ITHACAoutput/Reconstruct/")
        ///
        void solveOnline_Simple(scalar mu_now, int NmodesUproj, int NmodesPproj, int NmodesNut = 0,
                                int NmodesSup = 0,
                                word Folder = "./ITHACAoutput/Reconstruct/");

        ///
        /// @brief      Gets the Umatrix for the online solver: it just assembles the pure matrix.
        ///
        /// @param[in]      U     Velocity vector field
        /// @param[in]      p     Pressure scalar field
        ///
        /// @return     Velocity linear system.
        ///
        fvVectorMatrix get_Umatrix_Online(volVectorField& U, volScalarField& p);

        ///
        /// @brief      Gets the Pmatrix for the online solver: it just assembls the pure matrix.
        ///
        /// @param[in]      U     Velocity vector field
        /// @param[in]      p     Pressure scalar field
        ///
        /// @return     Pressure linear system.
        ///
        fvScalarMatrix get_Pmatrix_Online(volVectorField& U, volScalarField& p);

        ///
        /// @brief      It checks if the number of imposed boundary conditions is correct
        /// and set the inlet velocity equal to the given one.
        ///
        /// @param[in]  vel   Inlet velocity
        ///
        void setOnlineVelocity(Eigen::MatrixXd vel);

        // Variables

        /// Lifted velocity modes.
        Modes<vector> ULmodes;

        /// Full problem.
        SteadyNSSimple* problem;

        /// Imposed boundary conditions.
        Eigen::MatrixXd vel_now;

        /// Counter.
        int counter = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



#endif





