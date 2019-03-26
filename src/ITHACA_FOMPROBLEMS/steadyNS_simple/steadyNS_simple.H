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
    steadyNS
Description
    Reduction class for a Stationary Navier-Stokes problem.
SourceFiles
    steadyNS.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the steadyNS class.
/// \dir
/// Directory containing the header and source files for the steadyNS class.

#ifndef steadyNS_simple_H
#define steadyNS_simple_H
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "pisoControl.H"
#include "fvOptions.H"
#include "steadyNS.H"
#include "reductionProblem.H"
#include "ITHACAstream.H"
#include "ITHACAparameters.H"
#include "ITHACAforces.H"
#include "volFields.H"
#include <iostream>


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                        Class SteadyNS Declaration
\*---------------------------------------------------------------------------*/

/// Implementation of a parametrized full order <b> steady NS problem </b> and preparation of the the reduced matrices for the online solve.
/** In this class are implemented the methods for the offline solve of a steady NS problem
and the for the generation of the reduced matrices for subsequent online solve, this class is a son
of the reduction problem class */
class steadyNS_simple: public steadyNS
{

    public:
        // Constructors
        /// Null constructor
        steadyNS_simple();

        /// Construct with argc and argv
        steadyNS_simple(int argc, char* argv[]);
        ~steadyNS_simple() {};

        // Variables
        /// Initialization for the full velocity linear system.
        fvVectorMatrix* Ueqn_global;
        /// Initialization for the full pressure linear system.
        fvScalarMatrix* Peqn_global;

        //Functions
        ///
        /// @brief      Gets the offline velocity matrix and solve the system. Procedure of the SIMPLE algorithm implemented
        /// in openFOAM is followed with same corrections and iterations.
        ///
        /// @param      U     Velocity field at previous step.
        /// @param      p     Pressure field at previous step.
        ///
        /// @return     The velocity matrix.
        ///
        fvVectorMatrix get_Umatrix(volVectorField& U, volScalarField& p);
        ///
        /// @brief      Gets the offline pressure matrix, solve the system and updates the residual for pressure.
        /// Procedure of the SIMPLE algorithm implemented in openFOAM is followed with same corrections and iterations.
        ///
        /// @param      U     Velocity field at previous step.
        /// @param      p     Pressure field at previous step.
        ///
        /// @return     The velocity matrix.
        ///
        fvScalarMatrix get_Pmatrix(volVectorField& U, volScalarField& p,
                                   scalar& presidual);
        ///
        /// @brief      Offline solver for the whole Navier-Stokes problem.
        ///
        /// @param[in]  mu_now  Viscosity (parametrized) we want to solve the problem with.
        ///
        void truthSolve2(List<scalar> mu_now);
};

#endif










