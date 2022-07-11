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
    UnsteadyBB
Description
    Reduction class for a Stationary Navier-Stokes problem.
SourceFiles
    UnsteadyBB.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the UnsteadyBB class.
/// \dir
/// Directory containing the header and source files for the UnsteadyBB class.

#ifndef UnsteadyCompressibleNS_H
#define UnsteadyCompressibleNS_H
//#include "fvCFD.H"
//#include "singlePhaseTransportModel.H"
//#include "turbulentTransportModel.H"
//#include "fvOptions.H"
#include "unsteadyNS.H"
#include "directionInterpolate.H"
#include "psiThermo.H"
#include "turbulentFluidThermoModel.H"7
//#include <iostream>


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                        Class SteadyBB Declaration
\*---------------------------------------------------------------------------*/

/// Implementation of a parametrized full order <b> unsteady Boussinesq problem </b> and preparation of the the reduced matrices for the online solve.
/** In this class are implemented the methods for the offline solve of a unsteady Boussinesq problem (two way coupling between momentum and energy equation)
and the for the generation of the reduced matrices for subsequent online solve, this class is a son
of the reduction problem class */
class UnsteadyCompressibleNS: public unsteadyNS
{
    public:
        // Constructors
        /// Null constructor
        UnsteadyCompressibleNS();

        /// Construct with argc and argv
        UnsteadyCompressibleNS(int argc, char* argv[]);
        ~UnsteadyCompressibleNS() {};

        // Member Functions
        /// List of snapshots for the temperature solution
        PtrList<volScalarField> Tfield;

        /// List of POD temperature modes
        volScalarModes Tmodes;

        // Dummy variables to transform rhoCentralFoam into a class

        /// Temperature field
        autoPtr<volScalarField> _T;

        /// Density field
        autoPtr<volScalarField> _rho;

        autoPtr<surfaceScalarField> _pos;

        autoPtr<surfaceScalarField> _neg;

        autoPtr<volVectorField> _rhoU;

        autoPtr<volScalarField> _rhoE;

        autoPtr<psiThermo> _pThermo;

        autoPtr<compressible::turbulenceModel> turbulence;

        // Methods

        //--------------------------------------------------------------------------
        /// Perform a truthsolve
        void truthSolve(List<scalar> params, fileName folder = "./ITHACAoutput/Offline/");

        //--------------------------------------------------------------------------
};

#endif
