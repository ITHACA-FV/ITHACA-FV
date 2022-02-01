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
    fsiBasic
Description
    fluid structure interaction  problem with dynamic mesh
SourceFiles
    fsiBasic.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the fsiBasic class.
/// \dir
/// Directory containing the header and source files for the fsiBasic class.

#ifndef fsiBasic_H
#define fsiBasic_H
#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"

#include "unsteadyNS.H"
#include "ITHACAstream.H"
#include "ITHACAutilities.H"
#include "ITHACAPOD.H"
#include "forces.H"
#include "IOmanip.H"
#include "IOstreams.H"
#include "objectRegistry.H"
#include "localEulerDdtScheme.H"
#include "reductionProblem.H"
#include <iostream>


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                        Class fsiBasic Declaration
\*---------------------------------------------------------------------------*/

/// Implementation of a parametrized full order <b> unsteady NS problem </b> and preparation of the the reduced matrices for the online solve.
/** In this class are implemented the methods for the offline solve of a unsteady NS problem
and the for the generation of the reduced matrices for subsequent online solve, this class is a son
of the steadyNS class */
class fsiBasic: public unsteadyNS
{
    public:
        // Constructors
        /// Construct Null
        fsiBasic();

        /// Construct with argc and argv
        fsiBasic(int argc, char* argv[]);
        ~fsiBasic() {};

        // Dummy variables to transform pimplefoam into a class
        /// pimpleControl
        //autoPtr<pimpleControl> _pimple;

        // /// Turbulence model
        // autoPtr<incompressible::turbulenceModel> turbulence;
        /// Velocity field
        // volVectorField& U;
        // /// Pressure field
        // volScalarField& p;
        // ///
        // surfaceScalarField& phi;
        // /// Dynamic mesh field
        autoPtr<Foam::dynamicFvMesh> meshPtr;
        //autoPtr<dynamicFvMesh> _mesh;

        // Functions

        //--------------------------------------------------------------------------
        /// @brief      Perform a truthsolve
        /// @param[in]  mu_now  The actual value of the parameter for this truthSolve. Used only
        /// to construct mu_interp matrix which is written out in a specified folder, also for par
        /// file in the Parameters folder.
        ///
        void truthSolve3(List<scalar> mu_now,
                        fileName folder = "./ITHACAoutput/Offline/");
        void liftSolve3();
};        
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif










