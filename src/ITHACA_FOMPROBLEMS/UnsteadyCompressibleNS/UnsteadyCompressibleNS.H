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
/// Header file of the UnsteadyCompressibleNS class.
/// \dir
/// Directory containing the header and source files for the UnsteadyCompressibleNS class.

#ifndef UnsteadyCompressibleNS_H
#define UnsteadyCompressibleNS_H

#include "unsteadyNS.H"
#include "directionInterpolate.H"
#include "psiThermo.H"
#include "turbulentFluidThermoModel.H"
#include "localEulerDdtScheme.H"
//#include "dynamicFvMesh.H"
#include "fvcSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                Class UnsteadyCompressibleNS Declaration
\*---------------------------------------------------------------------------*/

/// Implementation of a parametrized full order <b> unsteady compressible Navier Stokes problem </b> and preparation of the the reduced matrices for the online solve.
/** In this class are implemented the methods for the offline solve of a unsteady compressible Navier Stokes problem
and for the generation of the reduced matrices for subsequent online solve. This class is a son of the unsteady NS problem class */
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

        /// List of snapshots for the density solution
        PtrList<volScalarField> rhofield;

        /// List of POD temperature modes
        volScalarModes Tmodes;

        /// List of POD density modes
        volScalarModes rhomodes;

        // Dummy variables to transform rhoCentralFoam into a class

        /// Density field
        autoPtr<volScalarField> _rho;

        /// Momentum 
        autoPtr<volVectorField> _rhoU;

        /// Total energy
        autoPtr<volScalarField> _rhoE;

        /// Volumetric flux
        autoPtr<surfaceScalarField> _phi;

        /// Volumetric flux associated with inward local propagation speed    
        autoPtr<surfaceScalarField> _pos;

        /// Volumetric flux associated with outward local propagation speed    
        autoPtr<surfaceScalarField> _neg;

        /// Thermodynamic model
        autoPtr<psiThermo> _pThermo;

        /// Turbulence model
        autoPtr<compressible::turbulenceModel> turbulence;

        // Methods

        //--------------------------------------------------------------------------
        /// Perform a truthsolve
        void truthSolve(fileName folder = "./ITHACAoutput/Offline/");

        //--------------------------------------------------------------------------
};

#endif
