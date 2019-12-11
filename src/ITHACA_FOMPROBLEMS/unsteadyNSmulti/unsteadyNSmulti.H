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
    unsteadyNS
Description
    Reduction class for a NON-Stationary NS multiphase problem
SourceFiles
    unsteadyNSmulti.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the unsteadyNSmulti class.
/// \dir
/// Directory containing the header and source files for the unsteadyNSmulti class.

#ifndef unsteadyNSmulti_H
#define unsteadyNSmulti_H
#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"
#include "unsteadyNS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                        Class SteadyNSmulti Declaration
\*---------------------------------------------------------------------------*/

/// Implementation of a parametrized full order <b> unsteady NS multiphase problem </b> and preparation of the the reduced matrices for the online solve.
/** In this class are implemented the methods for the offline solve of a unsteady NS multiphase problem
and the for the generation of the reduced matrices for subsequent online solve, this class is a son
of the unsteadyNS class */
class unsteadyNSmulti: public unsteadyNS
{
    public:
        // Constructors
        /// Construct Null
        unsteadyNSmulti();

        /// Construct with argc and argv
        unsteadyNSmulti(int argc, char* argv[]);

        // AutoPtr member variables
        autoPtr<volScalarField> _p_rgh;

        autoPtr<immiscibleIncompressibleTwoPhaseMixture> _mixture;

        autoPtr<volScalarField> _rho;

        autoPtr<surfaceScalarField> _rhoPhi;    

        autoPtr<surfaceScalarField> _alphaPhi10;

        autoPtr<tmp<surfaceScalarField>> _talphaPhi1Corr0;

        autoPtr<uniformDimensionedVectorField> _g;

        autoPtr<uniformDimensionedScalarField> _hRef;

        autoPtr<volScalarField> _gh;

        autoPtr<surfaceScalarField> _ghf;

        tmp<volScalarField> rAU;

        autoPtr<surfaceVectorField> Uf;
        
        bool correctPhi;

        bool checkMeshCourantNo;  

        bool moveMeshOuterCorrectors;

        bool alphaRestart;

        bool LTS;


        // Functions

        //--------------------------------------------------------------------------
        /// @brief      Perform a truthsolve
        /// @param[in]  mu_now  The actual value of the parameter for this truthSolve. Used only
        /// to construct mu_interp matrix which is written out in the Offline folder, also for par file in
        /// the Parameters folder.
        ///
        void truthSolve(List<scalar> mu_now);

};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif










