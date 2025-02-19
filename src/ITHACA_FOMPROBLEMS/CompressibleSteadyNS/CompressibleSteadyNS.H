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

#ifndef CompressibleSteadyNS_H
#define CompressibleSteadyNS_H
#include "fvCFD.H"
#include "fluidThermo.H"
#include "turbulentFluidThermoModel.H"
#include "simpleControl.H"
#include "pressureControl.H"
#include "fvOptions.H"
#include "steadyNS.H"
#include "SteadyNSTurb.H"
#include "SteadyNSSimple.H"
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
class CompressibleSteadyNS: public SteadyNSTurb
{
    public:
        // Constructors
        /// Null constructor
        CompressibleSteadyNS();

        /// Construct with argc and argv
        CompressibleSteadyNS(int argc, char* argv[]);
        ~CompressibleSteadyNS() {};

        // Functions

        //--------------------------------------------------------------------------
        /// @brief      Perform a truthsolve
        /// @param[in]  mu_now  The actual value of the parameter for this truthSolve. Used only
        /// to construct mu_interp matrix which is written out in the Offline folder, also for par file in
        /// the Parameters folder.
        ///
        //void truthSolve(List<scalar> mu_now);
        void truthSolve(word folder);

        // Momuntum equation terms
        fvVectorMatrix getUmatrix(volVectorField& U);//, Vector<double>& uresidual_v);
        fvVectorMatrix getNLTerm(volVectorField& U);
        fvVectorMatrix getViscTerm(volVectorField& U);
        volVectorField getGradP(volScalarField& p);

        // Energy equation terms
        fvScalarMatrix getEmatrix(volVectorField& U,
                                  volScalarField& p);//, scalar& eresidual);
        fvScalarMatrix getFluxTerm();
        volScalarField getKinEnTerm(volVectorField& U, volScalarField& p);
        fvScalarMatrix getDiffTerm();

        // Pressure equation terms
        fvScalarMatrix getPmatrix(fvVectorMatrix& Ueqn, volVectorField& U, volScalarField& p);
        surfaceScalarField getPhiHbyA(fvVectorMatrix& Ueqn, volVectorField& U,
                                      volScalarField& p);
        volScalarField getDivPhiHbyA(fvVectorMatrix& Ueqn, volVectorField& U, volScalarField& p);
        surfaceScalarField getRhorAUf(fvVectorMatrix& Ueqn);
        fvScalarMatrix getPoissonTerm(fvVectorMatrix& Ueqn, volScalarField& p);

        //--------------------------------------------------------------------------
        /// Function to change the viscosity
        ///
        /// @param[in]  mu    viscosity (scalar)
        ///
        void changeViscosity(double mu_new);

        //--------------------------------------------------------------------------
        /// set all variables back to the values into the 0 folder
        /// 
        void restart();

        autoPtr<fluidThermo> pThermo;

        autoPtr<volScalarField> _rho;
        
        autoPtr<volScalarField> _p;

        autoPtr<pressureControl> _pressureControl;

        autoPtr<compressible::turbulenceModel> turbulence;

        autoPtr<dimensionedScalar> _initialMass;

        autoPtr<volScalarField> _psi;

        autoPtr<volScalarField> _E;

        /// List of pointers to store initial solutions

        /// Initial Pressure field (for restart purposes)
        autoPtr<volScalarField> _p0;

        /// Initial Energy field (for restart purposes)
        autoPtr<volScalarField> _E0;

        /// Initial Velocity field (for restart purposes)
        autoPtr<volVectorField> _U0;

        /// Initial Density field (for restart purposes)
        autoPtr<volScalarField> _rho0;

        /// Initial Flux field (for restart purposes)
        autoPtr<surfaceScalarField> _phi0;


        /// List of pointers used to store the energy solutions
        PtrList<volScalarField> Efield;

        /// List of pointers used to form the energy modes
        volScalarModes Emodes;

        autoPtr<fvVectorMatrix> Ueqn_global;

        autoPtr<fvScalarMatrix> Peqn_global;

        autoPtr<fvScalarMatrix> Eeqn_global;

        autoPtr<surfaceScalarField> phiHbyA;

        autoPtr<surfaceScalarField> rhorAUf;

        autoPtr<volVectorField> HbyA;

        /// Tolerance for the residual of the stationary problems, there is the same tolerance for velocity and pressure
        scalar tolerance;

        /// Number of maximum iterations to be done for the computation of the truth solution
        scalar maxIter;

        /// Distancing between intermediate steps (for turbulent case only)
        label middleStep;
        //label middleStepInit;

        /// Export also intermediate fields
        bool middleExport;



};

#endif










