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
ReducedFsi
Description
Reduced classs for fluid structure interaction  problem with dynamic mesh
SourceFiles
ReducedFsi.C
\*---------------------------------------------------------------------------*/
#ifndef ReducedFsi_H
#define ReducedFsi_H

#include "fsiBasic.H"
#include "ITHACAPOD.H"
#include "ReducedProblem.H"
#include "pointConstraints.H"
#include "mathematicalConstants.H"
#include "ReducedSimpleSteadyNS.H"

class ReducedFsi: public reducedSimpleSteadyNS
{
    public:
        // Constructors
        /// Construct Null
        ReducedFsi();
        ///
        /// @param[in]  problem  a full order Fsi problem
        ///
        explicit ReducedFsi(fsiBasic& problem);
        ~ReducedFsi() {};

        /// Full problem.
        fsiBasic* problem;
        /// The respective modes
        volScalarModes Pmodes;
        volVectorModes Umodes;
        pointVectorModes Dmodes;
        /// PtrList to save the meshes
        PtrList<dynamicFvMesh> OnlineMeshes;

        ///////// time control variables
        scalar startTime = 0.0;
        scalar finalTime = 0.0;
        scalar timeStep = 0.0;
        scalar writeEvery = timeStep;
        scalar nextWrite = 0.0;
        /// List scalar for access the centerofmass
        List<scalar> centerofmassx;
        List<scalar> centerofmassy;
        List<scalar> centerofmassz;
        /// List to save ROM lift and drag forces
        List<scalar> romforcey;
        List<scalar> romforcex;
        // List of scaler for saving interpolation  coefficients using RBF
        List<scalar> pdcoeffrbf;
        /// List of POD coefficients
        List<Eigen::MatrixXd> CoeffU;
        List<Eigen::MatrixXd> CoeffP;

        PtrList<volScalarField> PredFields;
        PtrList<volVectorField> UredFields;
        /// List of Online nodes displacement
        PtrList<pointVectorField> Dfield;
        PtrList<volScalarField> Pfield;
        PtrList<volVectorField> Ufield;
        PtrList<pointField> ListOfpoints;
        label counter = problem->counter;
        //std::ofstream Errors;

        // Functions


        ///
        /// @param[in]  mu_now  is the parameter. It is defined in
        /// with an Eigen::MatrixXd and must have one col and as many rows depending on the number of parameters.
        ///
        /// @param[in]  NmodesUproj The number of modes required for velocity projection
        ///
        /// @param[in]  NmodesPproj The number of modes required for pressure projection
        ///
        /// @param[in]  NmodesDproj The number of modes required for  interpolation
        ///
        /// @param[in]  Folder The folder where we want to esport solutions (default is "./ITHACAoutput/Online/")
        ///
        void solveOnline_Pimple(int NmodesUproj,
                                int NmodesPproj,
                                int NmodesDproj,
                                fileName folder = "./ITHACAoutput/Online/");
        /// method to set all fields back to values in 0 folder
        void restart();
        /// Method to save snapshots at a given time
        bool checkWrite(Time& timeObject);
        ///  Method to construct PODI
        void PODI(Eigen::MatrixXd coeffL2, Eigen::MatrixXd muu, label NPdModes);
        /// Routine to save forces and center of mass
        void prepareRomData(const word& outputPath);
};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif



