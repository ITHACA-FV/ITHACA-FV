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

#ifndef HyperReducedCompressibleUnSteadyNS_H
#define HyperReducedCompressibleUnSteadyNS_H

#include "fvCFD.H"
#include "IOmanip.H"
#include "ReducedProblem.H"
#include "UnsteadyProblem.H"
#include "ReducedCompressibleSteadyNS.H"
#include "ITHACAutilities.H"
#include <Eigen/Eigen>
#include "Modes.H"
#include "CompressibleUnSteadyRhoPimple.H"

/*---------------------------------------------------------------------------*\
              Class HyperReducedCompressibleUnSteadyNS Declaration
\*---------------------------------------------------------------------------*/

/// Class where it is implemented a reduced problem for the <b> UnSteady RhoPimple </b> problem.
/** In this class are implemented the methods to solve a hyper reduced version of the unstationary Navier-stokes
equations. This class must be used together with the compressible unsteadyNS class */
class HyperReducedCompressibleUnSteadyNS: public ReducedCompressibleSteadyNS
{
    private:

    public:
        // Constructors
        /// Construct Null
        HyperReducedCompressibleUnSteadyNS();

        /// Construct Null
        ///
        /// @param[in]  problem  a full order steadyNS problem
        /// @param[in]  tipo     Type of pressure stabilisation method you want to use "SUP" for supremizer, "PPE" for pressure Poisson equation.
        ///
        explicit HyperReducedCompressibleUnSteadyNS(CompressibleUnSteadyRhoPimple&
                problem);

        ~HyperReducedCompressibleUnSteadyNS() {};

        // Functions

        /// Method to perform an online solve using a PPE stabilisation method
        ///
        /// @param[in]  mu_now  The online viscosity.
        /// @param[in]  NmodesUproj  The number of modes one wants to use for velocity projection.
        /// @param[in]  NmodesPproj  The number of modes one wants to use for pressure projection.
        /// @param[in]  NmodesEproj  The number of modes one wants to use for energy projection.
        ///
        void SolveHyperReducedSys(int NmodesUproj,
                                  int NmodesPproj,
                                  int NmodesEproj,
                                  fileName folder = "./ITHACAoutput/Online/");

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

        ///
        /// @brief      It assembles the reduced oeprators using the modes
        ///
        /// @param[in]  S  The SparseMatrix from which you want to extract the rows indices.
        /// @param[in]  s  The SparseMatrix from which you want to extract the rows indices.
        /// @param[in]  MagicPoints coming from DEIM
        /// @param[in]  Modes modes given by the DEIM
        ///
        std::tuple<Eigen::MatrixXd, Eigen::VectorXd> HyperReducedSys(
            Eigen::SparseMatrix<double>& S,
            Eigen::VectorXd& se,
            List<label>& MagicPoints,
            Eigen::MatrixXd& Modes);

        /// Method to check and write solutions
        bool checkWrite(Time& timeObject);
        // Variables

        /// velocity modes.
        volVectorModes Umodes;
        volScalarModes Emodes, Pmodes;
        ///////// time control variables
        scalar startTime = 0.0;
        scalar finalTime = 0.0;
        scalar timeStep = 0.0;
        scalar writeEvery = timeStep;
        scalar nextWrite = 0.0;
        /// List to save the solutions
        PtrList<volScalarField> Pfield, Efield;
        PtrList<volVectorField> Ufield;
        // /// List of snapshots for the solution for eddy viscosity
        // PtrList<volScalarField> nutFields;

        // /// List of POD modes for eddy viscosity
        volScalarModes nutModes;

        /// Full problem.
        CompressibleUnSteadyRhoPimple* problem;

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






