/*---------------------------------------------------------------------------*\
     ██╗████████╗██╗  ██╗ █████╗  ██████╗ █████╗       ███████╗██╗   ██╗
     ██║╚══██╔══╝██║  ██║██╔══██╗██╔════╝██╔══██╗      ██╔════╝██║   ██║
     ██║   ██║   ███████║███████║██║     ███████║█████╗█████╗  ██║   ██║
     ██║   ██║   ██╔══██║██╔══██║██║     ██╔══██║╚════╝██╔══╝  ╚██╗ ██╔╝
     ██║   ██║   ██║  ██║██║  ██║╚██████╗██║  ██║      ██║      ╚████╔╝
     ╚═╝   ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝      ╚═╝       ╚═══╝

 * In real Time Highly Advanced Computational Applications for Finite Volumes
 * Copyright (C) 2020 by the ITHACA-FV authors
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
    reducedProblem
Description
    A reduced problem for the unsteady NS equations
SourceFiles
    ReducedUnsteadyNSExplicit.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the ReducedUnsteadyNSExplicit class.
/// \dir
/// Directory containing the header and source files for the ReducedUnsteadyNSExplicit class

#ifndef ReducedUnsteadyNSExplicit_H
#define ReducedUnsteadyNSExplicit_H

#include "fvCFD.H"
#include "IOmanip.H"
#include "ReducedSteadyNS.H"
#include "ReducedUnsteadyNS.H"
#include "UnsteadyNSExplicit.H"
#include "unsteadyNS.H"
#include <Eigen/Dense>
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>

/*---------------------------------------------------------------------------*\
                        Class reducedProblem Declaration
\*---------------------------------------------------------------------------*/

/// Class where it is implemented a reduced problem for the <b> unsteady Navier-stokes </b> problem.
/** This class must be used together with the steadyNS class */
class ReducedUnsteadyNSExplicit: public reducedUnsteadyNS
{
    private:

    public:
        // Constructors
        /// Construct Null
        ReducedUnsteadyNSExplicit();
        /// Construct Null
        ///
        /// @param      problem  a full order UnsteadyNSExplicit problem
        ///
        explicit ReducedUnsteadyNSExplicit(UnsteadyNSExplicit& problem);

        ~ReducedUnsteadyNSExplicit() {};

        /// Pointer to the FOM problem
        UnsteadyNSExplicit* problem;

        // Functions

        /// Method to perform an online solve without a pressure stabilisation method
        ///
        /// @param[in]  vel_now   The vector of online velocity. It is defined in
        /// with an Eigen::MatrixXd and must have one col and as many rows as the number
        /// of parametrized boundary conditions.
        /// @param[in]  startSnap The first snapshot taken from the offline snapshots
        /// and used to get the reduced initial condition.
        ///
        void solveOnline(Eigen::MatrixXd vel_now, label startSnap = 0);

        /// Method to reconstruct the solutions from an online solve with a
        /// supremizer stabilisation technique. stabilisation method
        ///
        /// @param[in]  exportFields  A boolean variable which determines whether to export fields or not
        /// @param[in]  folder        The folder where to output the solutions in case on wants to
        ///
        void reconstruct(bool exportFields = false,
                         fileName folder = "./online_rec");


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



#endif

