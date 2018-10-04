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
    reducedProblem
Description
    A general reduced problem class is implemented
SourceFiles
    reducedProblem.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the reducedProblem class.
/// \dir
/// Directory containing the header and source files for the reducedProblem class.
/// \dir ../reducedProblems
/// Directory containing the implementation of all the reduced problems and of the methods to solve online reduced problems.

#ifndef reducedProblem_H
#define reducedProblem_H

#include "fvCFD.H"
#include "IOmanip.H"
#include "reductionProblem.H"
#include "../thirdparty/Eigen/Eigen/Eigen"
#include "newton_argument.H"


/*---------------------------------------------------------------------------*\
                        Class reducedProblem Declaration
\*---------------------------------------------------------------------------*/
/// Base class for the implementation of a reduced problem
class reducedProblem
{
    public:
        // Constructors
        /// Construct Null
        reducedProblem();
        virtual ~reducedProblem() {};

        /// Construct with reduced Problem
        ///
        /// @param  problem a reductionProblem or a child of it
        ///
        explicit reducedProblem(reductionProblem& problem);

        /// Viscosity
        scalar nu;

        // Function
        /// Virtual Method to perform and online Solve
        virtual void solveOnline();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

