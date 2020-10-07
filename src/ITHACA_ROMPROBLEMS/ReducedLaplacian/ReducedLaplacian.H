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
    reducedLaplacian
Description
    A reduced problem for the Laplace problem
SourceFiles
    reducedLaplacian.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the reducedLaplacian class.
/// \dir
/// Directory containing the header and source files for the reducedLaplacian class


#ifndef ReducedLaplacian_H
#define ReducedLaplacian_H

#include "fvCFD.H"
#include "IOmanip.H"
#include "laplacianProblem.H"
#include "ReducedProblem.H"
#include <Eigen/Eigen>

/*---------------------------------------------------------------------------*\
                        Class reducedLaplacian Declaration
\*---------------------------------------------------------------------------*/

/// Class to solve the online reduced problem associated with a Laplace problem
/** In this class are implemented the methods to solve a reduced version of the Laplace problem
    This class is used together with the laplacianProblem class */
class reducedLaplacian: public reducedProblem
{
    private:

    public:
        // Constructors
        /// Construct Null
        reducedLaplacian();

        /// Construct from a laplacian problem
        explicit reducedLaplacian(laplacianProblem& problem);

        ~reducedLaplacian() {};

        /// Problem object
        laplacianProblem* problem;

        /// Online solution
        Eigen::MatrixXd online_solution;

        /// Counter for online sol
        int count_online_solve = 1;

        /// Function to perform an online solve given a certain mu
        ///
        /// @param[in]  mu    Actual value of the parameters that are multiplying,
        /// the affine expansion of the operators.
        ///
        void solveOnline(Eigen::MatrixXd mu);

        /// Function to recover the solution given the online solution
        ///
        /// @param[in]  folder      The folder where you want to store the results (default is "./ITHACAOutput/online_rec")
        /// @param[in]  printevery  Variable to recover only every printevery the online solutions
        ///
        void reconstruct(fileName folder = "./ITHACAOutput/online_rec",
                         int printevery = 1);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



#endif






