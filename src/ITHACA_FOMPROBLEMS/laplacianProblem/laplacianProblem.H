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
    laplacianProblem
Description
    A general full order implementation of a laplacian problem
SourceFiles
    laplacianProblem.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the laplacianProblem class.
/// \dir
/// Directory containing the header and source files for the laplacianProblem class.


#ifndef laplacianProblem_H
#define laplacianProblem_H
#include "fvCFD.H"
#include "IOmanip.H"
#include "Time.H"
#include "reductionProblem.H"
#include <iostream>
#include "ITHACAPOD.H"
#include "ITHACAutilities.H"
#include "Modes.H"

/// Class to implement a full order laplacian parametrized problem
class laplacianProblem: public reductionProblem
{

    public:
        // Constructors
        laplacianProblem();
        /// Construct with argc and argv
        laplacianProblem(int argc, char* argv[]);
        ~laplacianProblem() {};

        // Member Functions
        /// List of snapshots for the solution
        PtrList<volScalarField> Tfield;

        /// List of snapshots for the solution
        PtrList<volScalarField> Tonline;

        /// List of POD modes
        volScalarModes Tmodes;

        /// List of operators
        PtrList<fvScalarMatrix> operator_list;

        /// Theta (coefficients of the affine expansion)
        List<scalar> theta;

        /// Nu (diffusivity)
        PtrList<volScalarField> nu_list;

        /// Number of modes reduced problem
        label NTmodes;

        // Reduced Matrices
        /// A matrices
        List<Eigen::MatrixXd> A_matrices;
        /// Source vector
        Eigen::MatrixXd source;

        /// Other Variables
        label counter = 1;

        // Dummy variables to transform laplacianFoam into a class
        /// Temperature field
        autoPtr<volScalarField> _T;
        /// Source Term
        autoPtr<volScalarField> _S;
        /// Diffusivity
        autoPtr<volScalarField> _nu;
        /// Mesh
        mutable autoPtr<fvMesh> _mesh;
        /// Time
        autoPtr<Time> _runTime;

        // Functions

        //--------------------------------------------------------------------------
        /// Perform a truthsolve
        void truthSolve(List<scalar> mu_now, word folder = "./ITHACAoutput/Offline/");

        //--------------------------------------------------------------------------
        /// Perform a projection onto the POD modes
        ///
        /// @param[in]  Nmodes  The number of modes used for the projection
        ///
        void project(label Nmodes);
};

#endif
