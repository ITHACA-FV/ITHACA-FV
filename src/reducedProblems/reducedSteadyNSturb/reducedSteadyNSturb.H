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
    reducedSteadyNSturb

Description
    A reduced problem for the stationary NS equations

SourceFiles
    reducedSteadyNSturb.C

\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the reducedSteadyNSturb class.
/// \dir
/// Directory containing the header and source files for the reducedSteadyNSturb class

#ifndef reducedSteadyNSturb_H
#define reducedSteadyNSturb_H

#include "fvCFD.H"
#include "IOmanip.H"
#include "reducedProblem.H"
#include "reducedSteadyNS.H"
#include "steadyNSturb.H"
#include "ITHACAutilities.H"
#include <Eigen/Dense>
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>

/// Structure to implement a newton object for a stationary NS problem
struct newton_steadyNSturb: public newton_argument<double>
{
public:
    newton_steadyNSturb() {}
    newton_steadyNSturb(int Nx, int Ny, steadyNSturb& problem): newton_argument<double>(Nx, Ny),
        problem(&problem),
        Nphi_u(problem.NUmodes + problem.liftfield.size() + problem.NSUPmodes),
        Nphi_nut(problem.Nnutmodes),
        Nphi_p(problem.NPmodes),
        N_BC(problem.inletIndex.rows())
    {}

    int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const;
    int df(const Eigen::VectorXd &x,  Eigen::MatrixXd &fjac) const;

    int Nphi_u;
    int Nphi_nut;
    int Nphi_p;
    int N_BC;
    scalar nu;
    Eigen::VectorXd BC;
    steadyNSturb* problem;
    Eigen::VectorXd nu_c;
    std::vector<SPLINTER::RBFSpline*> SPLINES;
};


/*---------------------------------------------------------------------------*\
                        Class reducedProblem Declaration
\*---------------------------------------------------------------------------*/


/// Class where it is implemented a reduced problem for the <b> steady turbulent Navier-stokes </b> problem.
/** In this class are implemented the methods to solve a reduced version of the stationary turbulent Navier-stokes
equations. This class must be used together with the steadyNS class */
class reducedSteadyNSturb: public reducedSteadyNS
{
private:

public:
    // Constructors
    /// Construct Null
    reducedSteadyNSturb();

    reducedSteadyNSturb(steadyNSturb& problem);

    ~reducedSteadyNSturb() {};


    /// List of snapshots for the solution for eddy viscosity
    PtrList<volScalarField> nutFields;

    /// List of POD modes for eddy viscosity
    PtrList<volScalarField> nuTmodes;

    /// Reconstructed eddy viscosity field
    PtrList<volScalarField> nutREC;

    /// Number of viscosity modes
    int Nphi_nut;

    /// Newton Object to solve the nonlinear problem
    newton_steadyNSturb newton_object;

    /// Pointer to the FOM problem
    steadyNSturb* problem;

    /// Method to reconstruct a solution from an online solve with a PPE stabilisation technique.
    /// stabilisation method
    ///
    /// @param      problem     a steadyNS full order problem.
    /// @param[in]  folder      The folder where you want to store the results (default is "./ITHACAOutput/online_rec").
    /// @param[in]  printevery  Variable to recover only every printevery online solutions default is 1.
    ///
    void solveOnline_sup(Eigen::MatrixXd vel_now);

    /// Method to reconstruct a solution from an online solve with a supremizer stabilisation technique.
    /// stabilisation method
    ///
    /// @param      problem     a steadyNS full order problem.
    /// @param[in]  folder      The folder where you want to store the results (default is "./ITHACAOutput/online_rec").
    /// @param[in]  printevery  Variable to recover only every printevery online solutions default is 1.
    ///
    void reconstruct_sup(fileName folder = "./ITHACAOutput/online_rec", int printevery = 1);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



#endif






