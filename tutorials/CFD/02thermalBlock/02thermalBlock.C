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
Description
    Example of a heat transfer Reduction Problem
SourceFiles
    02thermalBlock.C
\*---------------------------------------------------------------------------*/

#include <iostream>
#include "fvCFD.H"
#include "IOmanip.H"
#include "Time.H"
#include "laplacianProblem.H"
#include "ReducedLaplacian.H"
#include "ITHACAPOD.H"
#include "ITHACAutilities.H"
#include <Eigen/Dense>
#define _USE_MATH_DEFINES
#include <cmath>

/// \brief Class where the tutorial number 2 is implemented.
/// \details It is a child of the laplacianProblem class and some of its
/// functions are overridden to be adapted to the specific case.
class tutorial02: public laplacianProblem
{
    public:
        explicit tutorial02(int argc, char* argv[])
            :
            laplacianProblem(argc, argv),
            T(_T()),
            nu(_nu()),
            S(_S())
        {}
        //! [tutorial02]
        /// Temperature field
        volScalarField& T;
        /// Diffusivity field
        volScalarField& nu;
        /// Source term field
        volScalarField& S;

        /// It perform an offline Solve
        void offlineSolve(word folder = "./ITHACAoutput/Offline/")
        {
            if (offline)
            {
                ITHACAstream::read_fields(Tfield, "T", folder);
                mu_samples =
                    ITHACAstream::readMatrix(folder + "/mu_samples_mat.txt");
            }
            else
            {
                List<scalar> mu_now(9);
                scalar IF = 0;

                for (label i = 0; i < mu.rows(); i++)
                {
                    for (label j = 0; j < mu.cols() ; j++)
                    {
                        // mu_now[i] =
                        theta[j] = mu(i, j);
                    }

                    assignIF(T, IF);
                    Info << i << endl;
                    truthSolve(mu_now, folder);
                }
            }
        }

        /// Define the source term function
        void SetSource()
        {
            volScalarField yPos = T.mesh().C().component(vector::Y);
            volScalarField xPos = T.mesh().C().component(vector::X);
            forAll(S, counter)
            {
                S[counter] = Foam::sin(xPos[counter] / 0.9 * M_PI) + Foam::sin(
                                 yPos[counter] / 0.9 * M_PI);
            }
        }

        /// Compute the diffusivity in each subdomain
        void compute_nu()
        {
            nu_list.resize(9);
            volScalarField nu1(nu);
            volScalarField nu2(nu);
            volScalarField nu3(nu);
            volScalarField nu4(nu);
            volScalarField nu5(nu);
            volScalarField nu6(nu);
            volScalarField nu7(nu);
            volScalarField nu8(nu);
            volScalarField nu9(nu);
            Eigen::MatrixXd Box1(2, 3);
            Box1 << 0, 0, 0, 0.3, 0.3, 0.1;
            Eigen::MatrixXd Box2(2, 3);
            Box2 << 0.3, 0, 0, 0.6, 0.3, 0.1;
            Eigen::MatrixXd Box3(2, 3);
            Box3 << 0.6, 0, 0, 0.91, 0.3, 0.1;
            Eigen::MatrixXd Box4(2, 3);
            Box4 << 0, 0.3, 0, 0.3, 0.6, 0.1;
            Eigen::MatrixXd Box5(2, 3);
            Box5 << 0.3, 0.3, 0, 0.6, 0.6, 0.1;
            Eigen::MatrixXd Box6(2, 3);
            Box6 << 0.6, 0.3, 0, 0.91, 0.6, 0.1;
            Eigen::MatrixXd Box7(2, 3);
            Box7 << 0, 0.6, 0, 0.3, 0.91, 0.1;
            Eigen::MatrixXd Box8(2, 3);
            Box8 << 0.3, 0.61, 0, 0.6, 0.91, 0.1;
            Eigen::MatrixXd Box9(2, 3);
            Box9 << 0.6, 0.6, 0, 0.9, 0.91, 0.1;
            ITHACAutilities::setBoxToValue(nu1, Box1, 1.0);
            ITHACAutilities::setBoxToValue(nu2, Box2, 1.0);
            ITHACAutilities::setBoxToValue(nu3, Box3, 1.0);
            ITHACAutilities::setBoxToValue(nu4, Box4, 1.0);
            ITHACAutilities::setBoxToValue(nu5, Box5, 1.0);
            ITHACAutilities::setBoxToValue(nu6, Box6, 1.0);
            ITHACAutilities::setBoxToValue(nu7, Box7, 1.0);
            ITHACAutilities::setBoxToValue(nu8, Box8, 1.0);
            ITHACAutilities::setBoxToValue(nu9, Box9, 1.0);
            nu_list.set(0, nu1);
            nu_list.set(1, nu2);
            nu_list.set(2, nu3);
            nu_list.set(3, nu4);
            nu_list.set(4, nu5);
            nu_list.set(5, nu6);
            nu_list.set(6, nu7);
            nu_list.set(7, nu8);
            nu_list.set(8, nu9);
        }

        /// Construct the operator_list where each term of the affine decomposition is stored
        void assemble_operator()
        {
            for (int i = 0; i < nu_list.size(); i++)
            {
                operator_list.append(fvm::laplacian(nu_list[i], T));
            }
        }

};


int main(int argc, char* argv[])
{
    // Create the example object of the tutorial02 type
    tutorial02 example(argc, argv);
    // Read some parameters from file
    ITHACAparameters* para = ITHACAparameters::getInstance(example._mesh(),
                             example._runTime());
    int NmodesTout = para->ITHACAdict->lookupOrDefault<int>("NmodesTout", 15);
    int NmodesTproj = para->ITHACAdict->lookupOrDefault<int>("NmodesTproj", 10);
    // Set the number of parameters
    example.Pnumber = 9;
    example.Tnumber = 500;
    // Set the parameters
    example.setParameters();
    // Set the parameter ranges, in all the subdomains the diffusivity varies between
    // 0.001 and 0.1
    example.mu_range.col(0) = Eigen::MatrixXd::Ones(9, 1) * 0.001;
    example.mu_range.col(1) = Eigen::MatrixXd::Ones(9, 1) * 0.1;
    // Generate the Parameters
    example.genRandPar(500);
    // Set the size of the list of values that are multiplying the affine forms
    example.theta.resize(9);
    // Set the source term
    example.SetSource();
    // Compute the diffusivity field for each subdomain
    example.compute_nu();
    // Assemble all the operators of the affine decomposition
    example.assemble_operator();
    // Perform an Offline Solve
    example.offlineSolve();
    // Perform a POD decomposition and get the modes
    ITHACAPOD::getModes(example.Tfield, example.Tmodes, example._T().name(),
                        example.podex, 0, 0,
                        NmodesTout);
    /// Create a new instance of the FOM problem for testing purposes
    tutorial02 FOM_test(argc, argv);
    FOM_test.offline = false;
    FOM_test.Pnumber = 9;
    FOM_test.Tnumber = 50;
    // Set the parameters
    FOM_test.setParameters();
    // Set the parameter ranges, in all the subdomains the diffusivity varies between
    // 0.001 and 0.1
    FOM_test.mu_range.col(0) = Eigen::MatrixXd::Ones(9, 1) * 0.001;
    FOM_test.mu_range.col(1) = Eigen::MatrixXd::Ones(9, 1) * 0.1;
    // Generate the Parameters
    FOM_test.genRandPar(50);
    // Set the size of the list of values that are multiplying the affine forms
    FOM_test.theta.resize(9);
    // Set the source term
    FOM_test.SetSource();
    // Compute the diffusivity field for each subdomain
    FOM_test.compute_nu();
    // Assemble all the operators of the affine decomposition
    FOM_test.assemble_operator();
    // Perform an Offline Solve
    FOM_test.offlineSolve("./ITHACAoutput/FOMtest");
    // Perform the Galerkin projection onto the space spanned by the POD modes
    example.project(NmodesTproj);
    // Create a reduced object
    reducedLaplacian reduced(example);

    // Solve the online reduced problem some new values of the parameters
    for (int i = 0; i < FOM_test.mu.rows(); i++)
    {
        reduced.solveOnline(FOM_test.mu.row(i));
    }

    // Reconstruct the solution and store it into Reconstruction folder
    reduced.reconstruct("./ITHACAoutput/Reconstruction");
    // Compute the error on the testing set
    Eigen::MatrixXd error = ITHACAutilities::errorL2Rel(FOM_test.Tfield,
                            reduced.Trec);
}
//--------
/// \dir 02thermalBlock Folder of the turorial 2
/// \file
/// \brief Implementation of a tutorial of a steady heat transfer problem

/// \example 02thermalBlock.C
/// \section intro_thermal Introduction to tutorial 2
/// The problems consists of a thermal block with a source term.
/// The problem equations are:
/// \f[
/// \nabla \cdot (k \nabla T) = S
/// \f]
/// where \f$k\f$ is the diffusivity, \f$T\f$ is the temperature and \f$S\f$ is the source term.
/// The problem discretised and formalized in matrix equation reads:
/// \f[
/// AT = S
/// \f]
/// where \f$A\f$ is the matrix of interpolation coefficients, \f$T\f$ is the vector of unknowns
/// and \f$S\f$ is the vector representing the source term.
/// The domain is subdivided in 9 different parts and each part has parametrized diffusivity. See the image below for a clarification.
/// \image html drawing.png
/// Both the full order and the reduced order
/// problem are solved exploiting the parametric affine decomposition of the differential operators:
/// \f[
/// A = \sum_{i=1}^N \theta_i(\mu) A_i
///  \f]
/// For the operations performed by each command check the comments in the source 02thermalBlock.C file.
///
/// \section code A look under the code
///
/// In this section are explained the main steps necessary to construct the tutorial N°2
///
/// \subsection header The necessary header files
///
/// First of all let's have a look to the header files that needs to be included and what they are responsible for:
///
/// The standard C++ header for input/output stream objects:
///
/// \dontinclude 02thermalBlock.C
///
/// \skip iostream
/// \until >
///
/// The OpenFOAM header files:
///
/// \skipline fvCFD
/// \until Time.H
///
/// The header file of ITHACA-FV necessary for this tutorial
///
/// \skipline ITHACAPOD
/// \until ITHACAuti
///
/// The Eigen library for matrix manipulation and linear and non-linear algebra operations:
///
/// \line Eigen
///
/// And we define some mathematical constants and include the standard header for common math operations:
///
/// \skipline MATH_DEF
/// \until cmath
///
/// \subsection classtuto02 Implementation of the tutorial02 class
///
/// Then we can define the tutorial02 class as a child of the laplacianProblem class
///
/// \skipline tutorial02
/// \until {}
///
/// The members of the class are the fields that needs to be manipulated during the
/// resolution of the problem
///
/// Inside the class it is defined the offline solve method according to the
/// specific parametrized problem that needs to be solved.
///
/// \skipline void
/// \until {
///
/// If the offline solve has already been performed than read the existing snapshots
///
/// \skipline if
/// \until
/// }
///
/// else perform the offline solve where a loop over all the parameters is performed:
/// \skipline for
/// \until }
///
/// a 0 internal constant value is assigned before each solve command with the lines
///
/// \skipline assignIF
///
/// and the solve operation is performed, see also the laplacianProblem class for the definition of the methods
///
/// \skipline truthSolve
///
/// The we need also to implement a method to set/define the source term that may be problem dependent.
/// In this case the source term is defined with a hat function:
///
///
/// \skipline SetSource
/// \until }
/// \skipline }
///
/// Define by:
///
/// \f[ S = \sin(\frac{\pi}{L}\cdot x) + \sin(\frac{\pi}{L}\cdot y) \f]
///
/// where \f$L\f$ is the dimension of the thermal block which is equal to 0.9.
///
/// \image html hat.jpg
///
/// With the following is defined a method to set compute the parameter of the affine expansion:
///
/// \skipline compute_nu
/// \until {
///
/// The list of parameters is resized according to number of parametrized regions
/// \skipline nu_list
///
/// The nine different volScalarFields to identify the viscosity in each domain are initialized:
///
/// \skipline volScalarField
/// \until nu9
///
/// and the 9 different boxes are defined:
///
/// \skipline Box1
/// \until Box9
///
/// and for each of the defined boxes the relative diffusivity field is set to 1 inside the box and remain 0 elsewhere:
///
/// \skipline ITHACA
/// \until Box9
///
///  See also the ITHACAutilities::setBoxToValue for more details.
///
/// The list of diffusivity fields is set with:
///
/// \skipline nu_list
/// \until }
///
/// \subsection main Definition of the main function
///
/// Once the tutorial02 class is defined the main function is defined,
/// an example of type tutorial02 is constructed:
///
/// \skipline argv)
///
/// the number of parameter is set:
///
/// \skipline Pnumber
/// \skipline setParameters
///
/// the range of the parameters is defined:
///
/// \skipline mu_range
/// \skipline mu_range
///
/// and 500 random combinations of the parameters are generated:
///
/// \skipline genRandPar
///
/// the size of the list of values that are multiplying the affine forms is set:
///
/// \skipline theta.resize
///
/// the source term is defined, the compute_nu and assemble_operator functions are called
///
/// \skipline .SetSource
/// \skipline .compute_nu
/// \skipline .assemble_operator
///
/// then the Offline full order Solve is performed:
///
/// \skipline offlineSolve
///
/// Once the Offline solve is performed the modes ar obtained using the ITHACAPOD::getModes function:
///
/// \skipline ITHACAPOD
///
/// and the projection is performed onto the POD modes using 10 modes
///
/// \skipline project(10)
///
/// Once the projection is performed we can construct a reduced object:
///
/// \skipline reducedLaplacian
///
/// and solve the reduced problem for some values of the parameters:
///
/// \skipline for
/// \until }
///
/// Finally, once the online solve has been performed we can reconstruct the solution:
///
/// \skipline reconstruct
///
/// \section plaincode The plain program
/// Here there's the plain code
///
