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
    Example of a compressible flow Reduction Problem
SourceFiles
    24rhoCentralFoam.C
\*---------------------------------------------------------------------------*/

//#include <iostream>
//#include "fvCFD.H"
//#include "IOmanip.H"
//#include "Time.H"
#include "UnsteadyCompressibleNS.H"
//#include "ReducedUnsteadyNS.H"
#include "ITHACAPOD.H"
//#include "ITHACAutilities.H"
//#include <Eigen/Dense>
//#define _USE_MATH_DEFINES
//#include <cmath>

/// \brief Class where the tutorial number 2 is implemented.
/// \details It is a child of the unsteadyNS class and some of its
/// functions are overridden to be adapted to the specific case.
class tutorial24: public UnsteadyCompressibleNS
{
    public:
        explicit tutorial24(int argc, char* argv[])
            :
            UnsteadyCompressibleNS(argc, argv),
            p(_p()),
            T(_T()),
            U(_U())
        {}

        // Fields to perform
        volScalarField& p;
        volScalarField& T;
        volVectorField& U;

        void offlineSolve(word folder = "./ITHACAoutput/Offline/")
        {
            if (offline)
            {
                ITHACAstream::read_fields(Pfield, p, folder);
                ITHACAstream::read_fields(Tfield, T, folder);
                ITHACAstream::read_fields(Ufield, U, folder);
            }
            else
            {
                Info << "WIP: Online solve" << endl;
                exit(0);
            }
        }
};


int main(int argc, char* argv[])
{
    // Create the example object of the tutorial24 type
    tutorial24 example(argc, argv);
    // Read some parameters from file
    ITHACAparameters* para = ITHACAparameters::getInstance(example._mesh(),
                             example._runTime());
    int nModes = para->ITHACAdict->lookupOrDefault<int>("nModes", 10);
    // Perform an Offline Solve
    example.offlineSolve();

    // Get Modes
    ITHACAPOD::getModes(example.Pfield, example.Pmodes, example._p().name(),
                        example.podex, 0, 0,
                        nModes);
    ITHACAPOD::getModes(example.Tfield, example.Tmodes, example._T().name(),
                        example.podex, 0, 0,
                        nModes);
    ITHACAPOD::getModes(example.Ufield, example.Umodes, example._U().name(),
                        example.podex, 0, 0,
                        nModes);            

    // Project modes and reconstruct solution           
    PtrList<volScalarField> projectedSnapshotsP, projectedSnapshotsT;
    PtrList<volVectorField> projectedSnapshotsU;

    example.Pmodes.projectSnapshots(example.Pfield, projectedSnapshotsP, nModes);
    example.Tmodes.projectSnapshots(example.Tfield, projectedSnapshotsT, nModes);
    example.Umodes.projectSnapshots(example.Ufield, projectedSnapshotsU, nModes);

    ITHACAstream::exportFields(projectedSnapshotsP, "./ITHACAoutput/Reconstruction", "P");
    ITHACAstream::exportFields(projectedSnapshotsT, "./ITHACAoutput/Reconstruction", "T");
    ITHACAstream::exportFields(projectedSnapshotsU, "./ITHACAoutput/Reconstruction", "U");
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
/// \subsection classtuto02 Implementation of the tutorial24 class
///
/// Then we can define the tutorial24 class as a child of the unsteadyNS class
///
/// \skipline tutorial24
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
/// and the solve operation is performed, see also the unsteadyNS class for the definition of the methods
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
/// Once the tutorial24 class is defined the main function is defined,
/// an example of type tutorial24 is constructed:
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
/// \skipline .project(NmodesTproj)
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
