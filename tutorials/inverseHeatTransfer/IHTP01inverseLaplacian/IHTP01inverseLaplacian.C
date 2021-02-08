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
    example_paramBC.of a heat transfer Reduction Problem
SourceFiles
    IHTP01inverseLaplacian.C
\*---------------------------------------------------------------------------*/

#include <iostream>
#include <float.h>
#include "interpolation.H"
#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include "IOmanip.H"
#include "Time.H"
#include "inverseLaplacianProblem_CG.H"
#include "inverseLaplacianProblem_paramBC.H"
#include "ITHACAPOD.H"
#include "ITHACAutilities.H"
#include <Eigen/Dense>
#define _USE_MATH_DEFINES
#include <cmath>
#include "Foam2Eigen.H"
#include "mixedFvPatchFields.H"
#include "cellDistFuncs.H"

#include "IHTP01inverseLaplacian.H"


int main(int argc, char* argv[])
{
    solverPerformance::debug = 1; //No verbose output
    IHTP01inverseLaplacian_paramBC example_paramBC(argc, argv);
    IHTP01inverseLaplacian_CG example_CG(argc, argv);
    //Setting parameters for the analytical benchmark
    double a = 5;
    double b = 10;
    double c = 15;
    double d = 20;
    example_paramBC.a = a;
    example_paramBC.b = b;
    example_paramBC.c = c;
    example_paramBC.d = d;
    example_CG.a = a;
    example_CG.b = b;
    example_CG.c = c;
    example_CG.d = d;
    // Reading tests to perform
    ITHACAparameters* para = ITHACAparameters::getInstance(example_paramBC._mesh(),
                             example_paramBC._runTime());
    label CGtest = para->ITHACAdict->lookupOrDefault<int>("CGtest", 0);
    label CGnoiseTest = para->ITHACAdict->lookupOrDefault<int>("CGnoiseTest", 0);
    label CGnoiseLevelTest =
        para->ITHACAdict->lookupOrDefault<int>("CGnoiseLevelTest", 0);
    label parameterizedBCtest =
        para->ITHACAdict->lookupOrDefault<int>("parameterizedBCtest", 0);
    label parameterizedBCtest_RBFwidth =
        para->ITHACAdict->lookupOrDefault<int>("parameterizedBCtest_RBFwidth", 0);
    label thermocouplesLocationTest_CG =
        para->ITHACAdict->lookupOrDefault<int>("thermocouplesLocationTest_CG", 0);
    label thermocouplesLocationTest_paramBC =
        para->ITHACAdict->lookupOrDefault<int>("thermocouplesLocationTest_paramBC", 0);
    label thermocouplesNumberTest_CG =
        para->ITHACAdict->lookupOrDefault<int>("thermocouplesNumberTest_CG", 0);
    label thermocouplesNumberTest_paramBC =
        para->ITHACAdict->lookupOrDefault<int>("thermocouplesNumberTest_paramBC", 0);
    // Reading parameters from ITHACAdict
    example_CG.cgIterMax = para->ITHACAdict->lookupOrDefault<int>("cgIterMax", 100);
    example_CG.Jtol =  para->ITHACAdict->lookupOrDefault<double>("Jtolerance",
                       0.000001);
    example_CG.JtolRel =
        para->ITHACAdict->lookupOrDefault<double>("JrelativeTolerance",
                0.001);
    double rbfShapePar = para->ITHACAdict->lookupOrDefault<double>("rbfShapePar",
                         0);
    M_Assert(rbfShapePar > 0, "rbfShapePar not specified");
    example_paramBC.k =
        para->ITHACAdict->lookupOrDefault<double>("thermalConductivity", 0);
    M_Assert(example_paramBC.k > 0, "thermalConductivity, k, not specified");
    example_paramBC.H =
        para->ITHACAdict->lookupOrDefault<double>("heatTranferCoeff", 0);
    M_Assert(example_paramBC.H > 0, "Heat transfer coeff, H, not specified");
    example_CG.k = example_paramBC.k;
    example_CG.H = example_paramBC.H;
    int rbfWidthTest_size =
        para->ITHACAdict->lookupOrDefault<int>("rbfWidthTest_size", 0);
    // Setting analytical solution
    volScalarField T_true(example_paramBC._T());

    for (label i = 0; i < T_true.internalField().size(); i++)
    {
        auto cx = T_true.mesh().C()[i].component(vector::X);
        auto cy = T_true.mesh().C()[i].component(vector::Y);
        auto cz = T_true.mesh().C()[i].component(vector::Z);
        T_true.ref()[i] = a * cx * cx + b * cx * cy + c * cy - a * cz * cz + c;
    }

    // Perform true solution
    fvMesh& mesh = example_paramBC._mesh();
    example_paramBC.hotSide_ind = mesh.boundaryMesh().findPatchID("hotSide");
    label hotSideSize = T_true.boundaryField()[example_paramBC.hotSide_ind].size();
    example_paramBC.g.resize(hotSideSize);
    example_paramBC.gTrue.resize(hotSideSize);
    forAll(example_paramBC.g, faceI)
    {
        scalar faceX =
            mesh.boundaryMesh()[example_paramBC.hotSide_ind].faceCentres()[faceI].x();
        example_paramBC.g[faceI] = example_paramBC.k * (b * faceX + c) ;
    }
    example_paramBC.gTrue = example_paramBC.g;
    example_paramBC.solveTrue();
    example_CG.g = example_paramBC.g;
    example_CG.gTrue = example_paramBC.g;
    example_CG.solveTrue();
    volScalarField& T(example_paramBC._T());
    Info << "Exporting analytical solution" << endl;
    ITHACAstream::exportSolution(T_true, "1", "./ITHACAoutput/true/",
                                 "analyticalSol");
    volScalarField error = T_true - T;
    ITHACAstream::exportSolution(error, "1", "./ITHACAoutput/true/", "error");
    Info << "L2 norm of the relative error = " << ITHACAutilities::errorL2Rel(
             T_true, T) << endl;
    Info << "Linf norm of the relative error = " <<
         ITHACAutilities::errorLinfRel(T_true, T) << endl;
    Info << "L2 norm of the error = " << ITHACAutilities::errorL2Abs(T_true,
            T) << endl;
    Info << "Linf norm of the error = " << ITHACAutilities::LinfNorm(error) << endl;
    Info << "L2 norm of T = " << ITHACAutilities::L2Norm(T) << endl;
    // Setting up The thermocouples
    example_paramBC.readThermocouples();
    example_paramBC.Tmeas = example_paramBC.fieldValueAtThermocouples(T_true);
    example_CG.readThermocouples();
    example_CG.Tmeas = example_CG.fieldValueAtThermocouples(T_true);

    // Inverse problem tests

    // Alifanov's regularization test
    if (CGtest)
    {
#include "CGtest.H"
    }

    // Parameterized heat flux test
    if (parameterizedBCtest)
    {
#include"parameterizedBCtest.H"
    }

    // Test that solves the inverse problem using the parameterization of the heat flux with different RBF shape parameters
    if (parameterizedBCtest_RBFwidth)
    {
#include"parameterizedBCtest_RBFwidth.H"
    }

    // Test the Alifanov's regularization changing the distance of the thermocouples from the hotSide
    if (thermocouplesLocationTest_CG)
    {
#include"thermocouplesLocation_CG.H"
    }

    // Test the parameterization method changing the distance of the thermocouples from the hotSide
    if (thermocouplesLocationTest_paramBC)
    {
#include"thermocouplesLocation_paramBC.H"
    }

    // Test the Alifanov's regularization changing the number of thermocouples in the thermocouples plane
    if (thermocouplesNumberTest_CG)
    {
#include"thermocouplesNumberTest_CG.H"
    }

    // Test the parameterization method changing the number of thermocouples in the thermocouples plane
    if (thermocouplesNumberTest_paramBC)
    {
#include"thermocouplesNumberTest_paramBC.H"
    }

    return 0;
}

//--------
/// \dir  IHTP01inverseLaplacian Folder of the UQ tutorial 2
/// \file
/// \brief Implementation of an inverse heat transfer problem. Given a set of pointwise temeprature measurements in the interior of the domain, we want to estimate the boundary heat flux.

/// \example IHTP01inverseLaplacian.C
/// \section intro_IHTP01inverseLaplacian Boundary heat flux estimation
/// The tutorial consists in the estimation of the boundary heat flux \f$g\f$ given a set of pointwise temperature measurements in the interior of the domain.
///
/// We consider the domain
///
/// \image html IHTP01inverseLaplacian_schematicDomain.png
///
/// In this tutorial, the direct problem is
///
///   \f{eqnarray*}{
///        -k \Delta T  &=& 0 &\text{ in } \Omega\\
///        -k \nabla T \cdot \mathbf{n}  &=& g &\text{ on } \Gamma_{s_{in}},\\
///        -k \nabla T \cdot \mathbf{n}  &=& q_L &\text{ on } \Gamma_L, L\in\{I, II, III, IV\},\\
///        -k \nabla T \cdot \mathbf{n}  &=& H(T - T_f) &\text{ on } \Gamma_{sf},\\
///   \f}
/// where
///   \f{eqnarray*}{
///        q_{I}  (\mathbf{x})     &= 2 k a H,           \quad && q_{III}(\mathbf{x}) = 0,\\
///        q_{II} (\mathbf{x})     &= -k ( 2 a L + b y), \quad && q_{IV} (\mathbf{x}) = k b y, \\
///        T_f    (\mathbf{x})     &=\frac{k (b x + c)}{h} + a x^2 + c y - a z^2 + b x W + c, &q_{I}  (\mathbf{x})     = 2 k a H.           \quad && \\
///   \f}
///
///   The objective of this tutorial is then to reconstruct the boundary heat flux $g$, given point wise temperature measurements in the interior of the domain.
/// We state this problem as using an optimization framework.
///
/// Let \f$\Psi:=\{\mathbf{x}_1, \mathbf{x}_2 , \dots, \mathbf{x}_M \}\f$ be a collection of points in \f$\Omega\f$.
/// We define the application \f$\mathbf{x}_i \in \Psi\rightarrow \hat{T}(\mathbf{x}_i)\in \rm I\!R^+\f$, \f$\hat{T}(\mathbf{x}_i)\f$ being the experimentally measured temperature at \f$\mathbf{x}_i \in \Psi\f$.
/// Then, we state the inverse problem as:
///
/// Given \f$\{ \hat{T}(\mathbf{x}_i) \}_{i=1}^M\f$, find the heat flux \f$g \in L^2(\Gamma_{s_{in}}\f$ that minimizes the functional \f$J:L^2(\Gamma_{s_{in}}) \rightarrow \rm I\!R^+\f$,
///   \f{eqnarray*}{
///         J[g]:=\frac{1}{2}\sum^M_{i=1} [T[g](\mathbf{x}_i) - \hat{T}(\mathbf{x}_i)]^2,
///   \f}
///  where \f$T[g](\mathbf{x}_i)\f$ is the solution of the direct problem at points \f$\mathbf{x}_i\f$, for all \f$i=1,2,\dots,M\f$.
///
/// To define an academic benchmark, we notice that if we select the boundary heat flux
/// \f{eqnarray*}{
///     g = g_{an} (\mathbf{x}) = k (b x + c),
/// \f}
/// the direct problem has the analytical solution
/// \f{eqnarray*}{
///     T_{an}(\mathbf{x})= a x^2+ b x y + c y - a z^2 + c.
/// \f}
/// We can then abitrarily locate virtual thermocouples points \f$ \{\mathbf{x}_1, \mathbf{x}_2 , \dots, \mathbf{x}_M \} \f$ in the domain and obtain the measured temperatures as
/// \f{eqnarray*}{
///     \hat{T}(\mathbf{x}_i) = T_{an}(\mathbf{x}_i).
/// \f}
///
/// For the solution of this inverse problem, we use two methods.
/// Namely, Alifanov's regularization and the parameterization method.
/// Both methods are described in https://arxiv.org/abs/2101.11985
///
/// \section plaincode The plain program
/// Here there's the plain code
///

