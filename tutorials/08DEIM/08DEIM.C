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
    Example of the reconstruction of a non-linear function using the DEIM
SourceFiles
    08DEIM.C
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include "simpleControl.H"
#include "ITHACAutilities.H"
#include "Foam2Eigen.H"
#include <chrono>
#include "fvMeshSubset.H"
#include "ITHACAstream.H"
#include "ITHACAPOD.H"
#include "DEIM.H"
#include <chrono>

class DEIM_function : public DEIM<PtrList<volScalarField>, volScalarField >
{
public:
    using DEIM::DEIM;
    static volScalarField evaluate_expression(volScalarField& S, Eigen::MatrixXd mu)
    {
        volScalarField yPos = S.mesh().C().component(vector::Y);
        volScalarField xPos = S.mesh().C().component(vector::X);
        for (auto i = 0; i < S.size(); i++)
        {
            S[i] = std::exp( - 2 * std::pow(xPos[i] - mu(0) - 1, 2) - 2 * std::pow(yPos[i] - mu(1) - 0.5 , 2));
        }
        return S;
    }
    Eigen::VectorXd onlineCoeffs(Eigen::MatrixXd mu)
    {
        theta.resize(fields.size());
        for (int i = 0; i < fields.size(); i++)
        {
            double on_coeff = evaluate_expression(fields[i], mu)[localMagicPoints[i]];
            theta(i) = on_coeff;
        }
        return theta;
    }
};

int main(int argc, char *argv[])
{
#include "setRootCase.H"
#include "createTime.H"
#include "createMesh.H"
    simpleControl simple(mesh);
#include "createFields.H"

    // List of volScalarField where the snapshots are stored
    PtrList<volScalarField> Sp;

    // Non linear function
    volScalarField S
    (
        IOobject
        (
            "S",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        T.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 0, -1, 1, 0, 0, 0), 0)
    );

    // Parameters used to train the non-linear function
    Eigen::MatrixXd pars = ITHACAutilities::rand(100, 2, -0.5, 0.5);

    // Perform the offline phase
    for (int i = 0; i < 100; i++)
    {
        DEIM_function::evaluate_expression(S, pars.row(i));
        Sp.append(S);
        ITHACAutilities::exportSolution(S, "./ITHACAoutput/Offline/" , name(i + 1));
    }

    // Create DEIM object with given number of basis functions
    DEIM_function c(Sp, 30, "Gaussiana");

    // Generate the submeshes with the depth of the layer
    c.generateSubmeshes(2, mesh, S);

    // Define a new online parameter
    Eigen::MatrixXd par_new(2, 1);
    par_new(0, 0) = 0;
    par_new(1, 0) = 0;

    // Online evaluation of the non linear function
    Eigen::VectorXd aprfield = c.MatrixOnline * c.onlineCoeffs(par_new);

    // Transform to an OpenFOAM field and export
    volScalarField S2("S_online", Foam2Eigen::Eigen2field(S, aprfield));
    ITHACAutilities::exportSolution(S2, "./ITHACAoutput/Online/" , name(1));

    // Evaluate the full order function and export it
    DEIM_function::evaluate_expression(S, par_new);
    ITHACAutilities::exportSolution(S, "./ITHACAoutput/Online/" , name(1));

    return 0;
}

/// \dir 08DEIM Folder of the turorial 8
/// \file 
/// \brief Implementation of tutorial 8 for an unsteady Navier-Stokes problem

/// \example 08DEIM.C
/// \section intro_08DEIM Introduction to tutorial 8
/// In this tutorial we implement test
///
/// The following image illustrates blabla
/// \image html cylinder.png
///
/// \section code08 A detailed look into the code
///
/// In this section we explain the main steps necessary to construct the tutorial N°9
/// 
/// \subsection header ITHACA-FV header files
///
/// First of all let's have a look at the header files that need to be included and what they are responsible for.
/// 
/// The header files of ITHACA-FV necessary for this tutorial are: <unsteadyNS.H> for the full order unsteady NS problem,
/// <ITHACAPOD.H> for the POD decomposition, <reducedUnsteadyNS.H> for the construction of the reduced order problem,
/// and finally <ITHACAstream.H> for some ITHACA input-output operations.
/// 
/// \section plaincode The plain program
/// Here there's the plain code


