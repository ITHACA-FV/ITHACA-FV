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

\*---------------------------------------------------------------------------*/

/// \file
/// Source file of the inverseLaplacianProblemTotalHeatMeasure_CG class.


#include "inverseLaplacianProblemTotalHeatMeasure_CG.H"

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructors
inverseLaplacianProblemTotalHeatMeasure_CG::inverseLaplacianProblemTotalHeatMeasure_CG() {}

inverseLaplacianProblemTotalHeatMeasure_CG::inverseLaplacianProblemTotalHeatMeasure_CG(
    int argc, char* argv[])
    :
    inverseLaplacianProblem_CG::inverseLaplacianProblem_CG(argc, argv)
{
}

// * * * * * * * * * * * * * * Full Order Methods * * * * * * * * * * * * * * //

int inverseLaplacianProblemTotalHeatMeasure_CG::conjugateGradient()
{
    fvMesh& mesh = _mesh();
    set_g();
    set_valueFraction();
    cgIter = 0;
    J = 0;
    P = g;
    gradJ = g;       //Gradient of the cost function [W/m2]
    gamma = 0.0;
    gamma_den = 0.0;
    label sampleI = 1;
    gList.resize(0);
    Tfield.resize(0);
    lambdaField.resize(0);
    deltaTfield.resize(0);
    M_Assert(std::abs(gIntegral_meas) > 1e-16, "First set up gIntegral_meas");
    M_Assert(std::abs(gIntegralWeight) > 1e-16, "First set up gIntegralWeight");
    M_Assert(!interpolation, "Interpolation not implemented yet");

    while (cgIter < cgIterMax)
    {
        Info << "Iteration " << cgIter + 1 << endl;
        restart();
        solveDirect();
        gIntegral = ITHACAutilities::integralOnPatch(mesh, g, "hotSide");

        if (saveSolInLists && cgIter == 0)
        {
            gList.append(g.clone());
        }

        volScalarField& T = _T();
        ITHACAstream::exportSolution(T, std::to_string(sampleI),
                                     "./ITHACAoutput/CGtest/", T.name());
        differenceBetweenDirectAndMeasure();

        if (conjugateGradientConvergenceCheck())
        {
            Jlist.conservativeResize(cgIter + 1, 1);
            Jlist(cgIter) = J;
            ITHACAstream::exportMatrix(Jlist, "costFunctionFull", "eigen", "./");
            return (1);
        }

        Jlist.conservativeResize(cgIter + 1, 1);
        Jlist(cgIter) = J;
        solveAdjoint();
        volScalarField& lambda = _lambda();
        ITHACAstream::exportSolution(lambda, std::to_string(sampleI),
                                     "./ITHACAoutput/CGtest/", lambda.name());
        computeGradJ();
        searchDirection();
        solveSensitivity();
        volScalarField& deltaT = _deltaT();
        ITHACAstream::exportSolution(deltaT, std::to_string(sampleI),
                                     "./ITHACAoutput/CGtest/", deltaT.name());
        sensibilitySolAtThermocouplesLocations();
        computeSearchStep();
        updateHeatFlux();

        if (saveSolInLists)
        {
            volScalarField& T = _T();
            volScalarField& lambda = _lambda();
            volScalarField& deltaT = _deltaT();
            gList.append(g.clone());
            Tfield.append(T.clone());
            lambdaField.append(lambda.clone());
            deltaTfield.append(deltaT.clone());
            sampleI++;
        }

        cgIter++;
    }

    return (0);
}

void inverseLaplacianProblemTotalHeatMeasure_CG::computeGradJ()
{
    fvMesh& mesh = _mesh();
    volScalarField& lambda = _lambda();
    gradJ_L2norm = 0;
    forAll (lambda.boundaryField()[hotSide_ind], faceI)
    {
        gradJ [faceI] = - lambda.boundaryField()[hotSide_ind][faceI]  +
                        gIntegralWeight * ( gIntegral - gIntegral_meas );
        gradJ_L2norm += gradJ[faceI] * gradJ[faceI]  *
                        mesh.magSf().boundaryField()[hotSide_ind][faceI];
    }
    gradJ_L2norm = Foam::sqrt(gradJ_L2norm);
    Info << "gradJ L2norm = " << gradJ_L2norm << endl;
}

void inverseLaplacianProblemTotalHeatMeasure_CG::computeSearchStep()
{
    fvMesh& mesh = _mesh();
    double Pintegral = ITHACAutilities::integralOnPatch(mesh, P, "hotSide");
    beta = Tdiff.dot(Tsens) - gIntegralWeight * Pintegral *
           (gIntegral - gIntegral_meas);
    double betaDiv = Tsens.dot(Tsens) - gIntegralWeight * Pintegral * Pintegral;
    beta = beta / betaDiv;
    Info << "beta = " << beta << endl;
}

int inverseLaplacianProblemTotalHeatMeasure_CG::conjugateGradientConvergenceCheck()
{
    double Jold = J;
    J = 0.5 * Tdiff.dot(Tdiff) + 0.5 * gIntegralWeight * (gIntegral -
            gIntegral_meas) * (gIntegral - gIntegral_meas);
    //reduce(J, sumOp<double>());
    Info << "J = " << J << endl;

    if (J <= Jtol)
    {
        Info << "Convergence reached in " << cgIter << " iterations" << endl;
        return (1);
    }
    else if (Foam::mag((Jold - J) / J) <= JtolRel)
    {
        Info << "Relative tolerance criteria meet in " << cgIter << " iterations" <<
             endl;
        Info << "|Jold - J| / |J| = " << Foam::mag((Jold - J) / J) << endl;
        return (1);
    }
    else
    {
        return (0);
    }
}



