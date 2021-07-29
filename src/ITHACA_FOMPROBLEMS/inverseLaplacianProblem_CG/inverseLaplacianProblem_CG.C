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
/// Source file of the inverseLaplacianProblem_CG class.


#include "inverseLaplacianProblem_CG.H"

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructors
inverseLaplacianProblem_CG::inverseLaplacianProblem_CG() {}

inverseLaplacianProblem_CG::inverseLaplacianProblem_CG(int argc, char* argv[])
    :
    inverseLaplacianProblem::inverseLaplacianProblem(argc, argv)
{
    Foam::Time& runTime = _runTime();
    Foam::fvMesh& mesh = _mesh();
#include "createFields.H"
}

// * * * * * * * * * * * * * * Full Order Methods * * * * * * * * * * * * * * //

void inverseLaplacianProblem_CG::set_valueFraction()
{
    fvMesh& mesh = _mesh();
    valueFraction.resize(mesh.boundaryMesh()["coldSide"].size());
    homogeneousBCcoldSide.resize(mesh.boundaryMesh()["coldSide"].size());
    valueFractionAdj.resize(mesh.boundaryMesh()["coldSide"].size());
    Eigen::VectorXd faceCellDist =
        ITHACAutilities::boudaryFaceToCellDistance(mesh, coldSide_ind);
    forAll (valueFraction, faceI)
    {
        scalar faceDist = faceCellDist(faceI);
        valueFraction[faceI] =  1.0 / (1.0 + (k / H / faceDist));
        valueFractionAdj[faceI] =  1 / (1 + (1 / k / H / faceDist));
        homogeneousBCcoldSide[faceI] =  0;
    }
    refGrad = homogeneousBCcoldSide;
}

void inverseLaplacianProblem_CG::assignAdjointBC()
{
    fvMesh& mesh = _mesh();
    volScalarField& lambda = _lambda();
    set_valueFraction();
    forAll(mesh.boundaryMesh(), patchI)
    {
        if (patchI == mesh.boundaryMesh().findPatchID("coldSide"))
        {
            ITHACAutilities::assignMixedBC(lambda, patchI, homogeneousBCcoldSide, refGrad,
                                           valueFraction);
        }
        else if (patchI == mesh.boundaryMesh().findPatchID("hotSide"))
        {
            ITHACAutilities::assignBC(lambda, patchI, homogeneousBC);
        }
        else
        {
            ITHACAutilities::assignBC(lambda, patchI, homogeneousBC);
        }
    }
    ITHACAutilities::assignIF(lambda, homogeneousBC);
}

volScalarField inverseLaplacianProblem_CG::assignAdjointBCandSource()
{
    volScalarField& lambda = _lambda();
    fvMesh& mesh = _mesh();
    assignAdjointBC();
    dimensionedScalar sourceDim("sourceDim", dimensionSet(1, -1, -3, -1, 0, 0, 0),
                                1);
    autoPtr<volScalarField> f_
    (
        new volScalarField("f", lambda)
    );
    volScalarField& f = f_();

    if (interpolation)
    {
        forAll(interpolationPlane.cellID, cellI)
        {
            f.ref()[interpolationPlane.cellID [cellI]] =
                interpolationPlane.Tdiff[cellI] * k;
        }
    }
    else
    {
        for (int i = 0; i < thermocouplesCellID.size(); i++)
        {
            f.ref()[thermocouplesCellID [i]] = Tdiff(i) * k /
                                               mesh.V()[thermocouplesCellID [i]];
        }
    }

    return (f * sourceDim).ref();
}

void inverseLaplacianProblem_CG::assignSensitivityBC()
{
    fvMesh& mesh = _mesh();
    volScalarField& deltaT = _deltaT();
    set_valueFraction();
    forAll(mesh.boundaryMesh(), patchI)
    {
        if (patchI == mesh.boundaryMesh().findPatchID("coldSide"))
        {
            ITHACAutilities::assignMixedBC(deltaT, patchI, homogeneousBCcoldSide, refGrad,
                                           valueFraction);
        }
        else if (patchI == mesh.boundaryMesh().findPatchID("hotSide"))
        {
            ITHACAutilities::assignBC(deltaT, patchI, - P / k);
        }
        else
        {
            ITHACAutilities::assignBC(deltaT, patchI, homogeneousBC);
        }
    }
}

void inverseLaplacianProblem_CG::solveAdjoint()
{
    restart();
    volScalarField& lambda = _lambda();
    Foam::Time& runTime = _runTime();
    volScalarField f = assignAdjointBCandSource();
    simpleControl& simple = _simple();

#if defined(OFVER) && (OFVER == 6)
    while (simple.loop(runTime))
#else
    while (simple.loop())
#endif
    {
        //Info << "Time = " << runTime.timeName() << nl << endl;
        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix TEqn
            (
                fvm::laplacian(DT, lambda) ==  - f
            );
            TEqn.solve();
        }
    }
}


void inverseLaplacianProblem_CG::solveSensitivity()
{
    restart();
    assignSensitivityBC();
    solve("sensitivity");
}

void inverseLaplacianProblem_CG::solve(const char* problemID)
{
    volScalarField& T = _T();
    volScalarField& deltaT = _deltaT();
    simpleControl& simple = _simple();
    Foam::Time& runTime = _runTime();

    if (strcmp( problemID, "direct") == 0)
    {
        ITHACAutilities::assignIF(T, homogeneousBC);
    }
    else if (strcmp( problemID, "sensitivity") == 0)
    {
        ITHACAutilities::assignIF(deltaT, homogeneousBC);
    }
    else
    {
        Info << "Problem name should be direct or sensitivity" << endl;
        exit(10);
    }

#if defined(OFVER) && (OFVER == 6)
    while (simple.loop(runTime))
#else
    while (simple.loop())
#endif
    {
        while (simple.correctNonOrthogonal())
        {
            if (strcmp( problemID, "direct") == 0)
            {
                fvScalarMatrix TEqn
                (
                    fvm::laplacian(DT, T)
                );
                TEqn.solve();
            }
            else if (strcmp( problemID, "sensitivity") == 0)
            {
                fvScalarMatrix TEqn
                (
                    fvm::laplacian(DT, deltaT)
                );
                TEqn.solve();
            }
        }
    }
}

void inverseLaplacianProblem_CG::defineThermocouplesPlane()
{
    Info << "Defining the plane for measurements interpolation" << endl;
    fvMesh& mesh = _mesh();
    //Define thermocouples plane
    bool firstCell = 1;
    forAll(thermocouplesCellID, cellI)
    {
        if (thermocouplesCellID[cellI] != -1)
        {
            thermocouplesCellProc[cellI] = Pstream::myProcNo();

            if (firstCell)
            {
                interpolationPlane.minX = mesh.C()[thermocouplesCellID[cellI]].component(0);
                interpolationPlane.maxX = mesh.C()[thermocouplesCellID[cellI]].component(0);
                interpolationPlane.Y = mesh.C()[thermocouplesCellID[cellI]].component(1);
                interpolationPlane.minZ = mesh.C()[thermocouplesCellID[cellI]].component(2);
                interpolationPlane.maxZ = mesh.C()[thermocouplesCellID[cellI]].component(2);
                firstCell = 0;
            }

            if (mesh.C()[thermocouplesCellID[cellI]].component(0) < interpolationPlane.minX)
            {
                interpolationPlane.minX = mesh.C()[thermocouplesCellID[cellI]].component(0);
            }

            if (mesh.C()[thermocouplesCellID[cellI]].component(0) > interpolationPlane.maxX)
            {
                interpolationPlane.maxX = mesh.C()[thermocouplesCellID[cellI]].component(0) ;
            }

            if (mesh.C()[thermocouplesCellID[cellI]].component(2) < interpolationPlane.minZ)
            {
                interpolationPlane.minZ = mesh.C()[thermocouplesCellID[cellI]].component(2);
            }

            if (mesh.C()[thermocouplesCellID[cellI]].component(2) > interpolationPlane.maxZ)
            {
                interpolationPlane.maxZ = mesh.C()[thermocouplesCellID[cellI]].component(2) ;
            }
        }
        else
        {
            Tmeas (cellI) = 0;
            thermocouplesCellProc[cellI] = -1;
        }
    }
}

void inverseLaplacianProblem_CG::sensibilitySolAtThermocouplesLocations()
{
    volScalarField& deltaT = _deltaT();

    if (interpolation)
    {
        forAll(interpolationPlane.cellID, cellI)
        {
            interpolationPlane.Tsens[cellI] =
                deltaT.internalField()[interpolationPlane.cellID [cellI]];
        }
    }
    else
    {
        Tsens = fieldValueAtThermocouples(deltaT);
    }
}

int inverseLaplacianProblem_CG::conjugateGradient()
{
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

    while (cgIter < cgIterMax)
    {
        Info << "Iteration " << cgIter + 1 << endl;
        restart();
        solveDirect();

        if (saveSolInLists && cgIter == 0)
        {
            gList.append(g.clone());
        }

        volScalarField& T = _T();
        //ITHACAstream::exportSolution(T, std::to_string(sampleI),
        //                             "./ITHACAoutput/CGtest/", T.name());
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
        //ITHACAstream::exportSolution(lambda, std::to_string(sampleI),
        //                             "./ITHACAoutput/CGtest/", lambda.name());
        computeGradJ();
        searchDirection();
        solveSensitivity();
        volScalarField& deltaT = _deltaT();
        //ITHACAstream::exportSolution(deltaT, std::to_string(sampleI),
        //                             "./ITHACAoutput/CGtest/", deltaT.name());
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

void inverseLaplacianProblem_CG::computeGradJ()
{
    fvMesh& mesh = _mesh();
    volScalarField& lambda = _lambda();
    gradJ_L2norm = 0;
    forAll (lambda.boundaryField()[hotSide_ind], faceI)
    {
        gradJ [faceI] = - lambda.boundaryField()[hotSide_ind][faceI];
        gradJ_L2norm += gradJ[faceI] * gradJ[faceI]  *
                        mesh.magSf().boundaryField()[hotSide_ind][faceI];
    }
    gradJ_L2norm = Foam::sqrt(gradJ_L2norm);
    Info << "gradJ L2norm = " << gradJ_L2norm << endl;
}

void inverseLaplacianProblem_CG::searchDirection()
{
    fvMesh& mesh = _mesh();
    gamma = 0.0;
    scalar gammaNum = 0;
    gammaNum = gradJ_L2norm;

    if (cgIter > 0)
    {
        //reduce(gamma, sumOp<double>());
        gamma = gammaNum / gamma_den;
        Info << "gamma = " << gamma << endl;
    }

    P = gradJ + gamma * P; //Updating P
    gamma_den = gammaNum;
    //reduce(gamma_den, sumOp<double>());
}

void inverseLaplacianProblem_CG::computeSearchStep()
{
    if (interpolation)
    {
        List<scalar> temp = interpolationPlane.Tdiff *
                            interpolationPlane.Tsens;
        beta = 0.0;
        forAll(interpolationPlane.cellVol, cellI)
        {
            beta += interpolationPlane.cellVol [cellI] * temp [cellI];
        }
        //reduce(beta, sumOp<double>());
        temp = interpolationPlane.Tsens * interpolationPlane.Tsens;
        scalar betaDiv = 0.0;
        forAll(interpolationPlane.cellVol, cellI)
        {
            betaDiv += interpolationPlane.cellVol [cellI] * temp [cellI];
        }
        //reduce(betaDiv, sumOp<double>());
        beta = beta / betaDiv;
        temp.clear();
    }
    else
    {
        beta = Tdiff.dot(Tsens);
        //reduce(beta, sumOp<double>());
        double betaDiv = Tsens.dot(Tsens);
        //reduce(betaDiv, sumOp<double>());
        beta = beta / betaDiv;
    }

    Info << "beta = " << beta << endl;
}

void inverseLaplacianProblem_CG::updateHeatFlux()
{
    g = g - beta * P;
}


int inverseLaplacianProblem_CG::conjugateGradientConvergenceCheck()
{
    double Jold = J;

    if (interpolation)
    {
        List<scalar> sqTdiff;
        sqTdiff = interpolationPlane.Tdiff * interpolationPlane.Tdiff;
        J = 0.0;
        forAll(sqTdiff, cellI)
        {
            J += 0.5 * sqTdiff[cellI] * interpolationPlane.cellVol[cellI];
        }
        sqTdiff.clear();
    }
    else
    {
        J = 0.5 * Tdiff.dot(Tdiff);
    }

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



int inverseLaplacianProblem_CG::isInPlane(double cx, double cy, double cz,
        Foam::vector thermocoupleCellDim)
{
    return (cx >= interpolationPlane.minX - thermocoupleCellDim[0] / 4 &&
            cy >= interpolationPlane.Y - thermocoupleCellDim[1] / 4 &&
            cy <= interpolationPlane.Y + thermocoupleCellDim[1] / 4 &&
            cz >= interpolationPlane.minZ - thermocoupleCellDim[2] / 4 &&
            cx <= interpolationPlane.maxX + thermocoupleCellDim[0] / 4 &&
            cz <= interpolationPlane.maxZ + thermocoupleCellDim[2] / 4
           );
}

void inverseLaplacianProblem_CG::writeFields(label folderNumber,
        const char* folder)
{
    fvMesh& mesh = _mesh();
    volScalarField& T = _T();
    volScalarField& lambda = _lambda();
    volScalarField& deltaT = _deltaT();
    autoPtr<volScalarField> gVolField_
    (
        new volScalarField("g", T)
    );
    volScalarField& gVolField = gVolField_();
    ITHACAutilities::assignIF(gVolField, homogeneousBC);
    //Access the mesh information for the boundary
    const polyPatch& cPatch = mesh.boundaryMesh()[hotSide_ind];
    //List of cells close to a boundary
    const labelUList& faceCells = cPatch.faceCells();
    forAll(cPatch, faceI)
    {
        //id of the owner cell having the face
        label faceOwner = faceCells[faceI] ;
        gVolField[faceOwner] = g[faceI];
    }
    ITHACAstream::exportSolution(T, std::to_string(folderNumber + 1), folder,
                                 "T");
    ITHACAstream::exportSolution(lambda, std::to_string(folderNumber + 1), folder,
                                 "lambda");
    ITHACAstream::exportSolution(deltaT, std::to_string(folderNumber + 1), folder,
                                 "deltaT");
    ITHACAstream::exportSolution(gVolField, std::to_string(folderNumber + 1),
                                 folder, "g");
}

//Interpolates the values in Tmeas on the interpolation plane defined in readThermocouples()
void inverseLaplacianProblem_CG::thermocouplesInterpolation()
{
    fvMesh& mesh = _mesh();
    volScalarField& T = _T();
    DataTable thermocouplesSamples;
    DenseVector x(2);
    unsigned i = 0;

    if (!interpolationPlaneDefined)
    {
        interpolationPlane.thermocoupleX.resize(thermocouplesNum);
        interpolationPlane.thermocoupleZ.resize(thermocouplesNum);

        //find the dimensions of the cells containing the thermocouples
        //I am assuming all thermocouples are a the same y coordinate
        //I am assuming all cells have the same dimensions
        if ( Pstream::master() == true )
        {
            unsigned flag = 1;

            while (flag == 1)
            {
                if (thermocouplesCellID [i] != -1)
                {
                    labelList pLabels(mesh.cells()[thermocouplesCellID [i]].labels(mesh.faces()));
                    pointField pLocal(pLabels.size(), Foam::vector::zero);
                    interpolationPlane.thermocoupleCellDim = cellDim (mesh.faces(),
                            mesh.points(),
                            mesh.cells()[thermocouplesCellID [i]],
                            pLabels,
                            pLocal);
                    flag = 0;
                }

                i++;
            }
        }

        //reduce(interpolationPlane.thermocoupleCellDim, sumOp<vector>());
        forAll(thermocouplesCellID, thermocoupleI)
        {
            if (thermocouplesCellID[thermocoupleI] != -1)
            {
                interpolationPlane.thermocoupleX [thermocoupleI] =
                    mesh.C()[thermocouplesCellID [thermocoupleI]].component(0);
                interpolationPlane.thermocoupleZ [thermocoupleI] =
                    mesh.C()[thermocouplesCellID [thermocoupleI]].component(2);
            }
            else
            {
                interpolationPlane.thermocoupleX [thermocoupleI] = 0;
                interpolationPlane.thermocoupleZ [thermocoupleI] = 0;
            }
        }
        //reduce(interpolationPlane.thermocoupleX, sumOp<List<scalar>>());
        //reduce(interpolationPlane.thermocoupleZ, sumOp<List<scalar>>());
    }

    forAll(thermocouplesCellID, thermocoupleI)
    {
        x(0) = interpolationPlane.thermocoupleX[thermocoupleI];
        x(1) = interpolationPlane.thermocoupleZ[thermocoupleI];
        thermocouplesSamples.addSample(x, Tmeas(thermocoupleI));
    }
    std::cout << Tmeas << std::endl;
    RBFSpline rbfspline(thermocouplesSamples, RadialBasisFunctionType::GAUSSIAN);
    auto inPlaneCellID = 0;
    forAll(T.internalField(), cellI)
    {
        auto cx = mesh.C()[cellI].component(Foam::vector::X);
        auto cy = mesh.C()[cellI].component(Foam::vector::Y);
        auto cz = mesh.C()[cellI].component(Foam::vector::Z);

        if (!interpolationPlaneDefined)
        {
            if (isInPlane(cx, cy, cz, interpolationPlane.thermocoupleCellDim))
            {
                auto planeSize = interpolationPlane.cellID.size() + 1;
                interpolationPlane.cellID.resize (planeSize);
                interpolationPlane.Tmeas.resize  (planeSize);
                interpolationPlane.Tdirect.resize(planeSize);
                interpolationPlane.Tdiff.resize  (planeSize);
                interpolationPlane.Tsens.resize  (planeSize);
                interpolationPlane.cellVol.resize(planeSize);
                x(0) = cx;
                x(1) = cz;
                interpolationPlane.cellID [planeSize - 1] = cellI;
                interpolationPlane.Tmeas  [planeSize - 1] =
                    rbfspline.eval(x);
                interpolationPlane.cellVol[planeSize - 1] =
                    mesh.V()[cellI];
            }
        }
        else
        {
            if (isInPlane(cx, cy, cz, interpolationPlane.thermocoupleCellDim))
            {
                x(0) = cx;
                x(1) = cz;
                interpolationPlane.Tmeas  [inPlaneCellID] =
                    rbfspline.eval(x);
                inPlaneCellID++;
            }
        }
    }
    interpolationPlaneDefined = 1;
}

//Interpolates the values in Tmeas on the interpolation plane defined in readThermocouples()
void inverseLaplacianProblem_CG::thermocouplesInterpolation(
    DenseMatrix& RBFweights, DenseMatrix& RBFbasis)
{
    fvMesh& mesh = _mesh();
    volScalarField& T = _T();
    DataTable thermocouplesSamples;
    DenseVector x(2);
    unsigned i = 0;

    if (!interpolationPlaneDefined)
    {
        interpolationPlane.thermocoupleX.resize(thermocouplesNum);
        interpolationPlane.thermocoupleZ.resize(thermocouplesNum);

        //find the dimensions of the cells containing the thermocouples
        //I am assuming all thermocouples are a the same y coordinate
        //I am assuming all cells have the same dimensions
        if ( Pstream::master() == true )
        {
            unsigned flag = 1;

            while (flag == 1)
            {
                if (thermocouplesCellID [i] != -1)
                {
                    labelList pLabels(mesh.cells()[thermocouplesCellID [i]].labels(mesh.faces()));
                    pointField pLocal(pLabels.size(), Foam::vector::zero);
                    interpolationPlane.thermocoupleCellDim = cellDim (mesh.faces(),
                            mesh.points(),
                            mesh.cells()[thermocouplesCellID [i]],
                            pLabels,
                            pLocal);
                    flag = 0;
                }

                i++;
            }
        }

        //reduce(interpolationPlane.thermocoupleCellDim, sumOp<vector>());
        forAll(thermocouplesCellID, thermocoupleI)
        {
            if (thermocouplesCellID[thermocoupleI] != -1)
            {
                interpolationPlane.thermocoupleX [thermocoupleI] =
                    mesh.C()[thermocouplesCellID [thermocoupleI]].component(0);
                interpolationPlane.thermocoupleZ [thermocoupleI] =
                    mesh.C()[thermocouplesCellID [thermocoupleI]].component(2);
            }
            else
            {
                interpolationPlane.thermocoupleX [thermocoupleI] = 0;
                interpolationPlane.thermocoupleZ [thermocoupleI] = 0;
            }
        }
        //reduce(interpolationPlane.thermocoupleX, sumOp<List<scalar>>());
        //reduce(interpolationPlane.thermocoupleZ, sumOp<List<scalar>>());
    }

    forAll(thermocouplesCellID, thermocoupleI)
    {
        x(0) = interpolationPlane.thermocoupleX[thermocoupleI];
        x(1) = interpolationPlane.thermocoupleZ[thermocoupleI];
        thermocouplesSamples.addSample(x, Tmeas(thermocoupleI));
    }
    std::cout << Tmeas << std::endl;
    RBFSpline rbfspline(thermocouplesSamples, RadialBasisFunctionType::GAUSSIAN);
    auto inPlaneCellID = 0;
    forAll(T.internalField(), cellI)
    {
        auto cx = mesh.C()[cellI].component(Foam::vector::X);
        auto cy = mesh.C()[cellI].component(Foam::vector::Y);
        auto cz = mesh.C()[cellI].component(Foam::vector::Z);

        if (!interpolationPlaneDefined)
        {
            if (isInPlane(cx, cy, cz, interpolationPlane.thermocoupleCellDim))
            {
                auto planeSize = interpolationPlane.cellID.size() + 1;
                interpolationPlane.cellID.resize (planeSize);
                interpolationPlane.Tmeas.resize  (planeSize);
                interpolationPlane.Tdirect.resize(planeSize);
                interpolationPlane.Tdiff.resize  (planeSize);
                interpolationPlane.Tsens.resize  (planeSize);
                interpolationPlane.cellVol.resize(planeSize);
                x(0) = cx;
                x(1) = cz;
                interpolationPlane.cellID [planeSize - 1] = cellI;
                interpolationPlane.Tmeas  [planeSize - 1] =
                    rbfspline.eval(x);
                interpolationPlane.cellVol[planeSize - 1] =
                    mesh.V()[cellI];
            }
        }
        else
        {
            if (isInPlane(cx, cy, cz, interpolationPlane.thermocoupleCellDim))
            {
                x(0) = cx;
                x(1) = cz;
                interpolationPlane.Tmeas  [inPlaneCellID] =
                    rbfspline.eval(x);
                inPlaneCellID++;
            }
        }
    }
    interpolationPlaneDefined = 1;
}

void inverseLaplacianProblem_CG::differenceBetweenDirectAndMeasure()
{
    volScalarField& T = _T();

    if (interpolation)
    {
        forAll(interpolationPlane.cellID, cellI)
        {
            interpolationPlane.Tdirect[cellI] =
                T.internalField()[interpolationPlane.cellID [cellI]];
        }
        interpolationPlane.Tdiff = interpolationPlane.Tdirect -
                                   interpolationPlane.Tmeas;
    }
    else
    {
        Tdirect = fieldValueAtThermocouples(T);
        Tdiff = Tdirect - Tmeas;
    }
}

void inverseLaplacianProblem_CG::restart(word fieldName)
{
    Info << "\nResetting time, mesh and fields: " << fieldName << "\n" << endl;
    _simple.clear();

    if (fieldName == "T" || fieldName == "all")
    {
        _T.clear();
    }

    if (fieldName == "lambda" || fieldName == "all")
    {
        _lambda.clear();
    }

    if (fieldName == "deltaT" || fieldName == "all")
    {
        _deltaT.clear();
    }

    argList& args = _args();
    Time& runTime = _runTime();
    //Reinitializing runTime
    instantList Times = runTime.times();
    runTime.setTime(Times[1], 1);
    Foam::fvMesh& mesh = _mesh();
    _simple = autoPtr<simpleControl>
              (
                  new simpleControl
                  (
                      mesh
                  )
              );

    if (fieldName == "T" || fieldName == "all")
    {
        //Info << "ReReading field T\n" << endl;
        _T = autoPtr<volScalarField>
             (
                 new volScalarField
                 (
                     IOobject
                     (
                         "T",
                         runTime.timeName(),
                         mesh,
                         IOobject::MUST_READ,
                         IOobject::AUTO_WRITE
                     ),
                     mesh
                 )
             );
    }

    if (fieldName == "deltaT" || fieldName == "all")
    {
        //Info << "ReReading field deltaT\n" << endl;
        _deltaT = autoPtr<volScalarField>
                  (
                      new volScalarField
                      (
                          IOobject
                          (
                              "deltaT",
                              runTime.timeName(),
                              mesh,
                              IOobject::MUST_READ,
                              IOobject::AUTO_WRITE
                          ),
                          mesh
                      )
                  );
    }

    if (fieldName == "lambda" || fieldName == "all")
    {
        //Info << "ReReading field lambda\n" << endl;
        _lambda = autoPtr<volScalarField>
                  (
                      new volScalarField
                      (
                          IOobject
                          (
                              "lambda",
                              runTime.timeName(),
                              mesh,
                              IOobject::MUST_READ,
                              IOobject::AUTO_WRITE
                          ),
                          mesh
                      )
                  );
    }

    Info << "Ready for new computation" << endl;
}


//void inverseLaplacianProblem_CG::offlineSolve()
//{
//    volScalarField& T = _T();
//    volScalarField& lambda = _lambda();
//    volScalarField& deltaT = _deltaT();
//
//    if (offline)
//    {
//        ITHACAstream::read_fields(Tfield, T, "./ITHACAoutput/Offline/");
//        ITHACAstream::read_fields(lambdaField, lambda, "./ITHACAoutput/Offline/");
//        ITHACAstream::read_fields(deltaTfield, deltaT, "./ITHACAoutput/Offline/");
//    }
//    else
//    {
//        offlinePhase = 1;
//        set_g();
//        readThermocouples();
//        set_valueFraction();
//        saveSolInLists = 0;
//        offlineSolutionI = 1;
//
//        for (label i = 0; i < muSamples.cols(); i++)
//        {
//            Info << "Sample number " << i + 1 << endl;
//            Tmeas = muSamples.col(i);
//
//            if (interpolation)
//            {
//                Info << "Interpolating thermocouples measurements in the " <<
//                     "plane defined by them" << endl;
//                thermocouplesInterpolation();
//            }
//            else
//            {
//                Info << "NOT interpolating thermocouples measurements" << endl;
//            }
//
//            if (!conjugateGradient())
//            {
//                Info << "Conjugate gradient method did not converge" <<
//                     "Exiting" << endl;
//                exit(10);
//            }
//        }
//
//        saveSolInLists = 0;
//        offlinePhase = 0;
//    }
//}

