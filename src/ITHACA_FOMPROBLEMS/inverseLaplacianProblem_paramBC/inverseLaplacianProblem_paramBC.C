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
/// Source file of the inverseLaplacianProblem_paramBC class.

#include "inverseLaplacianProblem_paramBC.H"

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructors
inverseLaplacianProblem_paramBC::inverseLaplacianProblem_paramBC() {}

inverseLaplacianProblem_paramBC::inverseLaplacianProblem_paramBC(int argc,
        char* argv[])
    :
    inverseLaplacianProblem::inverseLaplacianProblem(argc, argv)
{
}

// * * * * * * * * * * * * * * Full Order Methods * * * * * * * * * * * * * * //

void inverseLaplacianProblem_paramBC::set_gBaseFunctions()
{
    volScalarField& T = _T();
    fvMesh& mesh = _mesh();

    if (baseFuncType == "rbf")
    {
        Info << "Radial Basis Functions are used." << endl;
        // The center of each function is the projection of each thermocouple
        // on the boundary hotSide

        if (!thermocouplesRead)
        {
            readThermocouples();
        }

        gBaseFunctions.resize(thermocouplesNum);
        gWeights.resize(thermocouplesNum);
        forAll(gBaseFunctions, funcI)
        {
            gBaseFunctions[funcI].resize(T.boundaryField()[hotSide_ind].size());
            scalar thermocoupleX = thermocouplesPos[funcI].x();
            scalar thermocoupleZ = thermocouplesPos[funcI].z();
            forAll (T.boundaryField()[hotSide_ind], faceI)
            {
                scalar faceX = mesh.boundaryMesh()[hotSide_ind].faceCentres()[faceI].x();
                scalar faceZ = mesh.boundaryMesh()[hotSide_ind].faceCentres()[faceI].z();
                scalar radius = Foam::sqrt((faceX - thermocoupleX) * (faceX - thermocoupleX) +
                                           (faceZ - thermocoupleZ) * (faceZ - thermocoupleZ));
                gBaseFunctions[funcI][faceI] = Foam::exp(- (shapeParameter *
                                               shapeParameter
                                               * radius * radius));
            }
        }
    }
    else if (baseFuncType == "pod")
    {
        printf("At line number %d in file %s\n", __LINE__, __FILE__);
        Info << "\nPod basis are not coded yet" << endl;
        Info << "Exiting" << endl;
        exit(101);
        //Eigen::MatrixXd temp =
        //    ITHACAstream::readMatrix("./ITHACAoutput/podMarquardt/gReducedBases_mat.txt");
        //gBaseFunctions.resize(temp.cols());
        //gWeights.resize(temp.cols());
        //forAll(gBaseFunctions, baseI)
        //{
        //    gBaseFunctions[baseI].resize(temp.rows());
        //    forAll(gBaseFunctions[baseI], faceI)
        //    {
        //        gBaseFunctions[baseI][faceI] = temp(faceI, baseI);
        //    }
        //}
    }
    else
    {
        printf("At line number %d in file %s\n", __LINE__, __FILE__);
        Info << "\nBasis function type not known. It can be rbf or pod" << endl;
        Info << "Exiting" << endl;
        exit(102);
    }
}

void inverseLaplacianProblem_paramBC::set_gBaseFunctionsPOD(label Nmodes)
{
    volScalarField& T = _T();
    fvMesh& mesh = _mesh();
    Eigen::MatrixXd gBaseFuncEigen;
    set_gParametrized("rbf");

    if (Nmodes == 0)
    {
        Info << "Selecting all modes." << endl;
        Nmodes = gBaseFunctions.size();
    }

    gBaseFuncEigen.resize(gBaseFunctions[0].size(), gBaseFunctions.size());
    Eigen::VectorXd faceAreaVect;
    faceAreaVect.resize(mesh.magSf().boundaryField()[hotSide_ind].size());
    forAll(gBaseFunctions, funcI)
    {
        forAll (T.boundaryField()[hotSide_ind], faceI)
        {
            if (funcI == 0)
            {
                faceAreaVect(faceI) = mesh.magSf().boundaryField()[hotSide_ind][faceI];
            }

            gBaseFuncEigen(faceI, funcI) = gBaseFunctions[funcI][faceI];
        }
    }
    Eigen::MatrixXd correlationMatrix = gBaseFuncEigen.transpose() *
                                        faceAreaVect.asDiagonal() * gBaseFuncEigen;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(correlationMatrix,
                                          Eigen::ComputeThinU | Eigen::ComputeThinV);
    gPODmodes = svd.matrixU().leftCols(Nmodes);
    Eigen::MatrixXd gBaseFuncEigen_new = gBaseFuncEigen * gPODmodes;
    Info << "gBaseFuncEigen_new size = " << gBaseFuncEigen_new.cols() << ", " <<
         gBaseFuncEigen_new.rows() << endl;
    gBaseFunctions.resize(Nmodes);
    gWeights.resize(Nmodes);
    forAll(gBaseFunctions, funcI)
    {
        forAll (T.boundaryField()[hotSide_ind], faceI)
        {
            gBaseFunctions[funcI][faceI] = gBaseFuncEigen_new(faceI, funcI);
        }
    }
}

void inverseLaplacianProblem_paramBC::set_gParametrized(word _baseFuncType,
        scalar _shapeParameter)
{
    shapeParameter = _shapeParameter;
    baseFuncType = _baseFuncType;
    volScalarField& T = _T();
    set_gBaseFunctions();
    g.resize(T.boundaryField()[hotSide_ind].size(), 0.0);
    forAll (gWeights, weigthI)
    {
        gWeights[weigthI] = 0; //-10000;
    }
    Info << "gWeights = " << gWeights << endl;
    forAll (T.boundaryField()[hotSide_ind], faceI)
    {
        g[faceI] = 0.0;
        forAll (gWeights, weigthI)
        {
            g[faceI] += gWeights[weigthI] * gBaseFunctions[weigthI][faceI];
        }
    }
}

void inverseLaplacianProblem_paramBC::update_gParametrized(List<scalar> weigths)
{
    M_Assert(weigths.size() == gBaseFunctions.size(),
             "weigths size different from basis functions size");
    volScalarField& T = _T();
    forAll (T.boundaryField()[hotSide_ind], faceI)
    {
        g[faceI] = 0.0;
        forAll (weigths, weigthI)
        {
            g[faceI] += weigths[weigthI] * gBaseFunctions[weigthI][faceI];
        }
    }
}

void inverseLaplacianProblem_paramBC::parameterizedBCoffline(bool force)
{
    fvMesh& mesh = _mesh();
    volScalarField& T(_T());
    Tbasis.resize(0);
    Tad_base.resize(0);
    char recomputeOffline;

    if (ITHACAutilities::check_file(folderOffline + "Theta_mat.txt") && force == 0)
    {
        do
        {
            metaData_offline metaData;
            std::ifstream fin(folderOffline + "metaData.txt");
            fin >> metaData.numberTC >> metaData.numberBasis >>
                metaData.basisType >> metaData.shapeParameter;
            fin.close();
            std::cout << "\nOffline FOUND with parameter:\n" <<
                      "Number of thermocouples = " << metaData.numberTC <<
                      "\nNumber of basis functions = " << metaData.numberBasis <<
                      "\nType of basis functions = " << metaData.basisType <<
                      "\nRBF shape parameters = " << metaData.shapeParameter <<
                      "\n\nShould I recompute it? [y/n]" << std::endl;
            std::cin >> recomputeOffline;
        }
        while ( !cin.fail() && recomputeOffline != 'y' && recomputeOffline != 'n' );
    }

    if (recomputeOffline == 'y')
    {
        force = 1;
    }

    if (ITHACAutilities::check_file(folderOffline + "Theta_mat.txt") && force == 0)
    {
        volScalarField Tad(_T());
        Tad.rename("Tad");
        Info << "\nOffline already computed." << endl;
        Info << "Check that the basis used for the parameterized BC are correct (RBF, POD, etc.)\n";
        Theta = ITHACAstream::readMatrix(folderOffline + "Theta_mat.txt");
        addSol = ITHACAstream::readMatrix(folderOffline + "addSol_mat.txt");
        ITHACAstream::read_fields(Tad_base, Tad, folderOffline, 0, 1);
        ITHACAstream::read_fields(Tbasis, T,
                                  folderOffline);
    }
    else
    {
        Info << "\nComputing offline" << endl;

        if (!ITHACAutilities::check_folder(folderOffline))
        {
            mkDir(folderOffline);
        }

        metaData_offline metaData(Tmeas.size(), gWeights.size(), baseFuncType,
                                  shapeParameter);
        std::ofstream fout(folderOffline + "metaData.txt");
        fout << metaData.numberTC << ' ' <<
             metaData.numberBasis << ' ' <<
             metaData.basisType << ' ' <<
             metaData.shapeParameter << ' ';
        fout.close();
        solveAdditional();
        ITHACAstream::exportMatrix(addSol, "addSol", "eigen", folderOffline);
        M_Assert(Tmeas.size() > 0, "Initialize Tmeas");
        M_Assert(gWeights.size() > 0, "Initialize gWeights");
        Theta.resize(Tmeas.size(), gWeights.size());

        for (label j = 0; j < Theta.cols(); j++)
        {
            gWeights = Foam::zero();
            gWeights[j] =  1;
            update_gParametrized(gWeights);
            Info << "Solving for j = " << j << endl;
            solveDirect();
            volScalarField& T = _T();
            Tbasis.append(T.clone());
            Tdirect = fieldValueAtThermocouples(T);

            for (label i = 0; i < Theta.rows(); i++)
            {
                Theta(i, j) = Tdirect(i) + addSol(i);
            }

            volScalarField gParametrizedField = list2Field(g);
            ITHACAstream::exportSolution(gParametrizedField, std::to_string(j + 1),
                                         folderOffline,
                                         "gParametrized");
        }

        ITHACAstream::exportFields(Tbasis, folderOffline, "T");
        ITHACAstream::exportMatrix(Theta, "Theta", "eigen", folderOffline);
        Info << "\nOffline part ENDED\n" << endl;
    }

    Eigen::MatrixXd A = Theta.transpose() * Theta;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A,
                                          Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::MatrixXd singularValues = svd.singularValues();

    if (singularValues.minCoeff() > 0)
    {
        double conditionNumber = singularValues.maxCoeff() / singularValues.minCoeff();
        Info << "Condition number = " << conditionNumber << endl;
    }

    ITHACAstream::exportMatrix(singularValues, "singularValues", "eigen",
                               folderOffline);
    offlineReady = 1;
}

Eigen::VectorXd inverseLaplacianProblem_paramBC::parameterizedBC(
    word linSys_solver,
    double regPar)
{
    M_Assert(offlineReady, "Compute offline before running parameterizedBC.");
    Info << endl << "Using quasilinearity of direct problem" << endl;
    //parameterizedBCoffline(folder, forceOffline);
    List<Eigen::MatrixXd> linSys;
    linSys.resize(2);
    Info << "Theta size = " << Theta.rows() << ", " << Theta.cols() << endl;
    linSys[0] = Theta.transpose() * Theta;
    linSys[1] = Theta.transpose() * (Tmeas + addSol);
    Eigen::VectorXd weigths;

    if (linSys_solver == "fullPivLU")
    {
        weigths = linSys[0].fullPivLu().solve(linSys[1]);
    }
    else if (linSys_solver == "jacobiSvd")
    {
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(linSys[0],
                                              Eigen::ComputeThinU | Eigen::ComputeThinV);
        weigths = svd.solve(linSys[1]);
    }
    else if (linSys_solver == "householderQr")
    {
        weigths = linSys[0].householderQr().solve(linSys[1]);
    }
    else if (linSys_solver == "ldlt")
    {
        weigths = linSys[0].ldlt().solve(linSys[1]);
    }
    else if (linSys_solver == "inverse")
    {
        weigths = linSys[0].inverse() * linSys[1];
    }
    else if (linSys_solver == "TSVD")
    {
        weigths = ITHACAregularization::TSVD(linSys[0], linSys[1], int(regPar));
    }
    else if (linSys_solver == "Tikhonov")
    {
        Info << "WARNING: Tikhonov regularization does not work properly" << endl;
        weigths = ITHACAregularization::Tikhonov(linSys[0], linSys[1], regPar);
    }
    else
    {
        printf("At line number %d in file %s\n", __LINE__, __FILE__);
        Info << "Select a linear system solver in this list:" << endl
             << "fullPivLU, jacobiSvd, householderQr, ldlt, inverse, TSVD" << endl;
        Info << "Exiting." << endl;
        exit(1);
    }

    parameterizedBCpostProcess(weigths);
    return weigths;
}

void inverseLaplacianProblem_paramBC::parameterizedBCpostProcess(
    Eigen::VectorXd weigths)
{
    //// Printing outputs at screen
    //std::cout << "Eigenvalues of Theta.transpose() * Theta are " << std::endl;
    //std::cout << linSys[0].eigenvalues() << std::endl;
    //Eigen::JacobiSVD<Eigen::MatrixXd> svd(linSys[0], Eigen::ComputeThinU | Eigen::ComputeThinV);
    //std::cout << "Singular values of Theta.transpose() * Theta are " << std::endl;
    //std::cout << svd.singularValues() << std::endl;
    //std::cout << "debug: weigths size = " << std::endl;
    //std::cout << weigths.size() << std::endl;
    //
    //residual =  linSys[0] * weigths - linSys[1];
    //std::cout << "Residual  = " << std::endl;
    //std::cout << residual << std::endl;
    //std::cout << "Residual 2-norm = " << std::endl;
    //std::cout << residual.squaredNorm() << std::endl;
    gWeights.resize(weigths.size());
    forAll(gWeights, weightI)
    {
        gWeights[weightI] = weigths(weightI);
    }
    update_gParametrized(gWeights);
    reconstructT();
    volScalarField& T = _T();
    Tdirect = fieldValueAtThermocouples(T);
    J = 0.5 * (Tdirect - Tmeas).dot(Tdirect - Tmeas);
    Info << "J = " << J << endl;
    Info << "End" << endl;
    Info << endl;
}

void inverseLaplacianProblem_paramBC::solveAdditional()
{
    restart();
    Tad_base.resize(0);
    fvMesh& mesh = _mesh();
    simpleControl& simple = _simple();
    volScalarField Tad(_T);
    Foam::Time& runTime = _runTime();
    set_valueFraction();
    List<scalar> RobinBC = - Tf;
    forAll(mesh.boundaryMesh(), patchI)
    {
        if (patchI == mesh.boundaryMesh().findPatchID("coldSide"))
        {
            ITHACAutilities::assignMixedBC(Tad, patchI, RobinBC, refGrad, valueFraction);
        }
        else if (patchI == mesh.boundaryMesh().findPatchID("hotSide"))
        {
            ITHACAutilities::assignBC(Tad, patchI, homogeneousBC);
        }
        else
        {
            ITHACAutilities::assignBC(Tad, patchI, homogeneousBC);
        }
    }
#if OFVER == 6

    while (simple.loop(runTime))
#else
    while (simple.loop())
#endif
    {
        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix TEqn
            (
                fvm::laplacian(DT, Tad)
            );
            TEqn.solve();
        }
    }

    addSol = fieldValueAtThermocouples(Tad);
    Tad_base.append(Tad.clone());
    ITHACAstream::exportSolution(Tad, "1",
                                 folderOffline,
                                 "Tad");
}

void inverseLaplacianProblem_paramBC::reconstructT()
{
    Info << "Reconstructing field T" << endl;
    restart();
    update_gParametrized(gWeights);
    volScalarField Trec = Tbasis[0];
    ITHACAutilities::assignIF(Trec, homogeneousBC);
    forAll(Tbasis, baseI)
    {
        Trec += gWeights[baseI] * (Tbasis[baseI] + Tad_base[0]);
    }
    Trec += - Tad_base[0];
    volScalarField& T = _T();
    T = Trec;
}
