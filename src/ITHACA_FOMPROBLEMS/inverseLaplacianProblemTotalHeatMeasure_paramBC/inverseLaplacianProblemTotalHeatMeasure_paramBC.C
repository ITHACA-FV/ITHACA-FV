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
/// Source file of the inverseLaplacianProblemTotalHeatMeasure_paramBC class.


#include "inverseLaplacianProblemTotalHeatMeasure_paramBC.H"

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructors
inverseLaplacianProblemTotalHeatMeasure_paramBC::inverseLaplacianProblemTotalHeatMeasure_paramBC() {}

inverseLaplacianProblemTotalHeatMeasure_paramBC::inverseLaplacianProblemTotalHeatMeasure_paramBC(
    int argc, char* argv[])
    :
    inverseLaplacianProblem_paramBC::inverseLaplacianProblem_paramBC(argc, argv)
{
}

// * * * * * * * * * * * * * * Full Order Methods * * * * * * * * * * * * * * //



void inverseLaplacianProblemTotalHeatMeasure_paramBC::parameterizedBCoffline(
    bool force)
{
    fvMesh& mesh = _mesh();
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
        Info << "\nOffline already computed." << endl;
        Info << "Check that the basis used for the parameterized BC are correct (RBF, POD, etc.)\n";
        Theta = ITHACAstream::readMatrix(folderOffline + "Theta_mat.txt");
        Phi = ITHACAstream::readMatrix(folderOffline + "Phi_mat.txt");
        phi = ITHACAstream::readMatrix(folderOffline + "phi_mat.txt");
        addSol = ITHACAstream::readMatrix(folderOffline + "addSol_mat.txt");
        volScalarField& T(_T());
        ITHACAstream::read_fields(Tad_base, "Tad", folderOffline, 0, 1);
        ITHACAstream::read_fields(Tbasis, "T",
                                  folderOffline);
    }
    else
    {
        Info << "\nComputing offline" << endl;
        solveAdditional();
        ITHACAstream::exportMatrix(addSol, "addSol", "eigen", folderOffline);
        M_Assert(Tmeas.size() > 0, "Initialize Tmeas");
        M_Assert(gWeights.size() > 0, "Initialize gWeights");
        Theta.resize(Tmeas.size(), gWeights.size());
        metaData_offline metaData(Tmeas.size(), gWeights.size(), baseFuncType,
                                  shapeParameter);
        std::ofstream fout(folderOffline + "metaData.txt");
        fout << metaData.numberTC << ' ' <<
             metaData.numberBasis << ' ' <<
             metaData.basisType << ' ' <<
             metaData.shapeParameter << ' ';
        fout.close();
        int Nbasis = Theta.cols();
        Phi.resize(gWeights.size(), gWeights.size());
        phi.resize(Nbasis);

        for (label j = 0; j < Nbasis; j++)
        {
            gWeights = Foam::zero();
            gWeights[j] =  1;
            update_gParametrized(gWeights);
            Info << "Solving for j = " << j << endl;
            solveDirect();
            phi(j) = ITHACAutilities::integralOnPatch(mesh, g, "hotSide");
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

        for (label baseI = 0; baseI < Nbasis; baseI++)
        {
            for (label baseJ = 0; baseJ < Nbasis; baseJ++)
            {
                Phi(baseI, baseJ) = phi(baseI) * phi(baseJ);
            }
        }

        ITHACAstream::exportFields(Tbasis, folderOffline, "T");
        ITHACAstream::exportMatrix(Theta, "Theta", "eigen", folderOffline);
        ITHACAstream::exportMatrix(Phi, "Phi", "eigen", folderOffline);
        ITHACAstream::exportMatrix(phi, "phi", "eigen", folderOffline);
        Info << "\nOffline part ENDED\n" << endl;
    }

    Eigen::MatrixXd A = Theta.transpose() * Theta + gIntegralWeight * Phi;
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
}

Eigen::VectorXd
inverseLaplacianProblemTotalHeatMeasure_paramBC::parameterizedBC(
    word linSys_solver,
    double regPar)
{
    Info << endl << "Using quasilinearity of direct problem" << endl;
    List<Eigen::MatrixXd> linSys;
    linSys.resize(2);
    Info << "debug: Theta size = " << Theta.rows() << ", " << Theta.cols() << endl;
    linSys[0] = Theta.transpose() * Theta + gIntegralWeight * Phi;
    linSys[1] = gIntegralWeight * gIntegral_meas * phi + Theta.transpose() *
                (Tmeas + addSol);
    Eigen::VectorXd weigths;
    M_Assert(std::abs(gIntegral_meas) > 1e-16, "First set up gIntegral_meas");
    M_Assert(std::abs(gIntegralWeight) > 1e-16, "First set up gIntegralWeight");

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
        weigths = ITHACAregularization::Tikhonov(linSys[0], linSys[1], regPar);
    }
    else
    {
        Info << "Select a linear system solver in this list:" << endl
             << "fullPivLU, jacobiSvd, householderQr, ldlt, inverse, TSVD" << endl;
        exit(1);
    }

    parameterizedBCpostProcess(weigths);
    return weigths;
}
