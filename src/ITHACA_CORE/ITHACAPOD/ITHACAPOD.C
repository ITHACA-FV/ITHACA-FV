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
/// source file for the ITHACAPOD class

#include "ITHACAPOD.H"
#include "EigenFunctions.H"

namespace ITHACAPOD
{

template<class Type, template<class> class PatchField, class GeoMesh>
void getNestedSnapshotMatrix(
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& snapshots,
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& ModesGlobal,
    word fieldName,
    label Npar, label NnestedOut)
{
    ITHACAparameters* para(ITHACAparameters::getInstance());
    List<PtrList<GeometricField<Type, PatchField, GeoMesh>>>
    SnapMatrixNested;
    label Nt = snapshots.size() / Npar;
    SnapMatrixNested.setSize(Nt);

    for (label i = 0; i < Npar; i++)
    {
        SnapMatrixNested[i].resize(Nt);

        for (label j = 0; j < Nt; j++)
        {
            SnapMatrixNested[i].set(j, snapshots[j + Nt * i].clone());
        }
    }

    List<PtrList<GeometricField<Type, PatchField, GeoMesh>>> UModesNested;
    PtrList<GeometricField<Type, PatchField, GeoMesh>>  y;
    UModesNested.setSize(Nt);

    for (label i = 0; i < Npar; i++)
    {
        getWeightedModes(SnapMatrixNested[i], UModesNested[i], 0, 0, 0,
                         NnestedOut);
    }

    for (label i = 0; i < Npar; i++)
    {
        y = UModesNested[i];

        for (label j = 0; j < y.size(); j++)
        {
            ModesGlobal.append(y[j].clone());
        }
    }
}

template void getNestedSnapshotMatrix(PtrList<volScalarField>&
                                      snapshots, PtrList<volScalarField>& ModesGlobal, word fieldName, label Npar,
                                      label NnestedOut);
template void getNestedSnapshotMatrix(PtrList<volVectorField>&
                                      snapshots, PtrList<volVectorField>& ModesGlobal, word fieldName, label Npar,
                                      label NnestedOut);

template<class Type, template<class> class PatchField, class GeoMesh>
void getModes(
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& snapshots,
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& modes,
    word fieldName, bool podex, bool supex, bool sup, label nmodes)
{
    ITHACAparameters* para(ITHACAparameters::getInstance());
    word PODkey = "POD_" + fieldName;
    word PODnorm = para->ITHACAdict->lookupOrDefault<word>(PODkey, "L2");
    M_Assert(PODnorm == "L2" ||
             PODnorm == "Frobenius", "The PODnorm can be only L2 or Frobenius");
    Info << "Performing POD for " << fieldName << " using the " << PODnorm <<
         " norm" << endl;

    if ((podex == 0 && sup == 0) || (supex == 0 && sup == 1))
    {
        if (para->eigensolver == "spectra" )
        {
            if (nmodes == 0)
            {
                nmodes = snapshots.size() - 2;
            }

            M_Assert(nmodes <= snapshots.size() - 2,
                     "The number of requested modes cannot be bigger than the number of Snapshots - 2");
        }
        else
        {
            if (nmodes == 0)
            {
                nmodes = snapshots.size();
            }

            M_Assert(nmodes <= snapshots.size(),
                     "The number of requested modes cannot be bigger than the number of Snapshots");
        }

        Eigen::MatrixXd SnapMatrix = Foam2Eigen::PtrList2Eigen(snapshots);
        List<Eigen::MatrixXd> SnapMatrixBC = Foam2Eigen::PtrList2EigenBC(snapshots);
        label NBC = snapshots[0].boundaryField().size();
        Eigen::MatrixXd _corMatrix;

        if (PODnorm == "L2")
        {
            _corMatrix = ITHACAutilities::getMassMatrix(snapshots);
        }
        else if (PODnorm == "Frobenius")
        {
            _corMatrix = ITHACAutilities::getMassMatrix(snapshots, 0, false);
        }

        if (Pstream::parRun())
        {
            List<double> vec(_corMatrix.data(), _corMatrix.data() + _corMatrix.size());
            reduce(vec, sumOp<List<double>>());
            std::memcpy(_corMatrix.data(), &vec[0], sizeof (double)*vec.size());
        }

        Eigen::VectorXd eigenValueseig;
        Eigen::MatrixXd eigenVectoreig;
        modes.resize(nmodes);
        Info << "####### Performing the POD using EigenDecomposition " <<
             fieldName << " #######" << endl;
        label ncv = snapshots.size();
        Spectra::DenseSymMatProd<double> op(_corMatrix);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esEg;

        if (para->eigensolver == "spectra")
        {
            Spectra::SymEigsSolver<double, Spectra::LARGEST_ALGE, Spectra::DenseSymMatProd<double>>
                    es(&op, nmodes, ncv);
            std::cout << "Using Spectra EigenSolver " << std::endl;
            es.init();
            es.compute(1000, 1e-10, Spectra::LARGEST_ALGE);
            M_Assert(es.info() == Spectra::SUCCESSFUL,
                     "The Eigenvalue Decomposition did not succeed");
            eigenVectoreig = es.eigenvectors().real();
            eigenValueseig = es.eigenvalues().real();
        }
        else if (para->eigensolver == "eigen")
        {
            std::cout << "Using Eigen EigenSolver " << std::endl;
            esEg.compute(_corMatrix);
            M_Assert(esEg.info() == Eigen::Success,
                     "The Eigenvalue Decomposition did not succeed");
            eigenVectoreig = esEg.eigenvectors().real().rowwise().reverse().leftCols(
                                 nmodes);
            eigenValueseig = esEg.eigenvalues().real().array().reverse();
        }

        if (eigenValueseig.array().minCoeff() < 0)
        {
            eigenValueseig = eigenValueseig.array() + 2 * abs(
                                 eigenValueseig.array().minCoeff());
        }

        Info << "####### End of the POD for " << snapshots[0].name() << " #######" <<
             endl;
        //Eigen::VectorXd eigenValueseigLam =
        //    eigenValueseig.real().array().abs().cwiseInverse().sqrt() ;
        //Eigen::MatrixXd modesEig = (SnapMatrix * eigenVectoreig) *
        //                           eigenValueseigLam.head(nmodes).asDiagonal();
        Eigen::MatrixXd modesEig = (SnapMatrix * eigenVectoreig);
        // Computing Normalization factors of the POD Modes
        Eigen::VectorXd V = ITHACAutilities::getMassMatrixFV(snapshots[0]);
        Eigen::VectorXd normFact(nmodes);

        for (label i = 0; i < nmodes; i++)
        {
            if (PODnorm == "L2")
            {
                normFact(i) = std::sqrt((modesEig.col(i).transpose() * V.asDiagonal() *
                                         modesEig.col(i))(0, 0));

                if (Pstream::parRun())
                {
                    normFact(i) = (modesEig.col(i).transpose() * V.asDiagonal() *
                                   modesEig.col(i))(0, 0);
                }
            }
            else if (PODnorm == "Frobenius")
            {
                normFact(i) = std::sqrt((modesEig.col(i).transpose() * modesEig.col(i))(0, 0));

                if (Pstream::parRun())
                {
                    normFact(i) = (modesEig.col(i).transpose() * modesEig.col(i))(0, 0);
                }
            }
        }

        if (Pstream::parRun())
        {
            List<double> vec(normFact.data(), normFact.data() + normFact.size());
            reduce(vec, sumOp<List<double>>());
            std::memcpy(normFact.data(), &vec[0], sizeof (double)*vec.size());
            normFact = normFact.cwiseSqrt();
        }

        List<Eigen::MatrixXd> modesEigBC;
        modesEigBC.resize(NBC);

        for (label i = 0; i < NBC; i++)
        {
            //modesEigBC[i] = (SnapMatrixBC[i] * eigenVectoreig) *
            //                eigenValueseigLam.head(nmodes).asDiagonal();
            modesEigBC[i] = (SnapMatrixBC[i] * eigenVectoreig);
        }

        for (label i = 0; i < nmodes; i++)
        {
            modesEig.col(i) = modesEig.col(i).array() / normFact(i);

            for (label j = 0; j < NBC; j++)
            {
                modesEigBC[j].col(i) = modesEigBC[j].col(i).array() / normFact(i);
            }
        }

        for (label i = 0; i < modes.size(); i++)
        {
            GeometricField<Type, PatchField, GeoMesh>  tmp2(snapshots[0].name(),
                    snapshots[0]);
            Eigen::VectorXd vec = modesEig.col(i);
            tmp2 = Foam2Eigen::Eigen2field(tmp2, vec);

            for (label k = 0; k < NBC; k++)
            {
                ITHACAutilities::assignBC(tmp2, k, modesEigBC[k].col(i));
            }

            modes.set(i, tmp2.clone());
        }

        eigenValueseig = eigenValueseig / eigenValueseig.sum();
        Eigen::VectorXd cumEigenValues(eigenValueseig);

        for (label j = 1; j < cumEigenValues.size(); ++j)
        {
            cumEigenValues(j) += cumEigenValues(j - 1);
        }

        Info << "####### Saving the POD bases for " << snapshots[0].name() <<
             " #######" << endl;

        if (sup)
        {
            ITHACAstream::exportFields(modes, "./ITHACAoutput/supremizer/",
                                       snapshots[0].name());
        }
        else
        {
            ITHACAstream::exportFields(modes, "./ITHACAoutput/POD/", snapshots[0].name());
        }

        Eigen::saveMarketVector(eigenValueseig,
                                "./ITHACAoutput/POD/Eigenvalues_" + snapshots[0].name(), para->precision,
                                para->outytpe);
        Eigen::saveMarketVector(cumEigenValues,
                                "./ITHACAoutput/POD/CumEigenvalues_" + snapshots[0].name(), para->precision,
                                para->outytpe);
    }
    else
    {
        Info << "Reading the existing modes" << endl;

        if (sup == 1)
        {
            ITHACAstream::read_fields(modes, fieldName + "sup",
                                      "./ITHACAoutput/supremizer/");
        }
        else
        {
            ITHACAstream::read_fields(modes, fieldName, "./ITHACAoutput/POD/");
        }
    }
}

template void getModes(PtrList<volVectorField>& snapshots,
                       PtrList<volVectorField>& modes, word fieldName, bool podex, bool supex,
                       bool sup, label nmodes);
template void getModes(PtrList<volScalarField>& snapshots,
                       PtrList<volScalarField>& modes, word fieldName, bool podex, bool supex,
                       bool sup, label nmodes);
template void getModes(PtrList<surfaceScalarField>& snapshots,
                       PtrList<surfaceScalarField>& modes, word fieldName, bool podex, bool supex,
                       bool sup, label nmodes);


template<class Type, template<class> class PatchField, class GeoMesh>
void getWeightedModes(
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& snapshots,
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& modes,
    word fieldName, bool podex,
    bool supex, bool sup, label nmodes)
{
    ITHACAparameters* para(ITHACAparameters::getInstance());

    if (nmodes == 0)
    {
        nmodes = snapshots.size() - 2;
    }

    M_Assert(nmodes <= snapshots.size() - 2,
             "The number of requested modes cannot be bigger than the number of Snapshots - 2");

    if ((podex == 0 && sup == 0) || (supex == 0 && sup == 1))
    {
        Eigen::MatrixXd SnapMatrix = Foam2Eigen::PtrList2Eigen(snapshots);
        List<Eigen::MatrixXd> SnapMatrixBC = Foam2Eigen::PtrList2EigenBC(snapshots);
        label NBC = snapshots[0].boundaryField().size();
        Eigen::MatrixXd _corMatrix = ITHACAutilities::getMassMatrix(snapshots);
        Eigen::VectorXd eigenValueseig;
        Eigen::MatrixXd eigenVectoreig;
        modes.resize(nmodes);
        Info << "####### Performing the POD using EigenDecomposition " <<
             snapshots[0].name() << " #######" << endl;
        label ncv = snapshots.size();
        Spectra::DenseSymMatProd<double> op(_corMatrix);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esEg;
        Spectra::SymEigsSolver<double, Spectra::LARGEST_ALGE, Spectra::DenseSymMatProd<double>>
                es(&op, nmodes, ncv);

        if (para->eigensolver == "spectra")
        {
            std::cout << "Using Spectra EigenSolver " << std::endl;
            es.init();
            es.compute(1000, 1e-10, Spectra::LARGEST_ALGE);
            M_Assert(es.info() == Spectra::SUCCESSFUL,
                     "The Eigenvalue Decomposition did not succeed");
            eigenVectoreig = es.eigenvectors().real();
            eigenValueseig = es.eigenvalues().real();
        }
        else if (para->eigensolver == "eigen")
        {
            std::cout << "Using Eigen EigenSolver " << std::endl;
            esEg.compute(_corMatrix);
            M_Assert(esEg.info() == Eigen::Success,
                     "The Eigenvalue Decomposition did not succeed");
            eigenVectoreig = esEg.eigenvectors().real().rowwise().reverse().leftCols(
                                 nmodes);
            eigenValueseig = esEg.eigenvalues().real().reverse();
        }

        Info << "####### End of the POD for " << snapshots[0].name() << " #######" <<
             endl;
        Eigen::VectorXd eigenValueseigLam =
            eigenValueseig.real().array().cwiseInverse().sqrt() ;
        Eigen::VectorXd eigenValueseigWeigted = eigenValueseig.head(
                nmodes).real().array() ;
        Eigen::MatrixXd modesEig = (SnapMatrix * eigenVectoreig) *
                                   eigenValueseigLam.head(nmodes).asDiagonal() *
                                   eigenValueseigWeigted.asDiagonal();
        List<Eigen::MatrixXd> modesEigBC;
        modesEigBC.resize(NBC);

        for (label i = 0; i < NBC; i++)
        {
            modesEigBC[i] = (SnapMatrixBC[i] * eigenVectoreig) *
                            eigenValueseigLam.asDiagonal() * eigenValueseigWeigted.asDiagonal();
        }

        for (label i = 0; i < modes.size(); i++)
        {
            GeometricField<Type, PatchField, GeoMesh> tmp2(snapshots[0].name(),
                    snapshots[0]);
            Eigen::VectorXd vec = modesEig.col(i);
            tmp2 = Foam2Eigen::Eigen2field(tmp2, vec);

            for (label k = 0; k < tmp2.boundaryField().size(); k++)
            {
                ITHACAutilities::assignBC(tmp2, k, modesEigBC[k].col(i));
            }

            modes.set(i, tmp2.clone());
        }

        eigenValueseig = eigenValueseig / eigenValueseig.sum();
        Eigen::VectorXd cumEigenValues(eigenValueseig);

        for (label j = 1; j < cumEigenValues.size(); ++j)
        {
            cumEigenValues(j) += cumEigenValues(j - 1);
        }

        Info << "####### Saving the POD bases for " << snapshots[0].name() <<
             " #######" << endl;

        //exportBases(modes, snapshots, sup);
        if (sup)
        {
            ITHACAstream::exportFields(modes, "./ITHACAoutput/supremizer/",
                                       snapshots[0].name());
        }
        else
        {
            ITHACAstream::exportFields(modes, "./ITHACAoutput/POD/", snapshots[0].name());
        }

        Eigen::saveMarketVector(eigenValueseig,
                                "./ITHACAoutput/POD/Eigenvalues_" + snapshots[0].name(), para->precision,
                                para->outytpe);
        Eigen::saveMarketVector(cumEigenValues,
                                "./ITHACAoutput/POD/CumEigenvalues_" + snapshots[0].name(), para->precision,
                                para->outytpe);
    }
    else
    {
        Info << "Reading the existing modes" << endl;

        if (sup == 1)
        {
            ITHACAstream::read_fields(modes, fieldName + "sup",
                                      "./ITHACAoutput/supremizer/");
        }
        else
        {
            ITHACAstream::read_fields(modes, fieldName, "./ITHACAoutput/POD/");
        }
    }
}

template void getWeightedModes(PtrList<volScalarField>& snapshots,
                               PtrList<volScalarField>& modes, word fieldName, bool podex, bool supex,
                               bool sup, label nmodes);

template void getWeightedModes(PtrList<volVectorField>& snapshots,
                               PtrList<volVectorField>& modes, word fieldName, bool podex, bool supex,
                               bool sup, label nmodes);

template<class Type, template<class> class PatchField, class GeoMesh>
void getModesSVD(
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& snapshots,
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& modes,
    word fieldName, bool podex,
    bool supex, bool sup, label nmodes)
{
    ITHACAparameters* para(ITHACAparameters::getInstance());

    if ((podex == 0 && sup == 0) || (supex == 0 && sup == 1))
    {
        PtrList<volVectorField> Bases;
        modes.resize(nmodes);
        Info << "####### Performing POD using Singular Value Decomposition for " <<
             snapshots[0].name() << " #######" << endl;
        Eigen::MatrixXd SnapMatrix = Foam2Eigen::PtrList2Eigen(snapshots);
        Eigen::VectorXd V = ITHACAutilities::getMassMatrixFV(snapshots[0]);
        Eigen::VectorXd V3dSqrt = V.array().sqrt();
        Eigen::VectorXd V3dInv = V3dSqrt.array().cwiseInverse();
        auto VMsqr = V3dSqrt.asDiagonal();
        auto VMsqrInv = V3dInv.asDiagonal();
        Eigen::MatrixXd SnapMatrix2 = VMsqr * SnapMatrix;
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(SnapMatrix2,
                                              Eigen::ComputeThinU | Eigen::ComputeThinV);
        Info << "####### End of the POD for " << snapshots[0].name() << " #######" <<
             endl;
        Eigen::VectorXd eigenValueseig;
        Eigen::MatrixXd eigenVectoreig;
        eigenValueseig = svd.singularValues().real();
        eigenVectoreig = svd.matrixU().real();
        Eigen::MatrixXd modesEig = VMsqrInv * eigenVectoreig;
        GeometricField<Type, PatchField, GeoMesh> tmb_bu(snapshots[0].name(),
                snapshots[0] * 0);

        for (label i = 0; i < nmodes; i++)
        {
            Eigen::VectorXd vec = modesEig.col(i);
            tmb_bu = Foam2Eigen::Eigen2field(tmb_bu, vec);
            modes.set(i, tmb_bu.clone());
        }

        eigenValueseig = eigenValueseig / eigenValueseig.sum();
        Eigen::VectorXd cumEigenValues(eigenValueseig);

        for (label j = 1; j < cumEigenValues.size(); ++j)
        {
            cumEigenValues(j) += cumEigenValues(j - 1);
        }

        Info << "####### Saving the POD bases for " << snapshots[0].name() <<
             " #######" << endl;

        if (sup)
        {
            ITHACAstream::exportFields(modes, "./ITHACAoutput/supremizer/",
                                       snapshots[0].name());
        }
        else
        {
            ITHACAstream::exportFields(modes, "./ITHACAoutput/POD/", snapshots[0].name());
        }

        //exportBases(modes, snapshots, sup);
        Eigen::saveMarketVector(eigenValueseig,
                                "./ITHACAoutput/POD/Eigenvalues_" + snapshots[0].name(), para->precision,
                                para->outytpe);
        Eigen::saveMarketVector(cumEigenValues,
                                "./ITHACAoutput/POD/CumEigenvalues_" + snapshots[0].name(), para->precision,
                                para->outytpe);
    }
    else
    {
        Info << "Reading the existing modes" << endl;

        if (sup == 1)
        {
            ITHACAstream::read_fields(modes, fieldName + "sup",
                                      "./ITHACAoutput/supremizer/");
        }
        else
        {
            ITHACAstream::read_fields(modes, fieldName, "./ITHACAoutput/POD/");
        }
    }
}

template void getModesSVD(PtrList<volScalarField>& snapshots,
                          PtrList<volScalarField>& modes, word fieldName, bool podex, bool supex,
                          bool sup, label nmodes);

template void getModesSVD(PtrList<volVectorField>& snapshots,
                          PtrList<volVectorField>& modes, word fieldName, bool podex, bool supex,
                          bool sup, label nmodes);


/// Construct the Correlation Matrix for Scalar Field
template<>
Eigen::MatrixXd corMatrix(PtrList<volScalarField>& snapshots)
{
    Info << "########## Filling the correlation matrix for " << snapshots[0].name()
         << "##########" << endl;
    Eigen::MatrixXd matrix( snapshots.size(), snapshots.size());

    for (label i = 0; i < snapshots.size(); i++)
    {
        for (label j = 0; j <= i; j++)
        {
            matrix(i, j) = fvc::domainIntegrate(snapshots[i] * snapshots[j]).value();
        }
    }

    for (label i = 1; i < snapshots.size(); i++)
    {
        for (label j = 0; j < i; j++)
        {
            matrix(j, i) = matrix(i, j);
        }
    }

    return matrix;
}


/// Construct the Correlation Matrix for Vector Field
template<>
Eigen::MatrixXd corMatrix(PtrList<volVectorField>& snapshots)
{
    Info << "########## Filling the correlation matrix for " << snapshots[0].name()
         << "##########" << endl;
    Eigen::MatrixXd matrix( snapshots.size(), snapshots.size());

    for (label i = 0; i < snapshots.size(); i++)
    {
        for (label j = 0; j <= i; j++)
        {
            matrix(i, j) = fvc::domainIntegrate(snapshots[i] & snapshots[j]).value();
        }
    }

    for (label i = 1; i < snapshots.size(); i++)
    {
        for (label j = 0; j < i; j++)
        {
            matrix(j, i) = matrix(i, j);
        }
    }

    return matrix;
}

/// Construct the Correlation Matrix for Vector Field
template<>
Eigen::MatrixXd corMatrix(List<Eigen::SparseMatrix<double>>&
                          snapshots)
{
    Info << "########## Filling the correlation matrix for the matrix list ##########"
         << endl;
    Eigen::MatrixXd matrix( snapshots.size(), snapshots.size());

    for (label i = 0; i < snapshots.size(); i++)
    {
        for (label j = 0; j <= i; j++)
        {
            double res = 0;

            for (label k = 0; k < snapshots[i].cols(); k++)
            {
                res += snapshots[i].col(k).dot(snapshots[j].col(k));
            }

            matrix(i, j) = res;
        }
    }

    for (label i = 1; i < snapshots.size(); i++)
    {
        for (label j = 0; j < i; j++)
        {
            matrix(j, i) = matrix(i, j);
        }
    }

    return matrix;
}

/// Construct the Correlation Matrix for Vector Field
template<>
Eigen::MatrixXd corMatrix(List<Eigen::VectorXd>& snapshots)
{
    Info << "########## Filling the correlation matrix for the matrix list ##########"
         << endl;
    Eigen::MatrixXd matrix( snapshots.size(), snapshots.size());

    for (label i = 0; i < snapshots.size(); i++)
    {
        for (label j = 0; j <= i; j++)
        {
            matrix(i, j) = (snapshots[i].transpose() * snapshots[j]).trace();
        }
    }

    for (label i = 1; i < snapshots.size(); i++)
    {
        for (label j = 0; j < i; j++)
        {
            matrix(j, i) = matrix(i, j);
        }
    }

    return matrix;
}



/// Export the Bases
template<class Type, template<class> class PatchField, class GeoMesh>
void exportBases(PtrList<GeometricField<Type, PatchField, GeoMesh>>& s,
                 PtrList<GeometricField<Type, PatchField, GeoMesh>>& bases,
                 word fieldName, bool sup)
{
    if (sup)
    {
        fileName fieldname;

        for (label i = 0; i < s.size(); i++)
        {
            mkDir("./ITHACAoutput/supremizer/" + name(i + 1));
            fieldname = "./ITHACAoutput/supremizer/" + name(i + 1) + "/" +
                        bases[i].name();
            OFstream os(fieldname);
            bases[i].writeHeader(os);
            os << s[i] << endl;
        }
    }
    else
    {
        fileName fieldname;

        for (label i = 0; i < s.size(); i++)
        {
            mkDir("./ITHACAoutput/POD/" + name(i + 1));
            fieldname = "./ITHACAoutput/POD/" + name(i + 1) + "/" + bases[i].name();
            OFstream os(fieldname);
            bases[i].writeHeader(os);
            os << s[i] << endl;
        }
    }
}
template void exportBases(PtrList<volVectorField>& s,
                          PtrList<volVectorField>& bases, word fieldName, bool sup);
template void exportBases(PtrList<volScalarField>& s,
                          PtrList<volScalarField>& bases, word fieldName, bool sup);

void exportEigenvalues(scalarField Eigenvalues, fileName name,
                       bool sup)
{
    if (sup)
    {
        fileName fieldname;
        mkDir("./ITHACAoutput/supremizer/");
        fieldname = "./ITHACAoutput/supremizer/Eigenvalues_" + name;
        OFstream os(fieldname);

        for (label i = 0; i < Eigenvalues.size(); i++)
        {
            os << Eigenvalues[i] << endl;
        }
    }
    else
    {
        fileName fieldname;
        mkDir("./ITHACAoutput/POD/");
        fieldname = "./ITHACAoutput/POD/Eigenvalues_" + name;
        OFstream os(fieldname);

        for (label i = 0; i < Eigenvalues.size(); i++)
        {
            os << Eigenvalues[i] << endl;
        }
    }
}

void exportcumEigenvalues(scalarField cumEigenvalues, fileName name,
                          bool sup)
{
    if (sup)
    {
        fileName fieldname;
        mkDir("./ITHACAoutput/supremizer/");
        fieldname = "./ITHACAoutput/supremizer/cumEigenvalues_" + name;
        OFstream os(fieldname);

        for (label i = 0; i < cumEigenvalues.size(); i++)
        {
            os << cumEigenvalues[i] << endl;
        }
    }
    else
    {
        fileName fieldname;
        mkDir("./ITHACAoutput/POD/");
        fieldname = "./ITHACAoutput/POD/cumEigenvalues_" + name;
        OFstream os(fieldname);

        for (label i = 0; i < cumEigenvalues.size(); i++)
        {
            os << cumEigenvalues[i] << endl;
        }
    }
}


std::tuple<List<Eigen::SparseMatrix<double>>, List<Eigen::VectorXd>>
        DEIMmodes(List<Eigen::SparseMatrix<double>>& A,
                  List<Eigen::VectorXd>& b, label nmodesA, label nmodesB, word MatrixName)
{
    ITHACAparameters* para(ITHACAparameters::getInstance());
    List<Eigen::SparseMatrix<double>> ModesA(nmodesA);
    List<Eigen::VectorXd> ModesB(nmodesB);

    if (!ITHACAutilities::check_folder("./ITHACAoutput/DEIM/" + MatrixName))
    {
        if (nmodesA > A.size() - 2 || nmodesB > A.size() - 2 )
        {
            std::cout <<
                      "The number of requested modes cannot be bigger than the number of Snapshots - 2"
                      << std::endl;
            exit(0);
        }

        scalarField eigenValuesA(nmodesA);
        scalarField cumEigenValuesA(nmodesA);
        scalarField eigenValuesB(nmodesB);
        scalarField cumEigenValuesB(nmodesB);
        List<scalarField> eigenVectorA(nmodesA);
        List<scalarField> eigenVectorB(nmodesB);

        for (label i = 0; i < nmodesA; i++)
        {
            eigenVectorA[i].setSize(A.size());
        }

        for (label i = 0; i < nmodesB; i++)
        {
            eigenVectorB[i].setSize(A.size());
        }

        Eigen::MatrixXd corMatrixA = corMatrix(A);
        Eigen::MatrixXd corMatrixB = corMatrix(b);
        Info << "####### Performing the POD for the Matrix List #######" << endl;
        Spectra::DenseSymMatProd<double> opA(corMatrixA);
        Spectra::DenseSymMatProd<double> opB(corMatrixB);
        Spectra::SymEigsSolver<double, Spectra::LARGEST_ALGE, Spectra:: DenseSymMatProd<double>>
                esA(&opA, nmodesA, nmodesA + 10);
        Spectra::SymEigsSolver<double, Spectra::LARGEST_ALGE, Spectra:: DenseSymMatProd<double>>
                esB(&opB, nmodesB, nmodesB + 10);
        esA.init();
        esB.init();
        esA.compute();

        if (nmodesB != 1)
        {
            esB.compute();
        }

        Info << "####### End of the POD for the Matrix List #######" << endl;
        Eigen::VectorXd eigenValueseigA;
        Eigen::MatrixXd eigenVectorseigA;
        Eigen::VectorXd eigenValueseigB;
        Eigen::MatrixXd eigenVectorseigB;

        if (esA.info() == Spectra::SUCCESSFUL)
        {
            eigenValueseigA = esA.eigenvalues().real();
            eigenVectorseigA = esA.eigenvectors().real();

            if (esB.info() == Spectra::SUCCESSFUL && nmodesB != 1)
            {
                eigenValueseigB = esB.eigenvalues().real();
                eigenVectorseigB = esB.eigenvectors().real();
            }
            else if (nmodesB == 1)
            {
                eigenValueseigB.resize(1);
                eigenVectorseigB.resize(A.size(), nmodesB);
                eigenValueseigB(0) = 1;
                eigenVectorseigB = eigenVectorseigB * 0;
                eigenVectorseigB(0, 0) = 1;
            }
            else
            {
                Info << "The Eigenvalue solver in ITHACAPOD.H did not converge, exiting the code"
                     << endl;
                exit(0);
            }
        }
        else
        {
            Info << "The Eigenvalue solver in ITHACAPOD.H did not converge, exiting the code"
                 << endl;
            exit(0);
        }

        for (label i = 0; i < nmodesA; i++)
        {
            eigenValuesA[i] = eigenValueseigA(i) / eigenValueseigA.sum();
        }

        for (label i = 0; i < nmodesB; i++)
        {
            eigenValuesB[i] = eigenValueseigB(i) / eigenValueseigB.sum();
        }

        for (label i = 0; i < nmodesA; i++)
        {
            for (label k = 0; k < A.size(); k++)
            {
                eigenVectorA[i][k] = eigenVectorseigA(k, i);
            }
        }

        for (label i = 0; i < nmodesB; i++)
        {
            for (label k = 0; k < A.size(); k++)
            {
                eigenVectorB[i][k] = eigenVectorseigB(k, i);
            }
        }

        cumEigenValuesA[0] = eigenValuesA[0];
        cumEigenValuesB[0] = eigenValuesB[0];

        for (label i = 1; i < nmodesA; i++)
        {
            cumEigenValuesA[i] = cumEigenValuesA[i - 1] + eigenValuesA[i];
        }

        for (label i = 1; i < nmodesB; i++)
        {
            cumEigenValuesB[i] = cumEigenValuesB[i - 1] + eigenValuesB[i];
        }

        Eigen::SparseMatrix<double> tmp_A;
        Eigen::VectorXd tmp_B;

        for (label i = 0; i < nmodesA; i++)
        {
            tmp_A = eigenVectorA[i][0] * A[0];

            for (label k = 1; k < A.size(); k++)
            {
                tmp_A += eigenVectorA[i][k] * A[k];
            }

            ModesA[i] = tmp_A;
        }

        for (label i = 0; i < nmodesB; i++)
        {
            tmp_B = eigenVectorB[i][0] * b[0];

            for (label k = 1; k < A.size(); k++)
            {
                tmp_B += eigenVectorB[i][k] * b[k];
            }

            ModesB[i] = tmp_B;
        }

        ITHACAstream::exportList(eigenValuesA,
                                 "./ITHACAoutput/DEIM/" + MatrixName + "/", "eigenValuesA_" + MatrixName);
        ITHACAstream::exportList(cumEigenValuesA,
                                 "./ITHACAoutput/DEIM/" + MatrixName + "/", "cumEigenValuesA_" + MatrixName);
        ITHACAstream::exportList(eigenValuesB,
                                 "./ITHACAoutput/DEIM/" + MatrixName + "/", "eigenValuesB_" + MatrixName);
        ITHACAstream::exportList(cumEigenValuesB,
                                 "./ITHACAoutput/DEIM/" + MatrixName + "/", "cumEigenValuesB_" + MatrixName);

        for (label i = 0; i < ModesA.size(); i++)
        {
            ITHACAstream::SaveSparseMatrix(ModesA[i],
                                           "./ITHACAoutput/DEIM/" + MatrixName + "/", "A_" + MatrixName + name(i));
        }

        for (label i = 0; i < ModesB.size(); i++)
        {
            ITHACAstream::SaveDenseMatrix(ModesB[i],
                                          "./ITHACAoutput/DEIM/" + MatrixName + "/", "B_" + MatrixName + name(i));
        }
    }
    else
    {
        for (label i = 0; i < nmodesA; i++)
        {
            ITHACAstream::ReadSparseMatrix(ModesA[i],
                                           "./ITHACAoutput/DEIM/" + MatrixName + "/", "A_" + MatrixName + name(i));
        }

        for (label i = 0; i < nmodesB; i++)
        {
            ITHACAstream::ReadDenseMatrix(ModesB[i],
                                          "./ITHACAoutput/DEIM/" + MatrixName + "/", "B_" + MatrixName + name(i));
        }
    }

    std::tuple <List<Eigen::SparseMatrix<double>>, List<Eigen::VectorXd>> tupla;
    tupla = std::make_tuple(ModesA, ModesB);
    return tupla;
}

void GrammSchmidt(Eigen::MatrixXd& Matrix)
{
    Eigen::MatrixXd Ortho = Matrix;
    Ortho = Matrix;

    for (label i = 0; i <  Matrix.cols(); i++)
    {
        for (label k = 0; k < i; k++)
        {
            double num = Ortho.col(k).transpose() * Matrix.col(i);
            double den = (Ortho.col(k).transpose() * Ortho.col(k));
            double fact = num / den;
            Ortho.col(i) -= fact * Ortho.col(k) ;
        }

        Ortho.col(i).normalize();
    }

    Matrix = Ortho;
}
template<class Type, template<class> class PatchField, class GeoMesh>
void getModes(PtrList<GeometricField<Type, PatchField, GeoMesh>>& snapshots,
              PtrList<GeometricField<Type, PatchField, GeoMesh>>& modes,
              PtrList<volScalarField>& Volumes,
              word fieldName, bool podex,
              bool supex, bool sup, label nmodes)
{
    ITHACAparameters* para(ITHACAparameters::getInstance());

    if (nmodes == 0 && para->eigensolver == "spectra")
    {
        nmodes = snapshots.size() - 2;
    }

    if (nmodes == 0 && para->eigensolver == "eigen")
    {
        nmodes = snapshots.size();
    }

    if (para->eigensolver == "spectra")
    {
        M_Assert(nmodes <= snapshots.size() - 2,
                 "The number of requested modes cannot be bigger than the number of Snapshots - 2");
    }

    if ((podex == 0 && sup == 0) || (supex == 0 && sup == 1))
    {
        Eigen::MatrixXd SnapMatrix = Foam2Eigen::PtrList2Eigen(snapshots);
        List<Eigen::MatrixXd> SnapMatrixBC = Foam2Eigen::PtrList2EigenBC(snapshots);
        label NBC = snapshots[0].boundaryField().size();
        Eigen::MatrixXd V = Foam2Eigen::PtrList2Eigen(Volumes);
        Eigen::MatrixXd _corMatrix(snapshots.size(), snapshots.size());
        Info << "Filling the correlation matrix for field " << snapshots[0].name() <<
             endl;

        for (label i = 0; i < snapshots.size(); i++)
        {
            for (label j = 0; j <= i; j++)
            {
                Eigen::VectorXd Mij = (V.col(i).array() * V.col(j).array());
                Mij = Mij.array().abs().sqrt();
                _corMatrix(i, j) = SnapMatrix.col(i).transpose() * Mij.asDiagonal() *
                                   SnapMatrix.col(j);
            }
        }

        std::cout << std::endl;

        for (label i = 1; i < snapshots.size(); i++)
        {
            for (label j = 0; j < i; j++)
            {
                _corMatrix(j, i) = _corMatrix(i, j);
            }
        }

        Eigen::VectorXd eigenValueseig;
        Eigen::MatrixXd eigenVectoreig;
        modes.resize(nmodes);
        Info << "####### Performing the POD using EigenDecomposition for " <<
             snapshots[0].name() << " #######" << endl;
        label ncv = snapshots.size();
        Spectra::DenseSymMatProd<double> op(_corMatrix);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esEg;

        if (para->eigensolver == "spectra")
        {
            Spectra::SymEigsSolver<double, Spectra::LARGEST_ALGE, Spectra::DenseSymMatProd<double>>
                    es(&op, nmodes, ncv);
            std::cout << "Using Spectra EigenSolver " << std::endl;
            es.init();
            es.compute(1000, 1e-10, Spectra::LARGEST_ALGE);
            M_Assert(es.info() == Spectra::SUCCESSFUL,
                     "The Eigenvalue Decomposition did not succeed");
            eigenVectoreig = es.eigenvectors().real();
            eigenValueseig = es.eigenvalues().real();
        }
        else if (para->eigensolver == "eigen")
        {
            std::cout << "Using Eigen EigenSolver " << std::endl;
            esEg.compute(_corMatrix);
            M_Assert(esEg.info() == Eigen::Success,
                     "The Eigenvalue Decomposition did not succeed");
            eigenVectoreig = esEg.eigenvectors().real().rowwise().reverse().leftCols(
                                 nmodes);
            eigenValueseig = esEg.eigenvalues().real().reverse().head(nmodes);
        }

        Info << "####### End of the POD for " << snapshots[0].name() << " #######" <<
             endl;
        Eigen::VectorXd eigenValueseigLam =
            eigenValueseig.real().array().cwiseInverse().sqrt() ;
        Eigen::MatrixXd modesEig = (SnapMatrix * eigenVectoreig) *
                                   eigenValueseigLam.asDiagonal();
        List<Eigen::MatrixXd> modesEigBC;
        modesEigBC.resize(NBC);

        for (label i = 0; i < NBC; i++)
        {
            modesEigBC[i] = (SnapMatrixBC[i] * eigenVectoreig) *
                            eigenValueseigLam.asDiagonal();
        }

        for (label i = 0; i < modes.size(); i++)
        {
            GeometricField<Type, PatchField, GeoMesh> tmp2(snapshots[0].name(),
                    snapshots[0] * 0);
            Eigen::VectorXd vec = modesEig.col(i);
            tmp2 = Foam2Eigen::Eigen2field(tmp2, vec);

            // Adjusting boundary conditions
            for (label k = 0; k < tmp2.boundaryField().size(); k++)
            {
                ITHACAutilities::assignBC(tmp2, k, modesEigBC[k].col(i));
            }

            modes.set(i, tmp2.clone());
        }

        eigenValueseig = eigenValueseig / eigenValueseig.sum();
        Eigen::VectorXd cumEigenValues(eigenValueseig);

        for (label j = 1; j < cumEigenValues.size(); ++j)
        {
            cumEigenValues(j) += cumEigenValues(j - 1);
        }

        Info << "####### Saving the POD bases for " << snapshots[0].name() << " #######"
             << endl;
        exportBases(modes, snapshots, fieldName, sup);
        Eigen::saveMarketVector(eigenValueseig,
                                "./ITHACAoutput/POD/Eigenvalues_" + snapshots[0].name(), para->precision,
                                para->outytpe);
        Eigen::saveMarketVector(cumEigenValues,
                                "./ITHACAoutput/POD/CumEigenvalues_" + snapshots[0].name(), para->precision,
                                para->outytpe);
    }
    else
    {
        Info << "Reading the existing modes" << endl;

        if (sup == 1)
        {
            ITHACAstream::read_fields(modes, "Usup", "./ITHACAoutput/supremizer/");
        }
        else
        {
            ITHACAstream::read_fields(modes, fieldName, "./ITHACAoutput/POD/");
        }
    }
}
template void getModes(PtrList<volScalarField>& snapshots,
                       PtrList<volScalarField>& modes, PtrList<volScalarField>& Volumes,
                       word fieldName, bool podex,
                       bool supex, bool sup, label nmodes);
template void getModes(PtrList<volVectorField>& snapshots,
                       PtrList<volVectorField>& modes, PtrList<volScalarField>& Volumes,
                       word fieldName, bool podex,
                       bool supex, bool sup, label nmodes);


template<typename type_matrix>
std::tuple<List<Eigen::SparseMatrix<double>>, List<Eigen::VectorXd>>
        DEIMmodes(PtrList<type_matrix>& MatrixList, label nmodesA, label nmodesB,
                  word MatrixName)
{
    ITHACAparameters* para(ITHACAparameters::getInstance());
    List<Eigen::SparseMatrix<double>> ModesA(nmodesA);
    List<Eigen::VectorXd> ModesB(nmodesB);

    if (!ITHACAutilities::check_folder("./ITHACAoutput/DEIM/" + MatrixName))
    {
        M_Assert(nmodesA <= MatrixList.size() - 2
                 && nmodesB <= MatrixList.size() - 2,
                 "The number of requested modes cannot be bigger than the number of Snapshots - 2");
        std::tuple<List<Eigen::SparseMatrix<double>>, List<Eigen::VectorXd>> snapshots =
                    Foam2Eigen::LFvMatrix2LSM(MatrixList);
        Eigen::MatrixXd corMatrixA = corMatrix(std::get<0>(snapshots));
        Eigen::MatrixXd corMatrixB = corMatrix(std::get<1>(snapshots));
        Eigen::VectorXd eigenValueseigA;
        Eigen::MatrixXd eigenVectorseigA;
        Eigen::VectorXd eigenValueseigB;
        Eigen::MatrixXd eigenVectorseigB;

        if (para->eigensolver == "spectra")
        {
            Info << "####### Performing the POD decomposition for the Matrix List using Spectra #######"
                 << endl;
            Spectra::DenseSymMatProd<double> opA(corMatrixA);
            Spectra::DenseSymMatProd<double> opB(corMatrixB);
            label ncvA = MatrixList.size();
            label ncvB = MatrixList.size();
            Spectra::SymEigsSolver<double, Spectra::LARGEST_ALGE, Spectra:: DenseSymMatProd<double>>
                    esA(&opA, nmodesA, ncvA);
            Spectra::SymEigsSolver<double, Spectra::LARGEST_ALGE, Spectra:: DenseSymMatProd<double>>
                    esB(&opB, nmodesB, ncvB);
            esA.init();
            esB.init();
            esA.compute(1000, 1e-10, Spectra::LARGEST_ALGE);
            M_Assert(esA.info() == Spectra::SUCCESSFUL,
                     "The Eigenvalue Decomposition did not succeed");

            if (nmodesB != 1)
            {
                esB.compute(1000, 1e-10, Spectra::LARGEST_ALGE);
                M_Assert(esB.info() == Spectra::SUCCESSFUL,
                         "The Eigenvalue Decomposition did not succeed");
            }

            Info << "####### End of the POD decomposition for the Matrix List #######" <<
                 endl;
            eigenValueseigA = esA.eigenvalues().real();
            eigenVectorseigA = esA.eigenvectors().real();

            if (nmodesB != 1)
            {
                eigenValueseigB = esB.eigenvalues().real();
                eigenVectorseigB = esB.eigenvectors().real();
            }
            else
            {
                eigenValueseigB.resize(1);
                eigenVectorseigB.resize(MatrixList.size(), nmodesB);
                eigenValueseigB(0) = 1;
                eigenVectorseigB = eigenVectorseigB * 0;
                eigenVectorseigB(0, 0) = 1;
            }
        }
        else if (para->eigensolver == "eigen")
        {
            Info << "####### Performing the POD decomposition for the Matrix List using Eigen #######"
                 << endl;
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esEgA;
            esEgA.compute(corMatrixA);
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esEgB;
            esEgB.compute(corMatrixB);
            M_Assert(esEgA.info() == Eigen::Success,
                     "The Eigenvalue Decomposition did not succeed");
            M_Assert(esEgB.info() == Eigen::Success,
                     "The Eigenvalue Decomposition did not succeed");
            eigenVectorseigA = esEgA.eigenvectors().real().rowwise().reverse().leftCols(
                                   nmodesA);
            eigenValueseigA = esEgA.eigenvalues().real().reverse().head(nmodesA);
            eigenVectorseigB = esEgB.eigenvectors().real().rowwise().reverse().leftCols(
                                   nmodesB);
            eigenValueseigB = esEgB.eigenvalues().real().reverse().head(nmodesB);
        }

        Eigen::SparseMatrix<double> tmp_A;
        Eigen::VectorXd tmp_B;

        for (label i = 0; i < nmodesA; i++)
        {
            tmp_A = eigenVectorseigA(0, i) * std::get<0>(snapshots)[0];

            for (label k = 1; k < MatrixList.size(); k++)
            {
                tmp_A += eigenVectorseigA(k, i) * std::get<0>(snapshots)[k];
            }

            ModesA[i] = tmp_A;
        }

        for (label i = 0; i < nmodesB; i++)
        {
            tmp_B = eigenVectorseigB(0, i) * std::get<1>(snapshots)[0];

            for (label k = 1; k < MatrixList.size(); k++)
            {
                tmp_B += eigenVectorseigB(k, i) * std::get<1>(snapshots)[k];
            }

            ModesB[i] = tmp_B;
        }

        eigenValueseigA = eigenValueseigA / eigenValueseigA.sum();
        eigenValueseigB = eigenValueseigB / eigenValueseigB.sum();
        Eigen::VectorXd cumEigenValuesA(eigenValueseigA);
        Eigen::VectorXd cumEigenValuesB(eigenValueseigB);

        for (label j = 1; j < cumEigenValuesA.size(); ++j)
        {
            cumEigenValuesA(j) += cumEigenValuesA(j - 1);
        }

        for (label j = 1; j < cumEigenValuesB.size(); ++j)
        {
            cumEigenValuesB(j) += cumEigenValuesB(j - 1);
        }

        for (label i = 0; i < ModesA.size(); i++)
        {
            ITHACAstream::SaveSparseMatrix(ModesA[i],
                                           "./ITHACAoutput/DEIM/" + MatrixName + "/", "A_" + MatrixName + name(i));
        }

        for (label i = 0; i < ModesB.size(); i++)
        {
            ITHACAstream::SaveDenseMatrix(ModesB[i],
                                          "./ITHACAoutput/DEIM/" + MatrixName + "/", "B_" + MatrixName + name(i));
        }

        Eigen::saveMarketVector(eigenValueseigA,
                                "./ITHACAoutput/DEIM/" + MatrixName + "/eigenValuesA", para->precision,
                                para->outytpe);
        Eigen::saveMarketVector(eigenValueseigB,
                                "./ITHACAoutput/DEIM/" + MatrixName + "/eigenValuesB", para->precision,
                                para->outytpe);
        Eigen::saveMarketVector(cumEigenValuesA,
                                "./ITHACAoutput/DEIM/" + MatrixName + "/cumEigenValuesA", para->precision,
                                para->outytpe);
        Eigen::saveMarketVector(cumEigenValuesB,
                                "./ITHACAoutput/DEIM/" + MatrixName + "/cumEigenValuesB", para->precision,
                                para->outytpe);
    }
    else
    {
        for (label i = 0; i < nmodesA; i++)
        {
            ITHACAstream::ReadSparseMatrix(ModesA[i],
                                           "./ITHACAoutput/DEIM/" + MatrixName + "/", "A_" + MatrixName + name(i));
        }

        for (label i = 0; i < nmodesB; i++)
        {
            ITHACAstream::ReadDenseMatrix(ModesB[i],
                                          "./ITHACAoutput/DEIM/" + MatrixName + "/", "B_" + MatrixName + name(i));
        }
    }

    std::tuple <List<Eigen::SparseMatrix<double>>, List<Eigen::VectorXd>> tupla;
    tupla = std::make_tuple(ModesA, ModesB);
    return tupla;
}

template std::tuple<List<Eigen::SparseMatrix<double>>, List<Eigen::VectorXd>>
DEIMmodes(PtrList<fvScalarMatrix>& MatrixList, label nmodesA,
          label nmodesB,
          word MatrixName);

template std::tuple<List<Eigen::SparseMatrix<double>>, List<Eigen::VectorXd>>
DEIMmodes(PtrList<fvVectorMatrix>& MatrixList, label nmodesA,
          label nmodesB,
          word MatrixName);

template<typename Type>
PtrList<GeometricField<Type, fvPatchField, volMesh>> DEIMmodes(
            PtrList<GeometricField<Type, fvPatchField, volMesh>>& SnapShotsMatrix,
            label nmodes,
            word FunctionName, word FieldName)
{
    ITHACAparameters* para(ITHACAparameters::getInstance());
    PtrList<GeometricField<Type, fvPatchField, volMesh>> modes;

    if (!ITHACAutilities::check_folder("./ITHACAoutput/DEIM/" + FunctionName))
    {
        if (nmodes == 0)
        {
            nmodes = SnapShotsMatrix.size() - 2;
        }

        M_Assert(nmodes <= SnapShotsMatrix.size() - 2,
                 "The number of requested modes cannot be bigger than the number of Snapshots - 2");

        if (nmodes > SnapShotsMatrix.size() - 2)
        {
            std::cout <<
                      "The number of requested modes cannot be bigger than the number of Snapshots - 2"
                      << std::endl;
            exit(0);
        }

        PtrList<GeometricField<Type, fvPatchField, volMesh>> Bases;
        Bases.resize(nmodes);
        modes.resize(nmodes);
        Eigen::MatrixXd _corMatrix;
        scalarField eigenValues(nmodes);
        scalarField cumEigenValues(nmodes);
        List<scalarField> eigenVector(nmodes);

        for (label i = 0; i < nmodes; i++)
        {
            eigenVector[i].setSize(SnapShotsMatrix.size());
        }

        _corMatrix = corMatrix(SnapShotsMatrix);
        Info << "####### Performing the POD decomposition for " <<
             SnapShotsMatrix[0].name() << " #######" << endl;
        label ncv = SnapShotsMatrix.size();
        Spectra::DenseGenMatProd<double> op(_corMatrix);
        Spectra::GenEigsSolver<double, Spectra::LARGEST_REAL, Spectra::DenseGenMatProd<double>>
                es(&op, nmodes, ncv);
        es.init();
        es.compute(1000, 1e-10, Spectra::LARGEST_REAL);
        Info << "####### End of the POD decomposition for " << SnapShotsMatrix[0].name()
             << " #######" << endl;
        Eigen::VectorXd eigenValueseig;
        Eigen::MatrixXd eigenVectoreig;

        if (es.info() == Spectra::SUCCESSFUL)
        {
            eigenValueseig = es.eigenvalues().real();
            eigenVectoreig = es.eigenvectors().real();
        }

        cumEigenValues[0] = eigenValueseig(0) / eigenValueseig.sum();
        eigenValues[0] = eigenValueseig(0) / eigenValueseig.sum();

        for (label i = 1; i < nmodes; i++)
        {
            cumEigenValues[i] = cumEigenValues[i - 1] + eigenValueseig(
                                    i) / eigenValueseig.sum();
            eigenValues[i] = eigenValueseig(i) / eigenValueseig.sum();
        }

        for (label i = 0; i < nmodes; i++)
        {
            for (label k = 0; k < SnapShotsMatrix.size(); k++)
            {
                eigenVector[i][k] = eigenVectoreig(k, i);
            }
        }

        GeometricField<Type, fvPatchField, volMesh> tmb_bu = SnapShotsMatrix[0];
        tmb_bu.rename(SnapShotsMatrix[0].name());

        for (label i = 0; i < nmodes; i++)
        {
            tmb_bu = eigenVector[i][0] * SnapShotsMatrix[0];

            for (label k = 1; k < SnapShotsMatrix.size(); k++)
            {
                tmb_bu += eigenVector[i][k] * SnapShotsMatrix[k];
            }

            Bases.set(i, tmb_bu.clone());
            Info << "creating the bases " << i << " for " << SnapShotsMatrix[0].name() <<
                 endl;
        }

        ITHACAutilities::normalizeFields(Bases);

        for (label i = 0; i < Bases.size(); i++)
        {
            auto tmp2(Bases[i]);
            tmp2.rename(SnapShotsMatrix[0].name());
            modes.set(i, tmp2.clone());
        }

        Info << "####### Saving the POD bases for " << SnapShotsMatrix[0].name() <<
             " #######" << endl;
        ITHACAutilities::createSymLink("./ITHACAoutput/DEIM");

        for (label i = 0; i < modes.size(); i++)
        {
            ITHACAstream::exportSolution(modes[i], name(i + 1), "./ITHACAoutput/DEIM",
                                         modes[i].name());
        }

        ITHACAstream::exportList(eigenValues, "./ITHACAoutput/DEIM/",
                                 "eigenValues_" + SnapShotsMatrix[0].name());
        ITHACAstream::exportList(cumEigenValues, "./ITHACAoutput/DEIM/",
                                 "cumEigenValues_" + SnapShotsMatrix[0].name());
    }
    else
    {
        Info << "Reading the existing modes" << endl;
        ITHACAstream::read_fields(modes, FieldName, "./ITHACAoutput/DEIM/");
    }

    return modes;
}

template<class Field_type, class Field_type_2>
void getModes(PtrList<Field_type>& snapshots,
              PtrList<Field_type>& modes, PtrList<Field_type_2>& fields2, word fieldName,
              bool podex,
              bool supex, bool sup, label nmodes)
{
    ITHACAparameters* para(ITHACAparameters::getInstance());

    if ((podex == 0 && sup == 0) || (supex == 0 && sup == 1))
    {
        if (para->eigensolver == "spectra" )
        {
            if (nmodes == 0)
            {
                nmodes = snapshots.size() - 2;
            }

            M_Assert(nmodes <= snapshots.size() - 2,
                     "The number of requested modes cannot be bigger than the number of Snapshots - 2");
        }
        else
        {
            if (nmodes == 0)
            {
                nmodes = snapshots.size();
            }

            M_Assert(nmodes <= snapshots.size(),
                     "The number of requested modes cannot be bigger than the number of Snapshots");
        }

        Eigen::MatrixXd SnapMatrix = Foam2Eigen::PtrList2Eigen(fields2);
        List<Eigen::MatrixXd> SnapMatrixBC = Foam2Eigen::PtrList2EigenBC(fields2);
        label NBC = fields2[0].boundaryField().size();
        auto VM = ITHACAutilities::getMassMatrixFV(fields2[0]);
        Eigen::MatrixXd _corMatrix = SnapMatrix.transpose() * VM.asDiagonal() *
                                     SnapMatrix;

        if (Pstream::parRun())
        {
            List<double> vec(_corMatrix.data(), _corMatrix.data() + _corMatrix.size());
            reduce(vec, sumOp<List<double>>());
            std::memcpy(_corMatrix.data(), &vec[0], sizeof (double)*vec.size());
        }

        Eigen::VectorXd eigenValueseig;
        Eigen::MatrixXd eigenVectoreig;
        modes.resize(nmodes);
        Info << "####### Performing the POD using EigenDecomposition " <<
             fields2[0].name() << " #######" << endl;
        label ncv = snapshots.size();
        Spectra::DenseSymMatProd<double> op(_corMatrix);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esEg;

        if (para->eigensolver == "spectra")
        {
            Spectra::SymEigsSolver<double, Spectra::LARGEST_ALGE, Spectra::DenseSymMatProd<double>>
                    es(&op, nmodes, ncv);
            std::cout << "Using Spectra EigenSolver " << std::endl;
            es.init();
            es.compute(1000, 1e-10, Spectra::LARGEST_ALGE);
            M_Assert(es.info() == Spectra::SUCCESSFUL,
                     "The Eigenvalue Decomposition did not succeed");
            eigenVectoreig = es.eigenvectors().real();
            eigenValueseig = es.eigenvalues().real();
        }
        else if (para->eigensolver == "eigen")
        {
            std::cout << "Using Eigen EigenSolver " << std::endl;
            esEg.compute(_corMatrix);
            M_Assert(esEg.info() == Eigen::Success,
                     "The Eigenvalue Decomposition did not succeed");
            eigenVectoreig = esEg.eigenvectors().real().rowwise().reverse().leftCols(
                                 nmodes);
            eigenValueseig = esEg.eigenvalues().real().reverse().head(nmodes);
        }

        Info << "####### End of the POD for " << snapshots[0].name() << " #######" <<
             endl;
        Eigen::VectorXd eigenValueseigLam =
            eigenValueseig.real().array().cwiseInverse().abs().sqrt() ;
        Eigen::MatrixXd modesEig = (SnapMatrix * eigenVectoreig) *
                                   eigenValueseigLam.asDiagonal();
        List<Eigen::MatrixXd> modesEigBC;
        modesEigBC.resize(NBC);

        for (label i = 0; i < modes.size(); i++)
        {
            Field_type tmp2(snapshots[0].name(), snapshots[0] * 0);

            for (label j = 0; j < snapshots.size(); j++)
            {
                tmp2 += eigenValueseigLam(i) * eigenVectoreig.col(i)(j) * snapshots[j];
            }

            modes.set(i, tmp2.clone());
        }

        eigenValueseig = eigenValueseig / eigenValueseig.sum();
        Eigen::VectorXd cumEigenValues(eigenValueseig);

        for (label j = 1; j < cumEigenValues.size(); ++j)
        {
            cumEigenValues(j) += cumEigenValues(j - 1);
        }

        Info << "####### Saving the POD bases for " << snapshots[0].name() <<
             " #######" << endl;

        if (sup)
        {
            ITHACAstream::exportFields(modes, "./ITHACAoutput/supremizer/",
                                       snapshots[0].name());
        }
        else
        {
            ITHACAstream::exportFields(modes, "./ITHACAoutput/POD/", snapshots[0].name());
        }

        Eigen::saveMarketVector(eigenValueseig,
                                "./ITHACAoutput/POD/Eigenvalues_" + snapshots[0].name(), para->precision,
                                para->outytpe);
        Eigen::saveMarketVector(cumEigenValues,
                                "./ITHACAoutput/POD/CumEigenvalues_" + snapshots[0].name(), para->precision,
                                para->outytpe);
    }
    else
    {
        Info << "Reading the existing modes" << endl;

        if (sup == 1)
        {
            ITHACAstream::read_fields(modes, fieldName, "./ITHACAoutput/supremizer/");
        }
        else
        {
            ITHACAstream::read_fields(modes, fieldName, "./ITHACAoutput/POD/");
        }
    }
}

template void getModes(PtrList<surfaceScalarField>& snapshots,
                       PtrList<surfaceScalarField>& modes,
                       PtrList<volVectorField>& fields2, word fieldName, bool podex, bool supex,
                       bool sup,
                       label nmodes);

template void getModes(PtrList<volScalarField>& snapshots,
                       PtrList<volScalarField>& modes,
                       PtrList<volVectorField>& fields2, word fieldName, bool podex, bool supex,
                       bool sup,
                       label nmodes);

template PtrList<GeometricField<scalar, fvPatchField, volMesh>>
DEIMmodes(
    PtrList<GeometricField<scalar, fvPatchField, volMesh>>& SnapShotsMatrix,
    label nmodes,
    word FunctionName, word FieldName);

template PtrList<GeometricField<vector, fvPatchField, volMesh>>
DEIMmodes(
    PtrList<GeometricField<vector, fvPatchField, volMesh>>& SnapShotsMatrix,
    label nmodes,
    word FunctionName, word FieldName);

}
