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


template<>
void ITHACAPOD::getNestedSnapshotMatrix(PtrList<volVectorField>& snapshotsU,
                                        PtrList<volVectorField>& UModesGlobal, int Npar, int NnestedOut)
{
    ITHACAparameters para;
    List<PtrList<volVectorField>> SnapMatrixNested;
    int Nt = snapshotsU.size() / Npar;
    SnapMatrixNested.setSize(Nt);

    for (int i = 0; i < Npar; i++)
    {
        SnapMatrixNested[i].resize(Nt);

        for (int j = 0; j < Nt; j++)
        {
            SnapMatrixNested[i].set(j, snapshotsU[j + Nt * i]);
        }
    }

    List<PtrList<volVectorField>> UModesNested;
    PtrList<volVectorField>  y;
    UModesNested.setSize(Nt);

    for (int i = 0; i < Npar; i++)
    {
        ITHACAPOD::getWeightedModes(SnapMatrixNested[i], UModesNested[i], 0, 0, 0,
                                    NnestedOut);
    }

    for (int i = 0; i < Npar; i++)
    {
        y = UModesNested[i];

        for (int j = 0; j < y.size(); j++)
        {
            UModesGlobal.append(y[j]);
        }
    }
}


template<>
void ITHACAPOD::getNestedSnapshotMatrix(PtrList<volScalarField>& snapshotsP,
                                        PtrList<volScalarField>& PModesGlobal, int Npar, int NnestedOut)
{
    ITHACAparameters para;
    List<PtrList<volScalarField>> SnapMatrixNested;
    int Nt = snapshotsP.size() / Npar;
    SnapMatrixNested.setSize(Nt);

    for (int i = 0; i < Npar; i++)
    {
        SnapMatrixNested[i].resize(Nt);

        for (int j = 0; j < Nt; j++)
        {
            SnapMatrixNested[i].set(j, snapshotsP[j + Nt * i]);
        }
    }

    List<PtrList<volScalarField>> PModesNested;
    PtrList<volScalarField>  y;
    PModesNested.setSize(Nt);

    for (int i = 0; i < Npar; i++)
    {
        ITHACAPOD::getWeightedModes(SnapMatrixNested[i], PModesNested[i], 0, 0, 0,
                                    NnestedOut);
    }

    for (int i = 0; i < Npar; i++)
    {
        y = PModesNested[i];

        for (int j = 0; j < y.size(); j++)
        {
            PModesGlobal.append(y[j]);
        }
    }
}

template<>
void ITHACAPOD::getModes(PtrList<volVectorField>& snapshotsU,
                         PtrList<volVectorField>& modes, bool podex, bool supex, bool sup, int nmodes)
{
    if (podex == 0)
    {
        ITHACAparameters para;

        if (para.eigensolver == "spectra" )
        {
            if (nmodes == 0)
            {
                nmodes = snapshotsU.size() - 2;
            }

            M_Assert(nmodes <= snapshotsU.size() - 2,
                     "The number of requested modes cannot be bigger than the number of Snapshots - 2");
        }

        else
        {
            if (nmodes == 0)
            {
                nmodes = snapshotsU.size();
            }

            M_Assert(nmodes <= snapshotsU.size(),
                     "The number of requested modes cannot be bigger than the number of Snapshots");
        }

        Eigen::MatrixXd SnapMatrix = Foam2Eigen::PtrList2Eigen(snapshotsU);
        List<Eigen::MatrixXd> SnapMatrixBC = Foam2Eigen::PtrList2EigenBC(snapshotsU);
        int NBC = snapshotsU[0].boundaryField().size();
        Eigen::VectorXd V = Foam2Eigen::field2Eigen(snapshotsU[0].mesh().V());
        Eigen::VectorXd V3d = (V.replicate(3, 1));
        auto VM = V3d.asDiagonal();
        Eigen::MatrixXd _corMatrix = SnapMatrix.transpose() * VM * SnapMatrix;

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
             snapshotsU[0].name() << " #######" << endl;
        int ncv = snapshotsU.size();
        Spectra::DenseSymMatProd<double> op(_corMatrix);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esEg;

        if (para.eigensolver == "spectra")
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

        else if (para.eigensolver == "eigen")
        {
            std::cout << "Using Eigen EigenSolver " << std::endl;
            esEg.compute(_corMatrix);
            M_Assert(esEg.info() == Eigen::Success,
                     "The Eigenvalue Decomposition did not succeed");
            eigenVectoreig = esEg.eigenvectors().real().rowwise().reverse().leftCols(
                                 nmodes);
            eigenValueseig = esEg.eigenvalues().real().reverse().head(nmodes);
        }

        Info << "####### End of the POD for " << snapshotsU[0].name() << " #######" <<
             endl;
        Eigen::VectorXd eigenValueseigLam =
            eigenValueseig.real().array().cwiseInverse().abs().sqrt() ;
        Eigen::MatrixXd modesEig = (SnapMatrix * eigenVectoreig) *
                                   eigenValueseigLam.asDiagonal();
        List<Eigen::MatrixXd> modesEigBC;
        modesEigBC.resize(NBC);

        for (int i = 0; i < NBC; i++)
        {
            modesEigBC[i] = (SnapMatrixBC[i] * eigenVectoreig) *
                            eigenValueseigLam.asDiagonal();
        }

        for (int i = 0; i < modes.size(); i++)
        {
            volVectorField tmp(snapshotsU[0].name(), snapshotsU[0] * 0);
            Eigen::VectorXd vec = modesEig.col(i);
            tmp = Foam2Eigen::Eigen2field(tmp, vec);

            for (int k = 0; k < NBC; k++)
            {
                ITHACAutilities::assignBC(tmp, k, modesEigBC[k].col(i));
            }

            modes.set(i, tmp);
        }

        eigenValueseig = eigenValueseig / eigenValueseig.sum();
        Eigen::VectorXd cumEigenValues(eigenValueseig);

        for (int j = 1; j < cumEigenValues.size(); ++j)
        {
            cumEigenValues(j) += cumEigenValues(j - 1);
        }

        Info << "####### Saving the POD bases for " << snapshotsU[0].name() <<
             " #######" << endl;

        //ITHACAPOD::exportBases(modes, snapshotsU, sup);
        if (sup)
        {
            ITHACAstream::exportFields(modes, "./ITHACAoutput/supremizer/",
                                       snapshotsU[0].name());
        }

        else
        {
            ITHACAstream::exportFields(modes, "./ITHACAoutput/POD/", snapshotsU[0].name());
        }

        Eigen::saveMarketVector(eigenValueseig,
                                "./ITHACAoutput/POD/Eigenvalues_" + snapshotsU[0].name(), para.precision,
                                para.outytpe);
        Eigen::saveMarketVector(cumEigenValues,
                                "./ITHACAoutput/POD/CumEigenvalues_" + snapshotsU[0].name(), para.precision,
                                para.outytpe);
    }

    else
    {
        Info << "Reading the existing modes" << endl;

        if (sup == 1)
        {
            ITHACAstream::read_fields(modes, snapshotsU[0], "./ITHACAoutput/supremizer/");
        }

        else
        {
            ITHACAstream::read_fields (modes, snapshotsU[0], "./ITHACAoutput/POD/");
        }
    }
}

template<>
void ITHACAPOD::getModes(PtrList<volScalarField>& snapshotsP,
                         PtrList<volScalarField>& modes, bool podex, bool supex, bool sup, int nmodes)
{
    if (podex == 0)
    {
        ITHACAparameters para;

        if (para.eigensolver == "spectra" )
        {
            if (nmodes == 0)
            {
                nmodes = snapshotsP.size() - 2;
            }

            M_Assert(nmodes <= snapshotsP.size() - 2,
                     "The number of requested modes cannot be bigger than the number of Snapshots - 2");
        }

        else
        {
            if (nmodes == 0)
            {
                nmodes = snapshotsP.size();
            }

            M_Assert(nmodes <= snapshotsP.size(),
                     "The number of requested modes cannot be bigger than the number of Snapshots");
        }

        Eigen::MatrixXd SnapMatrix = Foam2Eigen::PtrList2Eigen(snapshotsP);
        List<Eigen::MatrixXd> SnapMatrixBC = Foam2Eigen::PtrList2EigenBC(snapshotsP);
        int NBC = snapshotsP[0].boundaryField().size();
        // int NBC = 0;
        // if (Pstream::parRun())
        // {
        //     for (int i = 0; i < snapshotsP[0].boundaryField().size(); i++)
        //     {
        //         if (snapshotsP[0].boundaryField()[i].type() == "processor")
        //         {
        //             NBC++;
        //         }
        //     }
        // }
        // else
        // {
        //     NBC = snapshotsP[0].boundaryField().size();
        // }
        Eigen::VectorXd V = Foam2Eigen::field2Eigen(snapshotsP[0].mesh().V());
        auto VM = V.asDiagonal();
        Eigen::MatrixXd _corMatrix = SnapMatrix.transpose() * VM * SnapMatrix;

        if (Pstream::parRun())
        {
            List<double> vec(_corMatrix.data(), _corMatrix.data() + _corMatrix.size());
            reduce(vec, sumOp<List<double>>());
            std::memcpy(_corMatrix.data(), &vec[0], sizeof (double)*vec.size());
        }

        Eigen::VectorXd eigenValueseig;
        Eigen::MatrixXd eigenVectoreig;
        modes.resize(nmodes);
        Info << "####### Performing the POD using EigenDecomposition for " <<
             snapshotsP[0].name() << " #######" << endl;
        int ncv = snapshotsP.size();
        Spectra::DenseSymMatProd<double> op(_corMatrix);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esEg;

        if (para.eigensolver == "spectra")
        {
            std::cout << "Using Spectra EigenSolver " << std::endl;
            Spectra::SymEigsSolver<double, Spectra::LARGEST_ALGE, Spectra::DenseSymMatProd<double>>
                    es(&op, nmodes, ncv);
            es.init();
            es.compute(1000, 1e-10, Spectra::LARGEST_ALGE);
            M_Assert(es.info() == Spectra::SUCCESSFUL,
                     "The Eigenvalue Decomposition did not succeed");
            eigenVectoreig = es.eigenvectors().real();
            eigenValueseig = es.eigenvalues().real();
        }

        else if (para.eigensolver == "eigen")
        {
            std::cout << "Using Eigen EigenSolver " << std::endl;
            esEg.compute(_corMatrix);
            M_Assert(esEg.info() == Eigen::Success,
                     "The Eigenvalue Decomposition did not succeed");
            eigenVectoreig = esEg.eigenvectors().real().rowwise().reverse().leftCols(
                                 nmodes);
            eigenValueseig = esEg.eigenvalues().real().reverse().head(nmodes);
        }

        Info << "####### End of the POD for " << snapshotsP[0].name() << " #######" <<
             endl;
        Eigen::VectorXd eigenValueseigLam =
            eigenValueseig.real().array().cwiseInverse().abs().sqrt() ;
        Eigen::MatrixXd modesEig = (SnapMatrix * eigenVectoreig) *
                                   eigenValueseigLam.asDiagonal();
        List<Eigen::MatrixXd> modesEigBC;
        modesEigBC.resize(NBC);

        for (int i = 0; i < NBC; i++)
        {
            modesEigBC[i] = (SnapMatrixBC[i] * eigenVectoreig) *
                            eigenValueseigLam.asDiagonal();
        }

        for (label i = 0; i < modes.size(); i++)
        {
            volScalarField tmp(snapshotsP[0].name(), snapshotsP[0] * 0);
            Eigen::VectorXd vec = modesEig.col(i);
            tmp = Foam2Eigen::Eigen2field(tmp, vec);

            // Adjusting boundary conditions
            for (int k = 0; k < tmp.boundaryField().size(); k++)
            {
                ITHACAutilities::assignBC(tmp, k, modesEigBC[k].col(i));
            }

            // const fvPatchList & patches = tmp.mesh().boundary();
            // forAll(patches, iPatch)
            // {
            //     UList<scalar> myList(tmp.boundaryField()[iPatch].patchInternalField());
            //     forAll(tmp.boundaryFieldRef()[iPatch], iCell)
            //     {
            //         tmp.boundaryFieldRef()[iPatch][iCell] = myList[iCell];
            //     }
            // }
            modes.set(i, tmp);
        }

        eigenValueseig = eigenValueseig / eigenValueseig.sum();
        Eigen::VectorXd cumEigenValues(eigenValueseig);

        for (int j = 1; j < cumEigenValues.size(); ++j)
        {
            cumEigenValues(j) += cumEigenValues(j - 1);
        }

        Info << "####### Saving the POD bases for " << snapshotsP[0].name() <<
             " #######" << endl;

        //ITHACAPOD::exportBases(modes, snapshotsP, sup);
        if (sup)
        {
            ITHACAstream::exportFields(modes, "./ITHACAoutput/supremizer/",
                                       snapshotsP[0].name());
        }

        else
        {
            ITHACAstream::exportFields(modes, "./ITHACAoutput/POD/", snapshotsP[0].name());
        }

        Eigen::saveMarketVector(eigenValueseig,
                                "./ITHACAoutput/POD/Eigenvalues_" + snapshotsP[0].name(), para.precision,
                                para.outytpe);
        Eigen::saveMarketVector(cumEigenValues,
                                "./ITHACAoutput/POD/CumEigenvalues_" + snapshotsP[0].name(), para.precision,
                                para.outytpe);
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
            ITHACAstream::read_fields (modes, snapshotsP[0], "./ITHACAoutput/POD/");
        }
    }
}

template<>
void ITHACAPOD::getWeightedModes(PtrList<volVectorField>& snapshotsU,
                                 PtrList<volVectorField>& modes, bool podex, bool supex, bool sup, int nmodes)
{
    if (nmodes == 0)
    {
        nmodes = snapshotsU.size() - 2;
    }

    M_Assert(nmodes <= snapshotsU.size() - 2,
             "The number of requested modes cannot be bigger than the number of Snapshots - 2");

    if (podex == 0)
    {
        ITHACAparameters para;
        Eigen::MatrixXd SnapMatrix = Foam2Eigen::PtrList2Eigen(snapshotsU);
        List<Eigen::MatrixXd> SnapMatrixBC = Foam2Eigen::PtrList2EigenBC(snapshotsU);
        int NBC = snapshotsU[0].boundaryField().size();
        Eigen::VectorXd V = Foam2Eigen::field2Eigen(snapshotsU[0].mesh().V());
        Eigen::VectorXd V3d = (V.replicate(3, 1));
        auto VM = V3d.asDiagonal();
        Eigen::MatrixXd _corMatrix = SnapMatrix.transpose() * VM * SnapMatrix;
        Eigen::VectorXd eigenValueseig;
        Eigen::MatrixXd eigenVectoreig;
        modes.resize(nmodes);
        Info << "####### Performing the POD using EigenDecomposition " <<
             snapshotsU[0].name() << " #######" << endl;
        int ncv = snapshotsU.size();
        Spectra::DenseSymMatProd<double> op(_corMatrix);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esEg;
        Spectra::SymEigsSolver<double, Spectra::LARGEST_ALGE, Spectra::DenseSymMatProd<double>>
                es(&op, nmodes, ncv);

        if (para.eigensolver == "spectra")
        {
            std::cout << "Using Spectra EigenSolver " << std::endl;
            es.init();
            es.compute(1000, 1e-10, Spectra::LARGEST_ALGE);
            M_Assert(es.info() == Spectra::SUCCESSFUL,
                     "The Eigenvalue Decomposition did not succeed");
            eigenVectoreig = es.eigenvectors().real();
            eigenValueseig = es.eigenvalues().real();
        }

        else if (para.eigensolver == "eigen")
        {
            std::cout << "Using Eigen EigenSolver " << std::endl;
            esEg.compute(_corMatrix);
            M_Assert(esEg.info() == Eigen::Success,
                     "The Eigenvalue Decomposition did not succeed");
            eigenVectoreig = esEg.eigenvectors().real().rowwise().reverse().leftCols(
                                 nmodes);
            eigenValueseig = esEg.eigenvalues().real().reverse().head(nmodes);
        }

        Info << "####### End of the POD for " << snapshotsU[0].name() << " #######" <<
             endl;
        Eigen::VectorXd eigenValueseigLam =
            eigenValueseig.real().array().cwiseInverse().sqrt() ;
        Eigen::VectorXd eigenValueseigWeigted = eigenValueseig.real().array() ;
        Eigen::MatrixXd modesEig = (SnapMatrix * eigenVectoreig) *
                                   eigenValueseigLam.asDiagonal() * eigenValueseigWeigted.asDiagonal();
        List<Eigen::MatrixXd> modesEigBC;
        modesEigBC.resize(NBC);

        for (int i = 0; i < NBC; i++)
        {
            modesEigBC[i] = (SnapMatrixBC[i] * eigenVectoreig) *
                            eigenValueseigLam.asDiagonal() * eigenValueseigWeigted.asDiagonal();
        }

        for (int i = 0; i < modes.size(); i++)
        {
            volVectorField tmp(snapshotsU[0].name(), snapshotsU[0] * 0);
            Eigen::VectorXd vec = modesEig.col(i);
            tmp = Foam2Eigen::Eigen2field(tmp, vec);

            for (int k = 0; k < tmp.boundaryField().size(); k++)
            {
                ITHACAutilities::assignBC(tmp, k, modesEigBC[k].col(i));
            }

            // const fvPatchList & patches = tmp.mesh().boundary();
            // // Adjusting boundary conditions
            // forAll(patches, iPatch)
            // {
            //     UList<vector> myList(tmp.boundaryField()[iPatch].patchInternalField());
            //     forAll(tmp.boundaryFieldRef()[iPatch], iCell)
            //     {
            //         tmp.boundaryFieldRef()[iPatch][iCell] = myList[iCell];
            //     }
            // }
            modes.set(i, tmp);
        }

        eigenValueseig = eigenValueseig / eigenValueseig.sum();
        Eigen::VectorXd cumEigenValues(eigenValueseig);

        for (int j = 1; j < cumEigenValues.size(); ++j)
        {
            cumEigenValues(j) += cumEigenValues(j - 1);
        }

        Info << "####### Saving the POD bases for " << snapshotsU[0].name() <<
             " #######" << endl;

        //ITHACAPOD::exportBases(modes, snapshotsU, sup);
        if (sup)
        {
            ITHACAstream::exportFields(modes, "./ITHACAoutput/supremizer/",
                                       snapshotsU[0].name());
        }

        else
        {
            ITHACAstream::exportFields(modes, "./ITHACAoutput/POD/", snapshotsU[0].name());
        }

        Eigen::saveMarketVector(eigenValueseig,
                                "./ITHACAoutput/POD/Eigenvalues_" + snapshotsU[0].name(), para.precision,
                                para.outytpe);
        Eigen::saveMarketVector(cumEigenValues,
                                "./ITHACAoutput/POD/CumEigenvalues_" + snapshotsU[0].name(), para.precision,
                                para.outytpe);
    }

    else
    {
        Info << "Reading the existing modes" << endl;

        if (sup == 1)
        {
            ITHACAstream::read_fields(modes, snapshotsU[0], "./ITHACAoutput/supremizer/");
        }

        else
        {
            ITHACAstream::read_fields (modes, snapshotsU[0], "./ITHACAoutput/POD/");
        }
    }
}


template<>
void ITHACAPOD::getWeightedModes(PtrList<volScalarField>& snapshotsP,
                                 PtrList<volScalarField>& modes, bool podex, bool supex, bool sup, int nmodes)
{
    if (nmodes == 0)
    {
        nmodes = snapshotsP.size() - 2;
    }

    M_Assert(nmodes <= snapshotsP.size() - 2,
             "The number of requested modes cannot be bigger than the number of Snapshots - 2");

    if (podex == 0)
    {
        ITHACAparameters para;
        Eigen::MatrixXd SnapMatrix = Foam2Eigen::PtrList2Eigen(snapshotsP);
        List<Eigen::MatrixXd> SnapMatrixBC = Foam2Eigen::PtrList2EigenBC(snapshotsP);
        int NBC = snapshotsP[0].boundaryField().size();
        Eigen::VectorXd V = Foam2Eigen::field2Eigen(snapshotsP[0].mesh().V());
        auto VM = V.asDiagonal();
        Eigen::MatrixXd _corMatrix = SnapMatrix.transpose() * VM * SnapMatrix;
        Eigen::VectorXd eigenValueseig;
        Eigen::MatrixXd eigenVectoreig;
        modes.resize(nmodes);
        Info << "####### Performing the POD using EigenDecomposition for " <<
             snapshotsP[0].name() << " #######" << endl;
        int ncv = snapshotsP.size();
        Spectra::DenseSymMatProd<double> op(_corMatrix);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esEg;
        Spectra::SymEigsSolver<double, Spectra::LARGEST_ALGE, Spectra::DenseSymMatProd<double>>
                es(&op, nmodes, ncv);

        if (para.eigensolver == "spectra")
        {
            std::cout << "Using Spectra EigenSolver " << std::endl;
            es.init();
            es.compute(1000, 1e-10, Spectra::LARGEST_ALGE);
            M_Assert(es.info() == Spectra::SUCCESSFUL,
                     "The Eigenvalue Decomposition did not succeed");
            eigenVectoreig = es.eigenvectors().real();
            eigenValueseig = es.eigenvalues().real();
        }

        else if (para.eigensolver == "eigen")
        {
            std::cout << "Using Eigen EigenSolver " << std::endl;
            esEg.compute(_corMatrix);
            M_Assert(esEg.info() == Eigen::Success,
                     "The Eigenvalue Decomposition did not succeed");
            eigenVectoreig = esEg.eigenvectors().real().rowwise().reverse().leftCols(
                                 nmodes);
            eigenValueseig = esEg.eigenvalues().real().reverse().head(nmodes);
        }

        Info << "####### End of the POD for " << snapshotsP[0].name() << " #######" <<
             endl;
        Eigen::VectorXd eigenValueseigLam =
            eigenValueseig.real().array().cwiseInverse().sqrt() ;
        Eigen::VectorXd eigenValueseigWeigted = eigenValueseig.real().array() ;
        Eigen::MatrixXd modesEig = (SnapMatrix * eigenVectoreig) *
                                   eigenValueseigLam.asDiagonal() * eigenValueseigWeigted.asDiagonal();
        List<Eigen::MatrixXd> modesEigBC;
        modesEigBC.resize(NBC);

        for (int i = 0; i < NBC; i++)
        {
            modesEigBC[i] = (SnapMatrixBC[i] * eigenVectoreig) *
                            eigenValueseigLam.asDiagonal() * eigenValueseigWeigted.asDiagonal();
        }

        for (label i = 0; i < modes.size(); i++)
        {
            volScalarField tmp(snapshotsP[0].name(), snapshotsP[0] * 0);
            Eigen::VectorXd vec = modesEig.col(i);
            tmp = Foam2Eigen::Eigen2field(tmp, vec);

            // Adjusting boundary conditions
            for (int k = 0; k < tmp.boundaryField().size(); k++)
            {
                ITHACAutilities::assignBC(tmp, k, modesEigBC[k].col(i));
            }

            // const fvPatchList & patches = tmp.mesh().boundary();
            // forAll(patches, iPatch)
            // {
            //     UList<scalar> myList(tmp.boundaryField()[iPatch].patchInternalField());
            //     forAll(tmp.boundaryFieldRef()[iPatch], iCell)
            //     {
            //         tmp.boundaryFieldRef()[iPatch][iCell] = myList[iCell];
            //     }
            // }
            modes.set(i, tmp);
        }

        eigenValueseig = eigenValueseig / eigenValueseig.sum();
        Eigen::VectorXd cumEigenValues(eigenValueseig);

        for (int j = 1; j < cumEigenValues.size(); ++j)
        {
            cumEigenValues(j) += cumEigenValues(j - 1);
        }

        Info << "####### Saving the POD bases for " << snapshotsP[0].name() <<
             " #######" << endl;

        //ITHACAPOD::exportBases(modes, snapshotsP, sup);
        if (sup)
        {
            ITHACAstream::exportFields(modes, "./ITHACAoutput/supremizer/",
                                       snapshotsP[0].name());
        }

        else
        {
            ITHACAstream::exportFields(modes, "./ITHACAoutput/POD/", snapshotsP[0].name());
        }

        Eigen::saveMarketVector(eigenValueseig,
                                "./ITHACAoutput/POD/Eigenvalues_" + snapshotsP[0].name(), para.precision,
                                para.outytpe);
        Eigen::saveMarketVector(cumEigenValues,
                                "./ITHACAoutput/POD/CumEigenvalues_" + snapshotsP[0].name(), para.precision,
                                para.outytpe);
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
            ITHACAstream::read_fields (modes, snapshotsP[0], "./ITHACAoutput/POD/");
        }
    }
}

template<>
void ITHACAPOD::getModesSVD(PtrList<volVectorField>& snapshotsU,
                            PtrList<volVectorField>& modes, bool podex, bool supex, bool sup, int nmodes)
{
    if (podex == 0)
    {
        PtrList<volVectorField> Bases;
        modes.resize(nmodes);
        Info << "####### Performing POD using Singular Value Decomposition for " <<
             snapshotsU[0].name() << " #######" << endl;
        Eigen::MatrixXd SnapMatrix = Foam2Eigen::PtrList2Eigen(snapshotsU);
        Eigen::VectorXd V = Foam2Eigen::field2Eigen(snapshotsU[0].mesh().V());
        Eigen::VectorXd V3d = (V.replicate(3, 1));
        Eigen::VectorXd V3dSqrt = V3d.array().sqrt();
        Eigen::VectorXd V3dInv = V3dSqrt.array().cwiseInverse();
        auto VMsqr = V3dSqrt.asDiagonal();
        auto VMsqrInv = V3dInv.asDiagonal();
        Eigen::MatrixXd SnapMatrix2 = VMsqr * SnapMatrix;
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(SnapMatrix2,
                                              Eigen::ComputeThinU | Eigen::ComputeThinV);
        Info << "####### End of the POD for " << snapshotsU[0].name() << " #######" <<
             endl;
        Eigen::VectorXd eigenValueseig;
        Eigen::MatrixXd eigenVectoreig;
        eigenValueseig = svd.singularValues().real();
        eigenVectoreig = svd.matrixU().real();
        Eigen::MatrixXd modesEig = VMsqrInv * eigenVectoreig;
        volVectorField tmb_bu(snapshotsU[0].name(), snapshotsU[0] * 0);

        for (label i = 0; i < nmodes; i++)
        {
            Eigen::VectorXd vec = modesEig.col(i);
            tmb_bu = Foam2Eigen::Eigen2field(tmb_bu, vec);
            modes.set(i, tmb_bu);
        }

        eigenValueseig = eigenValueseig / eigenValueseig.sum();
        Eigen::VectorXd cumEigenValues(eigenValueseig);

        for (int j = 1; j < cumEigenValues.size(); ++j)
        {
            cumEigenValues(j) += cumEigenValues(j - 1);
        }

        Info << "####### Saving the POD bases for " << snapshotsU[0].name() <<
             " #######" << endl;
        ITHACAparameters para;

        if (sup)
        {
            ITHACAstream::exportFields(modes, "./ITHACAoutput/supremizer/",
                                       snapshotsU[0].name());
        }

        else
        {
            ITHACAstream::exportFields(modes, "./ITHACAoutput/POD/", snapshotsU[0].name());
        }

        //ITHACAPOD::exportBases(modes, snapshotsU, sup);
        Eigen::saveMarketVector(eigenValueseig,
                                "./ITHACAoutput/POD/Eigenvalues_" + snapshotsU[0].name(), para.precision,
                                para.outytpe);
        Eigen::saveMarketVector(cumEigenValues,
                                "./ITHACAoutput/POD/CumEigenvalues_" + snapshotsU[0].name(), para.precision,
                                para.outytpe);
    }

    else
    {
        Info << "Reading the existing modes" << endl;

        if (sup == 1)
        {
            ITHACAstream::read_fields(modes, snapshotsU[0], "./ITHACAoutput/supremizer/");
        }

        else
        {
            ITHACAstream::read_fields (modes, snapshotsU[0], "./ITHACAoutput/POD/");
        }
    }
}

template<>
void ITHACAPOD::getModesSVD(PtrList<volScalarField>& snapshotsP,
                            PtrList<volScalarField>& modes, bool podex, bool supex, bool sup, int nmodes)
{
    if (podex == 0)
    {
        PtrList<volVectorField> Bases;
        modes.resize(nmodes);
        Info << "####### Performing POD using Singular Value Decomposition for " <<
             snapshotsP[0].name() << " #######" << endl;
        Eigen::MatrixXd SnapMatrix = Foam2Eigen::PtrList2Eigen(snapshotsP);
        Eigen::VectorXd V = Foam2Eigen::field2Eigen(snapshotsP[0].mesh().V());
        Eigen::VectorXd VSqrt = V.array().sqrt();
        Eigen::VectorXd VInv = VSqrt.array().cwiseInverse();
        auto VMsqr = VSqrt.asDiagonal();
        auto VMsqrInv = VInv.asDiagonal();
        Eigen::MatrixXd SnapMatrix2 = VMsqr * SnapMatrix;
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(SnapMatrix2,
                                              Eigen::ComputeThinU | Eigen::ComputeThinV);
        Info << "####### End of the POD for " << snapshotsP[0].name() << " #######" <<
             endl;
        Eigen::VectorXd eigenValueseig;
        Eigen::MatrixXd eigenVectoreig;
        eigenValueseig = svd.singularValues().real();
        eigenVectoreig = svd.matrixU().real();
        Eigen::MatrixXd modesEig = VMsqrInv * eigenVectoreig;
        volScalarField tmb_bu(snapshotsP[0].name(), snapshotsP[0] * 0);

        for (label i = 0; i < nmodes; i++)
        {
            Eigen::VectorXd vec = modesEig.col(i);
            tmb_bu = Foam2Eigen::Eigen2field(tmb_bu, vec);
            modes.set(i, tmb_bu);
        }

        eigenValueseig = eigenValueseig / eigenValueseig.sum();
        Eigen::VectorXd cumEigenValues(eigenValueseig);

        for (int j = 1; j < cumEigenValues.size(); ++j)
        {
            cumEigenValues(j) += cumEigenValues(j - 1);
        }

        Info << "####### Saving the POD bases for " << snapshotsP[0].name() <<
             " #######" << endl;
        ITHACAparameters para;

        if (sup)
        {
            ITHACAstream::exportFields(modes, "./ITHACAoutput/supremizer",
                                       snapshotsP[0].name());
        }

        else
        {
            ITHACAstream::exportFields(modes, "./ITHACAoutput/POD", snapshotsP[0].name());
        }

        //ITHACAPOD::exportBases(modes, snapshotsP, sup);
        Eigen::saveMarketVector(eigenValueseig,
                                "./ITHACAoutput/POD/Eigenvalues_" + snapshotsP[0].name(), para.precision,
                                para.outytpe);
        Eigen::saveMarketVector(cumEigenValues,
                                "./ITHACAoutput/POD/CumEigenvalues_" + snapshotsP[0].name(), para.precision,
                                para.outytpe);
    }

    else
    {
        Info << "Reading the existing modes" << endl;

        if (sup == 1)
        {
            ITHACAstream::read_fields(modes, snapshotsP[0], "./ITHACAoutput/supremizer");
        }

        else
        {
            ITHACAstream::read_fields (modes, snapshotsP[0], "./ITHACAoutput/POD");
        }
    }
}

/// Normalize the bases
template<>
void ITHACAPOD::normalizeBases(PtrList<volScalarField>& Bases)
{
    scalar magSumSquare;

    for (label j = 0; j < Bases.size(); j++)
    {
        magSumSquare = Foam::sqrt(fvc::domainIntegrate(Bases[j] * Bases[j]).value());

        if (magSumSquare > SMALL)
        {
            Bases[j] /= magSumSquare;
        }

        Bases[j].correctBoundaryConditions();
    }
}

template<>
void ITHACAPOD::normalizeBases(PtrList<volVectorField>& Bases)
{
    scalar magSumSquare;

    for (label j = 0; j < Bases.size(); j++)
    {
        magSumSquare = Foam::sqrt(fvc::domainIntegrate(Bases[j] & Bases[j]).value());

        if (magSumSquare > SMALL)
        {
            Bases[j] /= magSumSquare;
        }

        //Bases[j].correctBoundaryConditions();
    }
}

/// Normalize the bases
void ITHACAPOD::normalizeBases(PtrList<volVectorField>& BasesU,
                               PtrList<volScalarField>& BasesP)
{
    scalar magSumSquare;

    for (label j = 0; j < BasesU.size(); j++)
    {
        magSumSquare = Foam::sqrt(fvc::domainIntegrate(BasesU[j] & BasesU[j]).value());

        if (magSumSquare > SMALL)
        {
            BasesP[j] /= magSumSquare;
        }

        BasesP[j].correctBoundaryConditions();
    }
}


/// Construct the Correlation Matrix for Scalar Field
template<>
Eigen::MatrixXd ITHACAPOD::corMatrix(PtrList<volScalarField>& snapshots)
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
Eigen::MatrixXd ITHACAPOD::corMatrix(PtrList<volVectorField>& snapshots)
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
Eigen::MatrixXd ITHACAPOD::corMatrix(List<Eigen::SparseMatrix<double>>&
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

            for (int k = 0; k < snapshots[i].cols(); k++)
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
Eigen::MatrixXd ITHACAPOD::corMatrix(List<Eigen::VectorXd>& snapshots)
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
template<>
void ITHACAPOD::exportBases(PtrList<volVectorField>& s,
                            PtrList<volVectorField>& bases, bool sup)
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


/// Export the Bases
template<>
void ITHACAPOD::exportBases(PtrList<volScalarField>& s,
                            PtrList<volScalarField>& bases, bool sup)
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

void ITHACAPOD::exportEigenvalues(scalarField Eigenvalues, fileName name,
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

void ITHACAPOD::exportcumEigenvalues(scalarField cumEigenvalues, fileName name,
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
        ITHACAPOD::DEIMmodes(List<Eigen::SparseMatrix<double>>& A,
                             List<Eigen::VectorXd>& b, int nmodesA, int nmodesB, word MatrixName)
{
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

        Eigen::MatrixXd corMatrixA = ITHACAPOD::corMatrix(A);
        Eigen::MatrixXd corMatrixB = ITHACAPOD::corMatrix(b);
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

        for (int i = 0; i < ModesA.size(); i++)
        {
            ITHACAstream::SaveSparseMatrix(ModesA[i],
                                           "./ITHACAoutput/DEIM/" + MatrixName + "/", "A_" + MatrixName + name(i));
        }

        for (int i = 0; i < ModesB.size(); i++)
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

void ITHACAPOD::GrammSchmidt(Eigen::MatrixXd& Matrix)
{
    Eigen::MatrixXd Ortho = Matrix;
    Ortho = Matrix;

    for (int i = 0; i <  Matrix.cols(); i++)
    {
        for (int k = 0; k < i; k++)
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

template<>
void ITHACAPOD::getModes(PtrList<volScalarField>& snapshots,
                         PtrList<volScalarField>& modes, PtrList<volScalarField>& Volumes, bool podex,
                         bool supex, bool sup, int nmodes)
{
    ITHACAparameters para;

    if (nmodes == 0 && para.eigensolver == "spectra")
    {
        nmodes = snapshots.size() - 2;
    }

    if (nmodes == 0 && para.eigensolver == "eigen")
    {
        nmodes = snapshots.size();
    }

    if (para.eigensolver == "spectra")
    {
        M_Assert(nmodes <= snapshots.size() - 2,
                 "The number of requested modes cannot be bigger than the number of Snapshots - 2");
    }

    if (podex == 0)
    {
        Eigen::MatrixXd SnapMatrix = Foam2Eigen::PtrList2Eigen(snapshots);
        List<Eigen::MatrixXd> SnapMatrixBC = Foam2Eigen::PtrList2EigenBC(snapshots);
        int NBC = snapshots[0].boundaryField().size();
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
        int ncv = snapshots.size();
        Spectra::DenseSymMatProd<double> op(_corMatrix);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esEg;

        if (para.eigensolver == "spectra")
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

        else if (para.eigensolver == "eigen")
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

        for (int i = 0; i < NBC; i++)
        {
            modesEigBC[i] = (SnapMatrixBC[i] * eigenVectoreig) *
                            eigenValueseigLam.asDiagonal();
        }

        for (label i = 0; i < modes.size(); i++)
        {
            volScalarField tmp(snapshots[0].name(), snapshots[0] * 0);
            Eigen::VectorXd vec = modesEig.col(i);
            tmp = Foam2Eigen::Eigen2field(tmp, vec);

            // Adjusting boundary conditions
            for (int k = 0; k < tmp.boundaryField().size(); k++)
            {
                ITHACAutilities::assignBC(tmp, k, modesEigBC[k].col(i));
            }

            modes.set(i, tmp);
        }

        eigenValueseig = eigenValueseig / eigenValueseig.sum();
        Eigen::VectorXd cumEigenValues(eigenValueseig);

        for (int j = 1; j < cumEigenValues.size(); ++j)
        {
            cumEigenValues(j) += cumEigenValues(j - 1);
        }

        Info << "####### Saving the POD bases for " << snapshots[0].name() << " #######"
             << endl;
        ITHACAPOD::exportBases(modes, snapshots, sup);
        Eigen::saveMarketVector(eigenValueseig,
                                "./ITHACAoutput/POD/Eigenvalues_" + snapshots[0].name(), para.precision,
                                para.outytpe);
        Eigen::saveMarketVector(cumEigenValues,
                                "./ITHACAoutput/POD/CumEigenvalues_" + snapshots[0].name(), para.precision,
                                para.outytpe);
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
            ITHACAstream::read_fields (modes, snapshots[0], "./ITHACAoutput/POD/");
        }
    }
}

template<>
void ITHACAPOD::getModes(PtrList<volVectorField>& snapshots,
                         PtrList<volVectorField>& modes, PtrList<volScalarField>& Volumes, bool podex,
                         bool supex, bool sup, int nmodes)
{
    ITHACAparameters para;

    if (nmodes == 0 && para.eigensolver == "spectra")
    {
        nmodes = snapshots.size() - 2;
    }

    if (nmodes == 0 && para.eigensolver == "eigen")
    {
        nmodes = snapshots.size();
    }

    if (para.eigensolver == "spectra")
    {
        M_Assert(nmodes <= snapshots.size() - 2,
                 "The number of requested modes cannot be bigger than the number of Snapshots - 2");
    }

    if (podex == 0)
    {
        Eigen::MatrixXd SnapMatrix = Foam2Eigen::PtrList2Eigen(snapshots);
        List<Eigen::MatrixXd> SnapMatrixBC = Foam2Eigen::PtrList2EigenBC(snapshots);
        int NBC = snapshots[0].boundaryField().size();
        Eigen::MatrixXd V = Foam2Eigen::PtrList2Eigen(Volumes);
        Eigen::MatrixXd V3d = (V.replicate(3, 1));
        Eigen::MatrixXd _corMatrix(snapshots.size(), snapshots.size());
        Info << "Filling the correlation matrix for field " << snapshots[0].name() <<
             endl;

        for (label i = 0; i < snapshots.size(); i++)
        {
            for (label j = 0; j <= i; j++)
            {
                Eigen::VectorXd Mij = (V3d.col(i).array() * V3d.col(j).array());
                Mij = Mij.array().abs().sqrt();
                _corMatrix(i, j) = SnapMatrix.col(i).transpose() * Mij.asDiagonal() *
                                   SnapMatrix.col(j);
            }
        }

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
        int ncv = snapshots.size();
        Spectra::DenseSymMatProd<double> op(_corMatrix);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esEg;

        if (para.eigensolver == "spectra")
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

        else if (para.eigensolver == "eigen")
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

        for (int i = 0; i < NBC; i++)
        {
            modesEigBC[i] = (SnapMatrixBC[i] * eigenVectoreig) *
                            eigenValueseigLam.asDiagonal();
        }

        for (label i = 0; i < modes.size(); i++)
        {
            volVectorField tmp(snapshots[0].name(), snapshots[0] * 0);
            Eigen::VectorXd vec = modesEig.col(i);
            tmp = Foam2Eigen::Eigen2field(tmp, vec);

            // Adjusting boundary conditions
            for (int k = 0; k < tmp.boundaryField().size(); k++)
            {
                ITHACAutilities::assignBC(tmp, k, modesEigBC[k].col(i));
            }

            modes.set(i, tmp);
        }

        eigenValueseig = eigenValueseig / eigenValueseig.sum();
        Eigen::VectorXd cumEigenValues(eigenValueseig);

        for (int j = 1; j < cumEigenValues.size(); ++j)
        {
            cumEigenValues(j) += cumEigenValues(j - 1);
        }

        Info << "####### Saving the POD bases for " << snapshots[0].name() << " #######"
             << endl;
        ITHACAPOD::exportBases(modes, snapshots, sup);
        Eigen::saveMarketVector(eigenValueseig,
                                "./ITHACAoutput/POD/Eigenvalues_" + snapshots[0].name(), para.precision,
                                para.outytpe);
        Eigen::saveMarketVector(cumEigenValues,
                                "./ITHACAoutput/POD/CumEigenvalues_" + snapshots[0].name(), para.precision,
                                para.outytpe);
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
            ITHACAstream::read_fields (modes, snapshots[0], "./ITHACAoutput/POD/");
        }
    }
}

