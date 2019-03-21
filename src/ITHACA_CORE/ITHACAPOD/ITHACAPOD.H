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
Class
    ITHACAPOD
Description
    Implementation of a POD decomposition of a general field
SourceFiles
    ITHACAPOD.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the ITHACAPOD class.
/// \dir
/// Directory containing the header and source files for the ITHACAPOD class.

#ifndef ITHACAPOD_H
#define ITHACAPOD_H

#include "fvCFD.H"
#include "ITHACAutilities.H"
#include "ITHACAstream.H"
#include "ITHACAparameters.H"
#include "Foam2Eigen.H"
#include "EigenFunctions.H"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wnon-virtual-dtor"
#include <Spectra/GenEigsSolver.h>
#include <Spectra/SymEigsSolver.h>
#include <Eigen/Eigen>
#include <unsupported/Eigen/SparseExtra>
#pragma GCC diagnostic pop

/*---------------------------------------------------------------------------*\
                        Class reductionProblem Declaration
\*---------------------------------------------------------------------------*/

/// Class of the computation of the POD, it exploits the method of snapshots.
class ITHACAPOD
{

    public:

        // Public Members

        //--------------------------------------------------------------------------
        /// Compute the bases or read them for a vector field
        ///
        /// @param[in]  snapshotsU  a PtrList of volVectorField snapshots.
        /// @param[out] modes      a PtrList of volVectorField where the modes are stored (it must be passed empty).
        /// @param[in]  podex       boolean variable 1 if the modes have been already computed and stored (in this case the function is reading them) 0 elsewhere.
        /// @param[in]  supex       boolean variable 1 if the supremizer modes have been already computed and stored (in this case the function is reading them) 0 elsewhere.
        /// @param[in]  sup         boolean variable 1 if you want to compute the supremizer modes 0 elsewhere.
        /// @param[in]  nmodes      int variable to set the number of modes to be stored, if set to 0 the maximum number of modes will computed.
        ///
        static void getModes(PtrList<volVectorField>& snapshotsU,
                             PtrList<volVectorField>& modes, bool podex, bool supex = 0, bool sup = 0,
                             int nmodes = 0);

        //--------------------------------------------------------------------------
        /// Compute the bases or read them for a scalar field
        ///
        /// @param[in]  snapshotsP  a PtrList of volScalarField snapshots.
        /// @param[out] modesP      a PtrList of volScalarField where the modes are stored (it must be passed empty).
        /// @param[in]  podex       boolean variable 1 if the modes have been already computed and stored (in this case the function is reading them) 0 elsewhere.
        /// @param[in]  supex       boolean variable 1 if the supremizer modes have been already computed and stored (in this case the function is reading them) 0 elsewhere.
        /// @param[in]  sup         boolean variable 1 if you want to compute the supremizer modes 0 elsewhere.
        /// @param[in]  nmodes      int variable to set the number of modes to be stored, if set to 0 the maximum number of modes will computed.
        ///
        static void getModes(PtrList<volScalarField>& snapshotsP,
                             PtrList<volScalarField>& modes, bool podex, bool supex = 0, bool sup = 0,
                             int nmodes = 0);

        //--------------------------------------------------------------------------
        /// @brief      Gets the bases for a scalar field using SVD instead of the method of snapshots
        ///
        /// @param[in]  snapshotsP  a PtrList of volScalarField snapshots.
        /// @param[out] modes      a PtrList of volScalarField where the modes are stored (it must be passed empty).
        /// @param[in]  podex       boolean variable 1 if the modes have been already computed and stored (in this case the function is reading them) 0 elsewhere.
        /// @param[in]  supex       boolean variable 1 if the supremizer modes have been already computed and stored (in this case the function is reading them) 0 elsewhere.
        /// @param[in]  sup         boolean variable 1 if you want to compute the supremizer modes 0 elsewhere.
        /// @param[in]  nmodes      int variable to set the number of modes to be stored, if set to 0 the maximum number of modes will computed.
        ///
        static void getModesSVD(PtrList<volScalarField>& snapshotsP,
                                PtrList<volScalarField>& modes, bool podex, bool supex = 0, bool sup = 0,
                                int nmodes = 0);
        //--------------------------------------------------------------------------
        /// @brief      Gets the bases for a vector field using SVD instead of the method of snapshots
        ///
        /// @param[in]  snapshotsP  a PtrList of volScalarField snapshots.
        /// @param[out] modes      a PtrList of volScalarField where the modes are stored (it must be passed empty).
        /// @param[in]  podex       boolean variable 1 if the modes have been already computed and stored (in this case the function is reading them) 0 elsewhere.
        /// @param[in]  supex       boolean variable 1 if the supremizer modes have been already computed and stored (in this case the function is reading them) 0 elsewhere.
        /// @param[in]  sup         boolean variable 1 if you want to compute the supremizer modes 0 elsewhere.
        /// @param[in]  nmodes      int variable to set the number of modes to be stored, if set to 0 the maximum number of modes will computed.
        ///
        static void getModesSVD(PtrList<volVectorField>& snapshotsU,
                                PtrList<volVectorField>& modes, bool podex, bool supex = 0, bool sup = 0,
                                int nmodes = 0);
        //--------------------------------------------------------------------------
        /// Nested-POD approach. Compute the nested snapshot matrix and weighted bases for a vector field
        ///
        /// @param[in]  snapshotsU      a PtrList of volVectorField snapshots.
        /// @param[out] UModesGlobal    a PtrList of volVectorField where the weighted modes of the nested snapshot matrix are stored (it must be passed empty).
        /// @param[in]  Npar            int variable which is equal with the number of parameters
        /// @param[in]  NnestedOut      int variable to set the number of weighted modes to be stored for the assembly of the nested matrix
        ///
        static void getNestedSnapshotMatrix(PtrList<volVectorField>& snapshotsU,
                                            PtrList<volVectorField>& UModesGlobal, int Npar, int NnestedOut);

        //--------------------------------------------------------------------------
        /// Nested-POD approach. Compute the nested snapshot matrix and weighted bases for a scalar field
        ///
        /// @param[in]  snapshotsP      a PtrList of volScalarField snapshots.
        /// @param[out] PModesGlobal    a PtrList of volScalarField where the weighted modes of the nested snapshot matrix are stored (it must be passed empty).
        /// @param[in]  Npar            int variable which is equal with the number of parameters
        /// @param[in]  NnestedOut      int variable to set the number of weighted modes to be stored for the assembly of the nested matrix
        ///
        static void getNestedSnapshotMatrix(PtrList<volScalarField>& snapshotsP,
                                            PtrList<volScalarField>& PModesGlobal, int Npar, int NnestedOut);

        //--------------------------------------------------------------------------

        /// Compute the weighted bases (using the nested-pod approach) or read them for a vector field
        ///
        /// @param[in]  snapshotsP  a PtrList of volVectorField snapshots.
        /// @param[out] modesP      a PtrList of volVectorField where the modes are stored (it must be passed empty).
        /// @param[in]  podex       boolean variable 1 if the modes have been already computed and stored (in this case the function is reading them) 0 elsewhere.
        /// @param[in]  supex       boolean variable 1 if the supremizer modes have been already computed and stored (in this case the function is reading them) 0 elsewhere.
        /// @param[in]  sup         boolean variable 1 if you want to compute the supremizer modes 0 elsewhere.
        /// @param[in]  nmodes      int variable to set the number of modes to be stored, if set to 0 the maximum number of modes will computed.
        ///
        static void getWeightedModes(PtrList<volVectorField>& snapshotsU,
                                     PtrList<volVectorField>& modesU, bool podex, bool supex = 0, bool sup = 0,
                                     int nmodes = 0);

        //--------------------------------------------------------------------------
        /// Compute the weighted bases (using the nested-pod approach) or read them for a scalar field
        ///
        /// @param[in]  snapshotsU  a PtrList of volScalarField snapshots.
        /// @param[out] modesU      a PtrList of volScalarField where the modes are stored (it must be passed empty).
        /// @param[in]  podex       boolean variable 1 if the modes have been already computed and stored (in this case the function is reading them) 0 elsewhere.
        /// @param[in]  supex       boolean variable 1 if the supremizer modes have been already computed and stored (in this case the function is reading them) 0 elsewhere.
        /// @param[in]  sup         boolean variable 1 if you want to compute the supremizer modes 0 elsewhere.
        ///
        static void getWeightedModes(PtrList<volScalarField>& snapshotsU,
                                     PtrList<volScalarField>& modesU, bool podex, bool supex = 0, bool sup = 0,
                                     int nmodes = 0);

        //--------------------------------------------------------------------------
        /// @brief      Perform GrammSchmidt orthonormalization on an Eigen Matrix
        ///
        /// @param      Matrix  The matrix
        ///
        static void GrammSchmidt(Eigen::MatrixXd& Matrix);

        //--------------------------------------------------------------------------
        /// Computes the correlation matrix given a vector field snapshot Matrix using the L2 norm
        ///
        /// @param[in] snapshots    a PtrList of volVectorField snapshots.
        ///
        /// @return    the Eigen::MatrixXd correlation matrix.
        ///
        static Eigen::MatrixXd corMatrix(PtrList<volVectorField>& snapshots);

        //--------------------------------------------------------------------------
        /// Computes the correlation matrix given a scalar field snapshot Matrix using the L2 norm.
        ///
        /// @param[in] snapshots    a PtrList of volScalarField snapshots.
        ///
        /// @return    the Eigen::MatrixXd correlation matrix.
        ///
        static Eigen::MatrixXd corMatrix(PtrList<volScalarField>& snapshots);

        //--------------------------------------------------------------------------
        /// Computes the correlation matrix given a list of matrices using the Frobenius inner product.
        ///
        /// @param[in]  snapshots  The snapshots matrices
        ///
        /// @return     the Eigen::MatrixXd correlation matrix.
        ///
        static Eigen::MatrixXd corMatrix(List<Eigen::SparseMatrix<double>>& snapshots);

        //--------------------------------------------------------------------------
        /// Computes the correlation matrix given a list of vectors using the Frobenius inner product.
        ///
        /// @param[in]  snapshots  The snapshots matrices
        ///
        /// @return     the Eigen::MatrixXd correlation matrix.
        ///
        static Eigen::MatrixXd corMatrix(List<Eigen::VectorXd>& snapshots);

        //--------------------------------------------------------------------------
        /// Export the basis for a vector field into the ITHACAOutput/POD or ITHACAOutput/supremizer
        ///
        /// @param[in] s        a PtrList of volVectorField where the snapshots associated with the bases are stored.
        /// @param[in] bases    a PtrList of volVectorField where the bases are stored.
        /// @param[in] sup      a boolean variable 1 if you want to export the supremizer bases (in ITHACAOutput/supremizer) 0 elsewhere (Default is 0).
        ///
        static void exportBases(PtrList<volVectorField>& s,
                                PtrList<volVectorField>& bases, bool sup = 0);

        //--------------------------------------------------------------------------
        /// Export the basis for a scalar field into the ITHACAOutput/POD or ITHACAOutput/supremizer
        ///
        /// @param[in] s        a PtrList of volScalarField where the snapshots associated with the bases are stored.
        /// @param[in] bases    a PtrList of volScalarField where the bases are stored.
        /// @param[in] sup      a boolean variable 1 if you want to export the supremizer bases (in ITHACAOutput/supremizer) 0 elsewhere (Default is 0).
        ///
        static void exportBases(PtrList<volScalarField>& s,
                                PtrList<volScalarField>& bases, bool sup = 0);

        //--------------------------------------------------------------------------
        /// Normalize the Bases for a vector field
        ///
        /// @param[in] s    a PtrList of volVectorField where the bases are stored.
        ///
        static void normalizeBases(PtrList<volVectorField>& s);

        //--------------------------------------------------------------------------
        /// Normalize the Bases for a scalar field
        ///
        /// @param[in] Bases    a PtrList of volScalarField where the bases are stored.
        ///
        static void normalizeBases(PtrList<volScalarField>& Bases);

        //--------------------------------------------------------------------------
        /// Normalize a volScalarField field using a volVectorField
        ///
        /// @param[in] Bases   a PtrList of volVectorField where the bases used for normalization are stored.
        /// @param[in] p       a PtrList of volScalarField  where the bases are stored.
        ///
        static void normalizeBases(PtrList<volVectorField>& Bases,
                                   PtrList<volScalarField>& p);

        //--------------------------------------------------------------------------
        /// Export the eigenvalues as a txt file
        ///
        /// @param[in] Eigenvalues  a scalarField of eigenvalues.
        /// @param[in] name         a fileName to indicate the name of the txt file.
        /// @param[in] sup          a boolean variable 1 if you want to export the eigenValues (in ITHACAOutput/supremizer) 0 elsewhere (Default is 0).
        ///
        static void exportEigenvalues(scalarField Eigenvalues, fileName name,
                                      bool sup = 0);

        //--------------------------------------------------------------------------
        /// Export the cumulative eigenvalues as a txt file
        ///
        /// @param[in] cumEigenvalues   a scalarField of cumEigenvalues.
        /// @param[in] name             a fileName to indicate the name of the txt file.
        /// @param[in] sup              a boolean variable 1 if you want to export the cumEigenValues (in ITHACAOutput/supremizer) 0 elsewhere (Default is 0).
        ///
        static void exportcumEigenvalues(scalarField cumEigenvalues,  fileName name,
                                         bool sup = 0);

        //--------------------------------------------------------------------------
        /// @brief      Get the DEIM modes for a generic non linear function
        ///
        /// @param[in]  SnapShotsMatrix  The snapshots matrix
        ///
        /// @tparam     T                templated object for the snapshots matrix type
        ///
        /// @return     T                the POD modes
        ///
        template<typename T>
        static PtrList<GeometricField<T, fvPatchField, volMesh>> DEIMmodes(
                    PtrList<GeometricField<T, fvPatchField, volMesh>>& SnapShotsMatrix, int nmodes,
                    word FunctionName);

        //--------------------------------------------------------------------------
        /// @brief      Get the DEIM modes for a generic non a parametrized matrix coming from a differential operator function
        ///
        /// @param[in]  MatrixList   The matrix list as a list of pointers
        /// @param[in]  nmodes       The number of modes
        /// @param[in]  MatrixName   The matrix name, used to save the matrix
        ///
        /// @tparam     the type of matrix, can be fvScalarMatrix or fvVectorMatrix
        ///
        /// @return     it returns a tuple containing the list of modes for the matrix and the list of modes for the source term
        ///
        template<typename type_matrix>
        static std::tuple<List<Eigen::SparseMatrix<double>>, List<Eigen::VectorXd>>
                DEIMmodes(PtrList<type_matrix>& MatrixList, int nmodesA, int nmodesB,
                          word MatrixName = "Matrix");

        //--------------------------------------------------------------------------
        /// @brief      Get the DEIM modes for a generic non a parametrized matrix coming from a differential operator function
        ///
        /// @param[in]  MatrixList   The matrix list as a list of pointers
        /// @param[in]  nmodes       The number of modes
        /// @param[in]  MatrixName   The matrix name, used to save the matrix
        ///
        /// @tparam     the type of matrix, can be fvScalarMatrix or fvVectorMatrix
        ///
        /// @return     it returns a tuple containing the list of modes for the matrix and the list of modes for the source term
        ///
        static std::tuple<List<Eigen::SparseMatrix<double>>, List<Eigen::VectorXd>>
                DEIMmodes(List<Eigen::SparseMatrix<double>>& A, List<Eigen::VectorXd>& b,
                          int nmodesA, int nmodesB, word MatrixName = "Matrix");

};

template<typename type_matrix>
std::tuple<List<Eigen::SparseMatrix<double>>, List<Eigen::VectorXd>>
        ITHACAPOD::DEIMmodes(PtrList<type_matrix>& MatrixList, int nmodesA, int nmodesB,
                             word MatrixName)
{
    List<Eigen::SparseMatrix<double>> ModesA(nmodesA);
    List<Eigen::VectorXd> ModesB(nmodesB);

    if (!ITHACAutilities::check_folder("./ITHACAoutput/DEIM/" + MatrixName))
    {
        M_Assert(nmodesA <= MatrixList.size() - 2
                 && nmodesB <= MatrixList.size() - 2,
                 "The number of requested modes cannot be bigger than the number of Snapshots - 2");
        std::tuple<List<Eigen::SparseMatrix<double>>, List<Eigen::VectorXd>> snapshots =
                    Foam2Eigen::LFvMatrix2LSM(MatrixList);
        Eigen::MatrixXd corMatrixA = ITHACAPOD::corMatrix(std::get<0>(snapshots));
        Eigen::MatrixXd corMatrixB = ITHACAPOD::corMatrix(std::get<1>(snapshots));
        Info << "####### Performing the POD decomposition for the Matrix List #######"
             << endl;
        Spectra::DenseSymMatProd<double> opA(corMatrixA);
        Spectra::DenseSymMatProd<double> opB(corMatrixB);
        int ncvA = MatrixList.size();
        int ncvB = MatrixList.size();
        Spectra::SymEigsSolver< double, Spectra::LARGEST_ALGE, Spectra:: DenseSymMatProd<double>>
                esA(&opA, nmodesA, ncvA);
        Spectra::SymEigsSolver< double, Spectra::LARGEST_ALGE, Spectra:: DenseSymMatProd<double>>
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
        Eigen::VectorXd eigenValueseigA;
        Eigen::MatrixXd eigenVectorseigA;
        Eigen::VectorXd eigenValueseigB;
        Eigen::MatrixXd eigenVectorseigB;
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

        for (int j = 1; j < cumEigenValuesA.size(); ++j)
        {
            cumEigenValuesA(j) += cumEigenValuesA(j - 1);
        }

        for (int j = 1; j < cumEigenValuesB.size(); ++j)
        {
            cumEigenValuesB(j) += cumEigenValuesB(j - 1);
        }

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

        ITHACAparameters para;
        Eigen::saveMarketVector(eigenValueseigA,
                                "./ITHACAoutput/DEIM/" + MatrixName + "/eigenValuesA", para.precision,
                                para.outytpe);
        Eigen::saveMarketVector(eigenValueseigB,
                                "./ITHACAoutput/DEIM/" + MatrixName + "/eigenValuesB", para.precision,
                                para.outytpe);
        Eigen::saveMarketVector(cumEigenValuesA,
                                "./ITHACAoutput/DEIM/" + MatrixName + "/cumEigenValuesA", para.precision,
                                para.outytpe);
        Eigen::saveMarketVector(cumEigenValuesB,
                                "./ITHACAoutput/DEIM/" + MatrixName + "/cumEigenValuesB", para.precision,
                                para.outytpe);
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


template<typename T>
PtrList<GeometricField<T, fvPatchField, volMesh>> ITHACAPOD::DEIMmodes(
            PtrList<GeometricField<T, fvPatchField, volMesh>>& SnapShotsMatrix, int nmodes,
            word FunctionName)
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

    PtrList<GeometricField<T, fvPatchField, volMesh>> Bases;
    PtrList<GeometricField<T, fvPatchField, volMesh>> modes;
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

    _corMatrix = ITHACAPOD::corMatrix(SnapShotsMatrix);
    Info << "####### Performing the POD decomposition for " <<
         SnapShotsMatrix[0].name() << " #######" << endl;
    int ncv = SnapShotsMatrix.size();
    Spectra::DenseGenMatProd<double> op(_corMatrix);
    Spectra::GenEigsSolver< double, Spectra::LARGEST_REAL, Spectra::DenseGenMatProd<double>>
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

    auto tmb_bu = SnapShotsMatrix[0];
    tmb_bu.rename(SnapShotsMatrix[0].name());

    for (label i = 0; i < nmodes; i++)
    {
        tmb_bu = eigenVector[i][0] * SnapShotsMatrix[0];

        for (label k = 1; k < SnapShotsMatrix.size(); k++)
        {
            tmb_bu += eigenVector[i][k] * SnapShotsMatrix[k];
        }

        Bases.set(i, tmb_bu);
        Info << "creating the bases " << i << " for " << SnapShotsMatrix[0].name() <<
             endl;
    }

    ITHACAPOD::normalizeBases(Bases);

    for (label i = 0; i < Bases.size(); i++)
    {
        auto tmp(Bases[i]);
        tmp.rename(SnapShotsMatrix[0].name());
        modes.set(i, tmp);
    }

    Info << "####### Saving the POD bases for " << SnapShotsMatrix[0].name() <<
         " #######" << endl;
    ITHACAutilities::createSymLink("./ITHACAoutput/DEIM");

    for (int i = 0; i < modes.size(); i++)
    {
        ITHACAstream::exportSolution(modes[i], name(i + 1), "./ITHACAoutput/DEIM",
                                     modes[i].name());
    }

    ITHACAstream::exportList(eigenValues, "./ITHACAoutput/DEIM/",
                             "eigenValues_" + SnapShotsMatrix[0].name());
    ITHACAstream::exportList(cumEigenValues, "./ITHACAoutput/DEIM/",
                             "cumEigenValues_" + SnapShotsMatrix[0].name());
    return modes;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif



