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
    DEIM
Description
    Implementation of the discrete empirical interpolation method
SourceFiles
    DEIM.C
\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "ITHACAPOD.H"
#include "Foam2Eigen.H"
#include "EigenFunctions.H"
#include "ITHACAutilities.H"
#include "fvMeshSubset.H"


template<typename T>
class DEIM
{
    public:


        ///
        /// @brief      Construct DEIM for non-linear function
        ///
        /// @param[in]  SnapShotsMatrix  The snapshots matrix
        /// @param[in]  MaxModes         The maximum number of modes
        /// @param[in]  FunctionName     The function name
        ///
        DEIM (PtrList<T>& SnapShotsMatrix, label MaxModes, word FunctionName,
              word FieldName);

        ///
        /// @brief      Construct DEIM for matrix with non-linear dependency
        ///
        /// @param      SnapShotsMatrix  The snapshots matrix
        /// @param[in]  MaxModesA        The maximum number of modes for the Matrix A
        /// @param[in]  MaxModesB        The maximum number of modes for the source term b
        /// @param[in]  MatrixName       The matrix name
        ///
        DEIM (PtrList<T>& SnapShotsMatrix, label MaxModesA, label MaxModesB,
              word MatrixName);

        /// The snapshots matrix containing the nonlinear function or operator
        PtrList<T> SnapShotsMatrix;

        /// The POD modes of the DEIM procedure that can be
        PtrList<T> modes;

        /// Matrix Modes
        std::tuple<List<Eigen::SparseMatrix<double>>, List<Eigen::VectorXd>>
                Matrix_Modes;

        /// Field used to evaluate the online non linear function
        //PtrList<S> fields;
        //PtrList<S> fieldsA;
        //PtrList<S> fieldsB;

        /// The maximum number of modes to be considered
        label MaxModes;
        label MaxModesA;
        label MaxModesB;

        /// The name of the non-linear function
        word FunctionName;

        /// A given tolerance
        double tolerance;

        /// Online vector
        Eigen::VectorXd theta;
        Eigen::VectorXd thetaA;
        Eigen::VectorXd thetaB;

        /// Int Number of Cells;
        label Ncells;

        /// source
        bool source;

        /// The name of the matrix
        word MatrixName;

        /// The matrix containing the modes
        Eigen::MatrixXd MatrixModes;

        /// Online Matrix
        Eigen::MatrixXd MatrixOnline;
        List<Eigen::SparseMatrix<double>> MatrixOnlineA;
        Eigen::MatrixXd MatrixOnlineB;

        /// The U matrix of the DEIM method
        Eigen::MatrixXd U;
        List<Eigen::SparseMatrix<double>> UA;
        Eigen::MatrixXd UB;

        /// Indices of the magic points
        List<label> magicPoints;
        List<Pair <label >> magicPointsA;
        List<label> magicPointsB;
        List<Pair<label>> xyz_A;
        List<label> xyz_B;

        /// Indices of the local magic points in the subMesh
        List<label> localMagicPoints;
        List<Pair <label>> localMagicPointsA;
        List<label> localMagicPointsB;

        /// The P matrix of the DEIM method
        Eigen::SparseMatrix<double> P;
        List<Eigen::SparseMatrix<double>> PA;
        Eigen::SparseMatrix<double> PB;

        /// List of submeshes
        PtrList<fvMeshSubset> submeshList;
        PtrList<fvMeshSubset> submeshListA;
        PtrList<fvMeshSubset> submeshListB;

        /// List of bool variables to check if the SubMeshes are available
        bool runSubMesh;
        bool runSubMeshA;
        bool runSubMeshB;

        /// @brief      Function to generate the submeshes
        ///
        /// @param[in]  layers  Number of layers used to generate each submesh
        /// @param      mesh    The mesh of the problem
        ///
        template <class S>
        PtrList<S> generateSubmeshes(label layers, fvMesh& mesh, S field,
                                     label secondTime = 0);

        ///
        /// @brief      { function_description }
        ///
        /// @param[in]  layers  The layers
        /// @param      mesh    The mesh
        /// @param[in]  field   The field
        ///
        template <class S>
        PtrList<S> generateSubmeshesMatrix(label layers, fvMesh& mesh, S field,
                                           label secondTime = 0);

        ///
        /// @brief      { function_description }
        ///
        /// @param[in]  layers  The layers
        /// @param      mesh    The mesh
        /// @param[in]  field   The field
        ///
        template <class S>
        PtrList<S> generateSubmeshesVector(label layers, fvMesh& mesh, S field,
                                           label secondTime = 0);

        ///
        /// @brief      Function to generate a PtrList of subfield in the location of the magic points computed for the Matrix (LHS)
        ///
        /// @param      field  For which you want the create the submeshes
        ///
        /// @tparam     F      type of field (can be volVectorField, volScalarField, surfaceScalar and surfaceVector)
        ///
        /// @return     a PtrList with the fields
        ///
        template <class F>
        PtrList<F> generateSubFieldsMatrix(F& field);

        ///
        /// @brief      Function to generate a PtrList of subfield in the location of the magic points computed for the Vector (RHS)
        ///
        /// @param      field  For which you want the create the submeshes
        ///
        /// @tparam     F      type of field (can be volVectorField, volScalarField, surfaceScalar and surfaceVector)
        ///
        /// @return     a PtrList with the fields
        ///
        template <class F>
        PtrList<F> generateSubFieldsVector(F& field);

        ///
        /// @brief      Function to get the onlineCoeffs of the DEIM approx. It is problem dependent so it must be overridden.
        ///
        void onlineCoeffs();

        ///
        /// @brief      Get local indices in the submeshes from indices in the global ones
        ///
        /// @param      points       The points
        /// @param      submeshList  The submesh list
        ///
        /// @return     The local indices
        ///
        List<label> global2local(List<label>& points,
                                 PtrList<fvMeshSubset>& submeshList);

        ///
        /// @brief      Get local indices in the submeshes from indices in the global ones
        ///
        /// @param      points       The points
        /// @param      submeshList  The submesh list
        ///
        /// @return     The local indices
        ///
        List<Pair <label >> global2local(List< Pair <label>>& points,
                                         PtrList<fvMeshSubset>& submeshList);

        ///
        /// @brief      check3DIndices in case of three dimensional fields
        ///
        /// @param      ind_rowA  The ind row a
        /// @param      ind_colA  The ind col a
        /// @param      xyz_rowA  The xyz row a
        /// @param      xyz_colA  The xyz col a
        ///
        void check3DIndices(label& ind_rowA, label&  ind_colA, label& xyz_rowA,
                            label& xyz_colA);

        ///
        /// @brief      check3DIndices in case of three dimensional fields
        ///
        /// @param      ind_rowA  The ind row a
        /// @param      xyz_rowA  The xyz row a
        ///
        void check3DIndices(label& ind_rowA, label& xyz_rowA);

        /// @brief      get the number of cells from the dimension of a LHS matrix
        ///
        /// @param      sizeM  The number of rows of the original Matrix
        ///
        /// @return     The number of cells
        ///
        label getNcells(label sizeM);

};

