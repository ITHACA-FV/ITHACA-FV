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
        DEIM (PtrList<T>& SnapShotsMatrix, int MaxModes, word FunctionName);

        ///
        /// @brief      Construct DEIM for matrix with non-linear dependency
        ///
        /// @param      SnapShotsMatrix  The snapshots matrix
        /// @param[in]  MaxModesA        The maximum number of modes for the Matrix A
        /// @param[in]  MaxModesB        The maximum number of modes for the source term b
        /// @param[in]  MatrixName       The matrix name
        ///
        DEIM (PtrList<T>& SnapShotsMatrix, int MaxModesA, int MaxModesB,
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
        int MaxModes;
        int MaxModesA;
        int MaxModesB;

        /// The name of the non-linear function
        word FunctionName;

        /// A given tolerance
        double tolerance;

        /// Online vector
        Eigen::VectorXd theta;
        Eigen::VectorXd thetaA;
        Eigen::VectorXd thetaB;

        /// Int size M;
        int sizeM;

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
        List<int> magicPoints;
        List<Pair <int >> magicPointsA;
        List<int> magicPointsB;
        List<Pair<int>> xyz_A;
        List<int> xyz_B;

        /// Indices of the local magic points in the subMesh
        List<int> localMagicPoints;
        List<Pair <int>> localMagicPointsA;
        List<int> localMagicPointsB;

        /// The P matrix of the DEIM method
        Eigen::SparseMatrix<double> P;
        List<Eigen::SparseMatrix<double>> PA;
        Eigen::SparseMatrix<double> PB;

        /// List of submeshes
        PtrList<fvMeshSubset> submeshList;
        PtrList<fvMeshSubset> submeshListA;
        PtrList<fvMeshSubset> submeshListB;

        /// @brief      Function to generate the submeshes
        ///
        /// @param[in]  layers  Number of layers used to generate each submesh
        /// @param      mesh    The mesh of the problem
        ///
        template <class S>
        PtrList<S> generateSubmeshes(int layers, fvMesh& mesh, S field,
                                     int secondTime = 0);

        ///
        /// @brief      { function_description }
        ///
        /// @param[in]  layers  The layers
        /// @param      mesh    The mesh
        /// @param[in]  field   The field
        ///
        template <class S>
        PtrList<S> generateSubmeshesMatrix(int layers, fvMesh& mesh, S field,
                                           int secondTime = 0);

        ///
        /// @brief      { function_description }
        ///
        /// @param[in]  layers  The layers
        /// @param      mesh    The mesh
        /// @param[in]  field   The field
        ///
        template <class S>
        PtrList<S> generateSubmeshesVector(int layers, fvMesh& mesh, S field,
                                           int secondTime = 0);

        ///
        /// @brief      Function to get the onlineCoeffs of the DEIM approx. It is problem dependent so it must be overridden.
        ///
        void onlineCoeffs();

        ///
        /// @brief      { function_description }
        ///
        /// @param      points       The points
        /// @param      submeshList  The submesh list
        ///
        /// @return     { description_of_the_return_value }
        ///
        List<int> global2local(List<int>& points, PtrList<fvMeshSubset>& submeshList);

        ///
        /// @brief      { function_description }
        ///
        /// @param      points       The points
        /// @param      submeshList  The submesh list
        ///
        /// @return     { description_of_the_return_value }
        ///
        List<Pair <int >> global2local(List< Pair <int>>& points,
                                       PtrList<fvMeshSubset>& submeshList);

        ///
        /// @brief      { function_description }
        ///
        /// @param      ind_rowA  The ind row a
        /// @param      ind_colA  The ind col a
        /// @param      xyz_rowA  The xyz row a
        /// @param      xyz_colA  The xyz col a
        ///
        void check3DIndices(int& ind_rowA, int&  ind_colA, int& xyz_rowA,
                            int& xyz_colA);

        ///
        /// @brief      { function_description }
        ///
        /// @param      ind_rowA  The ind row a
        /// @param      xyz_rowA  The xyz row a
        ///
        void check3DIndices(int& ind_rowA, int& xyz_rowA);

};

