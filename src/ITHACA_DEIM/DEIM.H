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

        //----------------------------------------------------------------------
        /// @brief      Construct DEIM for non-linear function
        ///
        /// @param[in]  SnapShotsMatrix  The snapshots matrix
        /// @param[in]  MaxModes         The maximum number of modes
        /// @param[in]  FunctionName     The function name, used to save the results
        /// @param[in]  FieldName        The field name
        ///
        DEIM (PtrList<T>& SnapShotsMatrix, label MaxModes, word FunctionName,
              word FieldName);

        //----------------------------------------------------------------------
        /// @brief      Construct DEIM for matrix with non-linear dependency
        ///
        /// @param      SnapShotsMatrix  The snapshots matrix
        /// @param[in]  MaxModesA        The maximum number of modes for the Matrix A
        /// @param[in]  MaxModesB        The maximum number of modes for the source term b
        /// @param[in]  MatrixName       The matrix name, used to save the results
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

        /// The maximum number of modes to be considered
        label MaxModes;
        label MaxModesA;
        label MaxModesB;

        /// The name of the non-linear function
        word FunctionName;

        /// A given tolerance
        double tolerance;

        /// Online vector
        ///@{
        Eigen::VectorXd theta;
        Eigen::VectorXd thetaA;
        Eigen::VectorXd thetaB;
        ///@}

        /// Int Number of Cells;
        label Ncells;

        /// source
        bool source;

        /// The name of the matrix
        word MatrixName;

        /// The matrix containing the modes
        Eigen::MatrixXd MatrixModes;

        /// Online Matrices
        ///@{
        Eigen::MatrixXd MatrixOnline;
        List<Eigen::SparseMatrix<double>> MatrixOnlineA;
        Eigen::MatrixXd MatrixOnlineB;
        ///@}

        /// The U matrix of the DEIM method
        ///@{
        Eigen::MatrixXd U;
        List<Eigen::SparseMatrix<double>> UA;
        Eigen::MatrixXd UB;
        ///@}

        /// Magic points in the case of the a nonlinear function
        autoPtr<IOList<label>> magicPoints;
        /// Magic points in the case of the a matrix function (rows indices)
        autoPtr<IOList<label>> magicPointsArow;
        /// Magic points in the case of the a matrix function (cols indices)
        autoPtr<IOList<label>> magicPointsAcol;
        /// Magic points in the case of the a matrix function, right hand side
        autoPtr<IOList<label>> magicPointsB;

        /// Magic points and indices of the surrounding layers
        ///@{
        autoPtr<IOList<labelList>> totalMagicPoints;
        autoPtr<IOList<labelList>> totalMagicPointsA;
        autoPtr<IOList<labelList>> totalMagicPointsB;
        ///@}

        /// List of unique indices that define the submesh
        ///@{
        autoPtr<IOList<label>> uniqueMagicPoints;
        autoPtr<IOList<label>> uniqueMagicPointsA;
        autoPtr<IOList<label>> uniqueMagicPointsB;
        ///@}

        /// Folder for nonlinear functions
        word Folder;
        ///Folder in the matrix case
        word FolderM;

        /// Definition of the x, y, z coordinate of the identified element in the matrix or source term
        /// 0 for x, 1 for y, 2 for z
        ///@{
        autoPtr<IOList<label>> xyz_Arow;
        autoPtr<IOList<label>> xyz_Acol;
        autoPtr<IOList<label>> xyz_B;
        ///@}

        /// Indices of the local magic points in the subMesh
        ///@{
        List<label> localMagicPoints;
        List<label> localMagicPointsArow;
        List<label> localMagicPointsAcol;
        List<label> localMagicPointsB;
        ///@}

        /// The P matrix of the DEIM method
        ///@{
        Eigen::SparseMatrix<double> P;
        List<Eigen::SparseMatrix<double>> PA;
        Eigen::SparseMatrix<double> PB;
        ///@}

        /// Submesh of the DEIM method
        ///@{
        autoPtr<fvMeshSubset> submesh;
        autoPtr<fvMeshSubset> submeshA;
        autoPtr<fvMeshSubset> submeshB;
        ///@}

        /// Bool variable to check if the SubMesh is available
        ///@{
        bool runSubMesh;
        bool runSubMeshA;
        bool runSubMeshB;
        ///@}

        //----------------------------------------------------------------------
        /// @brief      Function to generate the submesh for the nonlinear function case
        ///
        /// @param[in]  layers      Number of layers used to generate each submesh
        /// @param      mesh        The mesh of the problem
        /// @param[in]  field       The field used for the submesh
        /// @param[in]  secondTime  If the function is called for the second one must be set to 1
        ///
        /// @tparam     S           Type of object can be volscalarfield, volvectorfield etc...
        ///
        /// @return     return field in the submesh
        ///
        template <class S>
        S generateSubmesh(label layers, fvMesh& mesh, S field,
                          label secondTime = 0);

        //----------------------------------------------------------------------
        /// @brief      Function to generate the submesh for the nonlinear matrix function case (RHS)
        ///
        /// @param[in]  layers      Number of layers used to generate each submesh
        /// @param      mesh        The mesh of the problem
        /// @param[in]  field       The field used for the submesh
        /// @param[in]  secondTime  If the function is called for the second one must be set to 1
        ///
        /// @tparam     S           Type of object can be volscalarfield, volvectorfield etc...
        ///
        /// @return     return field in the submesh
        ///
        template <class S>
        S generateSubmeshMatrix(label layers, fvMesh& mesh, S field,
                                label secondTime = 0);

        //----------------------------------------------------------------------
        /// @brief      Function to generate the submesh for the nonlinear matrix function case (LHS)
        ///
        /// @param[in]  layers      Number of layers used to generate each submesh
        /// @param      mesh        The mesh of the problem
        /// @param[in]  field       The field used for the submesh
        /// @param[in]  secondTime  If the function is called for the second one must be set to 1
        ///
        /// @tparam     S           Type of object can be volscalarfield, volvectorfield etc...
        ///
        /// @return     return field in the submesh
        template <class S>
        S generateSubmeshVector(label layers, fvMesh& mesh, S field,
                                label secondTime = 0);

        //----------------------------------------------------------------------
        /// @brief      Function to generate a a subfield in the location of the magic points computed for the Matrix (LHS)
        ///
        /// @param      field  For which you want the create the submeshes
        ///
        /// @tparam     F      type of field (can be volVectorField, volScalarField, surfaceScalar and surfaceVector)
        ///
        /// @return     a the interpolated field in the submesh
        ///
        template <class F>
        F generateSubFieldMatrix(F& field);

        //----------------------------------------------------------------------
        /// @brief      Function to generate a a subfield in the location of the magic points computed for the Matrix (RHS)
        ///
        /// @param      field  For which you want the create the submeshes
        ///
        /// @tparam     F      type of field (can be volVectorField, volScalarField, surfaceScalar and surfaceVector)
        ///
        /// @return     a the interpolated field in the submesh
        ///
        template <class F>
        F generateSubFieldVector(F& field);

        //----------------------------------------------------------------------
        /// @brief      Function to get the onlineCoeffs of the DEIM approx. It is problem dependent so it must be overridden.
        ///
        void onlineCoeffs();

        //----------------------------------------------------------------------
        /// @brief      Get local indices in the submeshe from indices in the global ones
        ///
        /// @param      points       The points
        /// @param      submesh      The submesh
        ///
        /// @return     The local indices
        ///
        List<label> global2local(List<label>& points,
                                 fvMeshSubset& submesh);

        //----------------------------------------------------------------------
        /// @brief      check3DIndices in case of three dimensional fields
        ///
        /// @param      ind_rowA  The ind row a
        /// @param      ind_colA  The ind col a
        /// @param      xyz_rowA  The xyz row a
        /// @param      xyz_colA  The xyz col a
        ///
        void check3DIndices(label& ind_rowA, label&  ind_colA, label& xyz_rowA,
                            label& xyz_colA);

        //----------------------------------------------------------------------
        /// @brief      check3DIndices in case of three dimensional fields
        ///
        /// @param      ind_rowA  The ind row a
        /// @param      xyz_rowA  The xyz row a
        ///
        void check3DIndices(label& ind_rowA, label& xyz_rowA);

        //----------------------------------------------------------------------
        /// @brief      get the number of cells from the dimension of a LHS matrix
        ///
        /// @param      sizeM  The number of rows of the original Matrix
        ///
        /// @return     The number of cells
        ///
        label getNcells(label sizeM);

};

