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
    HyperReduction
Description
    Implementation of the hyper-reduction methods.
SourceFiles
    hyperReduction.C
\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "ITHACAPOD.H"
#include "Foam2Eigen.H"
#include "EigenFunctions.H"
#include "ITHACAutilities.H"
#include "fvMeshSubset.H"
#include <set>

template <typename... SnapshotsLists>
class HyperReduction
{
public:

    enum class HyperReductionMethod { GNAT };
    const word methodNames[1] = {"GNAT"};
    const word methodName;

    using SnapshotsTuple = std::tuple<std::decay_t<SnapshotsLists>...>;
    using FieldsTuple = std::tuple<std::decay_t<SnapshotsLists>...>;
    static constexpr auto n_fields = sizeof...(SnapshotsLists);
    unsigned int vectorial_dim;

    template <std::size_t N>
    using NthFieldListType = typename std::tuple_element<N, SnapshotsTuple>::type;

    template <std::size_t N>
    using NthFieldType = typename NthFieldListType<N>::value_type;

    //----------------------------------------------------------------------
    /// @brief      Construct HyperReduction.
    ///
    /// @param[in]  snapshotsLists  The snapshots matrix
    /// @param[in]  maxModes         The maximum number of modes
    ///
    HyperReduction(HyperReductionMethod hrMethod,
                   label maxModes,
                   label maxNodes,
                   Eigen::VectorXi initialSeeds,
                   word functionName,
                   SnapshotsLists&&... snapshotsLists);

    /// The maximum number of modes to be considered
    label maxModes;

    /// The maximum number of modes to be considered
    label maxNodes;

    label n_snapshots;
    
    // Initial nodes
    Eigen::VectorXi initialSeeds;

    /// The name of the non-linear function e.g. HR_method/residual
    word functionName;

    /// The snapshots matrix containing the nonlinear function or operator
    SnapshotsTuple snapshotsTuple;

    /// The POD modes of the HyperReduction procedure that can be
    SnapshotsTuple modes;

    /// Matrix Modes
    Eigen::MatrixXd matrixModes;

    // TODO
    Eigen::MatrixXd pinv_PU;

    /// TODO
    Eigen::MatrixXd maskAux;

    /// TODO
    std::set<label> nodePointsSet;

    /// Int Number of Cells;
    label n_cells;

    /// Int Number of Cells in submeshes;
    label n_cellsSubfields;

    /// Submeshes
    autoPtr<volVectorField> submesh_field;

    /// Online Matrices
    ///@{
    Eigen::MatrixXd MatrixOnline;

    /// The U matrix of the HyperReduction method
    Eigen::MatrixXd U;

    /// Magic points in the case of the a nonlinear function
    autoPtr<IOList<label>> nodePoints;
    /// Magic points in the case of the a nonlinear function
    autoPtr<IOList<label>> nodePointsFields;
    /// Magic points in the case of the a matrix function (rows indices)
    autoPtr<IOList<label>> nodePointsArow;
    /// Magic points in the case of the a matrix function (cols indices)
    autoPtr<IOList<label>> nodePointsAcol;
    /// Magic points in the case of the a matrix function, right hand side
    autoPtr<IOList<label>> nodePointsB;

    /// Magic points and indices of the surrounding layers
    ///@{
    autoPtr<IOList<labelList>> totalnodePoints;
    autoPtr<IOList<labelList>> totalnodePointsA;
    autoPtr<IOList<labelList>> totalnodePointsB;
    ///@}

    /// List of unique indices that define the submesh
    ///@{
    autoPtr<IOList<label>> uniquenodePoints;
    autoPtr<IOList<label>> uniquenodePointsA;
    autoPtr<IOList<label>> uniquenodePointsB;
    ///@}

    /// Folder for nonlinear functions
    word Folder;
    /// Folder in the matrix case
    word FolderM;

    /// Definition of the x, y, z coordinate of the identified element in the matrix or source term
    /// 0 for x, 1 for y, 2 for z
    ///@{
    autoPtr<IOList<label>> xyz;
    autoPtr<IOList<label>> xyz_Arow;
    autoPtr<IOList<label>> xyz_Acol;
    autoPtr<IOList<label>> xyz_B;
    ///@}

    /// Indices of the local magic points in the subMesh
    ///@{
    List<label> localnodePoints;
    List<label> localnodePointsArow;
    List<label> localnodePointsAcol;
    List<label> localnodePointsB;
    ///@}

    /// The P matrix of the HyperReduction method
    ///@{
    Eigen::SparseMatrix<double> P;
    List<Eigen::SparseMatrix<double>> PA;
    Eigen::SparseMatrix<double> PB;
    ///@}

    /// Submesh of the HyperReduction method
    ///@{
    autoPtr<fvMeshSubset> submesh;
    autoPtr<fvMeshSubset> submeshA;
    autoPtr<fvMeshSubset> submeshB;

    /// Non-repeated submeshes cells (only x coordinates)
    Eigen::VectorXi nonRepeatedSubmeshes;

    /// Non-repeated submeshes cells (only x, y coordinates)
    Eigen::VectorXi nonRepeatedSubmeshes_xy;

    void compute_basis();
                                      
    void evaluatePinv(Eigen::SparseMatrix<double> &Projector, Eigen::MatrixXd &Modes);
    
    voidgenerateSubmeshes(label layers, fvMesh &mesh, label secondTime);

    List<label> global2local(List<label> &points, fvMeshSubset &submesh);

    template<typename... SnapshotsLists>
    void createMasks()

    template<typename FieldType>
    void interpolateField(FieldType& field)

    template<typename SnapshotsList>
    void compute_matrix_modes(SnapshotsList sList);
    template<typename Field>
    inline constexpr unsigned int get_field_dim(){return std::is_same<typename Field::DimensionedField::FieldType, Foam::scalar>::value? 1: 3;}

    template <typename LastList>
    inline constexpr unsigned int compute_vectorial_dim(LastList x) {return get_field_dim<typename std::decay_t<LastList>::value_type>();}

    template<typename List, typename... RemainingLists>
    inline constexpr unsigned int compute_vectorial_dim(List&& head, RemainingLists&&... tail)  {
        return get_field_dim<typename std::decay_t<List>::value_type>() + compute_vectorial_dim<RemainingLists...>(tail...);}
};
