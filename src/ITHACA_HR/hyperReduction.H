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
#ifndef hyperReduction_H
#define hyperReduction_H

#include "fvCFD.H"
#include "ITHACAPOD.H"
#include "Foam2Eigen.H"
#include "EigenFunctions.H"
#include "ITHACAutilities.H"
#include "fvMeshSubset.H"
#include <set>

enum class HyperReductionMethod
{
    GappyDEIM
};

template <typename... SnapshotsLists>
class HyperReduction
{
public:
    const word methodNames[1] = {"GappyDEIM"};
    word methodName;

    using SnapshotsListTuple = std::tuple<std::decay_t<SnapshotsLists>...>;
    using FieldsTuple = std::tuple<typename std::decay_t<SnapshotsLists>::value_type ...>;
    static constexpr auto n_fields = sizeof...(SnapshotsLists);

    // TODO get vectorial dim from class template parameters
    unsigned int vectorial_dim;

    template <std::size_t N>
    using NthFieldListType = typename std::tuple_element<N, SnapshotsListTuple>::type;

    template <std::size_t N>
    using NthFieldType = typename NthFieldListType<N>::value_type;

    //----------------------------------------------------------------------
    /// @brief      Construct HyperReduction class, interpolation-based.`SnapshotsLists` is a variadic argument of PtrLists of fields e.g. for incompressible Navier-Stokes PtrList<volVectorField>, PtrList<volScalarField> for velocity and pressure.
    ///
    /// @param[in]  hrMethod          the chosen HR method
    /// @param[in]  n_modes           dimension of the HR basis
    /// @param[in]  n_nodes           number of nodes
    /// @param[in]  initialSeeds      initial nodes
    /// @param[in]  problemName       name of the function to be hyper-reduced
    /// @param[in]  snapshotsLists    lists of fields
    ///
    HyperReduction(label n_modes,
                   label n_nodes,
                   Eigen::VectorXi initialSeeds,
                   word problemName,
                   SnapshotsLists &&...snapshotsLists);
    
    // ITHACA parameters dict
    ITHACAparameters* para;

    /// The maximum number of modes to be considered
    label n_modes;

    /// The maximum number of modes to be considered
    label n_nodes;

    /// The length of the snapshots lists
    label n_snapshots;

    // Initial nodes
    Eigen::VectorXi initialSeeds;

    /// The name of the non-linear function e.g. HR_method/residual
    word problemName;

    /// The snapshots matrix containing the nonlinear function or operator
    SnapshotsListTuple snapshotsListTuple;

    /// Names of the fields
    List<word> fieldNames;

    /// Dimensions of the fields
    List<unsigned int> fieldDims;

    /// Int Number of Cells;
    label n_cells;

    /// Int Number of Cells in submeshes;
    label n_cellsSubfields;

    /// Nodes in the case of the a nonlinear function
    autoPtr<IOList<label>> nodePoints;

    /// List of label lists of the nodes and corresponding surrounding nodes
    autoPtr<IOList<labelList>> totalNodePoints;

    /// List of the unique indices of the nodes that define the submesh
    autoPtr<IOList<label>> uniqueNodePoints;

    /// Folder for nonlinear functions
    word Folder;

    /// Indices of the local node points in the subMesh
    List<label> localNodePoints;
    
    // Outputs
    /// Orthonormal basis of HR
    Eigen::MatrixXd basisMatrix;

    /// Renormalized basis of HR
    Eigen::MatrixXd renormalizedBasisMatrix;

    // TODO
    Eigen::MatrixXd pinvPU;

    /// The P matrix of the HyperReduction method
    Eigen::SparseMatrix<double> P;

    /// Normalizing weights of the offline snapshots
    Eigen::VectorXd normalizingWeights;

    /// Submesh of the HyperReduction method
    autoPtr<fvMeshSubset> submesh;
    
    /// Submeshes
    autoPtr<volVectorField> submesh_field;

    //----------------------------------------------------------------------
    /// @brief      TODO
    ///
    void offlineStage();
    
    //----------------------------------------------------------------------
    /// @brief      TODO
    ///
    void offlineGappyDEIM();

    //----------------------------------------------------------------------
    /// @brief      TODO
    ///
    void offlineEQP();

    //----------------------------------------------------------------------
    /// @brief      TODO
    ///
    /// @tparam     SnapshotsList  PtrList of volScalarField or volVectorField
    /// @param[in]  sList  The list of snapshots
    ///
    void getModesSVD(SnapshotsListTuple &SnapshotsListTuple, Eigen::MatrixXd &modesSVD, Eigen::VectorXd &fieldWeights);
    
    //----------------------------------------------------------------------
    /// @brief      TODO
    ///
    /// @tparam     SnapshotsList  PtrList of volScalarField or volVectorField
    /// @param[in]  sList  The list of snapshots
    ///
    void updateNodes(Eigen::SparseMatrix<double>& P, label& ind_max, Eigen::VectorXd& mp_not_mask, std::deque<label>& nodePointsDeque);

    //----------------------------------------------------------------------
    /// @brief      TODO
    ///
    /// @tparam     SnapshotsList  PtrList of volScalarField or volVectorField
    /// @param[in]  sList  The list of snapshots
    ///
    template <typename SnapshotsList>
    void stackSnapshots(SnapshotsList sList, Eigen::MatrixXd &snapshotsMatrix, Eigen::VectorXd &fieldWeights);

    //----------------------------------------------------------------------
    /// @brief      TODO
    ///
    /// @tparam     SnapshotsList  PtrList of volScalarField or volVectorField
    /// @param[in]  sList  The list of snapshots
    ///
    template <typename SnapshotsList>
    void stackNames(SnapshotsList sList){
        fieldNames.append(sList[0].name());
    }

    //----------------------------------------------------------------------
    /// @brief      TODO
    ///
    /// @tparam     SnapshotsList  PtrList of volScalarField or volVectorField
    /// @param[in]  sList  The list of snapshots
    ///
    template <typename SnapshotsList>
    inline void stackDimensions(SnapshotsList sList){
        fieldDims.append(get_field_dim<typename SnapshotsList::value_type>());}

    //----------------------------------------------------------------------
    /// @brief      Compute the pseudo-inverse of the matrix M restricted with the projector P
    ///
    /// @param[in]  Projector  projector into interpolation nodes
    /// @param[in]  Modes      matrix to restrict and invert
    ///
    void evaluatePinv(Eigen::SparseMatrix<double> &Projector, Eigen::MatrixXd &Modes, Eigen::VectorXd& fieldWeights);

    //----------------------------------------------------------------------
    /// @brief      Compute the submesh common to all fields in SnapshotsLists
    ///
    /// @param[in]  layers          projector into interpolation nodes
    /// @param[in]  mesh            matrix to restrict and invert
    ///
    void generateSubmesh(label layers, fvMesh &mesh);

    //----------------------------------------------------------------------
    /// @brief      Get local indices in the submesh from indices in the global ones
    ///
    /// @param      points       The points
    /// @param      submesh      The submesh
    ///
    /// @return     The local indices
    ///
    List<label> global2local(List<label> &points, fvMeshSubset &submesh);

    /// @brief      TODO
    ///
    /// @tparam     SnapshotsList  PtrList of volScalarField or volVectorField
    /// @param[in]  sList  The list of snapshots
    ///
    template <typename FieldType>
    inline autoPtr<FieldType> interpolateField(FieldType &field){
    return autoPtr<FieldType>(new FieldType(submesh->interpolate(field)));
}

    /// @brief      TODO
    ///
    /// @tparam     SnapshotsList  PtrList of volScalarField or volVectorField
    /// @param[in]  sList  The list of snapshots
    ///
    // TODO expand with sphericalTensor, symmTensor, tensor
    template <typename Field>
    inline constexpr unsigned int get_field_dim()
    {
        return std::is_same<Field, volScalarField>::value ? 1 : 3;
    }

    /// @brief      TODO
    ///
    /// @tparam     SnapshotsList  PtrList of volScalarField or volVectorField
    /// @param[in]  sList  The list of snapshots
    ///
    template <typename LastList>
    inline constexpr unsigned int compute_vectorial_dim(LastList x) { return get_field_dim<typename std::decay_t<LastList>::value_type>(); }

    /// @brief      TODO
    ///
    /// @tparam     SnapshotsList  PtrList of volScalarField or volVectorField
    /// @param[in]  sList  The list of snapshots
    ///
    template <typename List, typename... RemainingLists>
    inline constexpr unsigned int compute_vectorial_dim(List &&head, RemainingLists &&...tail)
    {
        return get_field_dim<typename std::decay_t<List>::value_type>() + compute_vectorial_dim<RemainingLists...>(tail...);
    }

    /// @brief      TODO
    ///
    /// @tparam     SnapshotsList  PtrList of volScalarField or volVectorField
    /// @param[in]  sList  The list of snapshots
    ///
    static inline double s_optimality(Eigen::MatrixXd& A)
    {   
        return std::pow(
            std::sqrt((A.transpose()*A).determinant())/\
            A.colwise().lpNorm<2>().prod(),
            1.0/A.cols());
    }

    /// @brief      TODO
    ///
    /// @tparam     SnapshotsList  PtrList of volScalarField or volVectorField
    /// @param[in]  sList  The list of snapshots
    ///
    template <typename... FieldsArgs>
    void eigen2fields(Eigen::VectorXd& eFields, FieldsArgs&&... oFields)
    {   
        unsigned int cumulativeSize{0}, ith_field{0};
        
        ([&](auto& oFieldsArg)
            {
                double ith_fieldSize = fieldDims[ith_field]*n_cells;
                Eigen::VectorXd vec = eFields.segment(cumulativeSize, ith_fieldSize);
                Foam2Eigen::Eigen2field(oFieldsArg, vec);
                cumulativeSize+=ith_fieldSize;
                ith_field++;
            }(oFields), ...);

    }
};

#endif
