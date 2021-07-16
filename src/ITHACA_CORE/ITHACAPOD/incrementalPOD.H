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
    incrementalPOD
Description
    Implementation of a incremental POD decomposition of a general field
SourceFiles
    incrementalPOD.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the incrementalPOD class.
/// \dir
/// Directory containing the header and source files for the incrementalPOD class.

#ifndef INCREMENTALPOD_H
#define INCREMENTALPOD_H

#include "fvCFD.H"
#include "ITHACAutilities.H"
#include "ITHACAstream.H"
#include "ITHACAparameters.H"
#include "Foam2Eigen.H"
#include "EigenFunctions.H"
#include "Modes.H"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wnon-virtual-dtor"
#include <Spectra/GenEigsSolver.h>
#include <Spectra/SymEigsSolver.h>
#include <Eigen/Eigen>
#include <unsupported/Eigen/SparseExtra>
#pragma GCC diagnostic pop

/*---------------------------------------------------------------------------*\
                        Class incrementalPOD Declaration
\*---------------------------------------------------------------------------*/

/// Implementation of a incremental POD algorithm according to
/// Oxberry et al. "Limited-memory adaptive snapshot selection for proper orthogonal
/// decomposition"

template<class Type, template<class> class PatchField, class GeoMesh>
class incrementalPOD : public Modes<Type, PatchField, GeoMesh>
{
    public:

        /// Tollerance for the incremental SVD
        double tolleranceSVD = 0;

        /// Rank of the reduced space
        int rank = 0;

        /// POD norm used for projection. It can be L2 or Frobenius
        word PODnorm;

        /// Vector of the cells volumes
        Eigen::VectorXd massVector;

        /// Vector of the singular values
        Eigen::VectorXd singularValues;

        /// Folder where the modes and singular values are saved
        word outputFolder = "./ITHACAoutput/incrementalPOD/";

        /// Constructors

        /// Construct Null
        incrementalPOD() {};

        //--------------------------------------------------------------------------
        /// @brief      Constructor with a first snapshot
        ///
        /// @param[in]  snapshot    A first snapshot
        /// @param[in]  _tol        Incremental SVD tollerance
        /// @param[in]  _maxRank   i Maximum POD rank
        ///
        incrementalPOD(GeometricField<Type, PatchField, GeoMesh>& snapshot,
                       double _tol, word _PODnorm);

        /// Destructor
        ~incrementalPOD() {};

        //--------------------------------------------------------------------------
        /// @brief      Initialize the incremental POD algorithm
        ///
        /// @param[in]  snapshot    A first snapshot
        ///
        void initialize(GeometricField<Type, PatchField, GeoMesh>& snapshot);

        //--------------------------------------------------------------------------
        /// @brief      Add a snapshot to the POD space
        ///
        /// @param[in]  snapshot    Snapshot to add
        ///
        void addSnapshot(GeometricField<Type, PatchField, GeoMesh>& snapshot);

        //--------------------------------------------------------------------------
        /// @brief      Fill the POD modes prtList from the Eigen matrix
        ///
        void fillPtrList();

        //--------------------------------------------------------------------------
        /// @brief      Project a field into the reduced space
        ///
        /// @param[in]  inputField  Field to project
        ///
        /// @return     The coefficient of the projected field into the reduced space
        ///
        Eigen::VectorXd project(
            GeometricField<Type, PatchField, GeoMesh>& inputField,
            label numberOfModes = 0);

        //--------------------------------------------------------------------------
        /// @brief      Reconstruct a field from the reduced coefficients
        ///
        /// @param[in/out]  inputField  Reconstructed field
        /// @param[in]      Coeff       Reduced coefficients
        /// @param[in]      Name        Name of the reconstructed field
        ///
        /// @return     Reconstructed field
        ///
        GeometricField<Type, PatchField, GeoMesh> reconstruct(
            GeometricField<Type, PatchField, GeoMesh>& inputField,
            Eigen::MatrixXd Coeff,
            word Name);

        //--------------------------------------------------------------------------
        /// @brief      Project a snapshot onto the POD space
        ///
        /// @param[in]  snapshot        Snapshot to project
        /// @param[in]  numberOfModes   Number of modes for projection (0 for all modes)
        ///
        /// @return     The projected snapshot
        ///
        GeometricField<Type, PatchField, GeoMesh> projectSnapshot(
            GeometricField<Type, PatchField, GeoMesh>& snapshot,
            label numberOfModes = 0);

        //--------------------------------------------------------------------------
        /// @brief     Project a list of snapshots onto the reduced space
        ///
        /// @param[in]      snapshots       Snapshots to be projected
        /// @param[in/out]  projSnapshots   Projected snapshots
        /// @param[in]  numberOfModes   Number of modes for projection (0 for all modes)
        ///
        void projectSnapshots(
            PtrList<GeometricField<Type, PatchField, GeoMesh>> snapshots,
            PtrList<GeometricField<Type, PatchField, GeoMesh>>& projSnapshots,
            label numberOfModes = 0);

        //--------------------------------------------------------------------------
        /// @brief      Write to modes to file
        ///
        void writeModes();
};

typedef incrementalPOD<scalar, fvPatchField, volMesh> scalarIncrementalPOD;
typedef incrementalPOD<vector, fvPatchField, volMesh> vectorIncrementalPOD;

namespace ITHACAPOD
{
//--------------------------------------------------------------------------
/// Weighted Gram-Schmidt orthogonalization
///
/// @param[in]  matrix  Matrix to be orthogonalized
/// @param[in]  weights Vector of weights
///
void weightedGramSchmidt(
    Eigen::MatrixXd& matrix,
    Eigen::VectorXd& weights)
{
    M_Assert(matrix.rows() == weights.rows(),
             "Matrix and weights must have the same number of rows");
    Eigen::MatrixXd Ortho = matrix;
    Ortho = matrix;

    for (label i = 0; i < matrix.cols(); i++)
    {
        for (label k = 0; k < i; k++)
        {
            double num = Ortho.col(k).transpose() * weights.asDiagonal()
                         * matrix.col(i);
            double den = (Ortho.col(k).transpose() * weights.asDiagonal()
                          * Ortho.col(k));
            double fact = num / den;
            Ortho.col(i) -= fact * Ortho.col(k) ;
        }

        Ortho.col(i).normalize();
    }

    matrix = Ortho;
};
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
