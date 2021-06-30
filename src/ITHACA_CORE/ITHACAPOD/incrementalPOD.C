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
/// source file for the incrementalPOD class

#include "incrementalPOD.H"
#include "EigenFunctions.H"

template<class Type, template<class> class PatchField, class GeoMesh>
incrementalPOD<Type, PatchField, GeoMesh>::incrementalPOD(
    GeometricField<Type, PatchField, GeoMesh>& snapshot,
    double _tol, word _PODnorm)
{
    Info << "WARNING: the projection of the BC has not been implemented yet!" <<
         endl;
    tolleranceSVD = _tol;
    PODnorm = _PODnorm;
    M_Assert(PODnorm == "L2" ||
             PODnorm == "Frobenius", "The PODnorm can be only L2 or Frobenius");
    initialize(snapshot);
}

template<class Type, template<class> class PatchField, class GeoMesh>
void incrementalPOD<Type, PatchField, GeoMesh>::initialize(
    GeometricField<Type, PatchField, GeoMesh>& snapshot)
{
    Info << "Initializing the incremental POD" << endl;
    M_Assert(tolleranceSVD > 0, "Set up the tollerance before initialization");
    M_Assert(rank == 0, "POD already initialized");
    double snapNorm = ITHACAutilities::L2Norm(snapshot);

    if (PODnorm == "L2")
    {
        Eigen::VectorXd snapshotEig = Foam2Eigen::field2Eigen(snapshot);
        massVector = ITHACAutilities::getMassMatrixFV(snapshot);
        singularValues.resize(1);
        singularValues[0] = Foam::sqrt((snapshotEig.transpose() *
                                        massVector.asDiagonal() *
                                        snapshotEig)(0, 0));
        snapshotEig = snapshotEig / singularValues[0];
        GeometricField<Type, PatchField, GeoMesh>  tmp(snapshot);
        tmp = Foam2Eigen::Eigen2field(tmp, snapshotEig);
        this->append(tmp.clone());
        this->toEigen();
        rank = 1;
        fillPtrList();
    }
    else
    {
        if (snapNorm > tolleranceSVD)
        {
            singularValues.resize(1);
            singularValues[0] = snapNorm;
            this->append(snapshot / snapNorm);
            this->toEigen();
            rank = 1;
        }
        else
        {
            rank = 0;
        }
    }

    Info << "Initialization ended" << endl;
}

template<class Type, template<class> class PatchField, class GeoMesh>
void incrementalPOD<Type, PatchField, GeoMesh>::addSnapshot(
    GeometricField<Type, PatchField, GeoMesh>& snapshot)
{
    Info << "********************************************************************"
         << endl;
    Info << "Adding a snapshot" << endl;
    Info << "Initial rank = " << rank << endl;
    Eigen::VectorXd snapshotEig = Foam2Eigen::field2Eigen(snapshot);
    double projectionError;
    double relProjectionError;

    if (rank == 0)
    {
        initialize(snapshot);
    }

    Eigen::VectorXd redCoeff = project(snapshot);
    GeometricField<Type, PatchField, GeoMesh> projectSnap = projectSnapshot(
                snapshot);
    GeometricField<Type, PatchField, GeoMesh> orthoSnap(snapshot - projectSnap);
    Eigen::VectorXd orthogonalSnap = Foam2Eigen::field2Eigen(orthoSnap);

    if (PODnorm == "L2")
    {
        projectionError = ITHACAutilities::L2Norm(orthoSnap);
        relProjectionError = projectionError / ITHACAutilities::L2Norm(snapshot);
    }
    else if (PODnorm == "Frobenius")
    {
        projectionError = snapshotEig.transpose() * snapshotEig -
                          (redCoeff.transpose() * redCoeff)(0, 0);
        projectionError = Foam::sqrt(std::abs(projectionError));
        relProjectionError = projectionError / (snapshotEig.transpose() * snapshotEig);
    }

    Info << "Relative projection error = " << relProjectionError << endl;
    Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(rank + 1, rank + 1);
    Q.topLeftCorner(rank, rank) = singularValues.asDiagonal();
    Q.topRightCorner(rank, 1) = redCoeff;
    Q(rank, rank) = projectionError;

    if (relProjectionError < tolleranceSVD)
    {
        Info << "Projection error is smaller than the tollerance" << endl;
        Q(rank, rank) = 0;
    }

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(Q,
                                          Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::MatrixXd newModes = svd.matrixU();
    Eigen::VectorXd newSingVal = svd.singularValues();

    if (relProjectionError < tolleranceSVD || rank >= snapshotEig.size())
    {
        Info << "Projection error is small or matrix is full rank." <<
             endl << "SVD rank held constant" << endl;
        this->EigenModes[0] = this->EigenModes[0] * newModes.topLeftCorner(rank, rank);
        singularValues = newSingVal.head(rank);
    }
    else
    {
        Info << "New snapshot is not in the span of the POD space." << endl
             << "SVD rank increases" << endl;
        orthogonalSnap = orthogonalSnap / projectionError;
        Eigen::MatrixXd temp(this->EigenModes[0].rows(),
                             this->EigenModes[0].cols() + 1);
        temp << this->EigenModes[0], orthogonalSnap;
        this->EigenModes[0] = temp * newModes;
        singularValues = newSingVal;
        rank++;
    }

    Info << "New POD rank = " << rank << endl;
    double orthogonalPar;

    if (PODnorm == "L2")
    {
        orthogonalPar = std::abs(
                            this->EigenModes[0].col(this->EigenModes[0].cols() - 1).transpose() *
                            massVector.asDiagonal() * this->EigenModes[0].col(0));
    }
    else if (PODnorm == "Frobenius")
    {
        orthogonalPar = std::abs(this->EigenModes[0].col(0).transpose() *
                                 this->EigenModes[0].col(this->EigenModes[0].cols() - 1));
    }

    double EPS = 2.2204e-16;
    Info << "Orthogonality = " << orthogonalPar << endl;

    if (orthogonalPar > std::min(tolleranceSVD, EPS * this->EigenModes[0].rows()))
    {
        Info << "Orthogonalization required" << endl;

        if (PODnorm == "L2")
        {
            ITHACAPOD::weightedGramSchmidt(this->EigenModes[0], massVector);
        }
        else if (PODnorm == "Frobenius")
        {
            Eigen::HouseholderQR<Eigen::MatrixXd> qr(this->EigenModes[0]);
            Eigen::MatrixXd thinQ(Eigen::MatrixXd::Identity(this->EigenModes[0].rows(),
                                  rank));
            this->EigenModes[0] = qr.householderQ() * thinQ;
        }
    }

    fillPtrList();
    Info << "********************************************************************"
         << endl;
}

template<class Type, template<class> class PatchField, class GeoMesh>
void incrementalPOD<Type, PatchField, GeoMesh>::fillPtrList()
{
    this->resize(rank);

    for (label i = 0; i < rank; i++)
    {
        GeometricField<Type, PatchField, GeoMesh>  tmp(this->first().name(),
                this->first());
        Eigen::VectorXd vec = this->EigenModes[0].col(i);
        tmp = Foam2Eigen::Eigen2field(tmp, vec);
        this->set(i, tmp.clone());
    }
}
template<class Type, template<class> class PatchField, class GeoMesh>
Eigen::VectorXd incrementalPOD<Type, PatchField, GeoMesh>::project(
    GeometricField<Type, PatchField, GeoMesh>& inputField,
    label numberOfModes)
{
    Eigen::VectorXd fieldEig = Foam2Eigen::field2Eigen(inputField);
    Eigen::VectorXd projField;

    if (numberOfModes == 0)
    {
        if (PODnorm == "L2")
        {
            projField = this->EigenModes[0].transpose() * fieldEig;
        }
        else if (PODnorm == "Frobenius")
        {
            projField = this->EigenModes[0].transpose() * fieldEig;
        }
    }
    else
    {
        M_Assert(numberOfModes <= this->EigenModes[0].cols(),
                 "Number of required modes for projection is higher then the number of available ones");

        if (PODnorm == "L2")
        {
            projField = (this->EigenModes[0].leftCols(numberOfModes)).transpose() *
                        massVector.asDiagonal() * fieldEig;
        }
        else if (PODnorm == "Frobenius")
        {
            projField = (this->EigenModes[0].leftCols(numberOfModes)).transpose() *
                        fieldEig;
        }
    }

    return projField;
}

template<class Type, template<class> class PatchField, class GeoMesh>
GeometricField<Type, PatchField, GeoMesh>
incrementalPOD<Type, PatchField, GeoMesh>::reconstruct(
    GeometricField<Type, PatchField, GeoMesh>& inputField,
    Eigen::MatrixXd Coeff,
    word Name)
{
    if (this->EigenModes.size() == 0)
    {
        this->toEigen();
    }

    label Nmodes = Coeff.rows();
    Eigen::VectorXd InField = this->EigenModes[0].leftCols(Nmodes) * Coeff;
    inputField = Foam2Eigen::Eigen2field(inputField, InField);
    inputField.rename(Name);
    Info << "WARNING: Boundary conditions are not reconstructed using incremental POD"
         << endl;
    return inputField;
}

template<class Type, template<class> class PatchField, class GeoMesh>
GeometricField<Type, PatchField, GeoMesh>
incrementalPOD<Type, PatchField, GeoMesh>::projectSnapshot(
    GeometricField<Type, PatchField, GeoMesh>& snapshot,
    label numberOfModes)
{
    Eigen::MatrixXd Modes;
    Eigen::MatrixXd M_vol;
    Eigen::MatrixXd M;
    Eigen::MatrixXd projSnapI;
    Eigen::MatrixXd projSnapCoeff;

    if (numberOfModes == 0)
    {
        Modes = this->EigenModes[0];
    }
    else
    {
        Modes = this->EigenModes[0].leftCols(numberOfModes);
    }

    GeometricField<Type, PatchField, GeoMesh> Fr = snapshot;
    Eigen::MatrixXd F_eigen = Foam2Eigen::field2Eigen(snapshot);

    if (PODnorm == "L2")
    {
        M = Modes.transpose() * massVector.asDiagonal() * Modes;
        projSnapI = Modes.transpose() * massVector.asDiagonal() * F_eigen;
    }
    else if (PODnorm == "Frobenius")
    {
        M = Modes.transpose() * Modes;
        projSnapI = Modes.transpose() * F_eigen;
    }
    else
    {
        std::cout << "Inner product not defined" << endl;
        exit(0);
    }

    projSnapCoeff = M.fullPivLu().solve(projSnapI);
    reconstruct(Fr, projSnapCoeff, "projSnap");
    return Fr;
}

template<class Type, template<class> class PatchField, class GeoMesh>
void incrementalPOD<Type, PatchField, GeoMesh>::projectSnapshots(
    PtrList<GeometricField<Type, PatchField, GeoMesh>> snapshots,
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& projSnapshots,
    label numberOfModes)
{
    M_Assert(numberOfModes <= this->size(),
             "The number of Modes used for the projection cannot be bigger than the number of available modes");
    projSnapshots.resize(snapshots.size());
    Eigen::MatrixXd Modes;

    if (numberOfModes == 0)
    {
        Modes = this->EigenModes[0];
    }
    else
    {
        Modes = this->EigenModes[0].leftCols(numberOfModes);
    }

    Eigen::MatrixXd M_vol;
    Eigen::MatrixXd M;
    Eigen::MatrixXd projSnapI;
    Eigen::MatrixXd projSnapCoeff;

    for (label i = 0; i < snapshots.size(); i++)
    {
        GeometricField<Type, PatchField, GeoMesh> Fr = snapshots[0];
        Eigen::MatrixXd F_eigen = Foam2Eigen::field2Eigen(snapshots[i]);

        if (PODnorm == "L2")
        {
            M = Modes.transpose() * massVector.asDiagonal() * Modes;
            projSnapI = Modes.transpose() * massVector.asDiagonal() * F_eigen;
        }
        else if (PODnorm == "Frobenius")
        {
            M = Modes.transpose() * Modes;
            projSnapI = Modes.transpose() * F_eigen;
        }
        else
        {
            std::cout << "Inner product not defined" << endl;
            exit(0);
        }

        projSnapCoeff = M.fullPivLu().solve(projSnapI);
        reconstruct(Fr, projSnapCoeff, "projSnap");
        projSnapshots.set(i, Fr.clone());
    }
}



template<class Type, template<class> class PatchField, class GeoMesh>
void incrementalPOD<Type, PatchField, GeoMesh>::writeModes()
{
    ITHACAstream::exportFields(this->toPtrList(),
                               outputFolder,
                               "base");
    ITHACAstream::exportMatrix(singularValues, "singularValues", "eigen",
                               outputFolder);
}


template class incrementalPOD<scalar, fvPatchField, volMesh>;
template class incrementalPOD<vector, fvPatchField, volMesh>;
