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
/// Source file of the Modes class.

#include "Modes.H"

template<class Type, template<class> class PatchField, class GeoMesh>
List<Eigen::MatrixXd> Modes<Type, PatchField, GeoMesh>::toEigen()
{
    NBC = 0;

    for (label i = 0; i < (this->first()).boundaryFieldRef().size(); i++)
    {
        if ((this->first()).boundaryFieldRef()[i].type() != "processor")
        {
            NBC++;
        }
    }

    EigenModes.resize(NBC + 1);
    EigenModes[0] = Foam2Eigen::PtrList2Eigen(this->toPtrList());
    List<Eigen::MatrixXd> BC = Foam2Eigen::PtrList2EigenBC(this->toPtrList());

    for (label i = 0; i < NBC; i++)
    {
        EigenModes[i + 1] = BC[i];
    }

    return EigenModes;
}

template<class Type, template<class> class PatchField, class GeoMesh>
List<Eigen::MatrixXd> Modes<Type, PatchField, GeoMesh>::project(
    fvMatrix<Type>& Af, label numberOfModes,
    word projType)
{
    M_Assert(projType == "G" || projType == "PG",
             "Projection type can be G for Galerking or PG for Petrov-Galerkin");
    List<Eigen::MatrixXd> LinSys;
    LinSys.resize(2);

    if (EigenModes.size() == 0)
    {
        toEigen();
    }

    Eigen::SparseMatrix<double> Ae;
    Eigen::VectorXd be;
    Foam2Eigen::fvMatrix2Eigen(Af, Ae, be);

    if (numberOfModes == 0)
    {
        if (projType == "G")
        {
            LinSys[0] = EigenModes[0].transpose() * Ae * EigenModes[0];
            LinSys[1] = EigenModes[0].transpose() * be;
        }

        if (projType == "PG")
        {
            LinSys[0] = (Ae * EigenModes[0]).transpose() * Ae * EigenModes[0];
            LinSys[1] = (Ae * EigenModes[0]).transpose() * be;
        }
    }
    else
    {
        M_Assert(numberOfModes <= EigenModes[0].cols(),
                 "Number of required modes for projection is higher then the number of available ones");

        if (projType == "G")
        {
            LinSys[0] = ((EigenModes[0]).leftCols(numberOfModes)).transpose() * Ae *
                        (EigenModes[0]).leftCols(numberOfModes);
            LinSys[1] = ((EigenModes[0]).leftCols(numberOfModes)).transpose() * be;
        }

        if (projType == "PG")
        {
            LinSys[0] = (Ae * ((EigenModes[0]).leftCols(numberOfModes))).transpose() * Ae *
                        (EigenModes[0]).leftCols(numberOfModes);
            LinSys[1] = (Ae * ((EigenModes[0]).leftCols(numberOfModes))).transpose() * be;
        }
    }

    return LinSys;
}

template<class Type, template<class> class PatchField, class GeoMesh>
Eigen::MatrixXd Modes<Type, PatchField, GeoMesh>::project(
    GeometricField<Type, PatchField, GeoMesh>&
    field, label numberOfModes, word projType, fvMatrix<Type>* Af)
{
    M_Assert(projType == "G" || projType == "PG",
             "Projection type can be G for Galerking or PG for Petrov-Galerkin");
    Eigen::MatrixXd fieldEig = Foam2Eigen::field2Eigen(field);
    auto vol = ITHACAutilities::getMassMatrixFV(field);
    Eigen::MatrixXd projField;

    if (EigenModes.size() == 0)
    {
        toEigen();
    }

    if (numberOfModes == 0)
    {
        if (projType == "G")
        {
            projField = EigenModes[0].transpose() * vol.asDiagonal() * fieldEig;
        }
        else if (projType == "PG")
        {
            M_Assert(Af != NULL,
                     "Using a Petrov-Galerkin projection you have to provide also the system matrix");
            Eigen::SparseMatrix<double> Ae;
            Eigen::VectorXd be;
            Foam2Eigen::fvMatrix2Eigen(*Af, Ae, be);
            projField = (Ae * EigenModes[0]).transpose() * vol.asDiagonal() * fieldEig;
        }
    }
    else
    {
        M_Assert(numberOfModes <= EigenModes[0].cols(),
                 "Number of required modes for projection is higher then the number of available ones");

        if (projType == "G")
        {
            projField = ((EigenModes[0]).leftCols(numberOfModes)).transpose() *
                        vol.asDiagonal() * fieldEig;
        }
        else if (projType == "PG")
        {
            M_Assert(Af != NULL,
                     "Using a Petrov-Galerkin projection you have to provide also the system matrix");
            Eigen::SparseMatrix<double> Ae;
            Eigen::VectorXd be;
            Foam2Eigen::fvMatrix2Eigen(*Af, Ae, be);
            projField = (Ae * ((EigenModes[0]).leftCols(numberOfModes))).transpose() *
                        vol.asDiagonal() * fieldEig;
        }
    }

    return projField;
}

template<class Type, template<class> class PatchField, class GeoMesh>
GeometricField<Type, PatchField, GeoMesh>
Modes<Type, PatchField, GeoMesh>::projectSnapshot(
    GeometricField<Type, PatchField, GeoMesh>&
    field, label numberOfModes, word projType, fvMatrix<Type>* Af)
{
    Eigen::MatrixXd proj = project(field, numberOfModes, projType, Af);
    GeometricField<Type, PatchField, GeoMesh> projSnap = field;
    reconstruct(projSnap, proj, projSnap.name());
    return projSnap;
}

template<class Type, template<class> class PatchField, class GeoMesh>
Eigen::MatrixXd Modes<Type, PatchField, GeoMesh>::project(
    PtrList<GeometricField<Type, PatchField, GeoMesh>>&
    fields,
    label numberOfModes, word projType, fvMatrix<Type>* Af)
{
    M_Assert(projType == "G" || projType == "PG",
             "Projection type can be G for Galerking or PG for Petrov-Galerkin");
    Eigen::MatrixXd fieldEig = Foam2Eigen::PtrList2Eigen(fields);
    auto vol = ITHACAutilities::getMassMatrixFV(fields[0]);
    Eigen::MatrixXd projField;

    if (EigenModes.size() == 0)
    {
        toEigen();
    }

    if (numberOfModes == 0)
    {
        if (projType == "G")
        {
            projField = EigenModes[0].transpose() * vol.asDiagonal() * fieldEig;
        }
        else if (projType == "PG")
        {
            M_Assert(Af != NULL,
                     "Using a Petrov-Galerkin projection you have to provide also the system matrix");
            Eigen::SparseMatrix<double> Ae;
            Eigen::VectorXd be;
            Foam2Eigen::fvMatrix2Eigen(*Af, Ae, be);
            projField = (Ae * EigenModes[0]).transpose() * vol.asDiagonal() * fieldEig;
        }
    }
    else
    {
        M_Assert(numberOfModes <= EigenModes[0].cols(),
                 "Number of required modes for projection is higher then the number of available ones");

        if (projType == "G")
        {
            projField = ((EigenModes[0]).leftCols(numberOfModes)).transpose() *
                        vol.asDiagonal() * fieldEig;
        }
        else if (projType == "PG")
        {
            M_Assert(Af != NULL,
                     "Using a Petrov-Galerkin projection you have to provide also the system matrix");
            Eigen::SparseMatrix<double> Ae;
            Eigen::VectorXd be;
            Foam2Eigen::fvMatrix2Eigen(*Af, Ae, be);
            projField = (Ae * ((EigenModes[0]).leftCols(numberOfModes))).transpose() *
                        vol.asDiagonal() * fieldEig;
        }
    }

    return projField;
}
template<class Type, template<class> class PatchField, class GeoMesh>
GeometricField<Type, PatchField, GeoMesh>
Modes<Type, PatchField, GeoMesh>::reconstruct(
    GeometricField<Type, PatchField, GeoMesh>& inputField,
    Eigen::MatrixXd Coeff,
    word Name)
{
    if (EigenModes.size() == 0)
    {
        toEigen();
    }

    label Nmodes = Coeff.rows();
    Eigen::VectorXd InField = EigenModes[0].leftCols(Nmodes) * Coeff;

    if (inputField.name() == "nut")
    {
        InField = (InField.array() < 0).select(0, InField);
    }

    inputField = Foam2Eigen::Eigen2field(inputField, InField);
    inputField.rename(Name);

    for (label i = 0; i < NBC; i++)
    {
        Eigen::VectorXd BF = EigenModes[i + 1].leftCols(Nmodes) * Coeff;

        if (inputField.name() == "nut")
        {
            BF = (BF.array() < 0).select(0, BF);
        }

        ITHACAutilities::assignBC(inputField, i, BF);
    }

    return inputField;
}
template<class Type, template<class> class PatchField, class GeoMesh>
PtrList<GeometricField<Type, PatchField, GeoMesh>>
        Modes<Type, PatchField, GeoMesh>::reconstruct(
            GeometricField<Type, PatchField, GeoMesh>& inputField,
            List < Eigen::MatrixXd> Coeff,
            word Name)
{
    PtrList<GeometricField<Type, PatchField, GeoMesh>> inputFields;
    inputFields.resize(0);

    for (label i = 0; i < Coeff.size(); i++)
    {
        inputField = reconstruct(inputField, Coeff[i], Name);
        inputFields.append(inputField.clone());
    }

    return inputFields;
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Modes<Type, PatchField, GeoMesh>::projectSnapshots(
    PtrList<GeometricField<Type, PatchField, GeoMesh>> snapshots,
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& projSnapshots,
    PtrList<volScalarField> Volumes,
    label numberOfModes,
    word innerProduct)
{
    if (EigenModes.size() == 0)
    {
        toEigen();
    }

    M_Assert(snapshots.size() == Volumes.size(),
             "The number of snapshots and the number of volumes vectors must be equal");
    M_Assert(numberOfModes <= this->size(),
             "The number of Modes used for the projection cannot be bigger than the number of available modes");
    M_Assert(innerProduct == "L2" || innerProduct == "Frobenius",
             "The chosen inner product is not implemented");
    projSnapshots.resize(snapshots.size());
    label dim = std::nearbyint(EigenModes[0].rows() /
                               Volumes[0].size()); //Checking if volumes and modes have the same size that means check if the problem is vector or scalar
    Eigen::MatrixXd totVolumes(Volumes[0].size()*dim, Volumes.size());

    for (label i = 0; i < Volumes.size(); i++)
    {
        totVolumes.col(i) = Foam2Eigen::field2Eigen(Volumes[i]);
    }

    totVolumes.replicate(dim, 1);
    Eigen::MatrixXd Modes;

    if (numberOfModes == 0)
    {
        Modes = EigenModes[0];
    }
    else
    {
        Modes = EigenModes[0].leftCols(numberOfModes);
    }

    Eigen::MatrixXd M;
    Eigen::MatrixXd projSnapI;
    Eigen::MatrixXd projSnapCoeff;

    for (label i = 0; i < snapshots.size(); i++)
    {
        GeometricField<Type, PatchField, GeoMesh> Fr = snapshots[0];
        Eigen::MatrixXd F_eigen = Foam2Eigen::field2Eigen(snapshots[i]);

        if (innerProduct == "L2")
        {
            M = Modes.transpose() * (totVolumes.col(i)).asDiagonal() * Modes;
            projSnapI = Modes.transpose() * (totVolumes.col(i)).asDiagonal() * F_eigen;
        }
        else //Frobenius
        {
            M = Modes.transpose() * Modes;
            projSnapI = Modes.transpose() * F_eigen;
        }

        projSnapCoeff = M.fullPivLu().solve(projSnapI);
        reconstruct(Fr, projSnapCoeff, "projSnap");
        projSnapshots.set(i, Fr.clone());
    }
}
template<class Type, template<class> class PatchField, class GeoMesh>
void Modes<Type, PatchField, GeoMesh>::projectSnapshots(
    PtrList<GeometricField<Type, PatchField, GeoMesh>> snapshots,
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& projSnapshots,
    PtrList<volScalarField> Volumes, word innerProduct)
{
    label numberOfModes = 0;
    projectSnapshots(snapshots, projSnapshots, Volumes, numberOfModes,
                     innerProduct);
}

template<class Type, template<class> class PatchField, class GeoMesh>
void Modes<Type, PatchField, GeoMesh>::projectSnapshots(
    PtrList<GeometricField<Type, PatchField, GeoMesh>> snapshots,
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& projSnapshots,
    label numberOfModes,
    word innerProduct)
{
    if (EigenModes.size() == 0)
    {
        toEigen();
    }

    M_Assert(numberOfModes <= this->size(),
             "The number of Modes used for the projection cannot be bigger than the number of available modes");
    M_Assert(innerProduct == "L2" || innerProduct == "Frobenius",
             "The chosen inner product is not implemented");
    projSnapshots.resize(snapshots.size());
    Eigen::MatrixXd Modes;

    if (numberOfModes == 0)
    {
        Modes = EigenModes[0];
    }
    else
    {
        Modes = EigenModes[0].leftCols(numberOfModes);
    }

    Eigen::MatrixXd M_vol;
    Eigen::MatrixXd M;
    Eigen::MatrixXd projSnapI;
    Eigen::MatrixXd projSnapCoeff;

    for (label i = 0; i < snapshots.size(); i++)
    {
        GeometricField<Type, PatchField, GeoMesh> Fr = snapshots[0];
        Eigen::MatrixXd F_eigen = Foam2Eigen::field2Eigen(snapshots[i]);

        if (innerProduct == "L2")
        {
            M_vol = ITHACAutilities::getMassMatrixFV(snapshots[i]);
        }
        else if (innerProduct == "Frobenius")
        {
            M_vol =  Eigen::VectorXd::Identity(F_eigen.rows(), 1);
        }
        else
        {
            std::cout << "Inner product not defined" << endl;
            exit(0);
        }

        M = Modes.transpose() * M_vol.asDiagonal() * Modes;
        projSnapI = Modes.transpose() * M_vol.asDiagonal() * F_eigen;
        projSnapCoeff = M.fullPivLu().solve(projSnapI);
        reconstruct(Fr, projSnapCoeff, "projSnap");
        projSnapshots.set(i, Fr.clone());
    }
}
template<class Type, template<class> class PatchField, class GeoMesh>
void Modes<Type, PatchField, GeoMesh>::projectSnapshots(
    PtrList<GeometricField<Type, PatchField, GeoMesh>> snapshots,
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& projSnapshots,
    word innerProduct)
{
    label numberOfModes = 0;
    projectSnapshots(snapshots, projSnapshots, numberOfModes, innerProduct);
}

template<class Type, template<class> class PatchField, class GeoMesh>
void Modes<Type, PatchField, GeoMesh>::operator=(const
        PtrList<GeometricField<Type, PatchField, GeoMesh>>& modes)
{
    this->resize(modes.size());

    for (label i = 0; i < modes.size(); i++)
    {
        (*this).set(i, modes[i].clone());
    }
}


template class Modes<scalar, fvPatchField, volMesh>;
template class Modes<vector, fvPatchField, volMesh>;
template class Modes<scalar, fvsPatchField, surfaceMesh>;
