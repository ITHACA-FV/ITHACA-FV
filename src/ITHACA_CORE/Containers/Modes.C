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
Eigen::MatrixXd Modes<Type, PatchField, GeoMesh>::project(
    GeometricField<Type, PatchField, GeoMesh>&
    field, label numberOfModes, word projType, fvMatrix<Type>* Af)
{
    M_Assert(projType == "F" || projType == "G" || projType == "PG",
             "Projection type can be F for Frobenius, G for Galerkin or PG for Petrov-Galerkin");
    Eigen::MatrixXd fieldEig = Foam2Eigen::field2Eigen(field);
    auto vol = ITHACAutilities::getMassMatrixFV(field);
    Eigen::MatrixXd projField;
    Eigen::MatrixXd M;
    Eigen::MatrixXd b;

    if (EigenModes.size() == 0)
    {
        toEigen();
    }

    if (numberOfModes == 0)
    {
        if (projType == "F")
        {
            vol = Eigen::VectorXd::Ones(vol.size());
            b = EigenModes[0].transpose() * vol.asDiagonal() * fieldEig;
            M = EigenModes[0].transpose() * EigenModes[0].transpose();
            projField = M.fullPivLu().solve(b);
        }
        else if (projType == "G")
        {
            projField = EigenModes[0].transpose() * vol.asDiagonal() * fieldEig;
        }
        else if (projType == "PG")
        {
            M_Assert(Af != NULL,
                     "Using a Petrov-Galerkin projection you have to provide also the system matrix");
            Eigen::SparseMatrix<double> Ae;
            Eigen::VectorXd be;
            Foam2Eigen::fvMatrix2Eigen(* Af, Ae, be);
            projField = (Ae * EigenModes[0]).transpose() * vol.asDiagonal() * fieldEig;
        }
    }
    else
    {
        M_Assert(numberOfModes <= EigenModes[0].cols(),
                 "Number of required modes for projection is higher then the number of available ones");

        if (projType == "F")
        {
            vol = Eigen::VectorXd::Ones(vol.size());
            M = EigenModes[0].leftCols(numberOfModes).transpose() * EigenModes[0].leftCols(
                    numberOfModes);
            b = ((EigenModes[0]).leftCols(numberOfModes)).transpose() *
                vol.asDiagonal() * fieldEig;
            projField = M.fullPivLu().solve(b);
        }
        else if (projType == "G")
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
            Foam2Eigen::fvMatrix2Eigen(* Af, Ae, be);
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
    PtrList<GeometricField<Type, PatchField, GeoMesh >>&
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
        if (projType == "F")
        {
            vol = Eigen::VectorXd::Ones(vol.size());
            projField = EigenModes[0].transpose() * vol.asDiagonal() * fieldEig;
        }
        else if (projType == "G")
        {
            projField = EigenModes[0].transpose() * vol.asDiagonal() * fieldEig;
        }
        else if (projType == "PG")
        {
            M_Assert(Af != NULL,
                     "Using a Petrov-Galerkin projection you have to provide also the system matrix");
            Eigen::SparseMatrix<double> Ae;
            Eigen::VectorXd be;
            Foam2Eigen::fvMatrix2Eigen(* Af, Ae, be);
            projField = (Ae * EigenModes[0]).transpose() * vol.asDiagonal() * fieldEig;
        }
    }
    else
    {
        M_Assert(numberOfModes <= EigenModes[0].cols(),
                 "Number of required modes for projection is higher then the number of available ones");

        if (projType == "F")
        {
            vol = Eigen::VectorXd::Ones(vol.size());
            projField = ((EigenModes[0]).leftCols(numberOfModes)).transpose() *
                        vol.asDiagonal() * fieldEig;
        }
        else if (projType == "G")
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
            Foam2Eigen::fvMatrix2Eigen(* Af, Ae, be);
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
PtrList<GeometricField<Type, PatchField, GeoMesh >>
Modes<Type, PatchField, GeoMesh>::reconstruct(
    GeometricField<Type, PatchField, GeoMesh>& inputField,
    List < Eigen::MatrixXd> Coeff,
    word Name)
{
    PtrList<GeometricField<Type, PatchField, GeoMesh >> inputFields;
    inputFields.resize(0);

    for (label i = 0; i < Coeff.size(); i++)
    {
        inputField = reconstruct(inputField, Coeff[i], Name);
        inputFields.append(inputField.clone());
    }

    return inputFields;
}


template<class Type, template<class> class PatchField, class GeoMesh >
void Modes<Type, PatchField, GeoMesh>::projectSnapshots(
    PtrList<GeometricField<Type, PatchField, GeoMesh >> snapshots,
    PtrList<GeometricField<Type, PatchField, GeoMesh >>& projSnapshots,
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
    Eigen::MatrixXd totVolumes(Volumes[0].size() * dim, Volumes.size());

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
    PtrList<GeometricField<Type, PatchField, GeoMesh >> snapshots,
    PtrList<GeometricField<Type, PatchField, GeoMesh >>& projSnapshots,
    PtrList<volScalarField> Volumes, word innerProduct)
{
    label numberOfModes = 0;
    projectSnapshots(snapshots, projSnapshots, Volumes, numberOfModes,
                     innerProduct);
}

template<class Type, template<class> class PatchField, class GeoMesh>
void Modes<Type, PatchField, GeoMesh>::projectSnapshots(
    PtrList<GeometricField<Type, PatchField, GeoMesh >> snapshots,
    PtrList<GeometricField<Type, PatchField, GeoMesh >>& projSnapshots,
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
            Foam::Info << "Inner product not defined" << Foam::endl;
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
    PtrList<GeometricField<Type, PatchField, GeoMesh >> snapshots,
    PtrList<GeometricField<Type, PatchField, GeoMesh >>& projSnapshots,
    word innerProduct)
{
    label numberOfModes = 0;
    projectSnapshots(snapshots, projSnapshots, numberOfModes, innerProduct);
}

template<class Type, template<class> class PatchField, class GeoMesh>
void Modes<Type, PatchField, GeoMesh>::operator=(const
    PtrList<GeometricField<Type, PatchField, GeoMesh >>& modes)
{
    this->resize(modes.size());

    for (label i = 0; i < modes.size(); i++)
    {
        (* this).set(i, modes[i].clone());
    }
}

template<class Type, template<class> class PatchField, class GeoMesh>
List<Eigen::MatrixXd>
Modes<Type, PatchField, GeoMesh>::project
(
    fvMatrix<Type>& Af,
    label numberOfModes,
    word projType
)
{
    FatalErrorInFunction
        << "fvMatrix projection is currently implemented only for volScalarField and volVectorField modes"
        << exit(FatalError);

    return List<Eigen::MatrixXd>();
}

template<>
List<Eigen::MatrixXd>
Modes<scalar, fvPatchField, volMesh>::project
(
    fvMatrix<scalar>& Af,
    label numberOfModes,
    word projType
)
{
    M_Assert
    (
        projType == "G" || projType == "PG",
        "Projection type can be G for Galerkin or PG for Petrov-Galerkin"
    );

    if (numberOfModes == 0)
    {
        numberOfModes = this->size();
    }

    M_Assert
    (
        numberOfModes <= this->size(),
        "Number of required modes is larger than number of available modes"
    );

    const label nModes = numberOfModes;

    List<Eigen::MatrixXd> LinSys(2);
    LinSys[0] = Eigen::MatrixXd::Zero(nModes, nModes);
    LinSys[1] = Eigen::MatrixXd::Zero(nModes, 1);

    class projectableMatrix
    :
        public fvScalarMatrix
    {
    public:

        projectableMatrix(const fvScalarMatrix& A)
        :
            fvScalarMatrix(A)
        {}

        void prepare()
        {
            addBoundaryDiag(diag(), 0);
        }

        void prepareSource(scalarField& source)
        {
            addBoundarySource(source, false);
        }
    };

    projectableMatrix Ap(Af);
    Ap.prepare();

    const lduInterfaceFieldPtrsList interfaces
    (
        Af.psi().boundaryField().scalarInterfaces()
    );

    // ------------------------------------------------------------
    // Compute A*phi_j for every mode once.
    // ------------------------------------------------------------

    List<scalarField> Amodes(nModes);

    for (label j = 0; j < nModes; ++j)
    {
        const scalarField& x =
            (*this)[j].primitiveField();

        Amodes[j].setSize(x.size());
        Amodes[j] = 0.0;

        Ap.Amul
        (
            Amodes[j],
            x,
            Ap.boundaryCoeffs(),
            interfaces,
            0
        );
    }

    // ------------------------------------------------------------
    // Reduced matrix
    //
    // G :  V^T A V
    // PG: (A V)^T A V
    // ------------------------------------------------------------

    for (label i = 0; i < nModes; ++i)
    {
        for (label j = 0; j < nModes; ++j)
        {
            scalar value = 0.0;

            if (projType == "G")
            {
                const scalarField& modeI =
                    (*this)[i].primitiveField();

                forAll(Amodes[j], celli)
                {
                    value +=
                        modeI[celli]
                       *Amodes[j][celli];
                }
            }
            else
            {
                forAll(Amodes[j], celli)
                {
                    value +=
                        Amodes[i][celli]
                       *Amodes[j][celli];
                }
            }

            reduce(value, sumOp<scalar>());

            LinSys[0](i,j) = value;
        }
    }

    // ------------------------------------------------------------
    // Reduced source
    //
    // G :  V^T b
    // PG: (A V)^T b
    // ------------------------------------------------------------

    scalarField rhs(Af.source());
    Ap.prepareSource(rhs);

    for (label i = 0; i < nModes; ++i)
    {
        scalar value = 0.0;

        if (projType == "G")
        {
            const scalarField& modeI =
                (*this)[i].primitiveField();

            forAll(rhs, celli)
            {
                value +=
                    modeI[celli]
                   *rhs[celli];
            }
        }
        else
        {
            forAll(rhs, celli)
            {
                value +=
                    Amodes[i][celli]
                   *rhs[celli];
            }
        }

        reduce(value, sumOp<scalar>());

        LinSys[1](i,0) = value;
    }

    return LinSys;
}

template<>
List<Eigen::MatrixXd>
Modes<vector, fvPatchField, volMesh>::project
(
    fvMatrix<vector>& Af,
    label numberOfModes,
    word projType
)
{
    M_Assert
    (
        projType == "G" || projType == "PG",
        "Projection type can be G for Galerkin or PG for Petrov-Galerkin"
    );

    if (numberOfModes == 0)
    {
        numberOfModes = this->size();
    }

    M_Assert
    (
        numberOfModes <= this->size(),
        "Number of required modes is larger than number of available modes"
    );

    const label nModes = numberOfModes;

    List<Eigen::MatrixXd> LinSys(2);
    LinSys[0] = Eigen::MatrixXd::Zero(nModes, nModes);
    LinSys[1] = Eigen::MatrixXd::Zero(nModes, 1);

    class projectableMatrix
    :
        public fvVectorMatrix
    {
    public:

        projectableMatrix(const fvVectorMatrix& A)
        :
            fvVectorMatrix(A)
        {}

        void prepareComponent
        (
            const scalarField& originalDiag,
            const direction cmpt
        )
        {
            diag() = originalDiag;
            addBoundaryDiag(diag(), cmpt);
        }

        void prepareSource(vectorField& source)
        {
            addBoundarySource(source, false);
        }
    };

    projectableMatrix Ap(Af);
    const scalarField originalDiag(Ap.diag());

    const lduInterfaceFieldPtrsList interfaces
    (
        Af.psi().boundaryField().scalarInterfaces()
    );

    // ------------------------------------------------------------
    // Compute A*phi_j for every vector mode once.
    // ------------------------------------------------------------

    List<vectorField> Amodes(nModes);

    for (label j = 0; j < nModes; ++j)
    {
        const vectorField& modeJ =
            (*this)[j].primitiveField();

        Amodes[j].setSize(modeJ.size());
        Amodes[j] = vector::zero;

        for
        (
            direction cmpt = 0;
            cmpt < vector::nComponents;
            ++cmpt
        )
        {
            Ap.prepareComponent(originalDiag, cmpt);

            scalarField x
            (
                modeJ.component(cmpt)
            );

            scalarField Ax
            (
                x.size(),
                0.0
            );

            FieldField<Field, scalar> bouCoeffsCmpt
            (
                Ap.boundaryCoeffs().component(cmpt)
            );

            Ap.Amul
            (
                Ax,
                x,
                bouCoeffsCmpt,
                interfaces,
                cmpt
            );

            forAll(Ax, celli)
            {
                Amodes[j][celli][cmpt] = Ax[celli];
            }
        }
    }

    // ------------------------------------------------------------
    // Reduced matrix
    //
    // G :  V^T A V
    // PG: (A V)^T A V
    // ------------------------------------------------------------

    for (label i = 0; i < nModes; ++i)
    {
        for (label j = 0; j < nModes; ++j)
        {
            scalar value = 0.0;

            if (projType == "G")
            {
                const vectorField& modeI =
                    (*this)[i].primitiveField();

                forAll(Amodes[j], celli)
                {
                    value +=
                        modeI[celli]
                      & Amodes[j][celli];
                }
            }
            else
            {
                forAll(Amodes[j], celli)
                {
                    value +=
                        Amodes[i][celli]
                      & Amodes[j][celli];
                }
            }

            reduce(value, sumOp<scalar>());

            LinSys[0](i,j) = value;
        }
    }

    // ------------------------------------------------------------
    // Reduced source
    //
    // G :  V^T b
    // PG: (A V)^T b
    // ------------------------------------------------------------

    vectorField rhs(Af.source());
    Ap.prepareSource(rhs);

    for (label i = 0; i < nModes; ++i)
    {
        scalar value = 0.0;

        if (projType == "G")
        {
            const vectorField& modeI =
                (*this)[i].primitiveField();

            forAll(rhs, celli)
            {
                value +=
                    modeI[celli]
                  & rhs[celli];
            }
        }
        else
        {
            forAll(rhs, celli)
            {
                value +=
                    Amodes[i][celli]
                  & rhs[celli];
            }
        }

        reduce(value, sumOp<scalar>());

        LinSys[1](i,0) = value;
    }

    Ap.diag() = originalDiag;

    return LinSys;
}

template class Modes<scalar, fvPatchField, volMesh>;
template class Modes<vector, fvPatchField, volMesh>;
template class Modes<scalar, fvsPatchField, surfaceMesh>;
template class Modes<vector, pointPatchField, pointMesh>;
