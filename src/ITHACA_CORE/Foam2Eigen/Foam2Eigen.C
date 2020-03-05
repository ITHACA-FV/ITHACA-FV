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

#include "Foam2Eigen.H"

/// \file
/// Source file of the foam2eigen class.

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<>
Eigen::VectorXd Foam2Eigen::field2Eigen(
    GeometricField<vector, fvPatchField, volMesh>& field)
{
    Eigen::VectorXd out;
    out.resize(int(field.size() * 3));

    for (int l = 0; l < field.size(); l++)
    {
        out(l) = field[l][0];
        out(field.size() +  l ) = field[l][1];
        out(2 * field.size() + l ) = field[l][2];
    }

    return out;
}

template<>
Eigen::VectorXd Foam2Eigen::field2Eigen(
    GeometricField<scalar, fvPatchField, volMesh>& field)
{
    Eigen::VectorXd out;
    out.resize(int(field.size()));

    for (int l = 0; l < field.size(); l++)
    {
        out(l) = field[l];
    }

    return out;
}

template<>
Eigen::VectorXd Foam2Eigen::field2Eigen(const Field<scalar>& field)
{
    Eigen::VectorXd out;
    out.resize(int(field.size()));

    for (int l = 0; l < field.size(); l++)
    {
        out(l) = field[l];
    }

    return out;
}

template<>
Eigen::VectorXd Foam2Eigen::field2Eigen(const Field<vector>& field)
{
    Eigen::VectorXd out;
    out.resize(int(field.size() * 3));

    for (int l = 0; l < field.size(); l++)
    {
        out(l) = field[l][0];
        out(field.size() +  l ) = field[l][1];
        out(2 * field.size() + l ) = field[l][2];
    }

    return out;
}


template<>
List<Eigen::VectorXd> Foam2Eigen::field2EigenBC(
    GeometricField<vector, fvPatchField, volMesh>& field)
{
    List<Eigen::VectorXd> Out;
    unsigned int size = field.boundaryField().size();
    Out.resize(size);

    for (unsigned int i = 0; i < size; i++ )
    {
        unsigned int sizei = field.boundaryField()[i].size();
        Out[i].resize(sizei * 3);

        for (unsigned int k = 0; k < sizei ; k++)
        {
            Out[i](k) = field.boundaryField()[i][k][0];
            Out[i](k + sizei) = field.boundaryField()[i][k][1];
            Out[i](k + 2 * sizei) = field.boundaryField()[i][k][2];
        }
    }

    return Out;
}

template<>
List<Eigen::VectorXd> Foam2Eigen::field2EigenBC(
    GeometricField<scalar, fvPatchField, volMesh>& field)
{
    List<Eigen::VectorXd> Out;
    unsigned int size = field.boundaryField().size();
    Out.resize(size);

    for (unsigned int i = 0; i < size; i++ )
    {
        unsigned int sizei = field.boundaryField()[i].size();
        Out[i].resize(sizei);

        for (unsigned int k = 0; k < sizei ; k++)
        {
            Out[i](k) = field.boundaryField()[i][k];
        }
    }

    return Out;
}

template<>
List<Eigen::MatrixXd> Foam2Eigen::PtrList2EigenBC(
    PtrList<GeometricField<scalar, fvPatchField, volMesh>>&
    fields, int Nfields)
{
    unsigned int Nf;
    M_Assert(Nfields <= fields.size(),
             "The Number of requested fields cannot be bigger than the number of requested entries.");

    if (Nfields == -1)
    {
        Nf = fields.size();
    }
    else
    {
        Nf = Nfields;
    }

    List<Eigen::MatrixXd> Out;
    unsigned int NBound = fields[0].boundaryField().size();
    Out.resize(NBound);

    for (unsigned int i = 0; i < NBound; i++)
    {
        int sizei = fields[0].boundaryField()[i].size();
        Out[i].resize(sizei, Nf);
    }

    for (unsigned int k = 0; k < Nf; k++)
    {
        List<Eigen::VectorXd> temp;
        temp = field2EigenBC(fields[k]);

        for (unsigned int i = 0; i < NBound; i++)
        {
            Out[i].col(k) = temp[i];
        }
    }

    return Out;
}

template<>
List<Eigen::MatrixXd> Foam2Eigen::PtrList2EigenBC(
    PtrList<GeometricField<vector, fvPatchField, volMesh>>&
    fields, int Nfields)
{
    unsigned int Nf;
    M_Assert(Nfields <= fields.size(),
             "The Number of requested fields cannot be bigger than the number of requested entries.");

    if (Nfields == -1)
    {
        Nf = fields.size();
    }
    else
    {
        Nf = Nfields;
    }

    List<Eigen::MatrixXd> Out;
    unsigned int NBound = fields[0].boundaryField().size();
    Out.resize(NBound);

    for (unsigned int i = 0; i < NBound; i++)
    {
        int sizei = fields[0].boundaryField()[i].size();
        Out[i].resize(sizei * 3, Nf);
    }

    for (unsigned int k = 0; k < Nf; k++)
    {
        List<Eigen::VectorXd> temp;
        temp = field2EigenBC(fields[k]);

        for (unsigned int i = 0; i < NBound; i++)
        {
            Out[i].col(k) = temp[i];
        }
    }

    return Out;
}


template<>
GeometricField<vector, fvPatchField, volMesh> Foam2Eigen::Eigen2field(
    GeometricField<vector, fvPatchField, volMesh>& field_in,
    Eigen::VectorXd& eigen_vector)
{
    GeometricField<vector, fvPatchField, volMesh> field_out(field_in);

    for (auto i = 0; i < field_out.size(); i++)
    {
        field_out.ref()[i][0] = eigen_vector(i);
        field_out.ref()[i][1] = eigen_vector(i + field_out.size());
        field_out.ref()[i][2] = eigen_vector(i + field_out.size() * 2);
    }

    field_out.correctBoundaryConditions();
    return field_out;
}

template<>
GeometricField<scalar, fvPatchField, volMesh> Foam2Eigen::Eigen2field(
    GeometricField<scalar, fvPatchField, volMesh>& field_in,
    Eigen::VectorXd& eigen_vector)
{
    GeometricField<scalar, fvPatchField, volMesh> field_out(field_in);

    // std::memcpy (&field_out.ref()[0], eigen_vector.data(), field_out.ref().size()*sizeof(double));

    for (auto i = 0; i < field_out.size(); i++)
    {
        field_out.ref()[i] = eigen_vector(i);
    }

    field_out.correctBoundaryConditions();
    return field_out;
}


template<>
Eigen::MatrixXd Foam2Eigen::PtrList2Eigen(
    PtrList<GeometricField<vector, fvPatchField, volMesh>>& fields,
    int Nfields)
{
    int Nf;
    M_Assert(Nfields <= fields.size(),
             "The Number of requested fields cannot be bigger than the number of requested entries.");

    if (Nfields == -1)
    {
        Nf = fields.size();
    }
    else
    {
        Nf = Nfields;
    }

    Eigen::MatrixXd out;
    out.resize(int(fields[0].size() * 3), Nf);

    for (int k = 0; k < Nf; k++)
    {
        out.col(k) = field2Eigen<vector>(fields[k]);
    }

    return out;
}

template <class type_f>
Eigen::MatrixXd Foam2Eigen::PtrList2Eigen(PtrList<Field<type_f>>& fields,
        int Nfields)
{
    int Nf;
    M_Assert(Nfields <= fields.size(),
             "The Number of requested fields cannot be bigger than the number of requested entries.");

    if (Nfields == -1)
    {
        Nf = fields.size();
    }
    else
    {
        Nf = Nfields;
    }

    Eigen::MatrixXd out;
    int nrows = (field2Eigen(fields[0])).rows();
    out.resize(nrows, Nf);

    for (int k = 0; k < Nf; k++)
    {
        out.col(k) = field2Eigen(fields[k]);
    }

    return out;
}

template Eigen::MatrixXd Foam2Eigen::PtrList2Eigen(PtrList<Field<scalar>>&
        fields, int Nfields);
template Eigen::MatrixXd Foam2Eigen::PtrList2Eigen(PtrList<Field<vector>>&
        fields, int Nfields);

template<>
Eigen::MatrixXd Foam2Eigen::PtrList2Eigen(
    PtrList<GeometricField<scalar, fvPatchField, volMesh>>& fields,
    int Nfields)
{
    int Nf;
    M_Assert(Nfields <= fields.size(),
             "The Number of requested fields cannot be bigger than the number of requested entries.");

    if (Nfields == -1)
    {
        Nf = fields.size();
    }
    else
    {
        Nf = Nfields;
    }

    Eigen::MatrixXd out;
    out.resize(int(fields[0].size()), Nf);

    for (int k = 0; k < Nf; k++)
    {
        out.col(k) = field2Eigen<scalar>(fields[k]);
    }

    return out;
}


template<>
void Foam2Eigen::fvMatrix2Eigen(const fvMatrix<scalar>& foam_matrix,
                                Eigen::MatrixXd& A,
                                Eigen::VectorXd& b)
{
    int sizeA = foam_matrix.diag().size();
    A.setZero(sizeA, sizeA);
    b.setZero(sizeA);

    for (auto i = 0; i < sizeA; i++)
    {
        A(i, i) = foam_matrix.diag()[i];
        b(i, 0) = foam_matrix.source()[i];
    }

    const lduAddressing& addr = foam_matrix.lduAddr();
    const labelList& lowerAddr = addr.lowerAddr();
    const labelList& upperAddr = addr.upperAddr();
    forAll(lowerAddr, i)
    {
        A(lowerAddr[i], upperAddr[i]) = foam_matrix.upper()[i];
        A(upperAddr[i], lowerAddr[i]) = foam_matrix.lower()[i];
    }
    forAll(foam_matrix.psi().boundaryField(), I)
    {
        const fvPatch& ptch = foam_matrix.psi().boundaryField()[I].patch();
        forAll(ptch, J)
        {
            int w = ptch.faceCells()[J];
            A(w, w) += foam_matrix.internalCoeffs()[I][J];
            b(w, 0)   += foam_matrix.boundaryCoeffs()[I][J];
        }
    }
}

template<>
void Foam2Eigen::fvMatrix2Eigen(const fvMatrix<vector>& foam_matrix,
                                Eigen::MatrixXd& A,
                                Eigen::VectorXd& b)
{
    int sizeA = foam_matrix.diag().size();
    A.resize(sizeA * 3, sizeA * 3);
    b.resize(sizeA * 3);

    for (auto i = 0; i < sizeA; i++)
    {
        A(i, i) = foam_matrix.diag()[i];
        A(sizeA + i, sizeA + i) = foam_matrix.diag()[i];
        A(2 * sizeA + i, 2 * sizeA + i) = foam_matrix.diag()[i];
        b(i) = foam_matrix.source()[i][0];
        b(sizeA + i) = foam_matrix.source()[i][1];
        b(2 * sizeA + i) = foam_matrix.source()[i][2];
    }

    const lduAddressing& addr = foam_matrix.lduAddr();
    const labelList& lowerAddr = addr.lowerAddr();
    const labelList& upperAddr = addr.upperAddr();
    forAll(lowerAddr, i)
    {
        A(lowerAddr[i], upperAddr[i]) = foam_matrix.upper()[i];
        A(lowerAddr[i] + sizeA, upperAddr[i] + sizeA) = foam_matrix.upper()[i];
        A(lowerAddr[i] + sizeA * 2, upperAddr[i] + sizeA * 2) = foam_matrix.upper()[i];
        A(upperAddr[i], lowerAddr[i]) = foam_matrix.lower()[i];
        A(upperAddr[i] + sizeA, lowerAddr[i] + sizeA) = foam_matrix.lower()[i];
        A(upperAddr[i] + sizeA * 2, lowerAddr[i] + sizeA * 2) = foam_matrix.lower()[i];
    }
    forAll(foam_matrix.psi().boundaryField(), I)
    {
        const fvPatch& ptch = foam_matrix.psi().boundaryField()[I].patch();
        forAll(ptch, J)
        {
            int w = ptch.faceCells()[J];
            A(w, w) += foam_matrix.internalCoeffs()[I][J][0];
            A(w + sizeA, w + sizeA) += foam_matrix.internalCoeffs()[I][J][1];
            A(w + sizeA * 2, w + sizeA * 2) += foam_matrix.internalCoeffs()[I][J][2];
            b(w)   += foam_matrix.boundaryCoeffs()[I][J][0];
            b(w + sizeA)   += foam_matrix.boundaryCoeffs()[I][J][1];
            b(w + sizeA * 2)   += foam_matrix.boundaryCoeffs()[I][J][2];
        }
    }
}

template<>
void Foam2Eigen::fvMatrix2Eigen(const fvMatrix<scalar>& foam_matrix,
                                Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b)
{
    int sizeA = foam_matrix.diag().size();
    int nel = foam_matrix.diag().size() + foam_matrix.upper().size() +
              foam_matrix.lower().size();
    A.resize(sizeA, sizeA);
    b.resize(sizeA);
    A.reserve(nel);
    typedef Eigen::Triplet<double> Trip;
    std::vector<Trip> tripletList;
    tripletList.reserve(nel);

    for (auto i = 0; i < sizeA; i++)
    {
        tripletList.push_back(Trip(i, i, foam_matrix.diag()[i]));
        b(i, 0) = foam_matrix.source()[i];
    }

    const lduAddressing& addr = foam_matrix.lduAddr();
    const labelList& lowerAddr = addr.lowerAddr();
    const labelList& upperAddr = addr.upperAddr();
    forAll(lowerAddr, i)
    {
        tripletList.push_back(Trip(lowerAddr[i], upperAddr[i], foam_matrix.upper()[i]));
        tripletList.push_back(Trip(upperAddr[i], lowerAddr[i], foam_matrix.lower()[i]));
    }
    forAll(foam_matrix.psi().boundaryField(), I)
    {
        const fvPatch& ptch = foam_matrix.psi().boundaryField()[I].patch();
        forAll(ptch, J)
        {
            int w = ptch.faceCells()[J];
            tripletList.push_back(Trip(w, w, foam_matrix.internalCoeffs()[I][J]));
            b(w, 0)   += foam_matrix.boundaryCoeffs()[I][J];
        }
    }
    A.setFromTriplets(tripletList.begin(), tripletList.end());
}

template<>
void Foam2Eigen::fvMatrix2Eigen(const fvMatrix<vector>& foam_matrix,
                                Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b)
{
    int sizeA = foam_matrix.diag().size();
    int nel = foam_matrix.diag().size() + foam_matrix.upper().size() +
              foam_matrix.lower().size();
    A.resize(sizeA * 3, sizeA * 3);
    A.reserve(nel * 3);
    b.resize(sizeA * 3);
    typedef Eigen::Triplet<double> Trip;
    std::vector<Trip> tripletList;
    tripletList.reserve(nel * 3);

    for (auto i = 0; i < sizeA; i++)
    {
        tripletList.push_back(Trip(i, i, foam_matrix.diag()[i]));
        tripletList.push_back(Trip(sizeA + i, sizeA + i, foam_matrix.diag()[i]));
        tripletList.push_back(Trip(2 * sizeA + i, 2 * sizeA + i,
                                   foam_matrix.diag()[i]));
        b(i) = foam_matrix.source()[i][0];
        b(sizeA + i) = foam_matrix.source()[i][1];
        b(2 * sizeA + i) = foam_matrix.source()[i][2];
    }

    const lduAddressing& addr = foam_matrix.lduAddr();
    const labelList& lowerAddr = addr.lowerAddr();
    const labelList& upperAddr = addr.upperAddr();
    forAll(lowerAddr, i)
    {
        tripletList.push_back(Trip(lowerAddr[i], upperAddr[i], foam_matrix.upper()[i]));
        tripletList.push_back(Trip(lowerAddr[i] + sizeA, upperAddr[i] + sizeA,
                                   foam_matrix.upper()[i]));
        tripletList.push_back(Trip(lowerAddr[i] + sizeA * 2, upperAddr[i] + sizeA * 2,
                                   foam_matrix.upper()[i]));
        tripletList.push_back(Trip(upperAddr[i], lowerAddr[i], foam_matrix.lower()[i]));
        tripletList.push_back(Trip(upperAddr[i] + sizeA, lowerAddr[i] + sizeA,
                                   foam_matrix.lower()[i]));
        tripletList.push_back(Trip(upperAddr[i] + sizeA * 2, lowerAddr[i] + sizeA * 2,
                                   foam_matrix.lower()[i]));
    }
    forAll(foam_matrix.psi().boundaryField(), I)
    {
        const fvPatch& ptch = foam_matrix.psi().boundaryField()[I].patch();
        forAll(ptch, J)
        {
            int w = ptch.faceCells()[J];
            tripletList.push_back(Trip(w, w, foam_matrix.internalCoeffs()[I][J][0]));
            tripletList.push_back(Trip(w + sizeA, w + sizeA,
                                       foam_matrix.internalCoeffs()[I][J][1]));
            tripletList.push_back(Trip(w + sizeA * 2, w + sizeA * 2,
                                       foam_matrix.internalCoeffs()[I][J][2]));
            b(w)   += foam_matrix.boundaryCoeffs()[I][J][0];
            b(w + sizeA)   += foam_matrix.boundaryCoeffs()[I][J][1];
            b(w + sizeA * 2)   += foam_matrix.boundaryCoeffs()[I][J][2];
        }
    }
    A.setFromTriplets(tripletList.begin(), tripletList.end());
}

template<>
void Foam2Eigen::fvMatrix2EigenM(fvMatrix<scalar>& foam_matrix,
                                 Eigen::MatrixXd& A)
{
    int sizeA = foam_matrix.diag().size();
    A.setZero(sizeA, sizeA);

    for (auto i = 0; i < sizeA; i++)
    {
        A(i, i) = foam_matrix.diag()[i];
    }

    const lduAddressing& addr = foam_matrix.lduAddr();
    const labelList& lowerAddr = addr.lowerAddr();
    const labelList& upperAddr = addr.upperAddr();
    forAll(lowerAddr, i)
    {
        A(lowerAddr[i], upperAddr[i]) = foam_matrix.upper()[i];
        A(upperAddr[i], lowerAddr[i]) = foam_matrix.lower()[i];
    }
    forAll(foam_matrix.psi().boundaryField(), I)
    {
        const fvPatch& ptch = foam_matrix.psi().boundaryField()[I].patch();
        forAll(ptch, J)
        {
            int w = ptch.faceCells()[J];
            A(w, w) += foam_matrix.internalCoeffs()[I][J];
        }
    }
}

template<>
void Foam2Eigen::fvMatrix2EigenM(fvMatrix<scalar>& foam_matrix,
                                 Eigen::SparseMatrix<double>& A)
{
    int sizeA = foam_matrix.diag().size();
    int nel = foam_matrix.diag().size() + foam_matrix.upper().size() +
              foam_matrix.lower().size();
    A.resize(sizeA, sizeA);
    A.reserve(nel);
    typedef Eigen::Triplet<double> Trip;
    std::vector<Trip> tripletList;
    tripletList.reserve(nel);

    for (auto i = 0; i < sizeA; i++)
    {
        tripletList.push_back(Trip(i, i, foam_matrix.diag()[i]));
    }

    const lduAddressing& addr = foam_matrix.lduAddr();
    const labelList& lowerAddr = addr.lowerAddr();
    const labelList& upperAddr = addr.upperAddr();
    forAll(lowerAddr, i)
    {
        tripletList.push_back(Trip(lowerAddr[i], upperAddr[i], foam_matrix.upper()[i]));
        tripletList.push_back(Trip(upperAddr[i], lowerAddr[i], foam_matrix.lower()[i]));
    }
    forAll(foam_matrix.psi().boundaryField(), I)
    {
        const fvPatch& ptch = foam_matrix.psi().boundaryField()[I].patch();
        forAll(ptch, J)
        {
            int w = ptch.faceCells()[J];
            tripletList.push_back(Trip(w, w, foam_matrix.internalCoeffs()[I][J]));
        }
    }
    A.setFromTriplets(tripletList.begin(), tripletList.end());
}

template<>
void Foam2Eigen::fvMatrix2EigenM(fvMatrix<vector>& foam_matrix,
                                 Eigen::MatrixXd& A)
{
    int sizeA = foam_matrix.diag().size();
    A.resize(sizeA * 3, sizeA * 3);

    for (auto i = 0; i < sizeA; i++)
    {
        A(i, i) = foam_matrix.diag()[i];
        A(sizeA + i, sizeA + i) = foam_matrix.diag()[i];
        A(2 * sizeA + i, 2 * sizeA + i) = foam_matrix.diag()[i];
    }

    const lduAddressing& addr = foam_matrix.lduAddr();
    const labelList& lowerAddr = addr.lowerAddr();
    const labelList& upperAddr = addr.upperAddr();
    forAll(lowerAddr, i)
    {
        A(lowerAddr[i], upperAddr[i]) = foam_matrix.upper()[i];
        A(lowerAddr[i] + sizeA, upperAddr[i] + sizeA) = foam_matrix.upper()[i];
        A(lowerAddr[i] + sizeA * 2, upperAddr[i] + sizeA * 2) = foam_matrix.upper()[i];
        A(upperAddr[i], lowerAddr[i]) = foam_matrix.lower()[i];
        A(upperAddr[i] + sizeA, lowerAddr[i] + sizeA) = foam_matrix.lower()[i];
        A(upperAddr[i] + sizeA * 2, lowerAddr[i] + sizeA * 2) = foam_matrix.lower()[i];
    }
    forAll(foam_matrix.psi().boundaryField(), I)
    {
        const fvPatch& ptch = foam_matrix.psi().boundaryField()[I].patch();
        forAll(ptch, J)
        {
            int w = ptch.faceCells()[J];
            A(w, w) += foam_matrix.internalCoeffs()[I][J][0];
            A(w + sizeA, w + sizeA) += foam_matrix.internalCoeffs()[I][J][1];
            A(w + sizeA * 2, w + sizeA * 2) += foam_matrix.internalCoeffs()[I][J][2];
        }
    }
}


template<>
void Foam2Eigen::fvMatrix2EigenM(fvMatrix<vector>& foam_matrix,
                                 Eigen::SparseMatrix<double>& A)
{
    int sizeA = foam_matrix.diag().size();
    int nel = foam_matrix.diag().size() + foam_matrix.upper().size() +
              foam_matrix.lower().size();
    A.resize(sizeA * 3, sizeA * 3);
    A.reserve(nel * 3);
    typedef Eigen::Triplet<double> Trip;
    std::vector<Trip> tripletList;
    tripletList.reserve(nel * 3);

    for (auto i = 0; i < sizeA; i++)
    {
        tripletList.push_back(Trip(i, i, foam_matrix.diag()[i]));
        tripletList.push_back(Trip(sizeA + i, sizeA + i, foam_matrix.diag()[i]));
        tripletList.push_back(Trip(2 * sizeA + i, 2 * sizeA + i,
                                   foam_matrix.diag()[i]));
    }

    const lduAddressing& addr = foam_matrix.lduAddr();
    const labelList& lowerAddr = addr.lowerAddr();
    const labelList& upperAddr = addr.upperAddr();
    forAll(lowerAddr, i)
    {
        tripletList.push_back(Trip(lowerAddr[i], upperAddr[i], foam_matrix.upper()[i]));
        tripletList.push_back(Trip(lowerAddr[i] + sizeA, upperAddr[i] + sizeA,
                                   foam_matrix.upper()[i]));
        tripletList.push_back(Trip(lowerAddr[i] + sizeA * 2, upperAddr[i] + sizeA * 2,
                                   foam_matrix.upper()[i]));
        tripletList.push_back(Trip(upperAddr[i], lowerAddr[i], foam_matrix.lower()[i]));
        tripletList.push_back(Trip(upperAddr[i] + sizeA, lowerAddr[i] + sizeA,
                                   foam_matrix.lower()[i]));
        tripletList.push_back(Trip(upperAddr[i] + sizeA * 2, lowerAddr[i] + sizeA * 2,
                                   foam_matrix.lower()[i]));
    }
    forAll(foam_matrix.psi().boundaryField(), I)
    {
        const fvPatch& ptch = foam_matrix.psi().boundaryField()[I].patch();
        forAll(ptch, J)
        {
            int w = ptch.faceCells()[J];
            tripletList.push_back(Trip(w, w, foam_matrix.internalCoeffs()[I][J][0]));
            tripletList.push_back(Trip(w + sizeA, w + sizeA,
                                       foam_matrix.internalCoeffs()[I][J][1]));
            tripletList.push_back(Trip(w + sizeA * 2, w + sizeA * 2,
                                       foam_matrix.internalCoeffs()[I][J][2]));
        }
    }
    A.setFromTriplets(tripletList.begin(), tripletList.end());
}

template<>
void Foam2Eigen::fvMatrix2EigenV(fvMatrix<scalar>& foam_matrix,
                                 Eigen::VectorXd& b)
{
    int sizeA = foam_matrix.diag().size();
    b.setZero(sizeA);

    for (auto i = 0; i < sizeA; i++)
    {
        b(i, 0) = foam_matrix.source()[i];
    }

    forAll(foam_matrix.psi().boundaryField(), I)
    {
        const fvPatch& ptch = foam_matrix.psi().boundaryField()[I].patch();
        forAll(ptch, J)
        {
            int w = ptch.faceCells()[J];
            b(w, 0)   += foam_matrix.boundaryCoeffs()[I][J];
        }
    }
}

template<>
void Foam2Eigen::fvMatrix2EigenV(fvMatrix<vector>& foam_matrix,
                                 Eigen::VectorXd& b)
{
    int sizeA = foam_matrix.diag().size();
    b.resize(sizeA * 3);

    for (auto i = 0; i < sizeA; i++)
    {
        b(i) = foam_matrix.source()[i][0];
        b(sizeA + i) = foam_matrix.source()[i][1];
        b(2 * sizeA + i) = foam_matrix.source()[i][2];
    }

    forAll(foam_matrix.psi().boundaryField(), I)
    {
        const fvPatch& ptch = foam_matrix.psi().boundaryField()[I].patch();
        forAll(ptch, J)
        {
            int w = ptch.faceCells()[J];
            b(w)   += foam_matrix.boundaryCoeffs()[I][J][0];
            b(w + sizeA)   += foam_matrix.boundaryCoeffs()[I][J][1];
            b(w + sizeA * 2)   += foam_matrix.boundaryCoeffs()[I][J][2];
        }
    }
}

template <>
Field<scalar> Foam2Eigen::Eigen2field(Field<scalar>& field,
                                      Eigen::MatrixXd& matrix)
{
    int sizeBC = field.size();
    M_Assert(matrix.cols() == 1,
             "The number of columns of the Input members is not correct, it should be 1");

    if (matrix.rows() == 1)
    {
        Eigen::MatrixXd new_matrix = matrix.replicate(sizeBC, 1);
        matrix.conservativeResize(sizeBC, 1);
        matrix = new_matrix;
    }

    std::string message = "The input Eigen::MatrixXd has size " + name(
                              matrix.rows()) +
                          ". It should have the same size of the Field, i.e. " +
                          name(sizeBC);
    M_Assert(matrix.rows() == sizeBC, message.c_str());

    for (auto i = 0; i < sizeBC; i++)
    {
        field[i] = matrix(i, 0);
    }

    return field;
}

template <>
Field<vector> Foam2Eigen::Eigen2field(Field<vector>& field,
                                      Eigen::MatrixXd& matrix)
{
    int sizeBC = field.size();
    M_Assert(matrix.cols() == 3,
             "The number of columns of the Input members is not correct, it should be 1");

    if (matrix.rows() == 1)
    {
        Eigen::MatrixXd new_matrix = matrix.replicate(sizeBC, 1);
        matrix.conservativeResize(sizeBC, 3);
        matrix = new_matrix;
    }

    // std::string message = "The size of the input Matrices " + name(
    //                           valueFrac.rows()) +
    //                       " must be equal to the dimension of the boundary condition you want to set.";
    M_Assert(matrix.rows() == sizeBC, "message.c_str()");

    for (auto i = 0; i < sizeBC; i++)
    {
        field[i][0] = matrix(i, 0);
        field[i][1] = matrix(i, 1);
        field[i][2] = matrix(i, 2);
    }

    return field;
}
