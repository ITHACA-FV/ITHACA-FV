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
template<template<class> class PatchField, class GeoMesh>
Eigen::VectorXd Foam2Eigen::field2Eigen(
    GeometricField<vector, PatchField, GeoMesh>& field)
{
    Eigen::VectorXd out;
    out.resize(label(field.size() * 3));

    for (label l = 0; l < field.size(); l++)
    {
        out(l) = field[l][0];
        out(field.size() +  l ) = field[l][1];
        out(2 * field.size() + l ) = field[l][2];
    }

    return out;
}

template Eigen::VectorXd Foam2Eigen::field2Eigen(
    volVectorField& field);

template<template<class> class PatchField, class GeoMesh>
Eigen::VectorXd Foam2Eigen::field2Eigen(
    GeometricField<scalar, PatchField, GeoMesh>& field)
{
    Eigen::VectorXd out;
    out.resize(label(field.size()));

    for (label l = 0; l < field.size(); l++)
    {
        out(l) = field[l];
    }

    return out;
}

template Eigen::VectorXd Foam2Eigen::field2Eigen(
    volScalarField& field);
template Eigen::VectorXd Foam2Eigen::field2Eigen(
    surfaceScalarField& field);

template<>
Eigen::VectorXd Foam2Eigen::field2Eigen(const Field<scalar>& field)
{
    Eigen::VectorXd out;
    out.resize(label(field.size()));

    for (label l = 0; l < field.size(); l++)
    {
        out(l) = field[l];
    }

    return out;
}

template<>
Eigen::VectorXd Foam2Eigen::field2Eigen(const Field<vector>& field)
{
    Eigen::VectorXd out;
    out.resize(label(field.size() * 3));

    for (label l = 0; l < field.size(); l++)
    {
        out(l) = field[l][0];
        out(field.size() +  l ) = field[l][1];
        out(2 * field.size() + l ) = field[l][2];
    }

    return out;
}

template<>
Eigen::VectorXd Foam2Eigen::field2Eigen(const
                                        DimensionedField<scalar, Foam::volMesh>& field)
{
    Eigen::VectorXd out;
    out.resize(label(field.size()));

    for (label l = 0; l < field.size(); l++)
    {
        out(l) = field[l];
    }

    return out;
}

template<template<class> class PatchField, class GeoMesh>
List<Eigen::VectorXd> Foam2Eigen::field2EigenBC(
    GeometricField<vector, PatchField, GeoMesh>& field)
{
    List<Eigen::VectorXd> Out;
    label size = field.boundaryField().size();
    Out.resize(size);

    for (label i = 0; i < size; i++ )
    {
        label sizei = field.boundaryField()[i].size();
        Out[i].resize(sizei * 3);

        for (label k = 0; k < sizei ; k++)
        {
            Out[i](k) = field.boundaryField()[i][k][0];
            Out[i](k + sizei) = field.boundaryField()[i][k][1];
            Out[i](k + 2 * sizei) = field.boundaryField()[i][k][2];
        }
    }

    return Out;
}

template List<Eigen::VectorXd> Foam2Eigen::field2EigenBC(
    volVectorField& field);


template<template<class> class PatchField, class GeoMesh>
List<Eigen::VectorXd> Foam2Eigen::field2EigenBC(
    GeometricField<scalar, PatchField, GeoMesh>& field)
{
    List<Eigen::VectorXd> Out;
    label size = field.boundaryField().size();
    Out.resize(size);

    for (label i = 0; i < size; i++ )
    {
        label sizei = field.boundaryField()[i].size();
        Out[i].resize(sizei);

        for (label k = 0; k < sizei ; k++)
        {
            Out[i](k) = field.boundaryField()[i][k];
        }
    }

    return Out;
}
template List<Eigen::VectorXd> Foam2Eigen::field2EigenBC(
    volScalarField& field);

template<template<class> class PatchField, class GeoMesh>
List<Eigen::MatrixXd> Foam2Eigen::PtrList2EigenBC(
    PtrList<GeometricField<scalar, PatchField, GeoMesh>>&
    fields, label Nfields)
{
    label Nf;
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
    label NBound = fields[0].boundaryField().size();
    Out.resize(NBound);

    for (label i = 0; i < NBound; i++)
    {
        label sizei = fields[0].boundaryField()[i].size();
        Out[i].resize(sizei, Nf);
    }

    for (label k = 0; k < Nf; k++)
    {
        List<Eigen::VectorXd> temp;
        temp = field2EigenBC(fields[k]);

        for (label i = 0; i < NBound; i++)
        {
            Out[i].col(k) = temp[i];
        }
    }

    return Out;
}

template List<Eigen::MatrixXd> Foam2Eigen::PtrList2EigenBC(
    PtrList<volScalarField>& fields, label Nfields);
template List<Eigen::MatrixXd> Foam2Eigen::PtrList2EigenBC(
    PtrList<surfaceScalarField>& fields, label Nfields);


template<template<class> class PatchField, class GeoMesh>
List<Eigen::MatrixXd> Foam2Eigen::PtrList2EigenBC(
    PtrList<GeometricField<vector, PatchField, GeoMesh>>&
    fields, label Nfields)
{
    label Nf;
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
    label NBound = fields[0].boundaryField().size();
    Out.resize(NBound);

    for (label i = 0; i < NBound; i++)
    {
        label sizei = fields[0].boundaryField()[i].size();
        Out[i].resize(sizei * 3, Nf);
    }

    for (label k = 0; k < Nf; k++)
    {
        List<Eigen::VectorXd> temp;
        temp = field2EigenBC(fields[k]);

        for (label i = 0; i < NBound; i++)
        {
            Out[i].col(k) = temp[i];
        }
    }

    return Out;
}

template List<Eigen::MatrixXd> Foam2Eigen::PtrList2EigenBC(
    PtrList<volVectorField>& fields, label Nfields);

template<template<class> class PatchField, class GeoMesh>
GeometricField<vector, PatchField, GeoMesh> Foam2Eigen::Eigen2field(
    GeometricField<vector, PatchField, GeoMesh>& field_in,
    Eigen::VectorXd& eigen_vector, bool correctBC)
{
    GeometricField<vector, PatchField, GeoMesh> field_out(field_in);

    for (auto i = 0; i < field_out.size(); i++)
    {
        field_out.ref()[i][0] = eigen_vector(i);
        field_out.ref()[i][1] = eigen_vector(i + field_out.size());
        field_out.ref()[i][2] = eigen_vector(i + field_out.size() * 2);
    }

    if (correctBC)
    {
        field_out.correctBoundaryConditions();
    }

    return field_out;
}

template volVectorField Foam2Eigen::Eigen2field(
    volVectorField& field_in, Eigen::VectorXd& eigen_vector, bool correctBC);

template<template<class> class PatchField, class GeoMesh>
GeometricField<scalar, PatchField, GeoMesh> Foam2Eigen::Eigen2field(
    GeometricField<scalar, PatchField, GeoMesh>& field_in,
    Eigen::VectorXd& eigen_vector, bool correctBC)
{
    GeometricField<scalar, PatchField, GeoMesh> field_out(field_in);

    for (auto i = 0; i < field_out.size(); i++)
    {
        field_out.ref()[i] = eigen_vector(i);
    }

    return field_out;
}

template surfaceScalarField Foam2Eigen::Eigen2field(
    surfaceScalarField& field_in,
    Eigen::VectorXd& eigen_vector,
    bool correctBC);

template<>
volScalarField Foam2Eigen::Eigen2field(
    volScalarField& field_in, Eigen::VectorXd& eigen_vector, bool correctBC)
{
    GeometricField<scalar, fvPatchField, volMesh> field_out(field_in);

    for (auto i = 0; i < field_out.size(); i++)
    {
        field_out.ref()[i] = eigen_vector(i);
    }

    if (correctBC)
    {
        field_out.correctBoundaryConditions();
    }

    return field_out;
}

template<>
Field<scalar> Foam2Eigen::Eigen2field(
    Field<scalar>& field, Eigen::MatrixXd& matrix, bool correctBC)
{
    label sizeBC = field.size();
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

template<>
Field<vector> Foam2Eigen::Eigen2field(
    Field<vector>& field, Eigen::MatrixXd& matrix, bool correctBC)
{
    label sizeBC = field.size();
    M_Assert(matrix.cols() == 3,
             "The number of columns of the Input members is not correct, it should be 1");

    if (matrix.rows() == 1)
    {
        Eigen::MatrixXd new_matrix = matrix.replicate(sizeBC, 1);
        matrix.conservativeResize(sizeBC, 3);
        matrix = new_matrix;
    }

    std::string message = "The input Eigen::MatrixXd has size " + name(
                              matrix.rows()) +
                          ". It should have the same size of the Field, i.e. " +
                          name(sizeBC);
    M_Assert(matrix.rows() == sizeBC, message.c_str());

    for (auto i = 0; i < sizeBC; i++)
    {
        field[i][0] = matrix(i, 0);
        field[i][1] = matrix(i, 1);
        field[i][2] = matrix(i, 2);
    }

    return field;
}

template<class Type, template<class> class PatchField, class GeoMesh>
Eigen::MatrixXd Foam2Eigen::PtrList2Eigen(
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& fields,
    label Nfields)
{
    label Nf;
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
    label nrows = (field2Eigen(fields[0])).rows();
    out.resize(nrows, Nf);

    for (label k = 0; k < Nf; k++)
    {
        out.col(k) = field2Eigen(fields[k]);
    }

    return out;
}

template Eigen::MatrixXd
Foam2Eigen::PtrList2Eigen<scalar, fvPatchField, volMesh>(PtrList<volScalarField>&
        fields, label Nfields);
template Eigen::MatrixXd Foam2Eigen::PtrList2Eigen(PtrList<surfaceScalarField>&
        fields, label Nfields);
template Eigen::MatrixXd
Foam2Eigen::PtrList2Eigen<vector, fvPatchField, volMesh>(PtrList<volVectorField>&
        fields, label Nfields);

template<>
void Foam2Eigen::fvMatrix2Eigen(fvMatrix<scalar> foam_matrix,
                                Eigen::MatrixXd& A,
                                Eigen::VectorXd& b)
{
    label sizeA = foam_matrix.diag().size();
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
            label w = ptch.faceCells()[J];
            const double intern = foam_matrix.internalCoeffs()[I][J];
            A(w, w) += intern;
            b(w, 0)   += foam_matrix.boundaryCoeffs()[I][J];
        }
    }
}

template<>
void Foam2Eigen::fvMatrix2Eigen(fvMatrix<vector> foam_matrix,
                                Eigen::MatrixXd& A,
                                Eigen::VectorXd& b)
{
    label sizeA = foam_matrix.diag().size();
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
            label w = ptch.faceCells()[J];
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
void Foam2Eigen::fvMatrix2Eigen(fvMatrix<scalar> foam_matrix,
                                Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b)
{
    label sizeA = foam_matrix.diag().size();
    label nel = foam_matrix.diag().size() + foam_matrix.upper().size() +
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
            label w = ptch.faceCells()[J];
            tripletList.push_back(Trip(w, w, foam_matrix.internalCoeffs()[I][J]));
            b(w, 0)   += foam_matrix.boundaryCoeffs()[I][J];
        }
    }
    A.setFromTriplets(tripletList.begin(), tripletList.end());
}

template<>
void Foam2Eigen::fvMatrix2Eigen(fvMatrix<vector> foam_matrix,
                                Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b)
{
    label sizeA = foam_matrix.diag().size();
    label nel = foam_matrix.diag().size() + foam_matrix.upper().size() +
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
            label w = ptch.faceCells()[J];
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
    label sizeA = foam_matrix.diag().size();
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
            label w = ptch.faceCells()[J];
            A(w, w) += foam_matrix.internalCoeffs()[I][J];
        }
    }
}

template<>
void Foam2Eigen::fvMatrix2EigenM(fvMatrix<scalar>& foam_matrix,
                                 Eigen::SparseMatrix<double>& A)
{
    label sizeA = foam_matrix.diag().size();
    label nel = foam_matrix.diag().size() + foam_matrix.upper().size() +
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
            label w = ptch.faceCells()[J];
            tripletList.push_back(Trip(w, w, foam_matrix.internalCoeffs()[I][J]));
        }
    }
    A.setFromTriplets(tripletList.begin(), tripletList.end());
}

template<>
void Foam2Eigen::fvMatrix2EigenM(fvMatrix<vector>& foam_matrix,
                                 Eigen::MatrixXd& A)
{
    label sizeA = foam_matrix.diag().size();
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
            label w = ptch.faceCells()[J];
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
    label sizeA = foam_matrix.diag().size();
    label nel = foam_matrix.diag().size() + foam_matrix.upper().size() +
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
            label w = ptch.faceCells()[J];
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
    label sizeA = foam_matrix.diag().size();
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
            label w = ptch.faceCells()[J];
            b(w, 0)   += foam_matrix.boundaryCoeffs()[I][J];
        }
    }
}

template<>
void Foam2Eigen::fvMatrix2EigenV(fvMatrix<vector>& foam_matrix,
                                 Eigen::VectorXd& b)
{
    label sizeA = foam_matrix.diag().size();
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
            label w = ptch.faceCells()[J];
            b(w)   += foam_matrix.boundaryCoeffs()[I][J][0];
            b(w + sizeA)   += foam_matrix.boundaryCoeffs()[I][J][1];
            b(w + sizeA * 2)   += foam_matrix.boundaryCoeffs()[I][J][2];
        }
    }
}

template<class Type, template<class> class PatchField, class GeoMesh>
Eigen::VectorXd Foam2Eigen::projectField(
    GeometricField<Type, PatchField, GeoMesh>& field,
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& modes,
    label Nmodes)
{
    Eigen::VectorXd fr;
    Eigen::MatrixXd Eig_Modes = PtrList2Eigen(modes, Nmodes);
    Eigen::VectorXd f = Foam2Eigen::field2Eigen(field);
    Eigen::VectorXd Volumes = field2Eigen(modes[0].mesh());
    Eigen::MatrixXd VolumesN(Volumes.rows(), 1);
    VolumesN = Volumes;

    if (Volumes.rows() != Eig_Modes.rows())
    {
        VolumesN.resize(Eig_Modes.rows(), 1);
        VolumesN.col(0).segment(0, Volumes.rows()) = Volumes;
        VolumesN.col(0).segment(Volumes.rows() + 1, Volumes.rows() * 2) = Volumes;
        VolumesN.col(0).segment(Volumes.rows() * 2 + 1, Volumes.rows() * 3) = Volumes;
    }

    fr = Eig_Modes.transpose() * (f.cwiseProduct(VolumesN));
    return fr;
}

template<class Type, template<class> class PatchField, class GeoMesh>
std::tuple<Eigen::MatrixXd, Eigen::VectorXd> Foam2Eigen::projectFvMatrix(
    fvMatrix<Type>& matrix,
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& modes, label Nmodes)
{
    Eigen::SparseMatrix<double> A;
    Eigen::MatrixXd Ar;
    Eigen::VectorXd b;
    Eigen::VectorXd br;
    Eigen::MatrixXd Eig_Modes = PtrList2Eigen(modes, Nmodes);
    Foam2Eigen::fvMatrix2Eigen(matrix, A, b);
    Eigen::VectorXd Volumes = field2Eigen(modes[0].mesh());
    Eigen::MatrixXd VolumesN(Volumes.rows(), Nmodes);

    if (Volumes.rows() != Eig_Modes.rows())
    {
        VolumesN.resize(Eig_Modes.rows(), Nmodes);
    }

    if (Volumes.rows() == Eig_Modes.rows())
    {
        for (auto i = 0; i < Nmodes; i++)
        {
            VolumesN.col(i) = Volumes;
        }
    }
    else
    {
        for (auto i = 0; i < Nmodes; i++)
        {
            VolumesN.col(i).segment(0, Volumes.rows()) = Volumes;
            VolumesN.col(i).segment(Volumes.rows() + 1, Volumes.rows() * 2) = Volumes;
            VolumesN.col(i).segment(Volumes.rows() * 2 + 1, Volumes.rows() * 3) = Volumes;
        }
    }

    Ar = Eig_Modes.transpose() * A * Eig_Modes;
    br = Eig_Modes.transpose() * b;
    std::tuple <Eigen::MatrixXd, Eigen::VectorXd> tupla;
    tupla = std::make_tuple(Ar, br);
    return tupla;
}

template<class Type, template<class> class PatchField, class GeoMesh>
Eigen::MatrixXd Foam2Eigen::MassMatrix(
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& modes, label Nmodes)
{
    Eigen::MatrixXd Mr;
    Eigen::MatrixXd Eig_Modes = PtrList2Eigen(modes, Nmodes);
    Eigen::VectorXd Volumes = field2Eigen(modes[0].mesh());
    Eigen::MatrixXd VolumesN(Volumes.rows(), Nmodes);

    if (Volumes.rows() != Eig_Modes.rows())
    {
        VolumesN.resize(Eig_Modes.rows(), Nmodes);
    }

    if (Volumes.rows() == Eig_Modes.rows())
    {
        for (auto i = 0; i < Nmodes; i++)
        {
            VolumesN.col(i) = Volumes;
        }
    }
    else
    {
        for (auto i = 0; i < Nmodes; i++)
        {
            VolumesN.col(i).segment(0, Volumes.rows()) = Volumes;
            VolumesN.col(i).segment(Volumes.rows() + 1, Volumes.rows() * 2) = Volumes;
            VolumesN.col(i).segment(Volumes.rows() * 2 + 1, Volumes.rows() * 3) = Volumes;
        }
    }

    Mr = Eig_Modes.transpose() * (Eig_Modes.cwiseProduct(VolumesN));
    return Mr;
}

template<class Type>
std::tuple<List<Eigen::SparseMatrix<double>>, List<Eigen::VectorXd>>
        Foam2Eigen::LFvMatrix2LSM(PtrList<fvMatrix<Type>>& MatrixList)
{
    List<Eigen::SparseMatrix<double>> SM_list;
    List<Eigen::VectorXd> V_list;
    label LSize =  MatrixList.size();
    SM_list.resize(LSize);
    V_list.resize(LSize);
    Eigen::SparseMatrix<double> A;
    Eigen::VectorXd b;

    for (label j = 0; j < LSize; j++)
    {
        fvMatrix2Eigen(MatrixList[j], A, b);
        SM_list[j] = A;
        V_list[j] = b;
    }

    std::tuple <List<Eigen::SparseMatrix<double>>, List<Eigen::VectorXd>> tupla;
    tupla = std::make_tuple(SM_list, V_list);
    return tupla;
}

template std::tuple<List<Eigen::SparseMatrix<double>>, List<Eigen::VectorXd>>
Foam2Eigen::LFvMatrix2LSM(PtrList<fvMatrix<scalar>>& MatrixList);
template std::tuple<List<Eigen::SparseMatrix<double>>, List<Eigen::VectorXd>>
Foam2Eigen::LFvMatrix2LSM(PtrList<fvMatrix<vector>>& MatrixList);

template<class type_matrix>
Eigen::Matrix<type_matrix, Eigen::Dynamic, Eigen::Dynamic>
Foam2Eigen::List2EigenMatrix ( List<type_matrix> list )
{
    Eigen::Matrix<type_matrix, Eigen::Dynamic, Eigen::Dynamic> matrix(list.size(),
            1);

    for (label i = 0; i < matrix.rows(); i++)
    {
        matrix(i, 0) = list[i];
    }

    return matrix;
}

template Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>
Foam2Eigen::List2EigenMatrix ( List<int> list );
template Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
Foam2Eigen::List2EigenMatrix ( List<double> list );

template<class type_matrix>
List<type_matrix> Foam2Eigen::EigenMatrix2List (
    Eigen::Matrix<type_matrix, Eigen::Dynamic, Eigen::Dynamic> matrix )
{
    if (matrix.cols() == 1)
    {
        List<type_matrix> list(matrix.rows());

        for (label i = 0; i < matrix.rows(); i++)
        {
            list[i] = matrix(i, 0);
        }

        return list;
    }
    else
    {
        Info << "Foam2Eigen::EigenMatrix2List only accepts matrices with 1 column, exiting"
             << endl;
        exit(11);
    }
}

template List<int> Foam2Eigen::EigenMatrix2List (
    Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> matrix );
template List<double> Foam2Eigen::EigenMatrix2List (
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix );
