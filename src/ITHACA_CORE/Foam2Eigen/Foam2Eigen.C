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

template <>
Eigen::MatrixXd Foam2Eigen::field2Eigen(
    volScalarField& field)
{
    Eigen::MatrixXd out = Eigen::Map<Eigen::MatrixXd>(const_cast<double*>(&
        (field[0])), field.size(), 1);
    return out;
};

template <>
Eigen::MatrixXd Foam2Eigen::field2Eigen(
    volVectorField& field)
{
    Eigen::MatrixXd out = Eigen::Map<Eigen::MatrixXd>(&field.ref()[0][0],
        field.size() * 3, 1);
    return out;
};

template <>
Eigen::MatrixXd Foam2Eigen::field2Eigen(
    volTensorField& field)
{
    Eigen::MatrixXd out = Eigen::Map<Eigen::MatrixXd>(&field.ref()[0][0],
        field.size() * 9, 1);
    return out;
};

template <>
Eigen::MatrixXd Foam2Eigen::field2Eigen(
    pointVectorField& field)
{
    Eigen::MatrixXd out = Eigen::Map<Eigen::MatrixXd>(&field.ref()[0][0],
        field.size() * 3, 1);
    return out;
};

template <>
Eigen::MatrixXd Foam2Eigen::field2Eigen(
    surfaceScalarField& field)
{
    Eigen::MatrixXd out = Eigen::Map<Eigen::MatrixXd>(&field.ref()[0], field.size(),
        1);
    return out;
};

template Eigen::MatrixXd Foam2Eigen::field2Eigen(
    volScalarField& field);

template Eigen::MatrixXd Foam2Eigen::field2Eigen(
    volTensorField& field);

template Eigen::MatrixXd Foam2Eigen::field2Eigen(
    volVectorField& field);

template Eigen::MatrixXd Foam2Eigen::field2Eigen(
    pointVectorField& field);

template Eigen::MatrixXd Foam2Eigen::field2Eigen(
    surfaceScalarField& field);

template <>
Eigen::Map<Eigen::MatrixXd> Foam2Eigen::field2EigenMap(
    volScalarField& field)
{
    Eigen::Map<Eigen::MatrixXd> output(field.ref().data(), field.size(), 1);
    return std::move(output);
}

template <>
Eigen::Map<Eigen::MatrixXd> Foam2Eigen::field2EigenMap(
    volVectorField& field)
{
    Eigen::Map<Eigen::MatrixXd> output(&field.ref()[0][0], field.size() * 3, 1);
    return std::move(output);
}

template <>
Eigen::Map<Eigen::MatrixXd> Foam2Eigen::field2EigenMap(
    volTensorField& field)
{
    Eigen::Map<Eigen::MatrixXd> output(&field.ref()[0][0], field.size() * 9, 1);
    return std::move(output);
}

template <>
Eigen::Map<Eigen::MatrixXd> Foam2Eigen::field2EigenMapBC(
    volScalarField& field, int BC_index)
{
    Eigen::Map<Eigen::MatrixXd> output(field.boundaryFieldRef()[BC_index].data(),
                                       field.boundaryField()[BC_index].size(), 1);
    return std::move(output);
}

template <>
Eigen::Map<Eigen::MatrixXd> Foam2Eigen::field2EigenMapBC(
    volVectorField& field, int BC_index)
{
    Eigen::Map<Eigen::MatrixXd> output(field.boundaryFieldRef()[BC_index][0].data(),
                                       field.boundaryField()[BC_index].size() * 3, 1);
    return std::move(output);
};


template Eigen::Map<Eigen::MatrixXd> Foam2Eigen::field2EigenMap(
    volScalarField& field);

template Eigen::Map<Eigen::MatrixXd> Foam2Eigen::field2EigenMap(
    volVectorField& field);

template Eigen::Map<Eigen::MatrixXd> Foam2Eigen::field2EigenMap(
    volTensorField& field);

template Eigen::Map<Eigen::MatrixXd> Foam2Eigen::field2EigenMapBC(
    volScalarField& field, int BC_index);

template Eigen::Map<Eigen::MatrixXd> Foam2Eigen::field2EigenMapBC(
    volVectorField& field, int BC_index);

template <>
Eigen::VectorXd Foam2Eigen::field2Eigen(const Field<scalar>& field)
{
    Eigen::VectorXd out = Eigen::Map<Eigen::MatrixXd>(const_cast<double*>(&
        (field[0])),
        field.size(), 1);
    return out;
}

template <>
Eigen::VectorXd Foam2Eigen::field2Eigen(const Field<vector>& field)
{
    Eigen::VectorXd out =  Eigen::Map<Eigen::MatrixXd>(const_cast<double*>(&
        (field[0][0])),
        field.size() * 3, 1);
    return out;
}

template <>
Eigen::VectorXd Foam2Eigen::field2Eigen(const Field<tensor>& field)
{
    Eigen::VectorXd out = Eigen::Map<Eigen::MatrixXd>(const_cast<double*>(&
        (field[0][0])),
        field.size() * 9, 1);
    return out;
}

template <>
Eigen::VectorXd Foam2Eigen::field2Eigen(const
                                        DimensionedField<scalar, Foam::volMesh>& field)
{
    Eigen::VectorXd out = Eigen::Map<Eigen::MatrixXd>(const_cast<double*>(&
        (field[0])),
        field.size(), 1);
    return out;
}

template <template <class> class PatchField, class GeoMesh>
List<Eigen::VectorXd> Foam2Eigen::field2EigenBC(
    GeometricField<tensor, PatchField, GeoMesh>& field)
{
    List<Eigen::VectorXd> Out;
    label size = field.boundaryField().size();
    Out.resize(size);

    for (label i = 0; i < size; i++)
    {
        Out[i] = Eigen::Map<Eigen::VectorXd>(const_cast<double*>(&
            (field.boundaryField()[i][0][0])),
                                             field.boundaryField()[i].size() * 9);
    }

    return Out;
}

template List<Eigen::VectorXd> Foam2Eigen::field2EigenBC(
    volTensorField& field);

template <template <class> class PatchField, class GeoMesh>
List<Eigen::VectorXd> Foam2Eigen::field2EigenBC(
    GeometricField<vector, PatchField, GeoMesh>& field)
{
    List<Eigen::VectorXd> Out;
    label size = field.boundaryField().size();
    Out.resize(size);
    constexpr bool check_vol = std::is_same<volMesh, GeoMesh>::value
                               || std::is_same<surfaceMesh, GeoMesh>::value;

    if  constexpr(check_vol)
    {
        for (label i = 0; i < size; i++ )
        {
            Out[i] = Eigen::Map<Eigen::VectorXd>(const_cast<double*>(&
                (field.boundaryField()[i][0][0])),
                                                 field.boundaryField()[i].size() * 3);
        }
    }
    else if  constexpr(std::is_same<pointMesh, GeoMesh>::value)
    {
        for (label i = 0; i < size;
                i++ ) //field.boundaryField()[i].patchInternalField()()[k][j];
        {
            Out[i] = Eigen::Map<Eigen::VectorXd>(const_cast<double*>(&
                (field.boundaryField()[i].patchInternalField()()[0][0])),
                                                 field.boundaryField()[i].patchInternalField()().size() * 3);
        }
    }

    return Out;
}

template List<Eigen::VectorXd> Foam2Eigen::field2EigenBC(
    volVectorField& field);

template <template <class> class PatchField, class GeoMesh>
List<Eigen::VectorXd> Foam2Eigen::field2EigenBC(
    GeometricField<scalar, PatchField, GeoMesh>& field)
{
    List<Eigen::VectorXd> Out;
    label size = field.boundaryField().size();
    Out.resize(size);

    for (label i = 0; i < size; i++)
    {
        Out[i] = Eigen::Map<Eigen::VectorXd>(const_cast<double*>(&
            (field.boundaryField()[i][0])), field.boundaryField()[i].size());
    }

    return Out;
}

template List<Eigen::VectorXd> Foam2Eigen::field2EigenBC(
    volScalarField& field);

template <template <class> class PatchField, class GeoMesh>
List<Eigen::MatrixXd> Foam2Eigen::PtrList2EigenBC(
    PtrList<GeometricField<scalar, PatchField, GeoMesh >>&
    fields,
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


template <template <class> class PatchField, class GeoMesh>
List<Eigen::MatrixXd> Foam2Eigen::PtrList2EigenBC(
    PtrList<GeometricField<vector, PatchField, GeoMesh >>&
    fields,
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

template List<Eigen::MatrixXd> Foam2Eigen::PtrList2EigenBC(
    PtrList<pointVectorField>& fields, label Nfields);

template <template <class> class PatchField, class GeoMesh>
List<Eigen::MatrixXd> Foam2Eigen::PtrList2EigenBC(
    PtrList<GeometricField<tensor, PatchField, GeoMesh >>&
    fields,
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

    List<Eigen::MatrixXd> Out;
    label NBound = fields[0].boundaryField().size();
    Out.resize(NBound);

    for (label i = 0; i < NBound; i++)
    {
        label sizei = fields[0].boundaryField()[i].size();
        Out[i].resize(sizei * 9, Nf);
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
    PtrList<volTensorField>& fields, label Nfields);
// Not consistent with the others, to be fixed, it is changing just the BC
template <template <class> class PatchField, class GeoMesh>
GeometricField<tensor, PatchField, GeoMesh> Foam2Eigen::Eigen2field(
    GeometricField<tensor, PatchField, GeoMesh>& field_in,
    Eigen::VectorXd& eigen_vector, bool correctBC)
{
    GeometricField<tensor, PatchField, GeoMesh> field_out(field_in);

    for (auto i = 0; i < field_out.size(); i++)
    {
        for (label j = 0; j < 9; j++)
        {
            field_out.ref()[i][j] = eigen_vector(i * 9 + j);
        }
    }

    if (correctBC)
    {
        field_out.correctBoundaryConditions();
    }

    return field_out;
}

template volTensorField Foam2Eigen::Eigen2field(
    volTensorField& field_in, Eigen::VectorXd& eigen_vector, bool correctBC);
// This is correct to assign also BCs
template <template <class> class PatchField, class GeoMesh>
GeometricField<vector, PatchField, GeoMesh> Foam2Eigen::Eigen2field(
    GeometricField<vector, PatchField, GeoMesh>& field_in,
    Eigen::VectorXd& eigen_vector, List<Eigen::VectorXd>& eigen_vector_boundary)
{
    GeometricField<vector, PatchField, GeoMesh> field_out(field_in);

    for (auto i = 0; i < field_out.size(); i++)
    {
        for (label j = 0; j < 3; j++)
        {
            field_out.ref()[i][j] = eigen_vector(i * 3 + j);
        }
    }

    for (unsigned int id = 0; id < field_out.boundaryField().size(); id++)
    {
        unsigned int idBSize = field_out.boundaryField()[id].size();

        for (unsigned int ith_bcell = 0; ith_bcell < idBSize; ith_bcell++)
        {
            ITHACAutilities::assignBC(field_out, id, eigen_vector_boundary[id]);
        }
    }

    return field_out;
}

template volVectorField Foam2Eigen::Eigen2field(
    volVectorField& field_in, Eigen::VectorXd& eigen_vector,
    List<Eigen::VectorXd>& eigen_vector_boundary);

template<template<class> class PatchField, class GeoMesh>
GeometricField<scalar, PatchField, GeoMesh> Foam2Eigen::Eigen2field(
    GeometricField<scalar, PatchField, GeoMesh>& field_in,
    Eigen::VectorXd& eigen_vector, List<Eigen::VectorXd>& eigen_vector_boundary)
{
    GeometricField<scalar, PatchField, GeoMesh> field_out(field_in);

    for (auto i = 0; i < field_out.size(); i++)
    {
        field_out.ref()[i] = eigen_vector(i);
    }

    for (unsigned int id = 0; id < field_out.boundaryField().size(); id++)
    {
        for (unsigned int ith_bcell = 0;
                ith_bcell < field_out.boundaryField()[id].size(); ith_bcell++)
        {
            ITHACAutilities::assignBC(field_out, id, eigen_vector_boundary[id]);
        }
    }

    return field_out;
}

template volScalarField Foam2Eigen::Eigen2field(
    volScalarField& field_in, Eigen::VectorXd& eigen_vector,
    List<Eigen::VectorXd>& eigen_vector_boundary);


// Now this is correct, the one for tensor is missing
template<template<class> class PatchField, class GeoMesh>
GeometricField<vector, PatchField, GeoMesh> Foam2Eigen::Eigen2field(
    GeometricField<vector, PatchField, GeoMesh>& field_in,
    Eigen::VectorXd& eigen_vector, bool correctBC)
{
    GeometricField<vector, PatchField, GeoMesh> field_out(field_in);

    for (auto i = 0; i < field_out.size(); i++)
    {
        for (label j = 0; j < 3; j++)
        {
            field_out.ref()[i][j] = eigen_vector(i * 3 + j);
        }
    }

    if (correctBC)
    {
        field_out.correctBoundaryConditions();
    }

    return field_out;
}

template volVectorField Foam2Eigen::Eigen2field(
    volVectorField& field_in, Eigen::VectorXd& eigen_vector, bool correctBC);
template pointVectorField Foam2Eigen::Eigen2field(
    pointVectorField& field_in, Eigen::VectorXd& eigen_vector, bool correctBC);

template <template <class> class PatchField, class GeoMesh>
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

template <>
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

template <>
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

template <>
Field<vector> Foam2Eigen::Eigen2field(
    Field<vector>& field, Eigen::MatrixXd& matrix, bool correctBC)
{
    label sizeBC = field.size();
    M_Assert(matrix.cols() == 1,
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
        for (label j = 0; j < 3; j++)
        {
            field[i][j] = matrix(i, j);
        }
    }

    return field;
}

template <>
Field<tensor> Foam2Eigen::Eigen2field(
    Field<tensor>& field, Eigen::MatrixXd& matrix, bool correctBC)
{
    label sizeBC = field.size();
    M_Assert(matrix.cols() == 1,
             "The number of columns of the Input members is not correct, it should be 1");

    if (matrix.rows() == 1)
    {
        Eigen::MatrixXd new_matrix = matrix.replicate(sizeBC, 1);
        matrix.conservativeResize(sizeBC, 9);
        matrix = new_matrix;
    }

    std::string message = "The input Eigen::MatrixXd has size " + name(
                              matrix.rows()) +
                          ". It should have the same size of the Field, i.e. " +
                          name(sizeBC);
    M_Assert(matrix.rows() == sizeBC, message.c_str());

    for (auto i = 0; i < sizeBC; i++)
    {
        for (label j = 0; j < 9; j++)
        {
            field[i][j] = matrix(i, j);
        }
    }

    return field;
}

template <class Type, template <class> class PatchField, class GeoMesh>
Eigen::MatrixXd Foam2Eigen::PtrList2Eigen(
    PtrList<GeometricField<Type, PatchField, GeoMesh >>& fields,
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
Foam2Eigen::PtrList2Eigen<scalar, fvPatchField, volMesh>
(PtrList<volScalarField>&
 fields,
 label Nfields);
template Eigen::MatrixXd Foam2Eigen::PtrList2Eigen(PtrList<surfaceScalarField>&
        fields,
        label Nfields);
template Eigen::MatrixXd
Foam2Eigen::PtrList2Eigen<vector, fvPatchField, volMesh>(PtrList<volVectorField>&
        fields,
        label Nfields);
template Eigen::MatrixXd
Foam2Eigen::PtrList2Eigen<vector, pointPatchField, pointMesh>
(PtrList<pointVectorField>&
 fields, label Nfields);

template Eigen::MatrixXd
Foam2Eigen::PtrList2Eigen<tensor, fvPatchField, volMesh>
(PtrList<volTensorField>&
 fields,
 label Nfields);


template <>
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
            b(w, 0) += foam_matrix.boundaryCoeffs()[I][J];
        }
    }
}

template <>
void Foam2Eigen::fvMatrix2Eigen
(
    fvMatrix<vector> foam_matrix,
    Eigen::MatrixXd& A,
    Eigen::VectorXd& b
)
{
    const label sizeA = foam_matrix.diag().size();
    const label nComp = 3;

    A.setZero(sizeA * nComp, sizeA * nComp);
    b.setZero(sizeA * nComp);

    auto idx = [nComp](label celli, direction cmpt)
    {
        return nComp * celli + cmpt;
    };


    // ------------------------------------------------------------
    // Diagonal + source
    // ------------------------------------------------------------

    for (label i = 0; i < sizeA; ++i)
    {
        for (direction cmpt = 0; cmpt < nComp; ++cmpt)
        {
            const label ii = idx(i, cmpt);

            A(ii, ii) = foam_matrix.diag()[i];

            b(ii) =
                foam_matrix.source()[i][cmpt];
        }
    }


    // ------------------------------------------------------------
    // Internal ldu coefficients
    // ------------------------------------------------------------

    const lduAddressing& addr =
        foam_matrix.lduAddr();

    const labelList& lowerAddr =
        addr.lowerAddr();

    const labelList& upperAddr =
        addr.upperAddr();


    forAll(lowerAddr, facei)
    {
        const label lowerCell =
            lowerAddr[facei];

        const label upperCell =
            upperAddr[facei];


        for (direction cmpt = 0; cmpt < nComp; ++cmpt)
        {
            const label lowerI =
                idx(lowerCell, cmpt);

            const label upperI =
                idx(upperCell, cmpt);


            A(lowerI, upperI) =
                foam_matrix.upper()[facei];

            A(upperI, lowerI) =
                foam_matrix.lower()[facei];
        }
    }


    // ------------------------------------------------------------
    // Boundary contributions
    // ------------------------------------------------------------

    forAll(foam_matrix.psi().boundaryField(), patchI)
    {
        const fvPatch& patch =
            foam_matrix.psi().boundaryField()[patchI].patch();


        forAll(patch, faceI)
        {
            const label celli =
                patch.faceCells()[faceI];


            for (direction cmpt = 0; cmpt < nComp; ++cmpt)
            {
                const label ii =
                    idx(celli, cmpt);


                A(ii, ii) +=
                    foam_matrix.internalCoeffs()
                    [patchI][faceI][cmpt];


                b(ii) +=
                    foam_matrix.boundaryCoeffs()
                    [patchI][faceI][cmpt];
            }
        }
    }
}


template <typename SparseMatType, typename VecType>
void Foam2Eigen::fvMat2Eigen(fvMatrix<scalar> foam_matrix,
                             SparseMatType& A,
                             VecType& b)
{
    //using Trip = Eigen::Triplet<typename SparseMatType::Scalar>;
    label sizeA = foam_matrix.diag().size();
    label nel = foam_matrix.diag().size() + foam_matrix.upper().size() +
                foam_matrix.lower().size();
    A.resize(sizeA, sizeA);
    b.resize(sizeA);
    A.reserve(nel);
    typedef Eigen::Triplet<double> Trip;
    std::vector<Trip> tripletList;
    tripletList.reserve(nel);

    for (label i = 0; i < sizeA; ++i)
    {
        tripletList.emplace_back(i, i, foam_matrix.diag()[i]);
        b(i) = foam_matrix.source()[i];
    }

    const lduAddressing& addr = foam_matrix.lduAddr();
    const labelList& lowerAddr = addr.lowerAddr();
    const labelList& upperAddr = addr.upperAddr();
    forAll(lowerAddr, i)
    {
        tripletList.emplace_back(lowerAddr[i], upperAddr[i], foam_matrix.upper()[i]);
        tripletList.emplace_back(upperAddr[i], lowerAddr[i], foam_matrix.lower()[i]);
    }
    forAll(foam_matrix.psi().boundaryField(), I)
    {
        const fvPatch& ptch = foam_matrix.psi().boundaryField()[I].patch();
        forAll(ptch, J)
        {
            label w = ptch.faceCells()[J];
            tripletList.emplace_back(w, w, foam_matrix.internalCoeffs()[I][J]);
            b(w) += foam_matrix.boundaryCoeffs()[I][J];
        }
    }
    A.setFromTriplets(tripletList.begin(), tripletList.end());
}

// Explicit instantiations
template void
Foam2Eigen::fvMat2Eigen<Eigen::SparseMatrix<double, Eigen::RowMajor>, Eigen::VectorXd>
(
    fvMatrix<scalar> foam_matrix,
    Eigen::SparseMatrix<double, Eigen::RowMajor>& A,
    Eigen::VectorXd& b);
template void
Foam2Eigen::fvMat2Eigen<Eigen::SparseMatrix<double, Eigen::ColMajor>, Eigen::VectorXd>
(
    fvMatrix<scalar> foam_matrix,
    Eigen::SparseMatrix<double, Eigen::ColMajor>& A,
    Eigen::VectorXd& b);


template<typename SparseMatType, typename VecType>
void Foam2Eigen::fvMat2Eigen
(
    fvMatrix<vector> foam_matrix,
    SparseMatType& A,
    VecType& b
)
{
    const label sizeA = foam_matrix.diag().size();
    const label nComp = 3;

    const label nel =
        foam_matrix.diag().size()
      + foam_matrix.upper().size()
      + foam_matrix.lower().size();

    A.resize(sizeA * nComp, sizeA * nComp);
    A.reserve(nel * nComp);

    b.resize(sizeA * nComp);
    b.setZero();

    typedef Eigen::Triplet<double> Trip;

    std::vector<Trip> tripletList;

    // Slightly more room because boundary diagonal contributions
    // are also inserted as triplets.
    tripletList.reserve
    (
        nel * nComp
      + foam_matrix.psi().boundaryField().size() * nComp
    );


    auto idx = [nComp](label celli, direction cmpt)
    {
        return nComp * celli + cmpt;
    };


    // ============================================================
    // Diagonal + source
    // ============================================================

    for (label celli = 0; celli < sizeA; ++celli)
    {
        for (direction cmpt = 0; cmpt < nComp; ++cmpt)
        {
            const label ii =
                idx(celli, cmpt);

            tripletList.push_back
            (
                Trip
                (
                    ii,
                    ii,
                    foam_matrix.diag()[celli]
                )
            );

            b(ii) =
                foam_matrix.source()[celli][cmpt];
        }
    }


    // ============================================================
    // Internal ldu coefficients
    // ============================================================

    const lduAddressing& addr =
        foam_matrix.lduAddr();

    const labelList& lowerAddr =
        addr.lowerAddr();

    const labelList& upperAddr =
        addr.upperAddr();


    forAll(lowerAddr, facei)
    {
        const label lowerCell =
            lowerAddr[facei];

        const label upperCell =
            upperAddr[facei];


        for (direction cmpt = 0; cmpt < nComp; ++cmpt)
        {
            const label lowerI =
                idx(lowerCell, cmpt);

            const label upperI =
                idx(upperCell, cmpt);


            tripletList.push_back
            (
                Trip
                (
                    lowerI,
                    upperI,
                    foam_matrix.upper()[facei]
                )
            );


            tripletList.push_back
            (
                Trip
                (
                    upperI,
                    lowerI,
                    foam_matrix.lower()[facei]
                )
            );
        }
    }


    // ============================================================
    // Boundary contributions
    // ============================================================

    forAll(foam_matrix.psi().boundaryField(), patchI)
    {
        const fvPatch& patch =
            foam_matrix.psi().boundaryField()[patchI].patch();


        forAll(patch, faceI)
        {
            const label celli =
                patch.faceCells()[faceI];


            for (direction cmpt = 0; cmpt < nComp; ++cmpt)
            {
                const label ii =
                    idx(celli, cmpt);


                // Diagonal boundary contribution
                tripletList.push_back
                (
                    Trip
                    (
                        ii,
                        ii,
                        foam_matrix.internalCoeffs()
                        [patchI][faceI][cmpt]
                    )
                );


                // RHS boundary contribution
                b(ii) +=
                    foam_matrix.boundaryCoeffs()
                    [patchI][faceI][cmpt];
            }
        }
    }


    A.setFromTriplets
    (
        tripletList.begin(),
        tripletList.end()
    );
}

template void
Foam2Eigen::fvMat2Eigen<Eigen::SparseMatrix<double, Eigen::RowMajor>, Eigen::VectorXd>
(
    fvMatrix<vector> foam_matrix,
    Eigen::SparseMatrix<double, Eigen::RowMajor>& A,
    Eigen::VectorXd& b);

template void
Foam2Eigen::fvMat2Eigen<Eigen::SparseMatrix<double, Eigen::ColMajor>, Eigen::VectorXd>
(
    fvMatrix<vector> foam_matrix,
    Eigen::SparseMatrix<double, Eigen::ColMajor>& A,
    Eigen::VectorXd& b);

/////////////////////////////////////////////////////////////////////////////////////////////
template <>
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
            b(w, 0) += foam_matrix.boundaryCoeffs()[I][J];
        }
    }
    A.setFromTriplets(tripletList.begin(), tripletList.end());
}

template <>
void Foam2Eigen::fvMatrix2Eigen
(
    fvMatrix<vector> foam_matrix,
    Eigen::SparseMatrix<double>& A,
    Eigen::VectorXd& b
)
{
    const label sizeA = foam_matrix.diag().size();
    const label nComp = vector::nComponents;
    const label nel =
        foam_matrix.diag().size()
      + foam_matrix.upper().size()
      + foam_matrix.lower().size();

    A.resize(sizeA * nComp, sizeA * nComp);
    A.reserve(nel * nComp);

    b.setZero(sizeA * nComp);

    typedef Eigen::Triplet<double> Trip;
    std::vector<Trip> tripletList;
    tripletList.reserve(nel * nComp);

    auto idx = [nComp](label celli, direction cmpt)
    {
        return nComp * celli + cmpt;
    };

    // Diagonal and source
    for (label celli = 0; celli < sizeA; ++celli)
    {
        for (direction cmpt = 0; cmpt < nComp; ++cmpt)
        {
            const label ii = idx(celli, cmpt);

            tripletList.emplace_back
            (
                ii,
                ii,
                foam_matrix.diag()[celli]
            );

            b(ii) = foam_matrix.source()[celli][cmpt];
        }
    }

    // Internal ldu coefficients
    const lduAddressing& addr = foam_matrix.lduAddr();
    const labelList& lowerAddr = addr.lowerAddr();
    const labelList& upperAddr = addr.upperAddr();

    forAll(lowerAddr, facei)
    {
        const label lowerCell = lowerAddr[facei];
        const label upperCell = upperAddr[facei];

        for (direction cmpt = 0; cmpt < nComp; ++cmpt)
        {
            const label lowerI = idx(lowerCell, cmpt);
            const label upperI = idx(upperCell, cmpt);

            tripletList.emplace_back
            (
                lowerI,
                upperI,
                foam_matrix.upper()[facei]
            );

            tripletList.emplace_back
            (
                upperI,
                lowerI,
                foam_matrix.lower()[facei]
            );
        }
    }

    // Boundary contributions
    forAll(foam_matrix.psi().boundaryField(), patchI)
    {
        const fvPatch& patch =
            foam_matrix.psi().boundaryField()[patchI].patch();

        forAll(patch, faceI)
        {
            const label celli = patch.faceCells()[faceI];

            for (direction cmpt = 0; cmpt < nComp; ++cmpt)
            {
                const label ii = idx(celli, cmpt);

                tripletList.emplace_back
                (
                    ii,
                    ii,
                    foam_matrix.internalCoeffs()[patchI][faceI][cmpt]
                );

                b(ii) +=
                    foam_matrix.boundaryCoeffs()[patchI][faceI][cmpt];
            }
        }
    }

    A.setFromTriplets(tripletList.begin(), tripletList.end());
}

template <>
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

template <>
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

template <>
void Foam2Eigen::fvMatrix2EigenM
(
    fvMatrix<vector>& foam_matrix,
    Eigen::MatrixXd& A
)
{
    const label sizeA = foam_matrix.diag().size();
    const label nComp = vector::nComponents;

    A.setZero(sizeA * nComp, sizeA * nComp);

    auto idx = [nComp](label celli, direction cmpt)
    {
        return nComp * celli + cmpt;
    };

    // Diagonal
    for (label celli = 0; celli < sizeA; ++celli)
    {
        for (direction cmpt = 0; cmpt < nComp; ++cmpt)
        {
            const label ii = idx(celli, cmpt);
            A(ii, ii) = foam_matrix.diag()[celli];
        }
    }

    // Internal ldu coefficients
    const lduAddressing& addr = foam_matrix.lduAddr();
    const labelList& lowerAddr = addr.lowerAddr();
    const labelList& upperAddr = addr.upperAddr();

    forAll(lowerAddr, facei)
    {
        const label lowerCell = lowerAddr[facei];
        const label upperCell = upperAddr[facei];

        for (direction cmpt = 0; cmpt < nComp; ++cmpt)
        {
            const label lowerI = idx(lowerCell, cmpt);
            const label upperI = idx(upperCell, cmpt);

            A(lowerI, upperI) = foam_matrix.upper()[facei];
            A(upperI, lowerI) = foam_matrix.lower()[facei];
        }
    }

    // Boundary diagonal contributions
    forAll(foam_matrix.psi().boundaryField(), patchI)
    {
        const fvPatch& patch =
            foam_matrix.psi().boundaryField()[patchI].patch();

        forAll(patch, faceI)
        {
            const label celli = patch.faceCells()[faceI];

            for (direction cmpt = 0; cmpt < nComp; ++cmpt)
            {
                const label ii = idx(celli, cmpt);

                A(ii, ii) +=
                    foam_matrix.internalCoeffs()[patchI][faceI][cmpt];
            }
        }
    }
}

template <>
void Foam2Eigen::fvMatrix2EigenM
(
    fvMatrix<vector>& foam_matrix,
    Eigen::SparseMatrix<double>& A
)
{
    const label sizeA = foam_matrix.diag().size();
    const label nComp = vector::nComponents;
    const label nel =
        foam_matrix.diag().size()
      + foam_matrix.upper().size()
      + foam_matrix.lower().size();

    A.resize(sizeA * nComp, sizeA * nComp);
    A.reserve(nel * nComp);

    typedef Eigen::Triplet<double> Trip;
    std::vector<Trip> tripletList;
    tripletList.reserve(nel * nComp);

    auto idx = [nComp](label celli, direction cmpt)
    {
        return nComp * celli + cmpt;
    };

    // Diagonal
    for (label celli = 0; celli < sizeA; ++celli)
    {
        for (direction cmpt = 0; cmpt < nComp; ++cmpt)
        {
            const label ii = idx(celli, cmpt);

            tripletList.emplace_back
            (
                ii,
                ii,
                foam_matrix.diag()[celli]
            );
        }
    }

    // Internal ldu coefficients
    const lduAddressing& addr = foam_matrix.lduAddr();
    const labelList& lowerAddr = addr.lowerAddr();
    const labelList& upperAddr = addr.upperAddr();

    forAll(lowerAddr, facei)
    {
        const label lowerCell = lowerAddr[facei];
        const label upperCell = upperAddr[facei];

        for (direction cmpt = 0; cmpt < nComp; ++cmpt)
        {
            const label lowerI = idx(lowerCell, cmpt);
            const label upperI = idx(upperCell, cmpt);

            tripletList.emplace_back
            (
                lowerI,
                upperI,
                foam_matrix.upper()[facei]
            );

            tripletList.emplace_back
            (
                upperI,
                lowerI,
                foam_matrix.lower()[facei]
            );
        }
    }

    // Boundary diagonal contributions
    forAll(foam_matrix.psi().boundaryField(), patchI)
    {
        const fvPatch& patch =
            foam_matrix.psi().boundaryField()[patchI].patch();

        forAll(patch, faceI)
        {
            const label celli = patch.faceCells()[faceI];

            for (direction cmpt = 0; cmpt < nComp; ++cmpt)
            {
                const label ii = idx(celli, cmpt);

                tripletList.emplace_back
                (
                    ii,
                    ii,
                    foam_matrix.internalCoeffs()[patchI][faceI][cmpt]
                );
            }
        }
    }

    A.setFromTriplets(tripletList.begin(), tripletList.end());
}

template <>
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
            b(w, 0) += foam_matrix.boundaryCoeffs()[I][J];
        }
    }
}

template <>
void Foam2Eigen::fvMatrix2EigenV
(
    fvMatrix<vector>& foam_matrix,
    Eigen::VectorXd& b
)
{
    const label sizeA = foam_matrix.diag().size();
    const label nComp = vector::nComponents;

    b.setZero(sizeA * nComp);

    auto idx = [nComp](label celli, direction cmpt)
    {
        return nComp * celli + cmpt;
    };

    // Source
    for (label celli = 0; celli < sizeA; ++celli)
    {
        for (direction cmpt = 0; cmpt < nComp; ++cmpt)
        {
            b(idx(celli, cmpt)) =
                foam_matrix.source()[celli][cmpt];
        }
    }

    // Boundary source contributions
    forAll(foam_matrix.psi().boundaryField(), patchI)
    {
        const fvPatch& patch =
            foam_matrix.psi().boundaryField()[patchI].patch();

        forAll(patch, faceI)
        {
            const label celli = patch.faceCells()[faceI];

            for (direction cmpt = 0; cmpt < nComp; ++cmpt)
            {
                b(idx(celli, cmpt)) +=
                    foam_matrix.boundaryCoeffs()[patchI][faceI][cmpt];
            }
        }
    }
}

template <class Type, template <class> class PatchField, class GeoMesh>
Eigen::VectorXd Foam2Eigen::projectField(
    GeometricField<Type, PatchField, GeoMesh>& field,
    PtrList<GeometricField<Type, PatchField, GeoMesh >>& modes,
    label Nmodes)
{
    Eigen::VectorXd fr;
    Eigen::MatrixXd Eig_Modes = PtrList2Eigen(modes, Nmodes);
    Eigen::VectorXd f = Foam2Eigen::field2Eigen(field);
    Eigen::VectorXd Volumes = field2Eigen(modes[0].mesh());
    Eigen::MatrixXd VolumesN(Eig_Modes.rows(), 1);

    M_Assert
    (
        Eig_Modes.rows() % Volumes.rows() == 0,
        "The number of Eigen field entries must be an integer multiple of the number of cells"
    );

    const label nComp = Eig_Modes.rows() / Volumes.rows();

    for (label celli = 0; celli < Volumes.rows(); ++celli)
    {
        for (label cmpt = 0; cmpt < nComp; ++cmpt)
        {
            VolumesN(nComp * celli + cmpt, 0) = Volumes(celli);
        }
    }

    fr = Eig_Modes.transpose() * (f.cwiseProduct(VolumesN));
    return fr;
}

template <class Type, template <class> class PatchField, class GeoMesh >
std::tuple<Eigen::MatrixXd, Eigen::VectorXd> Foam2Eigen::projectFvMatrix(
    fvMatrix<Type>& matrix,
    PtrList<GeometricField<Type, PatchField, GeoMesh >>& modes, label Nmodes)
{
    Eigen::SparseMatrix<double> A;
    Eigen::MatrixXd Ar;
    Eigen::VectorXd b;
    Eigen::VectorXd br;
    Eigen::MatrixXd Eig_Modes = PtrList2Eigen(modes, Nmodes);
    Foam2Eigen::fvMatrix2Eigen(matrix, A, b);
    Eigen::VectorXd Volumes = field2Eigen(modes[0].mesh());
    Eigen::MatrixXd VolumesN(Eig_Modes.rows(), Nmodes);

    M_Assert
    (
        Eig_Modes.rows() % Volumes.rows() == 0,
        "The number of Eigen mode entries must be an integer multiple of the number of cells"
    );

    const label nComp = Eig_Modes.rows() / Volumes.rows();

    for (label modeI = 0; modeI < Nmodes; ++modeI)
    {
        for (label celli = 0; celli < Volumes.rows(); ++celli)
        {
            for (label cmpt = 0; cmpt < nComp; ++cmpt)
            {
                VolumesN(nComp * celli + cmpt, modeI) = Volumes(celli);
            }
        }
    }

    Ar = Eig_Modes.transpose() * A * Eig_Modes;
    br = Eig_Modes.transpose() * b;
    std::tuple<Eigen::MatrixXd, Eigen::VectorXd> tupla;
    tupla = std::make_tuple(Ar, br);
    return tupla;
}

template <class Type, template <class> class PatchField, class GeoMesh>
Eigen::MatrixXd Foam2Eigen::MassMatrix(
    PtrList<GeometricField<Type, PatchField, GeoMesh >>& modes, label Nmodes)
{
    Eigen::MatrixXd Mr;
    Eigen::MatrixXd Eig_Modes = PtrList2Eigen(modes, Nmodes);
    Eigen::VectorXd Volumes = field2Eigen(modes[0].mesh());
    Eigen::MatrixXd VolumesN(Eig_Modes.rows(), Nmodes);

    M_Assert
    (
        Eig_Modes.rows() % Volumes.rows() == 0,
        "The number of Eigen mode entries must be an integer multiple of the number of cells"
    );

    const label nComp = Eig_Modes.rows() / Volumes.rows();

    for (label modeI = 0; modeI < Nmodes; ++modeI)
    {
        for (label celli = 0; celli < Volumes.rows(); ++celli)
        {
            for (label cmpt = 0; cmpt < nComp; ++cmpt)
            {
                VolumesN(nComp * celli + cmpt, modeI) = Volumes(celli);
            }
        }
    }

    Mr = Eig_Modes.transpose() * (Eig_Modes.cwiseProduct(VolumesN));
    return Mr;
}

template <class Type>
std::tuple<List<Eigen::SparseMatrix<double >>, List<Eigen::VectorXd >>
Foam2Eigen::LFvMatrix2LSM(PtrList<fvMatrix<Type >>& MatrixList)
{
    List<Eigen::SparseMatrix<double >> SM_list;
    List<Eigen::VectorXd> V_list;
    label LSize = MatrixList.size();
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

    std::tuple<List<Eigen::SparseMatrix<double >>, List<Eigen::VectorXd >> tupla;
    tupla = std::make_tuple(SM_list, V_list);
    return tupla;
}

template std::tuple<List<Eigen::SparseMatrix<double >>, List<Eigen::VectorXd >>
Foam2Eigen::LFvMatrix2LSM(PtrList<fvMatrix<scalar >>& MatrixList);
template std::tuple<List<Eigen::SparseMatrix<double >>, List<Eigen::VectorXd >>
Foam2Eigen::LFvMatrix2LSM(PtrList<fvMatrix<vector >>& MatrixList);

template <class type_matrix>
Eigen::Matrix<type_matrix, Eigen::Dynamic, Eigen::Dynamic>
Foam2Eigen::List2EigenMatrix(List<type_matrix> list)
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
Foam2Eigen::List2EigenMatrix(List<int> list);
template Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
Foam2Eigen::List2EigenMatrix(List<double> list);

template <class type_matrix>
List<type_matrix> Foam2Eigen::EigenMatrix2List(
    Eigen::Matrix<type_matrix, Eigen::Dynamic, Eigen::Dynamic> matrix)
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

template List<int> Foam2Eigen::EigenMatrix2List(
    Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> matrix);
template List<double> Foam2Eigen::EigenMatrix2List(
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix);

template <>
Eigen::MatrixXd Foam2Eigen::field2Eigen(const List<vector>& field)
{
    Eigen::MatrixXd out;
    out.resize(label(field.size() * 3), 1);

    for (label l = 0; l < field.size(); l++)
    {
        for (label j = 0; j < 3; j++)
        {
            out(j + l * 3, 0) = field[l][j];
        }
    }

    return out;
}
// There might repetions of functions (Matrix vs Vector). I wouls only use matrices
template <>
Eigen::MatrixXd Foam2Eigen::field2Eigen(const List<scalar>& field)
{
    Eigen::MatrixXd out;
    out.resize(label(field.size()), 1);

    for (label l = 0; l < field.size(); l++)
    {
        out(l, 0) = field[l];
    }

    return out;
}
