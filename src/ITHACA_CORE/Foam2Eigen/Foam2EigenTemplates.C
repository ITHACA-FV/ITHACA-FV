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
/// Source file of the foam2eigen class.

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<>
Eigen::VectorXd Foam2Eigen::field2Eigen(volVectorField& field)
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
Eigen::VectorXd Foam2Eigen::field2Eigen(volScalarField& field)
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
Eigen::VectorXd Foam2Eigen::field2Eigen(fvMesh const& field)
{
    Eigen::VectorXd out;
    out.resize(int(field.V().size()));

    for (int l = 0; l < field.V().size(); l++)
    {
        out(l) = field.V()[l];
    }

    return out;
}

template<>
volVectorField Foam2Eigen::Eigen2field(volVectorField& field_in,
                                       Eigen::VectorXd& eigen_vector)
{
    volVectorField field_out(field_in);

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
volScalarField Foam2Eigen::Eigen2field(volScalarField& field_in,
                                       Eigen::VectorXd& eigen_vector)
{
    volScalarField field_out(field_in);

    for (auto i = 0; i < field_out.size(); i++)
    {
        field_out.ref()[i] = eigen_vector(i);
    }

    field_out.correctBoundaryConditions();
    return field_out;
}




template<>
Eigen::MatrixXd Foam2Eigen::PtrList2Eigen(PtrList<volVectorField>& fields,
        int Nfields)
{
    int Nf;

    if (Nfields > fields.size())
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
        out.col(k) = field2Eigen<volVectorField>(fields[k]);
    }

    return out;
}

template<>
Eigen::MatrixXd Foam2Eigen::PtrList2Eigen(PtrList<volScalarField>& fields,
        int Nfields)
{
    int Nf;

    if (Nfields > fields.size())
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
        out.col(k) = field2Eigen<volScalarField>(fields[k]);
    }

    return out;
}


template<>
void Foam2Eigen::fvMatrix2Eigen(volScalarField& T, fvScalarMatrix& foam_matrix,
                                Eigen::MatrixXd& A, Eigen::VectorXd& b)
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
    forAll(T.boundaryField(), I)
    {
        const fvPatch& ptch = T.boundaryField()[I].patch();
        forAll(ptch, J)
        {
            int w = ptch.faceCells()[J];
            A(w, w) += foam_matrix.internalCoeffs()[I][J];
            b(w, 0)   += foam_matrix.boundaryCoeffs()[I][J];
        }
    }
}

template<>
void Foam2Eigen::fvMatrix2Eigen(volScalarField& T, fvScalarMatrix& foam_matrix,
                                Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b)
{
    int sizeA = foam_matrix.diag().size();
    int nel = foam_matrix.diag().size() + foam_matrix.upper().size() +
              foam_matrix.lower().size();;
    A.resize(sizeA, sizeA);
    b.resize(sizeA);
    A.reserve(nel);

    for (auto i = 0; i < sizeA; i++)
    {
        A.insert(i, i) = foam_matrix.diag()[i];
        b(i, 0) = foam_matrix.source()[i];
    }

    const lduAddressing& addr = foam_matrix.lduAddr();
    const labelList& lowerAddr = addr.lowerAddr();
    const labelList& upperAddr = addr.upperAddr();
    forAll(lowerAddr, i)
    {
        A.insert(lowerAddr[i], upperAddr[i]) = foam_matrix.upper()[i];
        A.insert(upperAddr[i], lowerAddr[i]) = foam_matrix.lower()[i];
    }
    forAll(T.boundaryField(), I)
    {
        const fvPatch& ptch = T.boundaryField()[I].patch();
        forAll(ptch, J)
        {
            int w = ptch.faceCells()[J];
            A.coeffRef(w, w) += foam_matrix.internalCoeffs()[I][J];
            b(w, 0)   += foam_matrix.boundaryCoeffs()[I][J];
        }
    }
}

template<>
void Foam2Eigen::fvMatrix2Eigen(volVectorField& T, fvVectorMatrix& foam_matrix,
                                Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b)
{
    int sizeA = foam_matrix.diag().size();
    int nel = foam_matrix.diag().size() + foam_matrix.upper().size() +
              foam_matrix.lower().size();
    A.resize(sizeA * 3, sizeA * 3);
    A.reserve(nel * 3);
    b.resize(sizeA * 3);

    for (auto i = 0; i < sizeA; i++)
    {
        A.insert(i, i) = foam_matrix.diag()[i];
        A.insert(sizeA + i, sizeA + i) = foam_matrix.diag()[i];
        A.insert(2 * sizeA + i, 2 * sizeA + i) = foam_matrix.diag()[i];
        b(i) = foam_matrix.source()[i][0];
        b(sizeA + i) = foam_matrix.source()[i][1];
        b(2 * sizeA + i) = foam_matrix.source()[i][2];
    }

    const lduAddressing& addr = foam_matrix.lduAddr();
    const labelList& lowerAddr = addr.lowerAddr();
    const labelList& upperAddr = addr.upperAddr();
    forAll(lowerAddr, i)
    {
        A.insert(lowerAddr[i], upperAddr[i]) = foam_matrix.upper()[i];
        A.insert(lowerAddr[i] + sizeA, upperAddr[i] + sizeA) = foam_matrix.upper()[i];
        A.insert(lowerAddr[i] + sizeA * 2,
                 upperAddr[i] + sizeA * 2) = foam_matrix.upper()[i];
        A.insert(upperAddr[i], lowerAddr[i]) = foam_matrix.lower()[i];
        A.insert(upperAddr[i] + sizeA, lowerAddr[i] + sizeA) = foam_matrix.lower()[i];
        A.insert(upperAddr[i] + sizeA * 2,
                 lowerAddr[i] + sizeA * 2) = foam_matrix.lower()[i];
    }
    forAll(T.boundaryField(), I)
    {
        const fvPatch& ptch = T.boundaryField()[I].patch();
        forAll(ptch, J)
        {
            int w = ptch.faceCells()[J];
            A.coeffRef(w, w) += foam_matrix.internalCoeffs()[I][J][0];
            A.coeffRef(w + sizeA, w + sizeA) += foam_matrix.internalCoeffs()[I][J][1];
            A.coeffRef(w + sizeA * 2,
                       w + sizeA * 2) += foam_matrix.internalCoeffs()[I][J][2];
            b(w)   += foam_matrix.boundaryCoeffs()[I][J][0];
            b(w + sizeA)   += foam_matrix.boundaryCoeffs()[I][J][1];
            b(w + sizeA * 2)   += foam_matrix.boundaryCoeffs()[I][J][2];
        }
    }
}

template <class type_m, class type_PtrList>
std::tuple<Eigen::MatrixXd, Eigen::VectorXd> Foam2Eigen::projectFvMatrix(
    type_m& matrix, type_PtrList& modes, int Nmodes)
{
    Eigen::SparseMatrix<double> A;
    Eigen::MatrixXd Ar;
    Eigen::VectorXd b;
    Eigen::VectorXd br;
    Eigen::MatrixXd Eig_Modes = PtrList2Eigen(modes, Nmodes);
    Foam2Eigen::fvMatrix2Eigen(modes[0], matrix, A, b);
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

    /*    Ar = Eig_Modes.transpose() * A * (Eig_Modes.cwiseProduct(VolumesN));
        br = Eig_Modes.transpose() * (b.cwiseProduct(Volumes));*/
    Ar = Eig_Modes.transpose() * A * Eig_Modes;
    br = Eig_Modes.transpose() * b;
    std::tuple <Eigen::MatrixXd, Eigen::VectorXd> tupla;
    tupla = std::make_tuple(Ar, br);
    return tupla;
}

template <class type_PtrList>
Eigen::MatrixXd Foam2Eigen::MassMatrix(type_PtrList& modes, int Nmodes)
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



template <class type_f, class type_PtrList>
Eigen::VectorXd Foam2Eigen::projectField(type_f& field, type_PtrList& modes,
        int Nmodes)
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




