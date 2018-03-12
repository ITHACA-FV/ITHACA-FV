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
    foam2eigen

Description
    Converts foam objects into Eigen and viceversa

SourceFiles
    foam2eigen.C

\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the ITHACAutilities class.
/// \dir
/// Directory containing the header and source files for the ITHACAutilities class.

#ifndef Foam2Eigen_H
#define Foam2Eigen_H

#include "fvCFD.H"
#include "IOmanip.H"
#include <tuple>
#include <sys/stat.h>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#include "../thirdparty/Eigen/Eigen/Eigen"
#include "../thirdparty/Eigen/Eigen/LU"
#pragma GCC diagnostic pop
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


/*---------------------------------------------------------------------------*\
  Class foam2eigen Declaration
  \*---------------------------------------------------------------------------*/

/// Class to converts foam objects into Eigen and viceversa
class Foam2Eigen
{
private:

public:

    //--------------------------------------------------------------------------
    /// @brief      Convert a ldu OpenFOAM matrix into a Eigen Matrix A and a source vector b
    ///
    /// @param[in]      foam_matrix       The foam matrix can be fvScalarMatrix or fvVectorMatrix
    /// @param[out]     A                 The matrix or the list of matrices, can be dense or sparse
    /// @param[out]     b                 The source term vector or the list of source term vectors, always dense
    ///
    /// @tparam         type_f            The type of field, can be a volvectorfield or a volscalarfield
    /// @tparam         type_foam_matrix  The type of foam matrix can be fvScalarMatrix or fvVectorMatrix
    /// @tparam         type_A            The type of matrix or the list of matrices, can be dense (Eigen::MatrixXd) or sparse (Eigen::SparseMatrix<double>) or a List of them
    /// @tparam         type_B            The type source term vector or the list can be Eigen::VectorXd or List<Eigen::VectorXd>
    ///
    template <class type_foam_matrix, class type_A, class type_B>
    static void fvMatrix2Eigen(type_foam_matrix& foam_matrix, type_A& A, type_B& b);


    //--------------------------------------------------------------------------
    /// @brief      Convert a PtrList List of snapshots to eigen matrix (only internalfield)
    ///
    /// @param      fields  The fields can be a PtrList<volScalarField> or PtrList<volVectorField>
    ///
    /// @tparam     type_f  Type of the fields can be PtrList<volScalarField> or PtrList<volVectorField>
    ///
    /// @return     A Eigen matrix containing as columns the snapshots
    ///
    ///
    ///
    template <class type_f>
    static Eigen::MatrixXd PtrList2Eigen(type_f& fields, int Nfields = 10000000);


    //--------------------------------------------------------------------------
    /// @brief      { function_description }
    ///
    /// @param      field   The field
    ///
    /// @tparam     type_f  { description }
    ///
    /// @return     { description_of_the_return_value }
    ///
    template <class type_f>
    static Eigen::VectorXd field2Eigen(type_f& field);

    template <class type_f>
    static type_f Eigen2field(type_f& field, Eigen::VectorXd& eigen_vector);

    template <class type_m, class type_PtrList>
    static std::tuple<Eigen::MatrixXd, Eigen::VectorXd> projectFvMatrix(type_m& matrix, type_PtrList& modes, int Nmodes);

    template <class type_f, class type_PtrList>
    static Eigen::VectorXd projectField(type_f& field, type_PtrList& modes, int Nmodes);

    template <class type_PtrList>
    static Eigen::MatrixXd MassMatrix(type_PtrList& modes, int Nmodes);

    template<class type_matrix>
    static std::tuple<List<Eigen::SparseMatrix<double> >, List<Eigen::VectorXd> > LFvMatrix2LSM(PtrList<type_matrix>& MatrixList);


};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//
//
//

template <class type_f>
Eigen::VectorXd Foam2Eigen::field2Eigen(type_f& field)
{
    Eigen::VectorXd out;
    out.resize(int(field.size()));

    for (int l = 0; l < field.size(); l++)
    {
        out(l) = field[l];
    }
    return out;
}


template<class type_matrix>
std::tuple<List<Eigen::SparseMatrix<double> >, List<Eigen::VectorXd> > Foam2Eigen::LFvMatrix2LSM(PtrList<type_matrix>& MatrixList)
{
    List<Eigen::SparseMatrix<double> > SM_list;
    List<Eigen::VectorXd> V_list;

    int LSize =  MatrixList.size();

    SM_list.resize(LSize);
    V_list.resize(LSize);
    Eigen::SparseMatrix<double> A;
    Eigen::VectorXd b;

    for (int j = 0; j < LSize; j++)
    {
        fvMatrix2Eigen(MatrixList[j], A, b);
        SM_list[j] = A;
        V_list[j] = b;
    }

    std::tuple <List<Eigen::SparseMatrix<double> >, List<Eigen::VectorXd> > tupla;
    tupla = std::make_tuple(SM_list, V_list);

    return tupla;
}

template <class type_m, class type_PtrList>
std::tuple<Eigen::MatrixXd, Eigen::VectorXd> Foam2Eigen::projectFvMatrix(type_m& matrix, type_PtrList& modes, int Nmodes)
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
Eigen::VectorXd Foam2Eigen::projectField(type_f& field, type_PtrList& modes, int Nmodes)
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

#endif


