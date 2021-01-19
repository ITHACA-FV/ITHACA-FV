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
    EigenFunctions
Description
    Set of routines to manipulate Eigen c++ objects
SourceFiles
    EigenFunctions.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the EigenFunctions class.
/// \dir
/// Directory containing the header and source files for the EigenFunctions class.

#ifndef EigenFunctions_H
#define EigenFunctions_H
#include <mutex>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#include <Eigen/Eigen>
#include <unsupported/Eigen/SparseExtra>
#include <unsupported/Eigen/CXX11/Tensor>
#pragma GCC diagnostic pop
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


/*---------------------------------------------------------------------------*\
  Class EigenFunctions Declaration
\*---------------------------------------------------------------------------*/

/// Class to perform operation on Eigen objects
class EigenFunctions
{
    public:

        //--------------------------------------------------------------------------
        /// @brief      Vectorize a Sparse Matrix
        ///
        /// @param[in]      Matrix  The matrix
        ///
        /// @tparam     T       type of object, i.e. double, float, ....
        ///
        /// @return     The Sparse Vector
        ///
        template <typename T>
        static Eigen::SparseMatrix<T> vectorize(Eigen::SparseMatrix<T>& Matrix);

        //--------------------------------------------------------------------------
        /// @brief      Find the maximum of a sparse Matrix (Useful for DEIM)
        ///
        /// @param[out]      mat      The matrix
        /// @param[out]      ind_row  The ind row
        /// @param[out]      ind_col  The ind col
        ///
        /// @tparam     T        type of object, i.e. double, float, ....
        ///
        /// @return     the maximum
        ///
        template <typename T>
        static T max(Eigen::SparseMatrix<T>& mat, label& ind_row, label& ind_col);

        //--------------------------------------------------------------------------
        /// @brief           Find the minimum of a sparse Matrix (Useful for DEIM)
        ///
        /// @param[out]      mat      The matrix
        /// @param[out]      ind_row  The index for the row
        /// @param[out]      ind_col  The ind col
        ///
        /// @tparam          T        type of object, i.e. double, float, ....
        ///
        /// @return     the minimum
        ///
        template <typename T>
        static T min(Eigen::SparseMatrix<T>& mat, label& ind_row, label& ind_col);

        //--------------------------------------------------------------------------
        /// @brief      sort eigenvalues
        ///
        /// @param[in,out]      eigenvalues   The eigenvalues
        /// @param[in,out]      eigenvectors  The eigenvectors
        ///
        static void sortEigenvalues(Eigen::VectorXd& eigenvalues,
                                    Eigen::MatrixXd& eigenvectors);

        //--------------------------------------------------------------------------
        /// @brief      Perform Frobenius inner Product between two list of sparse matrices A and B
        ///
        /// \f[ \mathbf{A} = [\mathbf{A_1}, \dots, \mathbf{A_N}] , \mathbf{B} = [\mathbf{B_1}, \dots, \mathbf{B_M}] , out_{ij} = tr(\mathbf{A_i}^T \mathbf{B_j} ) \f]
        ///
        /// @param[in]      A     List of Matrices A
        /// @param[in]      B     List of Matrices B
        ///
        /// @tparam     T     type of object, i.e. double, float, ....
        ///
        /// @return     out Matrix containing the inner products between the matrix lists
        ///
        template <typename T>
        static Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> innerProduct(
            List <Eigen::SparseMatrix<T>>& A, List <Eigen::SparseMatrix<T>>& B);

        //--------------------------------------------------------------------------
        /// @brief      Perform Frobenius inner Product between a list of sparse matrices A and a sparse matrix B
        ///
        ///   \f[\mathbf{A} = [\mathbf{A_1}, \dots, \mathbf{A_N}] , out_i = tr(\mathbf{A_i}^T \mathbf{B} ) \f]
        ///
        /// @param[in]      A     List of Matrices A
        /// @param[in]      B     Sparse Matrix B
        ///
        /// @tparam     T     type of object, i.e. double, float, ....
        ///
        /// @return     out Vector containing the inner products between the matrix lists
        ///
        template <typename T>
        static Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> innerProduct(
            List <Eigen::SparseMatrix<T>>& A, Eigen::SparseMatrix<T>& B);

        //--------------------------------------------------------------------------
        /// @brief      Perform Frobenius inner Product between two sparse matrices A and B
        ///
        ///   \f[out = tr(\mathbf{A}^T \mathbf{B} ) \f]
        ///
        ///
        /// @param[in]      A     Sparse Matrix A
        /// @param[in]      B     Sparse Matrix B
        ///
        /// @tparam     T     type of object, i.e. double, float, ....
        ///
        /// @return     out Value of the inner products between the two matrices
        ///
        template <typename T>
        static T innerProduct(Eigen::SparseMatrix<T>& A, Eigen::SparseMatrix<T>& B);

        //--------------------------------------------------------------------------
        /// @brief      Sparse Matrix-Vector product between a list of sparse matrices and a vector of coefficients
        ///
        ///   \f[ \mathbf{A} = [\mathbf{A_1}, \dots, \mathbf{A_N}] , \mathbf{c} = [c_1, ...., c_N] \f]
        ///
        ///   \f[ \mathbf{out} = \sum_{k=1}^N \mathbf{A_k} c_k  \f]
        ///
        /// @param[in]      A     Sparse Matrix A
        /// @param[in]      C     Dense Vector C
        ///
        /// @tparam     T     type of object, i.e. double, float, ....
        ///
        /// @return     out Sparse Matrix containing the sum of all the matrices multiplied by the vector coefficients
        ///
        template <typename T>
        static Eigen::SparseMatrix<T> MVproduct(List<Eigen::SparseMatrix<T>>& A,
                                                Eigen::DenseBase<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& C);

        //--------------------------------------------------------------------------
        /// @brief      Dense Matrix-Vector product between a list of sparse matrices and a vector of coefficients
        ///
        ///   \f[ \mathbf{A} = [\mathbf{A_1}, \dots, \mathbf{A_N}] , \mathbf{c} = [c_1, ...., c_N] \f]
        ///
        ///   \f[ \mathbf{out} = \sum_{k=1}^N \mathbf{A_k} c_k  \f]
        ///
        /// @param[in]  A     Dense Matrix A
        /// @param[in]  C     Dense Vector C
        ///
        /// @tparam     T     type of object, i.e. double, float, ....
        ///
        /// @return     Dense Matrix containing the sum of all the matrices multiplied by the vector coefficients
        ///
        template <typename T>
        static Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> MVproduct(
            const std::vector< Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& A,
            Eigen::DenseBase<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& C);

        //--------------------------------------------------------------------------
        /// @brief      Sparse Matrix - Dense Matrix product between a list of sparse matrices and a dense matrix
        ///
        ///   \f[ \mathbf{A} = [\mathbf{A_1}, \dots, \mathbf{A_N}] , \mathbf{C} = [\mathbf{c_{\bullet 1}}, ...., \mathbf{c_{\bullet N}}], \mathbf{out} = [\mathbf{O_1},\dots,\mathbf{O_N}] \f]
        ///
        ///   \f[ \mathbf{O_i} = \sum_{k=1}^N \mathbf{A_k} c_{ki}  \f]
        ///
        /// @param[in]      A     List of sparse Matrices
        /// @param[in]      C     Dense Matrix
        ///
        /// @tparam     T     type of object, i.e. double, float, ....
        ///
        /// @return     List of Sparse Matrices containing the sum of all the matrices multiplied by the matrix coefficients by columns
        ///
        template <typename T>
        static List<Eigen::SparseMatrix<T>> MMproduct(List<Eigen::SparseMatrix<T>>& A,
                                         Eigen::DenseBase<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& C);

        //--------------------------------------------------------------------------
        /// @brief      Conditioning number of a dense matrix
        ///
        /// @param[in]  A     input Matrix
        ///
        /// @tparam     T     type of object, i.e. double, float, ....
        ///
        /// @return     Conditioning number
        ///
        template <typename T>
        static T condNumber(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& A);



        //--------------------------------------------------------------------------
        /// @brief      A function that computes the product of  g.T c a, where c is a third dim tensor
        ///
        /// @param[in]  g     The first vector
        /// @param[in]  c     The three dim tensor
        /// @param[in]  a     The second vector
        ///
        /// @return     The vector which corresponds to the product g.T c a
        ///
        template <typename T>
        static Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> vectorTensorProduct(
            const Eigen::Matrix<T, Eigen::Dynamic, 1>& g,
            const Eigen::Tensor<T, 3 >& c,
            const Eigen::Matrix<T, Eigen::Dynamic, 1>& a);

};

template <typename T>
Eigen::SparseMatrix<T> EigenFunctions::vectorize(Eigen::SparseMatrix<T>& Matrix)
{
    typedef Eigen::Triplet<double> Trip;
    std::vector<Trip> tripletList;
    tripletList.reserve(Matrix.nonZeros());
    Eigen::SparseMatrix<T> vec;
    vec.resize(1, Matrix.size());
    vec.reserve(Matrix.nonZeros());

    for (label k = 0; k < Matrix.outerSize(); ++k)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(Matrix, k); it; ++it)
        {
            tripletList.push_back(Trip(0, it.row() + it.col()*Matrix.rows(), it.value()));
        }
    }

    vec.setFromTriplets(tripletList.begin(), tripletList.end());
}

template <typename T>
T EigenFunctions::max(Eigen::SparseMatrix<T>& mat, label& ind_row,
                      label& ind_col)
{
    label i = 0;
    T max = 0;

    for (label k = 0; k < mat.outerSize(); ++k)
    {
        for (typename Eigen::SparseMatrix<T>::InnerIterator it(mat, k); it; ++it)
        {
            if (max < it.value() || i == 0)
            {
                max = it.value();
                ind_row = it.row();
                ind_col = it.col();
                i++;
            }
        }
    }

    return max;
}

template <typename T>
T EigenFunctions::min(Eigen::SparseMatrix<T>& mat, label& ind_row,
                      label& ind_col)
{
    label i = 0;
    T min = 0;

    for (label k = 0; k < mat.outerSize(); ++k)
    {
        for (typename Eigen::SparseMatrix<T>::InnerIterator it(mat, k); it; ++it)
        {
            if (min > it.value() || i == 0)
            {
                min = it.value();
                ind_row = it.row();
                ind_col = it.col();
                i++;
            }
        }
    }

    return min;
}

template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> EigenFunctions::innerProduct(
    List <Eigen::SparseMatrix<T>>& A, List <Eigen::SparseMatrix<T>>& B)
{
    label rows = A.size();
    label cols = B.size();
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> out;
    out.resize(rows, cols);

    for (label i = 0; i < rows; i++)
    {
        for (label j = 0; j < cols; j++)
        {
            out(i, j) = innerProduct(A[i], B[j]);
        }
    }

    return out;
}

template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> EigenFunctions::innerProduct(
    List<Eigen::SparseMatrix<T>>& A, Eigen::SparseMatrix<T>& B)
{
    label rows = A.size();
    label cols = 1;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> out;
    out.resize(rows, cols);

    for (label i = 0; i < rows; i++)
    {
        out(i, 0) = innerProduct(A[i], B);
    }

    return out;
}

template <typename T>
T EigenFunctions::innerProduct(Eigen::SparseMatrix<T>& A,
                               Eigen::SparseMatrix<T>& B)
{
    T res = 0;

    for (label k = 0; k < A.cols(); k++)
    {
        res += A.col(k).dot(B.col(k));
    }

    return res;
}

template <typename T>
Eigen::SparseMatrix<T> EigenFunctions::MVproduct(List<Eigen::SparseMatrix<T>>& A
        , Eigen::DenseBase<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& C)
{
    Eigen::SparseMatrix<T> out;
    out = A[0] * C(0);

    for (label i = 1; i < A.size(); i++)
    {
        out += A[i] * C(i);
    }

    return out;
}

template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> EigenFunctions::MVproduct(
    const std::vector< Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& A,
    Eigen::DenseBase<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& C)
{
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> out;
    out = A[0] * C(0);

    for (label i = 1; i < A.size(); i++)
    {
        out += A[i] * C(i);
    }

    return out;
}



template <typename T>
List<Eigen::SparseMatrix<T>> EigenFunctions::MMproduct(
                              List<Eigen::SparseMatrix<T>>& A,
                              Eigen::DenseBase<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& C)
{
    List<Eigen::SparseMatrix<T>> out;
    out.resize(C.cols());
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> col;

    for (label i = 0; i < C.cols(); i++)
    {
        col = C.col(i);
        out[i] = MVproduct(A, col);
    }

    return out;
}

template <typename T>
T EigenFunctions::condNumber(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>&
                             A)
{
    Eigen::JacobiSVD<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> svd(A);
    T cond = svd.singularValues()(0) / svd.singularValues()(
                 svd.singularValues().size() - 1);
    return cond;
}

// Additional Eigen Functions
namespace Eigen
{
using namespace Eigen;

template<typename VectorType>
bool saveMarketVector (const VectorType& vec, const std::string& filename,
                       label prec, std::_Ios_Fmtflags outytpe = std::ios_base::scientific)
{
    typedef typename VectorType::Scalar Scalar;
    std::ofstream out(filename.c_str(), std::ios::out);

    if (!out)
    {
        return false;
    }

    out.flags(outytpe);
    out.precision(prec);

    if (internal::is_same<Scalar, std::complex<float>>::value
            || internal::is_same<Scalar, std::complex<double>>::value)
    {
        out << "%%MatrixMarket matrix array complex general\n";
    }
    else
    {
        out << "%%MatrixMarket matrix array real general\n";
    }

    out << vec.size() << " " << 1 << "\n";

    for (label i = 0; i < vec.size(); i++)
    {
        internal::putVectorElt(vec(i), out);
    }

    out.close();
    return true;
}

template<typename VectorType>
Matrix<VectorType, Dynamic, Dynamic> SliceFromTensor(
    Eigen::Tensor<VectorType, 3>& tensor, label dim, label index1)
{
    Eigen::Tensor<VectorType, 2> t2 = tensor.chip(index1, dim);
    Matrix<VectorType, Dynamic, Dynamic> m = Eigen::Map<Eigen::Matrix<
            VectorType,           /* scalar element type */
            Eigen::Dynamic,  /* num_rows is a run-time value */
            Eigen::Dynamic,  /* num_cols is a run-time value */
            Eigen::ColMajor  /* tensorflow::Tensor is always row-major */>>(
                t2.data(),  /* ptr to data */
                t2.dimensions()[0],           /* num_rows */
                t2.dimensions()[1]            /* num_cols */);
    return m;
}

template<typename VectorType>
Matrix<VectorType, Dynamic, Dynamic> SliceFromTensor(
    const Eigen::Tensor<VectorType, 3>& tensor, label dim, label index1)
{
    Eigen::Tensor<VectorType, 2> t2 = tensor.chip(index1, dim);
    Matrix<VectorType, Dynamic, Dynamic> m = Eigen::Map<Eigen::Matrix<
            VectorType,           /* scalar element type */
            Eigen::Dynamic,  /* num_rows is a run-time value */
            Eigen::Dynamic,  /* num_cols is a run-time value */
            Eigen::ColMajor  /* tensorflow::Tensor is always row-major */>>(
                t2.data(),  /* ptr to data */
                t2.dimensions()[0],           /* num_rows */
                t2.dimensions()[1]            /* num_cols */);
    return m;
}

}





#endif


