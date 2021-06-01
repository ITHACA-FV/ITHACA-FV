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
    ITHACAstream
Description
    Methods for input and outputs from and to files
SourceFiles
    ITHACAstream.C
    ITHACAstreamTemplates.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the ITHACAstream class, it contains the implementation of
/// several methods for input output operations.
/// \dir
/// Directory containing the header, source and template files for the ITHACAstream class, ITHACAparameters and Assert Functions.

#ifndef ITHACAstream_H
#define ITHACAstream_H

#include "fvCFD.H"
#include "IOmanip.H"
#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>
#include <algorithm>
#include <fstream>
#include <string>
#include <stdexcept>
#include <sstream>
#include <vector>
#include <cstdio>
#include <typeinfo>
#include <iostream>
#include <cassert>
#include <zlib.h>
#include <map>
#include "ITHACAassert.H"
#include "ITHACAparameters.H"
#include "ITHACAutilities.H"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wignored-attributes"
#include <Eigen/Eigen>
#include "EigenFunctions.H"

#include <unsupported/Eigen/CXX11/Tensor>
#pragma GCC diagnostic pop

#define MAXBUFSIZE (static_cast<int> (1e6))
#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


/*---------------------------------------------------------------------------*\
                    Namespace ITHACAstream Declaration
\*---------------------------------------------------------------------------*/

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#include "cnpy.H"
#pragma GCC diagnostic pop

/// Namespace for input-output manipulation
namespace ITHACAstream
{

//--------------------------------------------------------------------------
/// @brief      Export an Fv Matrix to folder together with its source term
///
/// @param      Matrix      The matrix
/// @param[in]  folder      The folder
/// @param[in]  MatrixName  The fieldname
///
/// @tparam     Type        can be scalar or vector
///
template<typename Type>
void exportFvMatrix(fvMatrix<Type>& Matrix, word folder,
                    word MatrixName);

//--------------------------------------------------------------------------
/// Export the reduced matrices in numpy (type=python), matlab (type=matlab) and txt (type=eigen) format
/// The function is implemented for a second order matrix
///
/// @param[in] matrix    Eigen::MatrixXd that you want to export.
/// @param[in] Name      Name of saved file.
/// @param[in] type      Format of the file: numpy (type="python"), matlab (type="matlab"), txt (type="eigen").
/// @param[in] folder    Directory in which the fali is saved.
///
/// @tparam    T         The matrix can be of integers, double or float
///
template<typename T, int dim>
void exportMatrix(Eigen::Matrix < T, -1, dim > & matrix, word Name,
                  word type = "python", word folder = "./Matrices");

//--------------------------------------------------------------------------
/// Export the reduced matrices in numpy (type=python), matlab (type=matlab) and txt (type=eigen) format
/// The function is implemented for a third order matrix
/// The eigen matrix is stored as a set of two dimensional matrices where
/// each matrix is in a different file
///
/// @param[in] matrix    Eigen::MatrixXd that you want to export.
/// @param[in] Name      Name of saved file.
/// @param[in] type      Format of the file: numpy (type="python"), matlab (type="matlab"), txt (type="eigen").
/// @param[in] folder    Directory in which the file is saved.
///
void exportMatrix(List <Eigen::MatrixXd>& matrix, word Name,
                  word type = "python", word folder = "./Matrices");

//--------------------------------------------------------------------------
/// Export the reduced tensor in numpy (tipo=python), matlab (tipo=matlab) and txt (tipo=eigen) format
///
/// @param[in]  tensor  The tensor that you want to export.
/// @param[in]  Name    Name of saved file.
/// @param[in]  type    Format of the file: numpy (type="python"), matlab (type="matlab"), txt (type="eigen").
/// @param[in]  folder  Directory in which the fali is saved.
/// @param[in]  folder  folder string to identify the folder where you want to save the file.
///
/// @tparam     T       The tensor can be of integers, double or float
///
template<typename T>
void exportTensor(Eigen::Tensor<T, 3 > tensor, word Name,
                  word type = "python", word folder = "./Matrices");


//----------------------------------------------------------------------
/// Function to read a list of fields from the name of the field and
/// casename
///
/// @param[in]  Lfield      a PtrList of OpenFOAM fields where you want
///                         to store the field.
/// @param[in]  Name        The name of the field you want to read.
/// @param[in]  casename    The folder where the field is stored.
/// @param[in]  first_snap  The first snapshots from which you want to
///                         start reading the field.
/// @param[in]  n_snap      The number of snapshots you want to read.
///
/// @tparam     Type        vector or scalar.
/// @tparam     PatchField  fvPatchField or fvsPatchField.
/// @tparam     GeoMesh     volMesh or surfaceMesh.
///
template<class Type, template<class> class PatchField, class GeoMesh>
void read_fields(PtrList<GeometricField<Type, PatchField, GeoMesh>>&
                 Lfield, word Name,
                 fileName casename, int first_snap = 0, int n_snap = 0);



//----------------------------------------------------------------------
/// Function to read a list of fields from the name of the field if it is already existing
///
/// @param[in]  Lfield      a PtrList of OpenFOAM fields where you want
///                         to store the field.
/// @param[in]  field       The field used as template to read other
///                         fields.
/// @param[in]  casename    The folder where the field is stored.
/// @param[in]  first_snap  The first snapshots from which you want to
///                         start reading the field.
/// @param[in]  n_snap      The number of snapshots you want to read.
///
/// @tparam     Type        vector or scalar.
/// @tparam     PatchField  fvPatchField or fvsPatchField.
/// @tparam     GeoMesh     volMesh or surfaceMesh.
///
template<class Type, template<class> class PatchField, class GeoMesh>
void read_fields(PtrList<GeometricField<Type, PatchField, GeoMesh>>&
                 Lfield,
                 GeometricField<Type, PatchField, GeoMesh>& field,
                 fileName casename, int first_snap = 0, int n_snap = 0);

//----------------------------------------------------------------------
/// Funtion to read a list of volVectorField from name of the field including all the intermediate snapshots
///
/// @param[in]  Lfield      a PtrList of OpenFOAM fields where you want
///                         to store the field.
/// @param[in]  field       The field used as template to read other
///                         fields.
/// @param[in]  casename    The folder where the field is stored.
///
/// @tparam     Type        vector or scalar.
/// @tparam     PatchField  fvPatchField or fvsPatchField.
/// @tparam     GeoMesh     volMesh or surfaceMesh.
///
template<class Type, template<class> class PatchField, class GeoMesh>
void readMiddleFields(PtrList<GeometricField<Type, PatchField, GeoMesh>>&
                      Lfield,
                      GeometricField<Type, PatchField, GeoMesh>& field,
                      fileName casename);

//----------------------------------------------------------------------
/// Function to read a list of volVectorField from name of the field including only converged snapshots
///
/// @param[in]  Lfield      a PtrList of OpenFOAM fields where you want
///                         to store the field.
/// @param[in]  field       The field used as template to read other
///                         fields.
/// @param[in]  casename    The folder where the field is stored.
///
/// @tparam     Type        vector or scalar.
/// @tparam     PatchField  fvPatchField or fvsPatchField.
/// @tparam     GeoMesh     volMesh or surfaceMesh.
///
template<class Type, template<class> class PatchField, class GeoMesh>
void readConvergedFields(
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& Lfield,
    GeometricField<Type, PatchField, GeoMesh>& field,
    fileName casename);

//--------------------------------------------------------------------------
/// Function to export a scalar of vector field
///
/// @param[in]  field      The field you want to export.
/// @param[in]  folder     The folder where you want to save the field.
/// @param[in]  fieldname  The name you want to give to the file.
///
/// @tparam     Type        vector or scalar.
/// @tparam     PatchField  fvPatchField or fvsPatchField.
/// @tparam     GeoMesh     volMesh or surfaceMesh.
///
template<class Type, template<class> class PatchField, class GeoMesh>
void exportFields(PtrList<GeometricField<Type, PatchField, GeoMesh>>&
                  field,
                  word folder, word fieldname);

//--------------------------------------------------------------------------
/// Read a two dimensional matrix from a txt file in Eigen format
/* One has to provide the complete filename with the absolute or relative path */
///
/// @param[in]  filename  The filename of the matrix.
///
/// @return     The Eigen::MatrixXd that contains the imported matrix.
///
Eigen::MatrixXd readMatrix(word filename);

//--------------------------------------------------------------------------
/// Read a three dimensional matrix from a txt file in Eigen format
/* One has to provide the folder containing the matrix files and the filename of the
 the matrix. Since it is stored as a List of matrices each matrix must be stored in a
 different file with the following format:
 matFileName(i)_mat.txt */
///
/// @param[in]  folder    The folder where the txt files are located
/// @param[in]  mat_name  The matrix name
///
/// @return     List <Eigen::MatrixXd> that contains the imported matrix.
///
List <Eigen::MatrixXd> readMatrix(word folder, word mat_name);

//--------------------------------------------------------------------------
/// Export a field to file in a certain folder and subfolder
///
/// @param[in] s          Field
/// @param[in] subfolder  Subfolder where the field is stored
/// @param[in] folder     Folder where the field is stored
/// @param[in] fieldName  Name of the field/file
///
/// @tparam     Type        vector or scalar.
/// @tparam     PatchField  fvPatchField or fvsPatchField.
/// @tparam     GeoMesh     volMesh or surfaceMesh.
///
template<class Type, template<class> class PatchField, class GeoMesh>
void exportSolution(GeometricField<Type, PatchField, GeoMesh>& s,
                    fileName subfolder, fileName folder,
                    word fieldName);

//--------------------------------------------------------------------------
/// Export a field to file in a certain folder and subfolder
///
/// @param[in] s          Field
/// @param[in] subfolder  Subfolder where the field is stored
/// @param[in] folder     Folder where the field is stored
///
/// @tparam     Type        vector or scalar.
/// @tparam     PatchField  fvPatchField or fvsPatchField.
/// @tparam     GeoMesh     volMesh or surfaceMesh.
///
template<class Type, template<class> class PatchField, class GeoMesh>
void exportSolution(GeometricField<Type, PatchField, GeoMesh>& s,
                    fileName subfolder, fileName folder);

//--------------------------------------------------------------------------
/// Export a list to file
///
/// @param[in]  list      The list
/// @param[in]  filename  Name of saved file
/// @param[in]  folder    Folder where the field is stored
///
/// @tparam     T         generic list type object
///
template<typename T>
void exportList(T& list, word folder, word filename);


//--------------------------------------------------------------------------
/// Export an Eigen sparse matrix into bynary format file
///
/// @param      m           Eigen::SparseMatrix
/// @param[in]  folder      The folder where you want to save the matrix
/// @param[in]  MatrixName  The matrix name for the output file
///
/// @tparam     T           type of the matrix, i.e. double, float, ...
/// @tparam     Nrows       number of matrix rows
/// @tparam     IND         number of matrix cols
///
template <typename T, int Nrows, typename IND>
void SaveSparseMatrix(Eigen::SparseMatrix<T, Nrows, IND>& m,
                      word folder, word MatrixName);

//--------------------------------------------------------------------------
/// Read an Eigen sparse matrix from a bynary format file
///
/// @param[in,out]  m           Eigen::SparseMatrix m
/// @param[in]      folder      The folder from where you want to read the matrix
/// @param[in]      MatrixName  The matrix name for the input file
///
/// @tparam     T           type of the matrix, i.e. double, float, ...
/// @tparam     Nrows       number of matrix rows
/// @tparam     IND         number of matrix cols
///
template <typename T, int Nrows, typename IND>
void ReadSparseMatrix(Eigen::SparseMatrix<T, Nrows, IND>& m,
                      word folder, word MatrixName);

//--------------------------------------------------------------------------
/// Save a dense matrix to a binary format file
///
/// @param[in]  Matrix      Eigen dense matrix
/// @param[in]  folder      Folder where you want to save the matrix
/// @param[in]  MatrixName  Name of the output file
///
/// @tparam     MatrixType  type of the matrix, i.e. double, float, ...
///
template <typename MatrixType>
void SaveDenseMatrix(MatrixType& Matrix, word folder, word MatrixName);

//--------------------------------------------------------------------------
/// Read a dense matrix from a binary format file
///
/// @param[in,out]  Matrix      Eigen dense matrix
/// @param[in]      folder      Folder from where you want to read the matrix
/// @param[in]      MatrixName  Name of the input file
///
/// @tparam     MatrixType  type of the matrix, i.e. double, float, ...
///
template <typename MatrixType>
void ReadDenseMatrix(MatrixType& Matrix, word folder, word MatrixName);


//----------------------------------------------------------------------
/// Save a dense tensor to file
///
/// @param[in]  Tensor      The tensor
/// @param[in]  folder      Folder where you want to save the tensor
/// @param[in]  MatrixName  The matrix name
///
/// @tparam     TensorType  type of the tensor, i.e. double, float, ...
///
template <typename TensorType>
void SaveDenseTensor(TensorType& Tensor, word folder, word MatrixName);


//----------------------------------------------------------------------
/// Read a dense tensor from file
///
/// @param[in/out] Tensor      The tensor
/// @param[in]     folder      Folder where you want to read the tensor
/// @param[in]     MatrixName  The matrix name
///
/// @tparam        TensorType   type of the tensor, i.e. double, float, ...
///
template <typename TensorType>
void ReadDenseTensor(TensorType& Tensor, word folder, word MatrixName);

//--------------------------------------------------------------------------
/// Save a sparse matrix list to bynary files
///
/// @param[in]  MatrixList  The matrix list
/// @param[in]  folder      The folder where you want to save files
/// @param[in]  MatrixName  The matrix name prefix, files are the named with an incremental number
///
/// @tparam     MatrixType   MatrixType  type of the matrix, i.e. double, float, ...
///
template <typename MatrixType>
void SaveSparseMatrixList(List<MatrixType>& MatrixList, word folder,
                          word MatrixName);

//--------------------------------------------------------------------------
/// Read a sparse matrix list to bynary files
///
/// @param[in,out]  MatrixList  The matrix list
/// @param[in]      folder      The folder from where you want to read the files
/// @param[in]      MatrixName  The matrix name prefix, files are the named with an incremental number
///
/// @tparam     MatrixType  MatrixType  type of the matrix, i.e. double, float, ...
/// @tparam     Nrows       number of matrix rows
/// @tparam     IND         number of matrix cols
///
template <typename T, int Nrows, typename IND>
void ReadSparseMatrixList(List<Eigen::SparseMatrix<T, Nrows, IND>>& m,
                          word folder, word MatrixName);

//--------------------------------------------------------------------------
/// Save a dense matrix list to bynary files
///
/// @param[in]  MatrixList  The matrix list
/// @param[in]  folder      The folder where you want to save files
/// @param[in]  MatrixName  The matrix name prefix, files are the named with an incremental number
///
/// @tparam     MatrixType   MatrixType  type of the matrix, i.e. double, float, ...
///
template <typename MatrixType>
void SaveDenseMatrixList(List<MatrixType>& MatrixList, word folder,
                         word MatrixName);


//--------------------------------------------------------------------------
/// Read a sparse matrix list to bynary files
///
/// @param[in,out]  MatrixList  The matrix list
/// @param[in]      folder      The folder from where you want to read the files
/// @param[in]      MatrixName  The matrix name prefix, files are the named with an incremental number
///
/// @tparam     MatrixType   MatrixType  type of the matrix, i.e. double, float, ...
///
template <typename MatrixType>
void ReadDenseMatrixList(List<MatrixType>& m, word folder,
                         word MatrixName);

//--------------------------------------------------------------------------
/// Count the number of files with a certain name prefix, used by reading functions of list of matrices
///
/// @param[in]  folder      The folder
/// @param[in]  MatrixName  The matrix name
///
/// @return     Number of files matching the request
///
int numberOfFiles(word folder, word MatrixName, word ext = "");

//--------------------------------------------------------------------------
/// Write points of a mesh to a file
///
/// @param[in]  points     The points you want to write
/// @param[in]  folder     The folder where you want to save the file
/// @param[in]  subfolder  The subfolder where you want to save the file
///
void writePoints(pointField points, fileName folder, fileName subfolder);

//--------------------------------------------------------------------------
/// Print progress bar given the percentage
///
/// @param[in]  percentage  The percentage
///
void printProgress(double percentage);

template<typename T>
void save(const List<Eigen::SparseMatrix<T>>& MatrixList, word folder,
          word MatrixName);

template<typename T>
void load(List<Eigen::SparseMatrix<T>>& MatrixList, word folder,
          word MatrixName);

};

namespace Foam
{
template<typename T>
Ostream& operator<< (Ostream& os, const Eigen::Matrix < T, -1, -1 > & mat)
{
    os << mat.rows() << mat.cols() << UList<T>(const_cast<Eigen::MatrixXd&>
            (mat).data(), mat.size());
    return os;
}

template<typename T>
Istream& operator>> (Istream& is, Eigen::Matrix < T, -1, -1 > & mat)
{
    int nrow, ncol;
    is >> nrow >> ncol;
    mat.resize(nrow, ncol);
    UList<T> list(mat.data(), nrow * ncol);
    is >> list;
    return is;
}

template<typename T>
Ostream& operator<< (Ostream& os, const Eigen::Matrix < T, -1, 1 > & mat)
{
    os << mat.rows() << mat.cols() << UList<T>(const_cast<Eigen::MatrixXd&>
            (mat).data(), mat.size());
    return os;
}

template<typename T>
Istream& operator>> (Istream& is, Eigen::Matrix < T, -1, 1 > & mat)
{
    int nrow, ncol;
    is >> nrow >> ncol;
    mat.resize(nrow, ncol);
    UList<T> list(mat.data(), nrow * ncol);
    is >> list;
    return is;
}

template<typename T>
Ostream& operator<< (Ostream& os, const Eigen::Tensor<T, 3>& tens)
{
    os << tens.dimension(0) << tens.dimension(1) << tens.dimension(
           2) << UList<T>(const_cast<Eigen::Tensor<T, 3>&>(tens).data(),
                          tens.size());
    return os;
}

template<typename T>
Istream& operator>> (Istream& is, Eigen::Tensor<T, 3>& tens)
{
    int d1, d2, d3;
    is >> d1 >> d2 >> d3;
    tens.resize(d1, d2, d3);
    UList<T> list(tens.data(), d1 * d2 * d3);
    is >> list;
    return is;
}

}

template <typename T, int whatever, typename IND>
void ITHACAstream::SaveSparseMatrix(Eigen::SparseMatrix<T, whatever, IND>& m,
                                    word folder, word MatrixName)
{
    typedef Eigen::Triplet<int> Trip;
    std::vector<Trip> res;
    m.makeCompressed();
    mkDir(folder);
    std::fstream writeFile;
    writeFile.open(folder + MatrixName, std::ios::binary | std::ios::out);

    if (writeFile.is_open())
    {
        IND rows, cols, nnzs, outS, innS;
        rows = m.rows()     ;
        cols = m.cols()     ;
        nnzs = m.nonZeros() ;
        outS = m.outerSize();
        innS = m.innerSize();
        writeFile.write(reinterpret_cast<const char*> (&rows), sizeof(IND));
        writeFile.write(reinterpret_cast<const char*> (&cols), sizeof(IND));
        writeFile.write(reinterpret_cast<const char*> (&nnzs), sizeof(IND));
        writeFile.write(reinterpret_cast<const char*> (&outS), sizeof(IND));
        writeFile.write(reinterpret_cast<const char*> (&innS), sizeof(IND));
        writeFile.write(reinterpret_cast<const char*>(m.valuePtr()),
                        sizeof(T  ) * m.nonZeros());
        writeFile.write(reinterpret_cast<const char*>(m.outerIndexPtr()),
                        sizeof(IND) * m.outerSize());
        writeFile.write(reinterpret_cast<const char*>(m.innerIndexPtr()),
                        sizeof(IND) * m.nonZeros());
        writeFile.close();
    }
}

template <typename T, int whatever, typename IND>
void ITHACAstream::ReadSparseMatrix(Eigen::SparseMatrix<T, whatever, IND>& m,
                                    word folder, word MatrixName)
{
    std::fstream readFile;
    readFile.open(folder + MatrixName, std::ios::binary | std::ios::in);

    if (!readFile.good())
    {
        std::cout << folder + MatrixName <<
                  " file does not exist, try to rerun the Offline Stage!" << std::endl;
        exit(EXIT_FAILURE);
    }
    else if (readFile.is_open())
    {
        IND rows, cols, nnz, inSz, outSz;
        readFile.read(reinterpret_cast<char*>(&rows ), sizeof(IND));
        readFile.read(reinterpret_cast<char*>(&cols ), sizeof(IND));
        readFile.read(reinterpret_cast<char*>(&nnz  ), sizeof(IND));
        readFile.read(reinterpret_cast<char*>(&inSz ), sizeof(IND));
        readFile.read(reinterpret_cast<char*>(&outSz), sizeof(IND));
        m.resize(rows, cols);
        m.makeCompressed();
        m.resizeNonZeros(nnz);
        readFile.read(reinterpret_cast<char*>(m.valuePtr()), sizeof(T  ) * nnz  );
        readFile.read(reinterpret_cast<char*>(m.outerIndexPtr()), sizeof(IND) * outSz);
        readFile.read(reinterpret_cast<char*>(m.innerIndexPtr()), sizeof(IND) * nnz );
        m.finalize();
        readFile.close();
    }
}

template <typename MatrixType>
void ITHACAstream::SaveDenseMatrix(MatrixType& Matrix, word folder,
                                   word MatrixName)
{
    mkDir(folder);
    std::ofstream out(folder + MatrixName,
                      std::ios::out | std::ios::binary | std::ios::trunc);
    typename MatrixType::Index rows = Matrix.rows(), cols = Matrix.cols();
    out.write(reinterpret_cast<char*> (&rows), sizeof(typename MatrixType::Index));
    out.write(reinterpret_cast<char*> (&cols), sizeof(typename MatrixType::Index));
    out.write(reinterpret_cast<char*> (Matrix.data()),
              rows * cols * sizeof(typename MatrixType::Scalar) );
    out.close();
}

template <typename MatrixType>
void ITHACAstream::ReadDenseMatrix(MatrixType& Matrix, word folder,
                                   word MatrixName)
{
    std::ifstream in;
    in.open(folder + MatrixName, std::ios::in | std::ios::binary);

    if (!in.good())
    {
        std::cout << folder + MatrixName <<
                  " file does not exist, try to rerun the Offline Stage!" << std::endl;
        exit(EXIT_FAILURE);
    }
    else if (in.is_open())
    {
        typename MatrixType::Index rows = 0, cols = 0;
        in.read(reinterpret_cast<char*> (&rows), sizeof(typename MatrixType::Index));
        in.read(reinterpret_cast<char*> (&cols), sizeof(typename MatrixType::Index));
        Matrix.resize(rows, cols);
        in.read( reinterpret_cast<char*>(Matrix.data()),
                 rows * cols * sizeof(typename MatrixType::Scalar) );
        in.close();
    }
}


template <typename TensorType>
void ITHACAstream::SaveDenseTensor(TensorType& Tensor, word folder,
                                   word MatrixName)
{
    std::ofstream out(folder + MatrixName,
                      std::ios::out | std::ios::binary | std::ios::trunc);
    typename TensorType::Dimensions dim = Tensor.dimensions();
    int tot = 1;

    for (int k = 0; k < dim.size(); k++)
    {
        tot *= dim[k];
    }

    out.write(reinterpret_cast<char*> (&dim),
              sizeof(typename TensorType::Dimensions));
    out.write(reinterpret_cast<char*> (Tensor.data()),
              tot * sizeof(typename TensorType::Scalar) );
    out.close();
}

template <typename TensorType>
void ITHACAstream::ReadDenseTensor(TensorType& Tensor, word folder,
                                   word MatrixName)
{
    std::ifstream in;
    in.open(folder + MatrixName, std::ios::in | std::ios::binary);
    typename TensorType::Dimensions dim;
    in.read(reinterpret_cast<char*> (&dim),
            sizeof(typename TensorType::Dimensions));
    auto dims = Tensor.dimensions();
    M_Assert(dims.size() == dim.size(),
             "The rank of the tensor you want to fill does not coincide with the rank of the tensor you are reading");
    int tot = 1;

    for (int k = 0; k < dim.size(); k++)
    {
        tot *= dim[k];
    }

    Tensor.resize(dim);
    in.read( reinterpret_cast<char*>(Tensor.data()),
             tot * sizeof(typename TensorType::Scalar) );
    in.close();
}

template <typename MatrixType>
void ITHACAstream::SaveSparseMatrixList(List<MatrixType>& MatrixList,
                                        word folder, word MatrixName)
{
    for (int i = 0; i < MatrixList.size(); i++)
    {
        SaveSparseMatrix(MatrixList[i], folder, MatrixName + name(i));
    }
}

template <typename MatrixType>
void ITHACAstream::SaveDenseMatrixList(List<MatrixType>& MatrixList,
                                       word folder, word MatrixName)
{
    for (int i = 0; i < MatrixList.size(); i++)
    {
        SaveDenseMatrix(MatrixList[i], folder, MatrixName + name(i));
    }
}

template <typename T, int whatever, typename IND>
void ITHACAstream::ReadSparseMatrixList(
    List<Eigen::SparseMatrix<T, whatever, IND>>& m, word folder, word MatrixName)
{
    int number_of_files = numberOfFiles(folder, MatrixName);
    std::cout << "Reading the Matrix " + folder + MatrixName << std::endl;
    M_Assert(number_of_files != 0,
             "Check if the file you are trying to read exists" );
    Eigen::SparseMatrix<T, whatever, IND> A;

    for (int i = 0; i < number_of_files; i++)
    {
        ReadSparseMatrix(A, folder, MatrixName + name(i));
        m.append(A);
    }
}


template <typename MatrixType>
void ITHACAstream::ReadDenseMatrixList(List<MatrixType>& m, word folder,
                                       word MatrixName)
{
    int number_of_files = numberOfFiles(folder, MatrixName);
    std::cout << "Reading the Matrix " + folder + MatrixName << std::endl;
    M_Assert(number_of_files != 0,
             "Check if the file you are trying to read exists" );
    MatrixType A;

    for (int i = 0; i < number_of_files; i++)
    {
        ReadDenseMatrix(A, folder, MatrixName + name(i));
        m.append(A);
    }
}





// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



#endif






