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
#include "ITHACAassert.H"
#include "ITHACAutilities.H"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#include "Eigen/Eigen"
#include "unsupported/Eigen/CXX11/Tensor"
#pragma GCC diagnostic pop

#define MAXBUFSIZE (static_cast<int> (1e6))


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


/*---------------------------------------------------------------------------*\
                        Class reductionProblem Declaration
\*---------------------------------------------------------------------------*/

/// Class for input-output manipulation
class ITHACAstream
{
    private:

    public:

        //--------------------------------------------------------------------------
        /// @brief      Export an Fv Matrix to folder together with its source term
        ///
        /// @param      Matrix      The matrix
        /// @param[in]  folder      The folder
        /// @param[in]  MatrixName  The fieldname
        ///
        /// @tparam     T           can be scalar or vector
        ///
        template<typename T>
        static void exportFvMatrix(fvMatrix<T>& Matrix, word folder, word MatrixName);

        //--------------------------------------------------------------------------
        /// Export the reduced matrices in numpy (tipo=python), matlab (tipo=matlab) format and txt (tipo=eigen) format
        /* In this case the function is implemented for a second order matrix  */
        ///
        /// @param[in] matrice Eigen::MatrixXd that you want to export.
        /// @param[in] Name string to identify the name you want to use to save the file.
        /// @param[in] tipo string to identify format to export the matrix if numpy (tipo="python"), if matlab (tipo="matlab") if txt (tipo="eigen").
        /// @param[in] folder string to identify the folder where you want to save the file.
        ///
        static void exportMatrix(Eigen::MatrixXd& matrice, word Name,
                                 word tipo = "python", word folder = "./Matrices");

        //--------------------------------------------------------------------------
        /// Export the reduced matrices in numpy (tipo=python), matlab (tipo=matlab) format and txt (tipo=eigen) format
        /* In this case the function is implemented for a third order matrix
        the eigen matrix is stored as a set of two dimensional matrices where
        each matrix is in a different file */
        ///
        /// @param[in] matrice List < Eigen::MatrixXd > that you want to export.
        /// @param[in] Name string to identify the name you want to use to save the file.
        /// @param[in] tipo format to export the matrix if numpy (tipo="python"), if matlab (tipo="matlab") if txt (tipo="eigen").
        /// @param[in] folder string to identify the folder where you want to save the file.
        ///
        static void exportMatrix(List < Eigen::MatrixXd >& matrice, word Name,
                                 word tipo = "python", word folder = "./Matrices");

        //--------------------------------------------------------------------------
        /// Funtion to read a list of volVectorField from name of the field and casename
        ///
        /// @param[in]  Lfield      a PtrList of volVectorField where you want to store the field.
        /// @param[in]  Name        The name of the field you want to read.
        /// @param[in]  casename    The folder where the field is stored.
        /// @param[in]  first_snap  The first snapshots from which you want to start reading the field.
        /// @param[in]  n_snap      The number of snapshots you want to read.
        ///
        static void read_fields(PtrList<volVectorField>& Lfield, word Name,
                                fileName casename, label first_snap = 0, label n_snap = 0);

        //--------------------------------------------------------------------------
        /// Funtion to read a list of volScalarField from name of the field and folder name
        ///
        /// @param[in]  Lfield      a PtrList of volScalarField where you want to store the field.
        /// @param[in]  Name        The name of the field you want to read.
        /// @param[in]  casename    The folder where the field is stored.
        /// @param[in]  first_snap  The first snapshots from which you want to start reading the field.
        /// @param[in]  n_snap      The number of snapshots you want to read.
        ///
        static void read_fields(PtrList<volScalarField>& Lfield, word Name,
                                fileName casename, label first_snap = 0, label n_snap = 0);

        //--------------------------------------------------------------------------
        /// Funtion to read a list of volVectorField from name of the field a volVectorField if it is already existing
        ///
        /// @param[in]  Lfield      a PtrList of volVectorField where you want to store the field.
        /// @param[in]  field       The field used as template to read other fields.
        /// @param[in]  casename    The folder where the field is stored.
        /// @param[in]  first_snap  The first snapshots from which you want to start reading the field.
        /// @param[in]  n_snap      The number of snapshots you want to read.
        ///
        static void read_fields(PtrList<volVectorField>& Lfield, volVectorField& field,
                                fileName casename, label first_snap = 0, label n_snap = 0);

        //--------------------------------------------------------------------------
        /// Funtion to read a list of volScalarField from casename and a list of volScalarField with the same name
        ///
        /// @param[in]  Lfield      a PtrList of volScalarField where you want to store the field.
        /// @param[in]  field       The field used as template to read other fields.
        /// @param[in]  casename    The folder where the field is stored.
        /// @param[in]  first_snap  The first snapshots from which you want to start reading the field.
        /// @param[in]  n_snap      The number of snapshots you want to read.
        static void read_fields(PtrList<volScalarField>& Lfield, volScalarField& field,
                                fileName casename, label first_snap = 0, label n_snap = 0);

        //--------------------------------------------------------------------------
        /// Function to export volVectorFields
        ///
        /// @param[in]  field      The field you want to export.
        /// @param[in]  folder     The folder where you want to save the field.
        /// @param[in]  fieldname  The name you want to give to the file.
        ///
        static void exportFields(PtrList<volVectorField>& field, word folder,
                                 word fieldname);

        //--------------------------------------------------------------------------
        /// Function to export volScalarFields
        ///
        /// @param[in]  field      The field you want to export.
        /// @param[in]  folder     The folder where you want to save the field.
        /// @param[in]  fieldname  The name you want to give to the file.
        ///
        static void exportFields(PtrList<volScalarField>& field, word folder,
                                 word fieldname);

        //--------------------------------------------------------------------------
        /// Read a two dimensional matrix from e a txt fle in Eigen format
        /*One has to provide the complete filename with the absolute or relative path*/
        ///
        /// @param[in]  filename  The filename of the matrix.
        ///
        /// @return     The Eigen::MatrixXd that contains the imported matrix.
        ///
        static Eigen::MatrixXd readMatrix(word filename);

        //--------------------------------------------------------------------------
        /// Read a three dimensional matrix from a txt file in Eigen format
        /*One has to provide the folder containing the matrix files and the filename of the
        the matrix. Since it is stored as a List of matrices each matrix must be stored in a
        different file with the following format:
        matFileName(i)_mat.txt*/
        ///
        /// @param[in]  folder    The folder where the txt files are located
        /// @param[in]  mat_name  The matrix name
        ///
        /// @return     List <Eigen::MatrixXd> that contains the imported matrix.
        ///
        static List <Eigen::MatrixXd> readMatrix(word folder, word mat_name);

        //--------------------------------------------------------------------------
        /// Export a generic field to file in a certain folder and subfolder
        /* This is a function to export a generic field into a certain subfolder */
        ///
        /// @param[in] s volVectorField or volScalarField.
        /// @param[in] subfolder string to indicated the subfolder where the field is stored.
        /// @param[in] folder string to indicated the folder where the field is stored.
        /// @param[in] fieldName string to indicate the name of the field/file.
        ///
        template<typename T>
        static void exportSolution(T& s, fileName subfolder, fileName folder,
                                   word fieldName);

        //--------------------------------------------------------------------------
        /// Export a generic field to file in a certain folder and subfolder
        /* This is a function to export a generic field into a certain subfolder */
        ///
        /// @param[in] s volVectorField or volScalarField.
        /// @param[in] subfolder string to indicated the subfolder where the field is stored.
        /// @param[in] folder string to indicated the folder where the field is stored.
        ///
        template<typename T>
        static void exportSolution(T& s, fileName subfolder, fileName folder);

        //--------------------------------------------------------------------------
        /// @brief      Export a list type object to file
        ///
        /// @param[in]  list      The list
        /// @param[in]  filename  The filename
        /// @param[in]  folder    The folder
        ///
        /// @tparam     T         generic list type object
        ///
        template<typename T>
        static void exportList(T& list, word folder, word filename);


        //--------------------------------------------------------------------------
        /// @brief      Saves an Eigen sparse matrix into bynary format file
        ///
        /// @param      m           Eigen::SparseMatrix m
        /// @param[in]  folder      The folder where you want to save the matrix
        /// @param[in]  MatrixName  The matrix name for the output file
        ///
        /// @tparam     T           type of the matrix, i.e. double, float, ...
        /// @tparam     whatever    number of matrix rows
        /// @tparam     IND         number of matrix cols
        ///
        template <typename T, int whatever, typename IND>
        static void SaveSparseMatrix(Eigen::SparseMatrix<T, whatever, IND>& m,
                                     word folder, word MatrixName);

        //--------------------------------------------------------------------------
        /// @brief      Reads an Eigen sparse matrix from a bynary format file
        ///
        /// @param[in,out]  m           Eigen::SparseMatrix m
        /// @param[in]      folder      The folder from where you want to read the matrix
        /// @param[in]      MatrixName  The matrix name for the input file
        ///
        /// @tparam     T           type of the matrix, i.e. double, float, ...
        /// @tparam     whatever    number of matrix rows
        /// @tparam     IND         number of matrix cols
        ///
        template <typename T, int whatever, typename IND>
        static void ReadSparseMatrix(Eigen::SparseMatrix<T, whatever, IND>& m,
                                     word folder, word MatrixName);

        //--------------------------------------------------------------------------
        /// @brief      Saves a dense matrix to a binary format file
        ///
        /// @param[in]  Matrix      The  Eigen dense matrix
        /// @param[in]  folder      he folder where you want to save the matrix
        /// @param[in]  MatrixName  The matrix name for the output file
        ///
        /// @tparam     MatrixType  type of the matrix, i.e. double, float, ...
        ///
        template <typename MatrixType>
        static void SaveDenseMatrix(MatrixType& Matrix, word folder, word MatrixName);

        //--------------------------------------------------------------------------
        /// @brief      Reads a dense matrix from a binary format file
        ///
        /// @param[in,out]  Matrix      The Eigen dense matrix
        /// @param[in]      folder      The folder from where you want to read the matrix
        /// @param[in]      MatrixName  The matrix name of the input file
        ///
        /// @tparam     MatrixType  type of the matrix, i.e. double, float, ...
        ///
        template <typename MatrixType>
        static void ReadDenseMatrix(MatrixType& Matrix, word folder, word MatrixName);


        //----------------------------------------------------------------------
        /// @brief      Saves a dense tensor.
        ///
        /// @param      Tensor      The tensor
        /// @param[in]  folder      The folder
        /// @param[in]  MatrixName  The matrix name
        ///
        /// @tparam     TensorType  type of the tensor, i.e. double, float, ...
        ///
        template <typename TensorType>
        static void SaveDenseTensor(TensorType& Tensor, word folder, word MatrixName);


        //----------------------------------------------------------------------
        /// @brief      Reads a dense tensor.
        ///
        /// @param      Tensor      The tensor
        /// @param[in]  folder      The folder
        /// @param[in]  MatrixName  The matrix name
        ///
        /// @tparam     TensorType   type of the tensor, i.e. double, float, ...
        ///
        template <typename TensorType>
        static void ReadDenseTensor(TensorType& Tensor, word folder, word MatrixName);

        //--------------------------------------------------------------------------
        /// @brief      Saves a sparse matrix list to bynary files
        ///
        /// @param[in]  MatrixList  The matrix list
        /// @param[in]  folder      The folder where you want to save files
        /// @param[in]  MatrixName  The matrix name prefix, files are the named with an incremental number
        ///
        /// @tparam     MatrixType   MatrixType  type of the matrix, i.e. double, float, ...
        ///
        template <typename MatrixType>
        static void SaveSparseMatrixList(List<MatrixType>& MatrixList, word folder,
                                         word MatrixName);

        //--------------------------------------------------------------------------
        /// @brief      Reads a sparse matrix list to bynary files
        ///
        /// @param[in,out]  MatrixList  The matrix list
        /// @param[in]      folder      The folder from where you want to read the files
        /// @param[in]      MatrixName  The matrix name prefix, files are the named with an incremental number
        ///
        /// @tparam     MatrixType   MatrixType  type of the matrix, i.e. double, float, ...
        ///
        template <typename T, int whatever, typename IND>
        static void ReadSparseMatrixList(List<Eigen::SparseMatrix<T, whatever, IND>>& m,
                                         word folder, word MatrixName);

        //--------------------------------------------------------------------------
        /// @brief      Saves a dense matrix list to bynary files
        ///
        /// @param[in]  MatrixList  The matrix list
        /// @param[in]  folder      The folder where you want to save files
        /// @param[in]  MatrixName  The matrix name prefix, files are the named with an incremental number
        ///
        /// @tparam     MatrixType   MatrixType  type of the matrix, i.e. double, float, ...
        ///
        template <typename MatrixType>
        static void SaveDenseMatrixList(List<MatrixType>& MatrixList, word folder,
                                        word MatrixName);


        //--------------------------------------------------------------------------
        /// @brief      Reads a sparse matrix list to bynary files
        ///
        /// @param[in,out]  MatrixList  The matrix list
        /// @param[in]      folder      The folder from where you want to read the files
        /// @param[in]      MatrixName  The matrix name prefix, files are the named with an incremental number
        ///
        /// @tparam     MatrixType   MatrixType  type of the matrix, i.e. double, float, ...
        ///
        template <typename MatrixType>
        static void ReadDenseMatrixList(List<MatrixType>& m, word folder,
                                        word MatrixName);

        //--------------------------------------------------------------------------
        /// @brief      Count the number of files with a certain name prefix, used by reading functions of list of matrices
        ///
        /// @param[in]  folder      The folder
        /// @param[in]  MatrixName  The matrix name
        ///
        /// @return     the number of files matching the request
        ///
        static int numberOfFiles(word folder, word MatrixName);
};

template<typename T>
void ITHACAstream::exportFvMatrix(fvMatrix<T>& Matrix, word folder,
                                  word MatrixName)
{
    Eigen::SparseMatrix<double> A;
    Eigen::VectorXd b;
    Foam2Eigen::fvMatrix2Eigen(Matrix, A, b);
    SaveSparseMatrix(A, folder + "/", "A_" + MatrixName);
    SaveDenseMatrix(b, folder + "/", "B_" + MatrixName);
}

template<typename T>
void ITHACAstream::exportSolution(T& s, fileName subfolder, fileName folder,
                                  word fieldName)
{
    if (!Pstream::parRun())
    {
        mkDir(folder + "/" + subfolder);
        ITHACAutilities::createSymLink(folder);
        T act(fieldName, s);
        fileName fieldname = folder + "/" + subfolder + "/" + fieldName;
        OFstream os(fieldname);
        act.writeHeader(os);
        os << act << endl;
    }
    else
    {
        mkDir(folder + "/processor" + name(Pstream::myProcNo()) + "/" + subfolder);
        ITHACAutilities::createSymLink(folder);
        T act(fieldName, s);
        fileName fieldname = folder + "/processor" + name(Pstream::myProcNo()) + "/" +
                             subfolder + "/" + fieldName;
        std::cout << fieldname << std::endl;
        OFstream os(fieldname);
        act.writeHeader(os);
        os << act << endl;
    }
}

template<typename T>
void ITHACAstream::exportSolution(T& s, fileName subfolder, fileName folder)
{
    if (!Pstream::parRun())
    {
        mkDir(folder + "/" + subfolder);
        ITHACAutilities::createSymLink(folder);
        fileName fieldname = folder + "/" + subfolder + "/" + s.name();
        OFstream os(fieldname);
        s.writeHeader(os);
        os << s << endl;
    }
    else
    {
        mkDir(folder + "/processor" + name(Pstream::myProcNo()) + "/" + subfolder);
        ITHACAutilities::createSymLink(folder);
        fileName fieldname = folder + "/processor" + name(Pstream::myProcNo()) + "/" +
                             subfolder + "/" + s.name();
        OFstream os(fieldname);
        s.writeHeader(os);
        os << s << endl;
    }
}


template<typename T>
void ITHACAstream::exportList(T& list, word folder, word filename)
{
    mkDir(folder);
    word fieldname = folder + filename;
    OFstream os(fieldname);

    for (label i = 0; i < list.size(); i++)
    {
        os << list[i] << endl;
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

    for (unsigned int k = 0; k < dim.size(); k++)
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

    for (unsigned int k = 0; k < dim.size(); k++)
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






