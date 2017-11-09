/*---------------------------------------------------------------------------*\
Copyright (C) 2017 by the ITHACA-FV authors

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
/// Directory containing the header, source and template files for the ITHACAstream class.

#ifndef ITHACAstream_H
#define ITHACAstream_H

#include "fvCFD.H"
#include "IOmanip.H"
#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>
#include <algorithm>  
#include "../thirdparty/Eigen/Eigen/Eigen"

#define MAXBUFSIZE  ((int) 1e6)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


/*---------------------------------------------------------------------------*\
                        Class reductionProblem Declaration
\*---------------------------------------------------------------------------*/

/// Class for input-output manipulation
class ITHACAstream
{
private:

public:
    /// Export a generic field to file in a certain folder and subfolder
    /* This is a function to export a generic field into a certain subfolder */
    ///
    /// @param[in] s volVectorField or volScalarField.
    /// @param[in] subfolder string to indicated the subfolder where the field is stored.
    /// @param[in] folder string to indicated the folder where the field is stored.
    /// @param[in] fieldName string to indicate the name of the field/file.
    /// 
    template<typename T>
    static void exportSolution(T& s, fileName subfolder, fileName folder, word fieldName);

    /// Export the reduced matrices in numpy (tipo=python), matlab (tipo=matlab) format and txt (tipo=eigen) format
    /* In this case the function is implemented for a second order matrix  */
    ///
    /// @param[in] matrice Eigen::MatrixXd that you want to export.
    /// @param[in] name string to identify the name you want to use to save the file.
    /// @param[in] tipo string to identify format to export the matrix if numpy (tipo="python"), if matlab (tipo="matlab") if txt (tipo="eigen").
    /// @param[in] folder string to identify the folder where you want to save the file.
    /// 
    static void exportMatrix(Eigen::MatrixXd& matrice, word name, word tipo = "python", word folder = "./Matrices");

    /// Export the reduced matrices in numpy (tipo=python), matlab (tipo=matlab) format and txt (tipo=eigen) format
    /* In this case the function is implemented for a third order matrix
    the eigen matrix is stored as a set of two dimensional matrices where
    each matrix is in a different file */
    ///
    /// @param[in] matrice List < Eigen::MatrixXd > that you want to export.
    /// @param[in] name string to identify the name you want to use to save the file.
    /// @param[in] tipo format to export the matrix if numpy (tipo="python"), if matlab (tipo="matlab") if txt (tipo="eigen").
    /// @param[in] folder string to identify the folder where you want to save the file.
    /// 
    static void exportMatrix(List < Eigen::MatrixXd >& matrice, word name, word tipo = "python", word folder = "./Matrices");

    /// Funtion to read a list of volVectorField from name of the field and casename
    /// 
    /// @param[in]  Lfield      a PtrList of volVectorField where you want to store the field.
    /// @param[in]  Name        The name of the field you want to read.
    /// @param[in]  caseName    The folder where the field is stored.
    /// @param[in]  first_snap  The first snapshots from which you want to start reading the field.
    /// @param[in]  n_snap      The number of snapshots you want to read.
    /// 
    static void read_fields(PtrList<volVectorField>& Lfield, word Name, fileName caseName, label first_snap=0, label n_snap=0);
    
    /// Funtion to read a list of volScalarField from name of the field and folder name
    ///
    /// @param[in]  Lfield      a PtrList of volScalarField where you want to store the field.
    /// @param[in]  Name        The name of the field you want to read.
    /// @param[in]  caseName    The folder where the field is stored.
    /// @param[in]  first_snap  The first snapshots from which you want to start reading the field.
    /// @param[in]  n_snap      The number of snapshots you want to read.
    ///
    static void read_fields(PtrList<volScalarField>& Lfield, word Name, fileName caseName, label first_snap=0, label n_snap=0); 

    /// Funtion to read a list of volVectorField from name of the field a volVectorField if it is already existing
    ///
    /// @param[in]  Lfield      a PtrList of volVectorField where you want to store the field.
    /// @param[in]  field       The field used as template to read other fields.
    /// @param[in]  caseName    The folder where the field is stored.
    /// @param[in]  first_snap  The first snapshots from which you want to start reading the field.
    /// @param[in]  n_snap      The number of snapshots you want to read.
    ///
    static void read_fields(PtrList<volVectorField>& Lfield, volVectorField& field, fileName caseName, label first_snap=0, label n_snap=0);

    /// Funtion to read a list of volScalarField from casename and a list of volScalarField with the same name
    ///
    /// @param[in]  Lfield      a PtrList of volScalarField where you want to store the field.
    /// @param[in]  field       The field used as template to read other fields.
    /// @param[in]  caseName    The folder where the field is stored.
    /// @param[in]  first_snap  The first snapshots from which you want to start reading the field.
    /// @param[in]  n_snap      The number of snapshots you want to read.
    static void read_fields(PtrList<volScalarField>& Lfield, volScalarField& field, fileName caseName, label first_snap=0, label n_snap=0);

	/// Function to export volVectorFields
    ///
    /// @param[in]  field      The field you want to export.
    /// @param[in]  folder     The folder where you want to save the field.
    /// @param[in]  fieldname  The name you want to give to the file. 
    ///
	static void exportFields(PtrList<volVectorField>& field, word folder, word fieldname);
	
	/// Function to export volScalarFields
    ///
    /// @param[in]  field      The field you want to export.
    /// @param[in]  folder     The folder where you want to save the field.
    /// @param[in]  fieldname  The name you want to give to the file. 
    ///
	static void exportFields(PtrList<volScalarField>& field, word folder, word fieldname);
	
    /// Read a two dimensional matrix from e a txt fle in Eigen format
    /*One has to provide the complete filename with the absolute or relative path*/
    ///
    /// @param[in]  filename  The filename of the matrix.
    ///
    /// @return     The Eigen::MatrixXd that contains the imported matrix.
    ///
    static Eigen::MatrixXd readMatrix(word filename);
    
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
};

#ifdef NoRepository
#   include "ITHACAstreamTemplates.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



#endif






