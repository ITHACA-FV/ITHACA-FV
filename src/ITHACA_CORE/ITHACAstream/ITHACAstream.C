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
/// Source file of the ITHACAstream class, it contains the implementation of
/// several methods for input output operations.

#include "ITHACAstream.H"


// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<typename Type>
void ITHACAstream::exportFvMatrix(fvMatrix<Type>& Matrix, word folder,
                                  word MatrixName)
{
    Eigen::SparseMatrix<double> A;
    Eigen::VectorXd b;
    Foam2Eigen::fvMatrix2Eigen(Matrix, A, b);
    SaveSparseMatrix(A, folder + "/", "A_" + MatrixName);
    SaveDenseMatrix(b, folder + "/", "B_" + MatrixName);
}

template <typename T>
void ITHACAstream::exportMatrix(Eigen::Matrix < T, -1, -1 > & matrix,
                                word Name, word type,
                                word folder)
{
    std::string message = "The extension \"" +  type +
                          "\" was not implemented. Check the list of possible extensions.";
    M_Assert(type == "python" || type == "matlab"
             || type == "eigen", message.c_str()
            );
    mkDir(folder);
    word est;

    if (type == "python")
    {
        est = ".py";
        OFstream str(folder + "/" + Name + "_mat" + est);
        str << Name << "=np.array([";

        for (label i = 0; i < matrix.rows(); i++)
        {
            for (label j = 0; j < matrix.cols(); j++)
            {
                if (j == 0)
                {
                    str << "[" << setprecision(10) << matrix(i, j);
                }
                else
                {
                    str << "," << setprecision(10) << matrix(i, j);
                }
            }

            if (i != (matrix.rows() - 1))
            {
                str << "]," << endl;
            }
        }

        str << "]])" << endl;
    }

    if (type == "matlab")
    {
        est = ".m";
        OFstream str(folder + "/" + Name + "_mat" + est);
        str << Name << "=[";

        for (label i = 0; i < matrix.rows(); i++)
        {
            for (label j = 0; j < matrix.cols(); j++)
            {
                str << " " << setprecision(10) << matrix(i, j);
            }

            if (i != (matrix.rows() - 1))
            {
                str << ";" << endl;
            }
        }

        str << "];" << endl;
    }

    if (type == "eigen")
    {
        std::ofstream ofs;
        ofs.open (folder + "/" + Name + "_mat.txt");
        ofs.precision(20);
        ofs << matrix << std::endl;
        ofs.close();
    }
}

template void ITHACAstream::exportMatrix(Eigen::Matrix < double, -1,
        -1 > & matrix, word Name, word type,
        word folder);

template void ITHACAstream::exportMatrix(Eigen::Matrix < int, -1,
        -1 > & matrix, word Name, word type,
        word folder);

template void ITHACAstream::exportMatrix(Eigen::Matrix < float, -1,
        -1 > & matrix, word Name, word type,
        word folder);

void ITHACAstream::exportMatrix(List <Eigen::MatrixXd>& matrix, word Name,
                                word type, word folder)
{
    std::string message = "The extension \"" +  type +
                          "\" was not implemented. Check the list of possible extensions.";
    M_Assert(type == "python" || type == "matlab"
             || type == "eigen", message.c_str()
            );
    mkDir(folder);
    word est;

    // Python Case
    if (type == "python")
    {
        est = ".py";
        OFstream str(folder + "/" + Name + "_mat" + est);
        str << Name <<  "=np.zeros([" << matrix.size() << "," << matrix[0].rows() <<
            "," << matrix[0].cols() << "])\n";

        for (label i = 0; i < matrix.size(); i++)
        {
            str << Name << "[" << i << ",:,:]=np.array([";

            for (label j = 0; j < matrix[0].rows(); j++)
            {
                for (label k = 0; k < matrix[0].cols(); k++)
                {
                    if ( k == 0)
                    {
                        str << "[" << setprecision(10) << matrix[i](j, k);
                    }
                    else
                    {
                        str << "," << setprecision(10) << matrix[i](j, k);
                    }
                }

                if (j != (matrix[0].rows() - 1))
                {
                    str << "]," << endl;
                }
            }

            str << "]])\n" << endl;
        }
    }
    // Matlab case
    else if (type == "matlab")
    {
        est = ".m";
        OFstream str(folder + "/" + Name + "_mat" + est);

        for (label i = 0; i < matrix.size(); i++)
        {
            str << Name << "(" << i + 1 << ",:,:)=[";

            for (label j = 0; j < matrix[0].rows(); j++)
            {
                for (label k = 0; k < matrix[0].cols(); k++)
                {
                    str << " " << setprecision(10) << matrix[i](j, k);
                }

                if (j != (matrix[0].rows() - 1))
                {
                    str << ";" << endl;
                }
            }

            str << "];" << endl;
        }
    }
    else if (type == "eigen")
    {
        for (label i = 0; i < matrix.size(); i++)
        {
            word Namei = Name + name(i);
            ITHACAstream::exportMatrix(matrix[i], Namei, "eigen", folder);
        }
    }
}

template<typename T>
void ITHACAstream::exportTensor(Eigen::Tensor<T, 3> tensor, word Name,
                                word type, word folder)
{
    std::string message = "The extension \"" +  type +
                          "\" was not implemented. Check the list of possible extensions.";
    M_Assert(type == "python" || type == "matlab"
             || type == "eigen", message.c_str()
            );
    mkDir(folder);
    word est;

    // Python Case
    if (type == "python")
    {
        est = ".py";
        OFstream str(folder + "/" + Name + "_mat" + est);
        str << Name <<  "=np.zeros([" << tensor.dimension(0) << "," <<
            Eigen::SliceFromTensor(
                tensor, 0,
                0).rows() <<
            "," << Eigen::SliceFromTensor(tensor, 0,
                                          0).cols() << "])\n";

        for (label i = 0; i < tensor.dimension(0); i++)
        {
            str << Name << "[" << i << ",:,:]=np.array([";

            for (label j = 0; j < Eigen::SliceFromTensor(tensor, 0, 0).rows(); j++)
            {
                for (label k = 0; k < Eigen::SliceFromTensor(tensor, 0, 0).cols(); k++)
                {
                    if ( k == 0)
                    {
                        str << "[" << setprecision(10) << Eigen::SliceFromTensor(tensor, 0,
                                i)(j, k);
                    }
                    else
                    {
                        str << "," << setprecision(10) << Eigen::SliceFromTensor(tensor, 0,
                                i)(j, k);
                    }
                }

                if (j != (Eigen::SliceFromTensor(tensor, 0,
                                                 0).rows() - 1))
                {
                    str << "]," << endl;
                }
            }

            str << "]])\n" << endl;
        }
    }
    // Matlab case
    else if (type == "matlab")
    {
        est = ".m";
        OFstream str(folder + "/" + Name + "_mat" + est);

        for (label i = 0; i < tensor.dimension(0); i++)
        {
            str << Name << "(" << i + 1 << ",:,:)=[";

            for (label j = 0; j < Eigen::SliceFromTensor(tensor, 0,
                    0).rows(); j++)
            {
                for (label k = 0; k < Eigen::SliceFromTensor(tensor, 0,
                        0).cols(); k++)
                {
                    str << " " << setprecision(10) << Eigen::SliceFromTensor(tensor, 0,
                            i)(j, k);
                }

                if (j != (Eigen::SliceFromTensor(tensor, 0,
                                                 0).rows() - 1))
                {
                    str << ";" << endl;
                }
            }

            str << "];" << endl;
        }
    }
    else if (type == "eigen")
    {
        for (label i = 0; i < tensor.dimension(0); i++)
        {
            Eigen::Matrix < T, -1, -1 > matrixAux = Eigen::SliceFromTensor(tensor, 0, i);
            word Namei = Name + name(i);
            exportMatrix(matrixAux, Namei, "eigen", folder);
        }
    }
}



template void ITHACAstream::exportTensor(Eigen::Tensor<double, 3> tensor,
        word Name,
        word type, word folder);

template void ITHACAstream::exportTensor(Eigen::Tensor<int, 3> tensor,
        word Name,
        word type, word folder);

template void ITHACAstream::exportTensor(Eigen::Tensor<float, 3> tensor,
        word Name,
        word type, word folder);

List<Eigen::MatrixXd> ITHACAstream::readMatrix(word folder, word mat_name)
{
    int file_count = 0;
    DIR* dirp;
    struct dirent* entry;
    dirp = opendir(folder.c_str());
    List <Eigen::MatrixXd> result;

    while ((entry = readdir(dirp)) != NULL)
    {
        if (entry->d_type == DT_REG)
        {
            file_count++;
        }
    }

    closedir (dirp);

    for (label j = 0; j < file_count ; j++)
    {
        word matname = folder + "/" + mat_name + name(j) + "_mat.txt";
        Eigen::MatrixXd temp = readMatrix(matname);
        result.append(temp);
    }

    return result;
}

Eigen::MatrixXd ITHACAstream::readMatrix(word filename)
{
    int cols = 0, rows = 0;
    double buff[MAXBUFSIZE];
    // Read numbers from file into buffer.
    std::ifstream infile;
    infile.open(filename.c_str());
    std::string message = "The matrix file \"" +  filename +
                          "\" does not exist. Check the existence of the file or the way it is named.";
    M_Assert(infile.good() != 0, message.c_str()
            );

    while (! infile.eof())
    {
        string line;
        getline(infile, line);
        int temp_cols = 0;
        std::stringstream stream(line);

        while (! stream.eof())
        {
            stream >> buff[cols * rows + temp_cols++];
        }

        if (temp_cols == 0)
        {
            continue;
        }

        if (cols == 0)
        {
            cols = temp_cols;
        }

        rows++;
    }

    infile.close();
    rows--;
    // Populate matrix with numbers.
    Eigen::MatrixXd result(rows, cols);

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            result(i, j) = buff[ cols * i + j ];
        }
    }

    return result;
}

template<class Type, template<class> class PatchField, class GeoMesh>
void ITHACAstream::read_fields(PtrList<GeometricField<Type, PatchField, GeoMesh>>& Lfield, word Name,
                               fileName casename, label first_snap, label n_snap)
{
    ITHACAparameters* para(ITHACAparameters::getInstance());
    fvMesh& mesh = para->mesh;

    if (!Pstream::parRun())
    {
        Info << "######### Reading the Data for " << Name << " #########" << endl;
        fileName rootpath(".");
        label last_s;
        Foam::Time runTime2(Foam::Time::controlDictName, rootpath, casename);

        if (first_snap >= runTime2.times().size())
        {
            Info << "Error the index of the first snapshot must be smaller than the number of snapshots"
                 << endl;
            exit(0);
        }

        if (n_snap == 0)
        {
            last_s = runTime2.times().size();
        }
        else
        {
            last_s = min(runTime2.times().size(), n_snap + 2);
        }

        for (label i = 2 + first_snap; i < last_s + first_snap; i++)
        {
            GeometricField<Type, PatchField, GeoMesh> tmp_field(
                IOobject
                (
                    Name,
                    casename + runTime2.times()[i].name(),
                    mesh,
                    IOobject::MUST_READ
                ),
                mesh
            );
            Lfield.append(tmp_field);
            printProgress(double(i + 1) / (last_s + first_snap));
        }

        std::cout << std::endl;
    }
    else
    {
        std::cerr << "File: ITHACAstream.C, Line: 343" << std::endl;
    }
}

template<class Type, template<class> class PatchField, class GeoMesh>
void ITHACAstream::read_fields(PtrList<GeometricField<Type, PatchField, GeoMesh>>& Lfield,
                               GeometricField<Type, PatchField, GeoMesh>& field, 
			       fileName casename, label first_snap, label n_snap)
{
    if (!Pstream::parRun())
    {
        Info << "######### Reading the Data for " << field.name() << " #########" <<
             endl;
        fileName rootpath(".");
        Foam::Time runTime2(Foam::Time::controlDictName, rootpath, casename);
        label last_s;

        if (first_snap >= runTime2.times().size())
        {
            Info << "Error the index of the first snapshot must be smaller than the number of snapshots"
                 << endl;
            exit(0);
        }

        if (n_snap == 0)
        {
            last_s = runTime2.times().size();
        }
        else
        {
            last_s = min(runTime2.times().size(), n_snap + 2);
        }

        for (label i = 2 + first_snap; i < last_s + first_snap; i++)
        {
            GeometricField<Type, PatchField, GeoMesh> tmp_field(
                IOobject
                (
                    field.name(),
                    casename + runTime2.times()[i].name(),
                    field.mesh(),
                    IOobject::MUST_READ
                ),
                field.mesh()
            );
            Lfield.append(tmp_field);
            printProgress(double(i + 1) / (last_s + first_snap));
        }

        std::cout << std::endl;
    }
    else
    {
        Info << "######### Reading the Data for " << field.name() << " #########" <<
             endl;
        word timename(field.mesh().time().rootPath() + "/" +
                      field.mesh().time().caseName() );
        timename = timename.substr(0, timename.find_last_of("\\/"));
        timename = timename + "/" + casename + "processor" + name(Pstream::myProcNo());
        label last_s = numberOfFiles(casename,
                                     "processor" + name(Pstream::myProcNo()) + "/");

        if (first_snap > last_s)
        {
            Info << "Error the index of the first snapshot must be smaller than the number of snapshots"
                 << endl;
            exit(0);
        }

        if (n_snap == 0)
        {
        }
        else
        {
            last_s = min(last_s, n_snap + 2);
        }

        for (label i = 2 + first_snap; i < last_s + first_snap; i++)
        {
            GeometricField<Type, PatchField, GeoMesh> tmp_field(
                IOobject
                (
                    field.name(),
                    timename + "/" + name(i),
                    field.mesh(),
                    IOobject::MUST_READ
                ),
                field.mesh()
            );
            Lfield.append(tmp_field);
            printProgress(double(i + 1) / (last_s + first_snap));
        }

        Info << endl;
    }
}

template<class Type, template<class> class PatchField, class GeoMesh>
void ITHACAstream::readMiddleFields(PtrList<GeometricField<Type, PatchField, GeoMesh>>& Lfield,
                                    GeometricField<Type, PatchField, GeoMesh>& field, fileName casename)
{
    int par = 1;
    M_Assert(ITHACAutilities::check_folder(casename + name(par)) != 0,
             "No parameter dependent solutions stored into Offline folder");

    while (ITHACAutilities::check_folder(casename + name(par)))
    {
        ITHACAstream::read_fields(Lfield, field, casename + name(par) + "/");
        par++;
    }
}

template<class Type, template<class> class PatchField, class GeoMesh>
void ITHACAstream::readConvergedFields(PtrList<GeometricField<Type, PatchField, GeoMesh>>& Lfield,
                                       GeometricField<Type, PatchField, GeoMesh>& field,
                                       fileName casename)
{
    int par = 1;
    M_Assert(ITHACAutilities::check_folder(casename + name(par)) != 0,
             "No parameter dependent solutions stored into Offline folder");
    std::cout << "######### Reading the Data for " << field.name() << " #########"
              << std::endl;

    while (ITHACAutilities::check_folder(casename + name(par)))
    {
        int last = 1;

        while (ITHACAutilities::check_folder(casename + name(par) + "/" + name(last)))
        {
            last++;
        }

        GeometricField<Type, PatchField, GeoMesh> tmpField(
            IOobject
            (
                field.name(),
                casename + name(par) + "/" + name(last - 1),
                field.mesh(),
                IOobject::MUST_READ
            ),
            field.mesh()
        );
        Lfield.append(tmpField);
        par++;
    }
}

int ITHACAstream::numberOfFiles(word folder, word MatrixName)
{
    int number_of_files = 0;
    std::ifstream in;
    in.open(folder + MatrixName + name(0), std::ios::in | std::ios::binary);

    while (in.good())
    {
        in.close();
        number_of_files++;
        in.open(folder + MatrixName + name(number_of_files),
                std::ios::in | std::ios::binary);
    }

    in.close();
    return number_of_files;
}

template<class Type, template<class> class PatchField, class GeoMesh>
void ITHACAstream::exportFields(PtrList<GeometricField<Type, PatchField, GeoMesh>>& field, 
				word folder, word fieldname)
{
    ITHACAutilities::createSymLink(folder);
    Info << "######### Exporting the Data for " << fieldname << " #########" <<
         endl;

    for (label j = 0; j < field.size() ; j++)
    {
        exportSolution(field[j], name(j + 1), folder, fieldname);
        printProgress(double(j + 1) / field.size());
    }

    std::cout << std::endl;
}

template<class Type, template<class> class PatchField, class GeoMesh>
void ITHACAstream::exportSolution(GeometricField<Type, PatchField, GeoMesh>& s, 
				  fileName subfolder, fileName folder,
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
        GeometricField<Type, PatchField, GeoMesh> act(fieldName, s);
        fileName fieldname = folder + "/processor" + name(Pstream::myProcNo()) + "/" +
                             subfolder + "/" + fieldName;
        std::cout << fieldname << std::endl;
        OFstream os(fieldname);
        act.writeHeader(os);
        os << act << endl;
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void ITHACAstream::exportSolution(GeometricField<Type, PatchField, GeoMesh>& s, 
				  fileName subfolder, fileName folder)
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


void ITHACAstream::writePoints(pointField points, fileName folder,
                               fileName subfolder)
{
    if (!Pstream::parRun())
    {
        mkDir(folder + "/" + subfolder);
        ITHACAutilities::createSymLink(folder);
        fileName fieldname = folder + "/" + subfolder + "/" + "points";
        OFstream os(fieldname);
        os << "FoamFile \n { \n version     2.0; \n format      ascii; \n class       vectorField; \n location    ""1 / polyMesh""; \n object      points; \n }"
           << endl;
        os << points << endl;
    }
    else
    {
        mkDir(folder + "/processor" + name(Pstream::myProcNo()) + "/" + subfolder);
        ITHACAutilities::createSymLink(folder);
        fileName fieldname = folder + "/processor" + name(Pstream::myProcNo()) + "/" +
                             subfolder + "/" + "points";
        OFstream os(fieldname);
        os << "FoamFile \n { \n version     2.0; \n format      ascii; \n class       vectorField; \n location    ""1 / polyMesh""; \n object      points; \n }"
           << endl;
        os << points << endl;
    }
}

void ITHACAstream::printProgress(double percentage)
{
    int val = static_cast<int>(percentage * 100);
    int lpad = static_cast<int> (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;

    if (Pstream::master())
    {
        printf ("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
        fflush (stdout);
    }
}

template void ITHACAstream::read_fields(PtrList<volScalarField>& Lfield,
                                        word Name,
                                        fileName casename, label first_snap, label n_snap);
template void ITHACAstream::read_fields(PtrList<volVectorField>& Lfield,
                                        word Name,
                                        fileName casename, label first_snap, label n_snap);
template void ITHACAstream::read_fields(PtrList<surfaceScalarField>& Lfield,
                                        word Name,
                                        fileName casename, label first_snap, label n_snap);
template void ITHACAstream::read_fields(PtrList<surfaceVectorField>& Lfield,
                                        word Name,
                                        fileName casename, label first_snap, label n_snap);
template void ITHACAstream::read_fields(PtrList<volScalarField>& Lfield,
                                        volScalarField& field, fileName casename, label first_snap, label n_snap);
template void ITHACAstream::read_fields(PtrList<volVectorField>& Lfield,
                                        volVectorField& field, fileName casename, label first_snap, label n_snap);
template void ITHACAstream::read_fields(PtrList<surfaceScalarField>& Lfield,
                                        surfaceScalarField& field, fileName casename, label first_snap, label n_snap);
template void ITHACAstream::read_fields(PtrList<surfaceVectorField>& Lfield,
                                        surfaceVectorField& field, fileName casename, label first_snap, label n_snap);
template void ITHACAstream::readMiddleFields(PtrList<volScalarField>& Lfield,
        volScalarField& field, fileName casename);
template void ITHACAstream::readMiddleFields(PtrList<volVectorField>& Lfield,
        volVectorField& field, fileName casename);
template void ITHACAstream::readMiddleFields(PtrList<surfaceScalarField>&
        Lfield, surfaceScalarField& field, fileName casename);
template void ITHACAstream::readMiddleFields(PtrList<surfaceVectorField>&
        Lfield, surfaceVectorField& field, fileName casename);
template void ITHACAstream::readConvergedFields(PtrList<volScalarField>& Lfield,
        volScalarField& field, fileName casename);
template void ITHACAstream::readConvergedFields(PtrList<volVectorField>& Lfield,
        volVectorField& field, fileName casename);
template void ITHACAstream::readConvergedFields(PtrList<surfaceScalarField>&
        Lfield, surfaceScalarField& field, fileName casename);
template void ITHACAstream::readConvergedFields(PtrList<surfaceVectorField>&
        Lfield, surfaceVectorField& field, fileName casename);

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

template <typename T, int Nrows, typename IND>
void ITHACAstream::SaveSparseMatrix(Eigen::SparseMatrix<T, Nrows, IND>& m,
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

template <typename T, int Nrows, typename IND>
void ITHACAstream::ReadSparseMatrix(Eigen::SparseMatrix<T, Nrows, IND>& m,
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

template <typename T, int Nrows, typename IND>
void ITHACAstream::ReadSparseMatrixList(
    List<Eigen::SparseMatrix<T, Nrows, IND>>& m, word folder, word MatrixName)
{
    int number_of_files = numberOfFiles(folder, MatrixName);
    std::cout << "Reading the Matrix " + folder + MatrixName << std::endl;
    M_Assert(number_of_files != 0,
             "Check if the file you are trying to read exists" );
    Eigen::SparseMatrix<T, Nrows, IND> A;

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
