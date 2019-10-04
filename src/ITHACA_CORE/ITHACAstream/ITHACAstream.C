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
void ITHACAstream::exportTensor(Eigen::Tensor<T, 3 > tensor, word Name,
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

template void ITHACAstream::exportTensor(Eigen::Tensor<double, 3 > tensor,
        word Name,
        word type, word folder);

template void ITHACAstream::exportTensor(Eigen::Tensor<int, 3 > tensor,
        word Name,
        word type, word folder);

template void ITHACAstream::exportTensor(Eigen::Tensor<float, 3 > tensor,
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

void ITHACAstream::writePoints(pointField points, fileName folder,
                               fileName subfolder)
{
    mkDir(folder + "/" + subfolder);
    ITHACAutilities::createSymLink(folder);
    fileName fieldname = folder + "/" + subfolder + "/" + "points";
    OFstream os(fieldname);
    os << "FoamFile \n { \n version     2.0; \n format      ascii; \n class       vectorField; \n location    ""1 / polyMesh""; \n object      points; \n }"
       << endl;
    os << points << endl;
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
