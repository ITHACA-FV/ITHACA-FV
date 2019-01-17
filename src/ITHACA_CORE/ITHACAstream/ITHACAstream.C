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


void ITHACAstream::exportFields(PtrList<volVectorField>& field, word folder,
                                word fieldname)
{
    for (label j = 0; j < field.size() ; j++)
    {
        exportSolution(field[j], name(j + 1), folder, fieldname);
    }

    ITHACAutilities::createSymLink(folder);
}

void ITHACAstream::exportFields(PtrList<volScalarField>& field, word folder,
                                word fieldname)
{
    for (label j = 0; j < field.size() ; j++)
    {
        exportSolution(field[j], name(j + 1), folder, fieldname);
    }

    ITHACAutilities::createSymLink(folder);
}

void ITHACAstream::exportMatrix(Eigen::MatrixXd& matrice, word Name, word tipo,
                                word folder)
{
    mkDir(folder);
    word est;

    if (tipo == "python")
    {
        est = ".py";
        OFstream str(folder + "/" + Name + "_mat" + est);
        str << Name << "=np.array([";

        for (label i = 0; i < matrice.rows(); i++)
        {
            for (label j = 0; j < matrice.cols(); j++)
            {
                if (j == 0)
                {
                    str << "[" << setprecision(10) << matrice(i, j);
                }
                else
                {
                    str << "," << setprecision(10) << matrice(i, j);
                }
            }

            if (i != (matrice.rows() - 1))
            {
                str << "]," << endl;
            }
        }

        str << "]])" << endl;
    }

    if (tipo == "matlab")
    {
        est = ".m";
        OFstream str(folder + "/" + Name + "_mat" + est);
        str << Name << "=[";

        for (label i = 0; i < matrice.rows(); i++)
        {
            for (label j = 0; j < matrice.cols(); j++)
            {
                str << " " << setprecision(10) << matrice(i, j);
            }

            if (i != (matrice.rows() - 1))
            {
                str << ";" << endl;
            }
        }

        str << "];" << endl;
    }

    if (tipo == "eigen")
    {
        std::ofstream ofs;
        ofs.open (folder + "/" + Name + "_mat.txt");
        ofs.precision(20);
        ofs << matrice << std::endl;
        ofs.close();
    }
}

void ITHACAstream::exportMatrix(List < Eigen::MatrixXd >& matrice, word Name,
                                word tipo, word folder)
{
    mkDir(folder);
    word est;

    // Python Case
    if (tipo == "python")
    {
        est = ".py";
        OFstream str(folder + "/" + Name + "_mat" + est);
        str << Name <<  "=np.zeros([" << matrice.size() << "," << matrice[0].rows() <<
            "," << matrice[0].cols() << "])\n";

        for (label i = 0; i < matrice.size(); i++)
        {
            str << Name << "[" << i << ",:,:]=np.array([";

            for (label j = 0; j < matrice[0].rows(); j++)
            {
                for (label k = 0; k < matrice[0].cols(); k++)
                {
                    if ( k == 0)
                    {
                        str << "[" << setprecision(10) << matrice[i](j, k);
                    }
                    else
                    {
                        str << "," << setprecision(10) << matrice[i](j, k);
                    }
                }

                if (j != (matrice[0].rows() - 1))
                {
                    str << "]," << endl;
                }
            }

            str << "]])\n" << endl;
        }
    }
    // Matlab case
    else if (tipo == "matlab")
    {
        est = ".m";
        OFstream str(folder + "/" + Name + "_mat" + est);

        for (label i = 0; i < matrice.size(); i++)
        {
            str << Name << "(" << i + 1 << ",:,:)=[";

            for (label j = 0; j < matrice[0].rows(); j++)
            {
                for (label k = 0; k < matrice[0].cols(); k++)
                {
                    str << " " << setprecision(10) << matrice[i](j, k);
                }

                if (j != (matrice[0].rows() - 1))
                {
                    str << ";" << endl;
                }
            }

            str << "];" << endl;
        }
    }
    else if (tipo == "eigen")
    {
        for (label i = 0; i < matrice.size(); i++)
        {
            word Namei = Name + name(i);
            ITHACAstream::exportMatrix(matrice[i], Namei, "eigen", folder);
        }
    }
}

List< Eigen::MatrixXd > ITHACAstream::readMatrix(word folder, word mat_name)
{
    int file_count = 0;
    DIR* dirp;
    struct dirent* entry;
    dirp = opendir(folder.c_str());
    List <Eigen::MatrixXd > result;

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
    ifstream infile;
    infile.open(filename.c_str());

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

void ITHACAstream::read_fields(PtrList<volVectorField>& Lfield, word Name,
                               fileName casename, label first_snap, label n_snap)
{
    Info << "######### Reading the Data for " << Name << " #########" << endl;
    fileName rootpath(".");
    label last_s;
    Foam::Time runTime2(Foam::Time::controlDictName, rootpath, casename);
    fvMesh mesh
    (
        Foam::IOobject
        (
            Foam::fvMesh::defaultRegion,
            casename + runTime2.timeName(),
            runTime2,
            Foam::IOobject::MUST_READ
        )
    );

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

    for (label i = 2 + first_snap; i < last_s; i++)
    {
        Info << "Reading " << Name << " number " << i - 1 << endl;
        volVectorField tmp_field(
            IOobject
            (
                Name,
                runTime2.times()[i].name(),
                mesh,
                IOobject::MUST_READ
            ),
            mesh
        );
        Lfield.append(tmp_field);
    }

    std::cout << std::endl;
}

void ITHACAstream::read_fields(PtrList<volScalarField>& Lfield, word Name,
                               fileName casename, label first_snap, label n_snap)
{
    Info << " ######### Reading the Data for " << Name << " #########" << endl;
    fileName rootpath(".");
    Foam::Time runTime2(Foam::Time::controlDictName, rootpath, casename);
    label last_s;
    fvMesh mesh
    (
        Foam::IOobject
        (
            Foam::fvMesh::defaultRegion,
            casename + runTime2.timeName(),
            runTime2,
            Foam::IOobject::MUST_READ
        )
    );

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

    for (label i = 2 + first_snap; i < last_s; i++)
    {
        Info << "Reading " << Name << " number " << i - 1 << endl;
        volScalarField tmp_field(
            IOobject
            (
                Name,
                runTime2.times()[i].name(),
                mesh,
                IOobject::MUST_READ
            ),
            mesh
        );
        Lfield.append(tmp_field);
    }

    std::cout << std::endl;
}

void ITHACAstream::read_fields(PtrList<volScalarField>& Lfield,
                               volScalarField& field, fileName casename, label first_snap, label n_snap)
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
        last_s = min(runTime2.times().size(), n_snap + 1);
    }

    for (label i = 2 + first_snap; i < last_s; i++)
    {
        Info << "Reading " << field.name() << " number " << i - 1 << endl;
        volScalarField tmp_field(
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
    }

    std::cout << std::endl;
}

void ITHACAstream::read_fields(PtrList<volVectorField>& Lfield,
                               volVectorField& field, fileName casename, label first_snap, label n_snap)
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
        last_s = min(runTime2.times().size(), n_snap + 1);
    }

    for (label i = 2 + first_snap; i < last_s; i++)
    {
        Info << "Reading " << field.name() << " number " << i - 1 << endl;
        volVectorField tmp_field(
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
    }

    std::cout << std::endl;
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



// * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //

