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
#include "ITHACAsystem.H"

namespace ITHACAutilities
{

void createSymLink(word folder)
{
    if (!Pstream::parRun())
    {
        mkDir(folder);
        word command1("ln -s  $(readlink -f constant/) " + folder + "/" +
                      " >/dev/null 2>&1");
        word command2("ln -s  $(readlink -f system/) " + folder + "/" +
                      " >/dev/null 2>&1");
        word command3("ln -s  $(readlink -f 0/) " + folder + "/" + " >/dev/null 2>&1");
        std::cout.setstate(std::ios_base::failbit);
        system(command1);
        system(command2);
        system(command3);
        std::cout.clear();
    }
    else
    {
        // Serial Part
        mkDir(folder);
        word command1s("ln -s  $(readlink -f constant/) " + folder + "/" +
                       " >/dev/null 2>&1");
        word command2s("ln -s  $(readlink -f system/) " + folder + "/" +
                       " >/dev/null 2>&1");
        word command3s("ln -s  $(readlink -f 0/) " + folder + "/" + " >/dev/null 2>&1");
        std::cout.setstate(std::ios_base::failbit);
        system(command1s);
        system(command2s);
        system(command3s);
        std::cout.clear();
        // Parallel Part
        mkDir(folder + "/processor" + name(Pstream::myProcNo()));
        word command1("ln -s  $(readlink -f processor" + name(Pstream::myProcNo()) +
                      "/constant/) " + folder + "/processor" +
                      name(Pstream::myProcNo()) + "/ >/dev/null 2>&1");
        word command2("ln -s  $(readlink -f system/) " + folder + "/processor" +
                      name(Pstream::myProcNo()) + "/ >/dev/null 2>&1");
        word command3("ln -s  $(readlink -f processor" + name(Pstream::myProcNo()) +
                      "/0/) " + folder + "/processor" +
                      name(Pstream::myProcNo()) + "/ >/dev/null 2>&1");
        std::cout.setstate(std::ios_base::failbit);
        system(command1);
        system(command2);
        system(command3);
        std::cout.clear();
    }
}

void createSymLink(word linkFolder, word destFolder)
{
    if (!check_folder(destFolder))
    {
        mkDir(destFolder);
    }

    system("ln -s  $(readlink -f " + linkFolder + ") " + destFolder +
           " >/dev/null 2>&1");
}

bool check_folder(word folder)
{
    struct stat sb;
    bool exist;

    if (stat(folder.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode))
    {
        exist = true;
    }
    else
    {
        exist = false;
    }

    return exist;
}

// Check if the offline data exist
bool check_off()
{
    bool off_exist = 0;

    if (Pstream::master())
    {
        if (check_folder("./ITHACAoutput/Offline"))
        {
            off_exist = true;
            Info << "Offline data already exist, reading existing data" << endl;
        }
        else
        {
            off_exist = false;
            Info << "Offline don't exist, performing the Offline Solve" << endl;
            mkDir("./ITHACAoutput/Offline");
        }
    }

    reduce(off_exist, sumOp<label>());

    if (!off_exist)
    {
        createSymLink("./ITHACAoutput/Offline");
    }

    return off_exist;
}

bool check_file(std::string fileName)
{
    std::ifstream infile(fileName);
    return infile.good();
}

// Check if the modes exists
bool check_pod()
{
    bool pod_exist = 0;

    if (Pstream::master())
    {
        if (check_folder("./ITHACAoutput/POD"))
        {
            pod_exist = true;
            Info << "POD data already exist, reading existing modes" << endl;
        }
        else if (!Pstream::parRun())
        {
            pod_exist = false;
            Info << "POD don't exist, performing a POD decomposition" << endl;
            mkDir("./ITHACAoutput/POD");
        }
    }

    reduce(pod_exist, sumOp<label>());

    if (!pod_exist)
    {
        createSymLink("./ITHACAoutput/POD");
    }

    return pod_exist;
}

// Check if the supremizer data exist
bool check_sup()
{
    bool sup_exist = 0;

    if (Pstream::master())
    {
        if (check_folder("./ITHACAoutput/supremizer"))
        {
            sup_exist = true;
            Info << "Supremizer data already exist, reading existing data" << endl;
        }
        else
        {
            sup_exist = false;
            Info << "Supremizers don't exist, performing a POD decomposition" << endl;
            mkDir("./ITHACAoutput/supremizer");
        }
    }

    reduce(sup_exist, sumOp<label>());

    if (!sup_exist)
    {
        createSymLink("./ITHACAoutput/supremizer");
    }

    return sup_exist;
}


}
