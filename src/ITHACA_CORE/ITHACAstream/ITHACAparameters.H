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
    IHACAparameters
Description
    Some parameters for input output informations
SourceFiles
    IHACAparameters.H
\*---------------------------------------------------------------------------*/
#ifndef Parameters_H
#define Parameters_H

#include <iostream>
#include "fvCFD.H"
#include "ITHACAassert.H"

/// Class for the definition of some general parameters, the parameters must be defined from the file ITHACAdict inside the
/// system folder.
class ITHACAparameters
{
    public:

        ITHACAparameters(fvMesh& mesh, Time& localTime);

        /// type of eigensolver used in the eigenvalue decomposition can be either be eigen or spectra
        word eigensolver;

        /// precision of the output Market Matrix objects (i.e. reduced matrices, eigenvalues, ...)
        label precision;

        ///
        bool exportPython;
        bool exportMatlab;
        bool exportTxt;
        bool debug;
        bool warnings;

        /// type of output format can be fixed or scientific
        std::_Ios_Fmtflags outytpe;

        /// Dictionary for input objects from file
        IOdictionary* ITHACAdict;

        /// runTime defined locally
        Time& runTime;

        /// Mesh object defined locally
        fvMesh& mesh;

        /// Pointer to the ITHACAparameters class
        static ITHACAparameters* instance;

    public:
        ///
        /// @brief      Gets an instance of ITHACAparameters, to be used if the instance is not existing.
        ///
        /// @param      mesh       The mesh
        /// @param      localTime  The runtime object
        ///
        /// @return     The instance.
        ///
        static ITHACAparameters* getInstance(fvMesh& mesh, Time& localTime);

        ///
        /// @brief      Gets an instance of ITHACAparameters, to be used if the instance is already existing.
        ///
        /// @return     The instance.
        ///
        static ITHACAparameters* getInstance();


        /// Delete empty constructor
        ~ITHACAparameters() = delete;


    private:

};

#endif
