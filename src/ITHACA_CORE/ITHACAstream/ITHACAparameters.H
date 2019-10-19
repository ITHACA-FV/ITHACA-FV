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

/// Class for the definition of some general parameters, the parameters must be defined from the file ITHACAdict inside the
/// system folder.
class ITHACAparameters
{
    public:
        ITHACAparameters()
        {
            Foam::Time runTime(Foam::Time::controlDictName, "./", "");
            Foam::fvMesh mesh
            (
                Foam::IOobject
                (
                    Foam::fvMesh::defaultRegion,
                    runTime.timeName(),
                    runTime,
                    Foam::IOobject::MUST_READ
                )
            );
            ITHACAdict = new IOdictionary
            (
                IOobject
                (
                    "ITHACAdict",
                    "./system",
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );
            precision = ITHACAdict->lookupOrDefault<int>("OutPrecision", 10);
            word typeout = ITHACAdict->lookupOrDefault<word>("OutType", "fixed");

            if (typeout == "fixed")
            {
                outytpe = std::ios_base::fixed;
            }

            if (typeout == "scientific")
            {
                outytpe = std::ios_base::scientific;
            }

            eigensolver = ITHACAdict->lookupOrDefault<word>("EigenSolver", "spectra");
            exportPython = ITHACAdict->lookupOrDefault<bool>("exportPython", 0);
            exportMatlab = ITHACAdict->lookupOrDefault<bool>("exportMatlab", 0);
            exportTxt = ITHACAdict->lookupOrDefault<bool>("exportTxt", 0);
            debug = ITHACAdict->lookupOrDefault<bool>("debug", 0);
        }
        /// type of eigensolver used in the eigenvalue decomposition can be either be eigen or spectra
        word eigensolver;

        /// precision of the output Market Matrix objects (i.e. reduced matrices, eigenvalues, ...)
        int precision;

        ///
        bool exportPython;
        bool exportMatlab;
        bool exportTxt;
        bool debug;

        /// type of output format can be fixed or scientific
        std::_Ios_Fmtflags outytpe;

        /// Dictionary for input objects from file
        IOdictionary* ITHACAdict;
    private:

};

#endif
