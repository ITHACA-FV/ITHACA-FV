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
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    lift_and_drag

Author
    Giovanni Stabile, SISSA MathLab (International School for Advanced Studies) gstabile@sissa.it

Description
    Application to recover the lift and the drag after the simulation is performed

\*---------------------------------------------------------------------------*/

/// \file
/// \brief Application to recover the lift and the drag after the simulation is performed
/// \details In order to use this file one needs to prepare a FORCESdict file, in order to 
/// check the syntax one needs to check the \ref FORCESdict file.

/// \file FORCESdict
/// \brief Example of a FORCESdict file

#include "fvCFD.H"
#include "IOmanip.H"
#include "IFstream.H"
#include "primitiveFields.H"
#include "FieldFields.H"
#include "scalarMatrices.H"
#include "SortableList.H"
#include "volFieldsFwd.H"
#include "forces.H"
#include "forceCoeffs.H"
#include "volFields.H"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdio.h>
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{

#include "setRootCase.H"
#include "createTime.H"
#include "createMesh.H"

    Info << argv[0] << endl;

    instantList Times = runTime.times();

    runTime.setTime(Times[2], 2);

    //Read FORCESdict
    IOdictionary FORCESdict
    (
        IOobject
        (
            "FORCESdict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
            )
        );

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
            )
        );

    word pName(FORCESdict.lookup("pName"));
    word UName(FORCESdict.lookup("UName"));

    dictionary forcesDict; 
    forcesDict.add("type", functionObjects::forces::typeName);
    forcesDict.add("patches", FORCESdict.lookup("patches"));
    forcesDict.add("origin", FORCESdict.lookup("pitchAxis"));
    forcesDict.add("pitchAxis", FORCESdict.lookup("pitchAxis"));
    forcesDict.add("CofR", FORCESdict.lookup("CofR"));
    forcesDict.add("liftDir", FORCESdict.lookup("liftDir"));
    forcesDict.add("dragDir", FORCESdict.lookup("dragDir"));
    forcesDict.add("magUInf", FORCESdict.lookup("magUInf"));
    forcesDict.add("lRef", FORCESdict.lookup("lRef"));
    forcesDict.add("Aref", FORCESdict.lookup("Aref"));
    forcesDict.add("rhoInf", FORCESdict.lookup("rhoInf"));
    forcesDict.add("rho", FORCESdict.lookup("rho"));

    functionObjects::forceCoeffs fc("FC", runTime, forcesDict);
    functionObjects::forces f("Forces", mesh, forcesDict);

    for (label i = 2; i < Times.size(); i++)
    {
        runTime.setTime(Times[i], i);
        mesh.readUpdate();

        volVectorField U
        (
            IOobject
            (
                UName,
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ
                ),
            mesh
            );

        U.rename("U");

        volScalarField P
        (
            IOobject
            (
                pName,
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ
                ),
            mesh
            );

        P.rename("p");

        fc.execute();
        fc.write();
        f.write();
    }
    Info << endl;
    Info << "End\n" << endl;
    return 0;
}

