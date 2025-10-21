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
    POD

Description
    Application to perform a POD decomposition of a general field on a given case

SourceFiles
    perform_POD.C

\*---------------------------------------------------------------------------*/

/// \file
/// \brief Application to perform POD on an already run case
/// \details In order to use this file one needs to prepare a ITHACAPODdict file, in order to 
/// check the syntax one needs to check the \ref ITHACAPODdict file.

/// \file ITHACAPODdict
/// \brief Example of a ITHACAPODdict file


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
#include "ITHACAPOD.H"
#include "ITHACAparameters.H"
#include "ITHACAstream.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<typename T>
void computeLift(PtrList<T>& Lfield,
                PtrList<T>& liftfield,
                PtrList<T>& omfield)
{
    for (label k = 0; k < liftfield.size(); k++)
    {
        for (label j = 0; j < Lfield.size(); j++)
        {
            if (k == 0)
            {
                autoPtr<T> p(new T("U", Lfield[j] - liftfield[k]));
                omfield.append(p);
            }
            else
            {
                autoPtr<T> p(new T("U", omfield[j] - liftfield[k]));
                omfield.set(j, p);
            }
        }
    }
}

int main(int argc, char *argv[])
{

#include "setRootCase.H"
#include "createTime.H"
#include "createMesh.H"

    PtrList<volVectorField> Vfield;
    PtrList<volScalarField> Sfield;
    PtrList<volVectorField> Vomfield;
    PtrList<volScalarField> Somfield;
    PtrList<volVectorField> Vmodes;
    PtrList<volScalarField> Smodes;

    ITHACAparameters* para = ITHACAparameters::getInstance(mesh,
                             runTime);

    bool pod_exist;
    struct stat sb;

    if (stat("./ITHACAoutput/POD", &sb) == 0 && S_ISDIR(sb.st_mode))
    {
        pod_exist = true;
    }
    else
    {
        pod_exist = false;
        Info << "POD don't exist, performing a POD decomposition" << endl;
        ITHACAutilities::createSymLink("./ITHACAoutput/POD");
    }
    if(pod_exist == 1)
    {
        Info << "The POD has already been performed, please delete the ITHACAoutput folder and try again." << endl;
        exit(0);
    }

    // Initialize Variables
    label nSnapshots;
    label startTime;
    label endTime;

    //Read FORCESdict
    IOdictionary ITHACAPODdict
    (
        IOobject
        (
            "ITHACAPODdict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Get times list from the case folder
    instantList Times = runTime.times();

    // Read Initial and last time from the POD dictionary
    const entry* existnsnap = ITHACAPODdict.lookupEntryPtr("Nsnapshots", false, true);
    const entry* existLT = ITHACAPODdict.lookupEntryPtr("FinalTime", false, true);

    // Check if the snapshots are lifted
    const bool lifted = ITHACAPODdict.lookupOrDefault<bool>("lifted", false);

    // Initiate variable from PODSolverDict
    if ((existnsnap) && (existLT))
    {
        Info << "Error you cannot define LatestTime and NSnapShots together" << endl;
        abort();
    }
    else if (existnsnap)
    {
        scalar InitialTime = ITHACAPODdict.lookupOrDefault<scalar>("InitialTime", 0);
        nSnapshots = readScalar(ITHACAPODdict.lookup("Nsnapshots"));
        startTime = Time::findClosestTimeIndex(runTime.times(), InitialTime);
        nSnapshots = min(nSnapshots , Times.size() - startTime);
        endTime = startTime + nSnapshots - 1;
        Info << nSnapshots << endl;
    }
    else
    {
        scalar InitialTime = ITHACAPODdict.lookupOrDefault<scalar>("InitialTime", 0);
        scalar FinalTime = ITHACAPODdict.lookupOrDefault<scalar>("FinalTime", 1000000);
        endTime = Time::findClosestTimeIndex(runTime.times(), FinalTime);
        startTime = Time::findClosestTimeIndex(runTime.times(), InitialTime);
        nSnapshots = endTime - startTime + 1;
        if (InitialTime > FinalTime)
        {
            Info << "FinalTime cannot be smaller than the InitialTime check your ITHACAPODdict file\n" << endl;
            abort();
        }
    }
    // Print out some Infos
    Info<< "startTime: " << Times[startTime].name() << "\n" 
        << "endTime: " << Times[endTime].name() << "\n" 
        << "nSnapshots: " << nSnapshots << "\n" << endl;

    // Set the initial time
    runTime.setTime(Times[startTime], startTime);

    wordList fieldlist
    (
        ITHACAPODdict.lookup("fields")
    );

    //word Name = ITHACAPODdict.lookup("fieldName");
    //word type = ITHACAPODdict.lookup("type");

    if (startTime == endTime)
    {
        Info << "The case has no snapshots to process, exiting the code" << endl;
        exit(0);
    }

    for (label k = 0; k < fieldlist.size(); k++)
    {
        dictionary& subDict = ITHACAPODdict.subDict(fieldlist[k]);
        scalar nmodes = readScalar(subDict.lookup("nmodes"));
        word field_name(subDict.lookup("field_name"));
        word field_type(subDict.lookup("field_type"));

        for (label i = startTime; i < endTime + 1; i++)
        {
            Info << "Reading snapshot " << Times[i].name() << " for field " << field_name << endl;
            runTime.setTime(Times[i], i);
            mesh.readUpdate();

            if (field_type == "vector")
            {
                autoPtr<volVectorField> p(new volVectorField
                (
                    IOobject
                    (
                        field_name,
                        runTime.timeName(),
                        mesh,
                        IOobject::MUST_READ
                    ),
                    mesh
                ));
                Vfield.append(p);
            }

            if (field_type == "scalar")
            {
                autoPtr<volScalarField> p(new volScalarField
                (
                    IOobject
                    (
                        field_name,
                        runTime.timeName(),
                        mesh,
                        IOobject::MUST_READ
                    ),
                    mesh
                ));
                Sfield.append(p);
            }
        }

        if (field_type == "vector")
        {
            if (lifted)
            {
                PtrList<volVectorField> liftFields;
                ITHACAstream::read_fields(liftFields, field_name, "./lift/");                
                computeLift<volVectorField>(Vfield, liftFields, Vomfield);

                ITHACAPOD::getModes(Vomfield, Vmodes, field_name, 0, 0, 0, nmodes, para->correctBC);
                Eigen::MatrixXd coeffs = ITHACAutilities::getCoeffs(Vomfield,
                                            Vmodes, nmodes);
                
                ITHACAstream::exportFields(Vomfield, "./ITHACAoutput/Offline", field_name+"omfield");
                ITHACAstream::exportMatrix(coeffs, field_name+"coeffs", "eigen",
                                       "./ITHACAoutput/Matrices/");
                Vfield.clear();
                Vmodes.clear();
                Vomfield.clear();
                liftFields.clear();
                Info << "Lifted POD modes computed for field " << field_name << endl;
            }
            else
            {
                ITHACAPOD::getModes(Vfield, Vmodes, field_name, 0, 0, 0, nmodes, para->correctBC);
                Eigen::MatrixXd coeffs = ITHACAutilities::getCoeffs(Vfield,
                                            Vmodes, nmodes);

                ITHACAstream::exportMatrix(coeffs, field_name+"coeffs", "eigen",
                                       "./ITHACAoutput/Matrices/");
                Vfield.clear();
                Vmodes.clear();
                Info << "POD modes computed for field " << field_name << endl;
            }
        }
        if (field_type == "scalar")
        {
            if (lifted)
            {
                PtrList<volScalarField> liftFields;
                ITHACAstream::read_fields(liftFields, field_name, "./lift/");
                computeLift<volScalarField>(Sfield, liftFields, Somfield);

                ITHACAPOD::getModes(Somfield, Smodes, field_name, 0, 0, 0, nmodes, para->correctBC);
                Eigen::MatrixXd coeffs = ITHACAutilities::getCoeffs(Somfield,
                                            Smodes, nmodes);

                ITHACAstream::exportFields(Somfield, "./ITHACAoutput/Offline", "Sofield");
                ITHACAstream::exportMatrix(coeffs, field_name+"coeffs", "eigen",
                                       "./ITHACAoutput/Matrices/");
                Sfield.clear();
                Smodes.clear();
                Somfield.clear();
                liftFields.clear();
                Info << "Lifted POD modes computed for field " << field_name << endl;
            }
            else
            {
                ITHACAPOD::getModes(Sfield, Smodes, field_name, 0, 0, 0, nmodes, para->correctBC);
                Eigen::MatrixXd coeffs = ITHACAutilities::getCoeffs(Sfield,
                                            Smodes, nmodes);
                ITHACAstream::exportMatrix(coeffs, field_name+"coeffs", "eigen",
                                       "./ITHACAoutput/Matrices/");
                Sfield.clear();
                Smodes.clear();
                Info << "POD modes computed for field " << field_name << endl;
            }
        }
    }
    Info << endl;
    Info << "End\n" << endl;
    return 0;
}


