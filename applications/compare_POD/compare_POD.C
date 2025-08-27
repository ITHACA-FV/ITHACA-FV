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
#include "PODTemplate.H"
#include "PODParameters.H"
#include <chrono>
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class FieldType>
auto loadSnapshotHelper(
    const word& fieldName,
    const fileName& timeDir,
    const typename FieldType::Mesh& mesh
) -> FieldType*
{
    return new FieldType(
        IOobject(
            fieldName,
            timeDir,
            mesh,
            IOobject::MUST_READ
        ),
        mesh
    );
}



int main(int argc, char *argv[])
{

#include "setRootCase.H"
#include "createTime.H"
#include "createMesh.H"

    PtrList<volVectorField> Vfield;
    PtrList<volScalarField> Sfield;
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
        if (Pstream::myProcNo()==0)
        {
        mkDir("./ITHACAoutput/POD");
        system("ln -s ../../constant ./ITHACAoutput/POD/constant");
        system("ln -s ../../0 ./ITHACAoutput/POD/0");
        system("ln -s ../../system ./ITHACAoutput/POD/system");
        }
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
        scalar FinalTime = ITHACAPODdict.lookupOrDefault<scalar>("FinalTime", 100000000000000);
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
    Info << "startTime: " << startTime << "\n" << "endTime: " << endTime << "\n" << "nSnapshots: " << nSnapshots << "\n" << endl;

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
            Info << "Reading snapshot " << i << " for field " << field_name << endl;
            runTime.setTime(Times[i], i);
            mesh.readUpdate();

            if (field_type == "vector")
            {
                Vfield.append(
                    loadSnapshotHelper<volVectorField>(
                        field_name,
                        runTime.timeName(),
                        mesh
                    )
                );
            }

            if (field_type == "scalar")
            {
                Sfield.append(
                    loadSnapshotHelper<volScalarField>(
                        field_name,
                        runTime.timeName(),
                        mesh
                    )
                );
            }
        }

        ITHACAPOD::PODParameters* podparam;
        podparam = new ITHACAPOD::PODParameters(argc,argv) ;


        if (field_type == "vector")
        {
            PtrList<volVectorField> f_spatialModesU;
            Eigen::MatrixXd m_temporalModesU;
            Eigen::MatrixXd m_temporalModesUSimulation;
            Eigen::MatrixXd covMatrixU;

            auto t1 = std::chrono::high_resolution_clock::now();
            ITHACAPOD::PODTemplate <volVectorField> ithacaFVPODVector(podparam,field_name);
            ithacaFVPODVector.getModes(f_spatialModesU, m_temporalModesU,
            m_temporalModesUSimulation, covMatrixU);
            auto t2 = std::chrono::high_resolution_clock::now();

            /* Getting number of milliseconds as a double. */
            std::chrono::duration<double, std::milli> ms_double = t2 - t1;
            Info << "*********************************************************" << endl;
            Info << "POD-RL for the speed took " << ms_double.count() / 1000 << " s" << endl;
            Info << "*********************************************************" << endl;


            ITHACAPOD::getModes(Vfield, Vmodes, field_name, 0, 0, 0, nmodes);
        }
        if (field_type == "scalar")
        {
            auto t3 = std::chrono::high_resolution_clock::now();
            PtrList<volScalarField> f_spatialModesp;
            Eigen::MatrixXd m_temporalModesp;
            Eigen::MatrixXd m_temporalModespSimulation;
            Eigen::MatrixXd covMatrixp;

            ITHACAPOD::PODTemplate <volScalarField> ithacaFVPODScalar(podparam,field_name);
            ithacaFVPODScalar.getModes(f_spatialModesp, m_temporalModesp,
            m_temporalModespSimulation, covMatrixp);
            auto t4 = std::chrono::high_resolution_clock::now();

            /* Getting number of milliseconds as a double. */
            std::chrono::duration<double, std::milli> ms_double = t4 - t3;
            Info << "*********************************************************" << endl;
            Info << "POD-RL for the pressure took " << ms_double.count() / 1000 << " s" << endl;
            Info << "*********************************************************" << endl;

            ITHACAPOD::getModes(Sfield, Smodes, field_name, 0, 0, 0, nmodes);
            
        }

        Vfield.clear();
        Sfield.clear();
    }
    Info << endl;
    Info << "End\n" << endl;
    return 0;
}


