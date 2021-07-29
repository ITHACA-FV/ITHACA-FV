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
/// Source file of the inverseLaplacianProblem class.


#include "inverseLaplacianProblem.H"

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructors
inverseLaplacianProblem::inverseLaplacianProblem() {}

inverseLaplacianProblem::inverseLaplacianProblem(int argc, char* argv[])
    :
    DT("DT", dimensionSet(1, 1, -3, -1, 0, 0, 0), 1.0)
{
    _args = autoPtr<argList>
            (
                new argList(argc, argv)
            );

    if (!_args->checkRootCase())
    {
        Foam::FatalError.exit();
    }

    argList& args = _args();
#include "createTime.H"
#include "createMesh.H"
    _simple = autoPtr<simpleControl>
              (
                  new simpleControl
                  (
                      mesh
                  )
              );
    simpleControl& simple = _simple();
#include "createFields.H"
#include "createThermocouples.H"
    thermocouplesPos = TCpos;
#include "createFvOptions.H"
    ITHACAdict = new IOdictionary
    (
        IOobject
        (
            "ITHACAdict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    para = ITHACAparameters::getInstance(mesh, runTime);
    offline = ITHACAutilities::check_off();
    podex = ITHACAutilities::check_pod();
    nProcs = Pstream::nProcs();
}

// * * * * * * * * * * * * * * Full Order Methods * * * * * * * * * * * * * * //

void inverseLaplacianProblem::set_g()
{
    volScalarField& T = _T();
    g.resize(T.boundaryField()[hotSide_ind].size(), 0.0);
    forAll (T.boundaryField()[hotSide_ind], faceI)
    {
        g[faceI] = 0.0;
    }
}

volScalarField inverseLaplacianProblem::list2Field(List<scalar> list,
        scalar innerField)
{
    volScalarField& T = _T();
    fvMesh& mesh = _mesh();
    volScalarField field(T);
    ITHACAutilities::assignIF(field, innerField);
    //Access the mesh information for the boundary
    const polyPatch& cPatch = mesh.boundaryMesh()[hotSide_ind];
    //List of cells close to a boundary
    const labelUList& faceCells = cPatch.faceCells();
    M_Assert(faceCells.size() == list.size(), "Inpust list has the wrong size.");
    forAll(cPatch, faceI)
    {
        //id of the owner cell having the face
        label faceOwner = faceCells[faceI] ;
        field[faceOwner] = list[faceI];
    }
    return field;
}

void inverseLaplacianProblem::set_valueFraction()
{
    fvMesh& mesh = _mesh();
    valueFraction.resize(mesh.boundaryMesh()["coldSide"].size());
    homogeneousBCcoldSide.resize(mesh.boundaryMesh()["coldSide"].size());
    coldSide_ind = mesh.boundaryMesh().findPatchID("coldSide");
    Eigen::VectorXd faceCellDist =
        ITHACAutilities::boudaryFaceToCellDistance(mesh, coldSide_ind);
    forAll (valueFraction, faceI)
    {
        scalar faceDist = faceCellDist(faceI);
        valueFraction[faceI] =  1.0 / (1.0 + (k / H / faceDist));
        homogeneousBCcoldSide[faceI] =  0;
    }
    refGrad = homogeneousBCcoldSide;
}


void inverseLaplacianProblem::assignDirectBC()
{
    fvMesh& mesh = _mesh();
    volScalarField& T = _T();
    set_valueFraction();
    forAll(mesh.boundaryMesh(), patchI)
    {
        if (patchI == mesh.boundaryMesh().findPatchID("coldSide"))
        {
            ITHACAutilities::assignMixedBC(T, patchI, Tf, refGrad, valueFraction);
        }
        else if (patchI == mesh.boundaryMesh().findPatchID("hotSide"))
        {
            ITHACAutilities::assignBC(T, patchI, - g / k);
        }
        else
        {
            ITHACAutilities::assignBC(T, patchI, homogeneousBC);
        }
    }
}

void inverseLaplacianProblem::solveDirect()
{
    restart();
    assignDirectBC();
    solve("direct");
}

void inverseLaplacianProblem::solve(const char* problemID)
{
    volScalarField& T = _T();
    simpleControl& simple = _simple();
    Foam::Time& runTime = _runTime();

    if (strcmp( problemID, "direct") == 0)
    {
        ITHACAutilities::assignIF(T, homogeneousBC);
    }
    else
    {
        Info << "Problem name should be direct or sensitivity" << endl;
        exit(10);
    }

#if defined(OFVER) && (OFVER == 6)

    while (simple.loop(runTime))
#else
    while (simple.loop())
#endif
    {
        while (simple.correctNonOrthogonal())
        {
            if (strcmp( problemID, "direct") == 0)
            {
                fvScalarMatrix TEqn
                (
                    fvm::laplacian(DT, T)
                );
                TEqn.solve();
            }
        }
    }
}

void inverseLaplacianProblem::readThermocouples()
{
    if (!thermocouplesRead)
    {
        word fileName = "./thermocouplesCellsID";

        if (ITHACAutilities::check_file(fileName + "_mat.txt"))
        {
            Info << "Reading thermocouples cells from file" << endl;
            Eigen::MatrixXi TCmatrix = ITHACAstream::readMatrix(fileName + "_mat.txt").cast
                                       <int> ();
            thermocouplesCellID = Foam2Eigen::EigenMatrix2List(TCmatrix);
        }
        else
        {
            Info << "Defining positions of thermocouples" << endl;
            fvMesh& mesh = _mesh();
            volScalarField& T = _T();
            thermocouplesCellID.resize(thermocouplesPos.size());
            forAll(thermocouplesPos, tcI)
            {
                thermocouplesCellID[tcI] = mesh.findCell(thermocouplesPos[tcI]);
            }
            volScalarField thermocouplesField(T);
            ITHACAutilities::assignIF(thermocouplesField, homogeneousBC);
            forAll(thermocouplesCellID, tcI)
            {
                thermocouplesField.ref()[thermocouplesCellID[tcI]] = 1;
            }
            ITHACAstream::exportSolution(thermocouplesField, "1", "./ITHACAoutput/debug/",
                                         "thermocouplesField,");
            Eigen::MatrixXi thermocouplesCellID_eigen = Foam2Eigen::List2EigenMatrix(
                        thermocouplesCellID);
            ITHACAstream::exportMatrix(thermocouplesCellID_eigen, fileName,
                                       "eigen", "./");
        }

        thermocouplesRead = 1;
        thermocouplesNum = thermocouplesPos.size();
    }
    else
    {
        WarningInFunction << "readThermocouples function called twice." << endl;
        WarningInFunction << "I am not doing the second reading." << endl;
    }
}

Eigen::VectorXd inverseLaplacianProblem::fieldValueAtThermocouples(
    volScalarField& field)
{
    if (!thermocouplesRead)
    {
        readThermocouples();
    }

    fvMesh& mesh = _mesh();
    dictionary interpolationDict =
        mesh.solutionDict().subDict("interpolationSchemes");
    autoPtr<Foam::interpolation<scalar>> fieldInterp =
                                          Foam::interpolation<scalar>::New(interpolationDict, field);
    Eigen::VectorXd fieldInt;
    fieldInt.resize(thermocouplesPos.size());
    forAll(thermocouplesPos, tcI)
    {
        fieldInt(tcI) = fieldInterp->interpolate(thermocouplesPos[tcI],
                        thermocouplesCellID[tcI]);
    }
    return fieldInt;
}

void inverseLaplacianProblem::differenceBetweenDirectAndMeasure()
{
    volScalarField& T = _T();
    Tdirect = fieldValueAtThermocouples(T);
    Tdiff = Tdirect - Tmeas;
}


Foam::vector inverseLaplacianProblem::cellDim(const faceList& ff,
        const pointField& pp,
        const cell& cc, labelList pLabels, pointField pLocal)
{
    forAll (pLabels, pointi)
    pLocal[pointi] = pp[pLabels[pointi]];
    double  xDim = Foam::max(pLocal & Foam::vector(1, 0, 0))
                   - Foam::min(pLocal & Foam::vector(1, 0, 0));
    double  yDim = Foam::max(pLocal & Foam::vector(0, 1, 0))
                   - Foam::min(pLocal & Foam::vector(0, 1, 0));
    double  zDim = Foam::max(pLocal & Foam::vector(0, 0, 1))
                   - Foam::min(pLocal & Foam::vector(0, 0, 1));
    Foam::vector dim (xDim, yDim, zDim);
    return dim;
}


void inverseLaplacianProblem::restart()
{
    _simple.clear();
    _T.clear();
    argList& args = _args();
    Time& runTime = _runTime();
    //Reinitializing runTime
    instantList Times = runTime.times();
    runTime.setTime(Times[1], 1);
    Foam::fvMesh& mesh = _mesh();
    _simple = autoPtr<simpleControl>
              (
                  new simpleControl
                  (
                      mesh
                  )
              );
    _T = autoPtr<volScalarField>
         (
             new volScalarField
             (
                 IOobject
                 (
                     "T",
                     runTime.timeName(),
                     mesh,
                     IOobject::MUST_READ,
                     IOobject::AUTO_WRITE
                 ),
                 mesh
             )
         );
    Info << "Ready for new computation" << endl;
}

