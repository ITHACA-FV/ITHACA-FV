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

#include "DEIM.H"

// Template function constructor
template<typename T>
DEIM<T>::DEIM (PtrList<T>& s, label MaxModes, word FunctionName, word FieldName)
    :
    SnapShotsMatrix(s),
    MaxModes(MaxModes),
    FunctionName(FunctionName)
{
    ITHACAparameters* para(ITHACAparameters::getInstance());
    Folder = "ITHACAoutput/DEIM/" + FunctionName;
    magicPoints = autoPtr<IOList<label>>
                  (
                      new IOList<label>
                      (
                          IOobject
                          (
                              "magicPoints",
                              para->runTime.time().constant(),
                              "../" + Folder,
                              para->mesh,
                              IOobject::READ_IF_PRESENT,
                              IOobject::NO_WRITE
                          )
                      )
                  );
    modes = ITHACAPOD::DEIMmodes(SnapShotsMatrix, MaxModes, FunctionName,
                                 FieldName);

    if (!magicPoints().headerOk())
    {
        Eigen::MatrixXd A;
        Eigen::VectorXd b;
        Eigen::VectorXd c;
        Eigen::VectorXd r;
        Eigen::VectorXd rho(1);
        MatrixModes = Foam2Eigen::PtrList2Eigen(modes);
        label ind_max, c1;
        double max = MatrixModes.cwiseAbs().col(0).maxCoeff(&ind_max, &c1);
        rho(0) = max;
        magicPoints().append(ind_max);
        U = MatrixModes.col(0);
        P.resize(MatrixModes.rows(), 1);
        P.insert(ind_max, 0) = 1;

        for (label i = 1; i < MaxModes; i++)
        {
            A = P.transpose() * U;
            b = P.transpose() * MatrixModes.col(i);
            c = A.fullPivLu().solve(b);
            r = MatrixModes.col(i) - U * c;
            max = r.cwiseAbs().maxCoeff(&ind_max, &c1);
            P.conservativeResize(MatrixModes.rows(), i + 1);
            P.insert(ind_max, i) = 1;
            U.conservativeResize(MatrixModes.rows(), i + 1);
            U.col(i) =  MatrixModes.col(i);
            rho.conservativeResize(i + 1);
            rho(i) = max;
            magicPoints().append(ind_max);
        }

        MatrixOnline = U * ((P.transpose() * U).fullPivLu().inverse());
        mkDir(Folder);
        cnpy::save(MatrixOnline, Folder + "/MatrixOnline.npy");
        magicPoints().write();
    }
    else
    {
        cnpy::load(MatrixOnline, Folder + "/MatrixOnline.npy");
    }
}


template<typename T>
DEIM<T>::DEIM (PtrList<T>& s, label MaxModesA, label MaxModesB, word MatrixName)
    :
    SnapShotsMatrix(s),
    MaxModesA(MaxModesA),
    MaxModesB(MaxModesB),
    MatrixName(MatrixName),
    runSubMesh(false),
    runSubMeshA(false),
    runSubMeshB(false)

{
    ITHACAparameters* para(ITHACAparameters::getInstance());
    FolderM = "ITHACAoutput/DEIM/" + MatrixName;
    magicPointsArow = autoPtr<IOList<label>>
                      (
                          new IOList<label>
                          (
                              IOobject
                              (
                                  "magicPointsArow",
                                  para->runTime.time().constant(),
                                  "../" + FolderM,
                                  para->mesh,
                                  IOobject::READ_IF_PRESENT,
                                  IOobject::NO_WRITE
                              )
                          )
                      );
    magicPointsAcol = autoPtr<IOList<label>>
                      (
                          new IOList<label>
                          (
                              IOobject
                              (
                                  "magicPointsAcol",
                                  para->runTime.time().constant(),
                                  "../" + FolderM,
                                  para->mesh,
                                  IOobject::READ_IF_PRESENT,
                                  IOobject::NO_WRITE
                              )
                          )
                      );
    magicPointsB = autoPtr<IOList<label>>
                   (
                       new IOList<label>
                       (
                           IOobject
                           (
                               "magicPointsB",
                               para->runTime.time().constant(),
                               "../" + FolderM,
                               para->mesh,
                               IOobject::READ_IF_PRESENT,
                               IOobject::NO_WRITE
                           )
                       )
                   );
    xyz_Arow = autoPtr<IOList<label>>
               (
                   new IOList<label>
                   (
                       IOobject
                       (
                           "xyz_Arow",
                           para->runTime.time().constant(),
                           "../" + FolderM,
                           para->mesh,
                           IOobject::READ_IF_PRESENT,
                           IOobject::NO_WRITE
                       )
                   )
               );
    xyz_Acol = autoPtr<IOList<label>>
               (
                   new IOList<label>
                   (
                       IOobject
                       (
                           "xyz_Acol",
                           para->runTime.time().constant(),
                           "../" + FolderM,
                           para->mesh,
                           IOobject::READ_IF_PRESENT,
                           IOobject::NO_WRITE
                       )
                   )
               );
    xyz_B = autoPtr<IOList<label>>
            (
                new IOList<label>
                (
                    IOobject
                    (
                        "xyz_B",
                        para->runTime.time().constant(),
                        "../" + FolderM,
                        para->mesh,
                        IOobject::READ_IF_PRESENT,
                        IOobject::NO_WRITE
                    )
                )
            );

    if (!(magicPointsArow().headerOk() && magicPointsAcol().headerOk() &&
            magicPointsB().headerOk() && xyz_Arow().headerOk() &&
            xyz_Acol().headerOk() && xyz_B().headerOk()))
    {
        Eigen::MatrixXd AA;
        Eigen::VectorXd bA;
        Eigen::MatrixXd cA;
        Eigen::SparseMatrix<double> rA;
        Eigen::VectorXd rhoA(1);
        Matrix_Modes = ITHACAPOD::DEIMmodes(SnapShotsMatrix, MaxModesA, MaxModesB,
                                            MatrixName);
        Ncells = getNcells(std::get<1>(Matrix_Modes)[0].rows());
        label ind_rowA, ind_colA, xyz_rowA, xyz_colA;
        ind_rowA = ind_colA = xyz_rowA = xyz_colA = 0;
        double maxA = EigenFunctions::max(std::get<0>(Matrix_Modes)[0], ind_rowA,
                                          ind_colA);
        label ind_rowAOF = ind_rowA;
        label ind_colAOF = ind_colA;
        check3DIndices(ind_rowAOF, ind_colAOF, xyz_rowA, xyz_colA);
        xyz_Arow().append(xyz_rowA);
        xyz_Acol().append(xyz_colA);
        rhoA(0) = maxA;
        magicPointsArow().append(ind_rowAOF);
        magicPointsAcol().append(ind_colAOF);
        UA.append(std::get<0>(Matrix_Modes)[0]);
        Eigen::SparseMatrix<double> Pnow(std::get<0>(Matrix_Modes)[0].rows(),
                                         std::get<0>(Matrix_Modes)[0].cols());
        Pnow.insert(ind_rowA, ind_colA) = 1;
        PA.append(Pnow);

        for (label i = 1; i < MaxModesA; i++)
        {
            AA = EigenFunctions::innerProduct(PA, UA);
            bA = EigenFunctions::innerProduct(PA, std::get<0>(Matrix_Modes)[i]);
            cA = AA.fullPivLu().solve(bA);
            rA = std::get<0>(Matrix_Modes)[i] - EigenFunctions::MVproduct(UA, cA);
            double maxA = EigenFunctions::max(rA, ind_rowA, ind_colA);
            rhoA.conservativeResize(i + 1);
            rhoA(i) = maxA;
            label ind_rowAOF = ind_rowA;
            label ind_colAOF = ind_colA;
            check3DIndices(ind_rowAOF, ind_colAOF, xyz_rowA, xyz_colA);
            xyz_Arow().append(xyz_rowA);
            xyz_Acol().append(xyz_colA);
            magicPointsArow().append(ind_rowAOF);
            magicPointsAcol().append(ind_colAOF);
            UA.append(std::get<0>(Matrix_Modes)[i]);
            Eigen::SparseMatrix<double> Pnow(std::get<0>(Matrix_Modes)[0].rows(),
                                             std::get<0>(Matrix_Modes)[0].cols());
            Pnow.insert(ind_rowA, ind_colA) = 1;
            PA.append(Pnow);
        }

        Eigen::MatrixXd Aaux = EigenFunctions::innerProduct(PA,
                               UA).fullPivLu().inverse();
        MatrixOnlineA = EigenFunctions::MMproduct(UA, Aaux);
        Eigen::MatrixXd AB;
        Eigen::VectorXd bB;
        Eigen::VectorXd cB;
        Eigen::VectorXd rB;
        Eigen::VectorXd rhoB(1);
        label ind_rowB, xyz_rowB, c1;
        double maxB = std::get<1>(Matrix_Modes)[0].cwiseAbs().maxCoeff(&ind_rowB,
                      &c1);
        label ind_rowBOF = ind_rowB;
        check3DIndices(ind_rowBOF, xyz_rowB);
        rhoB(0) = maxB;
        xyz_B().append(xyz_rowB);
        magicPointsB().append(ind_rowBOF);
        UB = std::get<1>(Matrix_Modes)[0];
        PB.resize(UB.rows(), 1);
        PB.insert(ind_rowB, 0) = 1;

        for (label i = 1; i < MaxModesB; i++)
        {
            AB = PB.transpose() * UB;
            bB = PB.transpose() * std::get<1>(Matrix_Modes)[i];
            cB = AB.fullPivLu().solve(bB);
            rB = std::get<1>(Matrix_Modes)[i] - UB * cB;
            maxB = rB.cwiseAbs().maxCoeff(&ind_rowB, &c1);
            ind_rowBOF = ind_rowB;
            check3DIndices(ind_rowBOF, xyz_rowB);
            xyz_B().append(xyz_rowB);
            PB.conservativeResize((std::get<1>(Matrix_Modes)[i]).size(), i + 1);
            PB.insert(ind_rowB, i) = 1;
            UB.conservativeResize((std::get<1>(Matrix_Modes)[i]).size(), i + 1);
            UB.col(i) =  std::get<1>(Matrix_Modes)[i];
            rhoB.conservativeResize(i + 1);
            rhoB(i) = maxB;
            magicPointsB().append(ind_rowBOF);
        }

        if (MaxModesB == 1 && std::get<1>(Matrix_Modes)[0].norm() < 1e-8)
        {
            MatrixOnlineB = Eigen::MatrixXd::Zero(std::get<1>(Matrix_Modes)[0].rows(), 1);
        }
        else if (MaxModesB != 1)
        {
            MatrixOnlineB = UB * ((PB.transpose() * UB).fullPivLu().inverse());
        }
        else
        {
            Eigen::MatrixXd aux = PB.transpose() * UB;
            MatrixOnlineB = UB * 1 / aux(0, 0);
        }

        mkDir(FolderM + "/lhs");
        mkDir(FolderM + "/rhs");
        ITHACAstream::save(MatrixOnlineA, FolderM, "/lhs/MatrixOnlineA");
        cnpy::save(MatrixOnlineB, FolderM + "/rhs/MatrixOnlineB.npy");
        magicPointsArow().write();
        magicPointsAcol().write();
        magicPointsB().write();
        xyz_Arow().write();
        xyz_Acol().write();
        xyz_B().write();
    }
    else
    {
        ITHACAstream::load(MatrixOnlineA, FolderM, "/lhs/MatrixOnlineA");
        cnpy::load(MatrixOnlineB, FolderM + "/rhs/MatrixOnlineB.npy");
    }
}


template<typename T>
template<typename S>
S DEIM<T>::generateSubmesh(label layers, fvMesh& mesh, S field,
                           label secondTime)
{
    ITHACAparameters* para(ITHACAparameters::getInstance());
    totalMagicPoints = autoPtr<IOList<labelList>>
                       (
                           new IOList<labelList>
                           (
                               IOobject
                               (
                                   "totalMagicPoints",
                                   para->runTime.time().constant(),
                                   "../" + Folder,
                                   para->mesh,
                                   IOobject::READ_IF_PRESENT,
                                   IOobject::NO_WRITE
                               )
                           )
                       );
    uniqueMagicPoints = autoPtr<IOList<label>>
                        (
                            new IOList<label>
                            (
                                IOobject
                                (
                                    "uniqueMagicPoints",
                                    para->runTime.time().constant(),
                                    "../" + Folder,
                                    para->mesh,
                                    IOobject::READ_IF_PRESENT,
                                    IOobject::NO_WRITE
                                )
                            )
                        );
    volScalarField Indici
    (
        IOobject
        (
            FunctionName + "_indices",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(FunctionName + "_indices", dimensionSet(0, 0, 0, 0, 0, 0, 0),
                          0.0)
    );
    submesh = autoPtr<fvMeshSubset>(new fvMeshSubset(mesh));

    if (!totalMagicPoints().headerOk())
    {
        List<label> indices;

        for (label i = 0; i < magicPoints().size(); i++)
        {
            indices = ITHACAutilities::getIndices(mesh, magicPoints()[i], layers);
            totalMagicPoints().append(indices);
        }

        uniqueMagicPoints() = ITHACAutilities::combineList(totalMagicPoints());
    }

#if OPENFOAM >= 1812
    submesh->setCellSubset(uniqueMagicPoints());
#else
    submesh->setLargeCellSubset(uniqueMagicPoints());
#endif
    submesh->subMesh().fvSchemes::readOpt() = mesh.fvSchemes::readOpt();
    submesh->subMesh().fvSolution::readOpt() = mesh.fvSolution::readOpt();
    submesh->subMesh().fvSchemes::read();
    submesh->subMesh().fvSolution::read();
    std::cout.clear();
    S f = submesh->interpolate(field);
    scalar zerodot25 = 0.25;
    ITHACAutilities::assignIF(Indici, zerodot25,
                              uniqueMagicPoints().List<label>::clone()());
    ITHACAutilities::assignONE(Indici, magicPoints());

    if (!secondTime)
    {
        localMagicPoints = global2local(magicPoints(), submesh());
        ITHACAstream::exportSolution(Indici, "1", "./ITHACAoutput/DEIM/" + FunctionName
                                    );
    }

    totalMagicPoints().write();
    uniqueMagicPoints().write();
    return f;
}

template<typename T>
template<typename S>
S DEIM<T>::generateSubmeshMatrix(label layers, fvMesh& mesh, S field,
                                 label secondTime)
{
    ITHACAparameters* para(ITHACAparameters::getInstance());
    totalMagicPointsA = autoPtr<IOList<labelList>>
                        (
                            new IOList<labelList>
                            (
                                IOobject
                                (
                                    "totalMagicPointsA",
                                    para->runTime.time().constant(),
                                    "../" + FolderM,
                                    para->mesh,
                                    IOobject::READ_IF_PRESENT,
                                    IOobject::NO_WRITE
                                )
                            )
                        );
    uniqueMagicPointsA = autoPtr<IOList<label>>
                         (
                             new IOList<label>
                             (
                                 IOobject
                                 (
                                     "uniqueMagicPointsA",
                                     para->runTime.time().constant(),
                                     "../" + FolderM,
                                     para->mesh,
                                     IOobject::READ_IF_PRESENT,
                                     IOobject::NO_WRITE
                                 )
                             )
                         );
    List<label> indices;
    volScalarField Indici
    (
        IOobject
        (
            MatrixName + "_A_indices",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(MatrixName + "_A_indices", dimensionSet(0, 0, 0, 0, 0, 0, 0),
                          Foam::scalar(0))
    );
    submeshA = autoPtr<fvMeshSubset>(new fvMeshSubset(mesh));

    for (label i = 0; i < magicPointsArow().size(); i++)
    {
        indices = ITHACAutilities::getIndices(mesh, magicPointsArow()[i], layers);
        indices.append(ITHACAutilities::getIndices(mesh, magicPointsAcol()[i],
                       layers));
        totalMagicPointsA().append(indices);
    }

    uniqueMagicPointsA() = ITHACAutilities::combineList(totalMagicPointsA());
#if OPENFOAM >= 1812
    submeshA->setCellSubset(uniqueMagicPointsA());
#else
    submeshA->setLargeCellSubset(uniqueMagicPointsA());
#endif
    submeshA->subMesh().fvSchemes::readOpt() = mesh.fvSchemes::readOpt();
    submeshA->subMesh().fvSolution::readOpt() = mesh.fvSolution::readOpt();
    submeshA->subMesh().fvSchemes::read();
    submeshA->subMesh().fvSolution::read();
    std::cout.clear();
    S f = submeshA->interpolate(field);
    scalar zerodot25 = 0.25;
    ITHACAutilities::assignIF(Indici, zerodot25,
                              uniqueMagicPointsA().List<label>::clone()());
    ITHACAutilities::assignONE(Indici, magicPointsArow());
    ITHACAutilities::assignONE(Indici, magicPointsAcol());

    if (!secondTime)
    {
        localMagicPointsArow = global2local(magicPointsArow(), submeshA());
        localMagicPointsAcol = global2local(magicPointsAcol(), submeshA());
        ITHACAstream::exportSolution(Indici, "1", "./ITHACAoutput/DEIM/" + MatrixName
                                    );
        totalMagicPointsA().write();
        uniqueMagicPointsA().write();
    }

    runSubMeshA = true;
    return f;
}

template<typename T>
template<typename S>
S DEIM<T>::generateSubmeshVector(label layers, fvMesh& mesh, S field,
                                 label secondTime)
{
    ITHACAparameters* para(ITHACAparameters::getInstance());
    totalMagicPointsB = autoPtr<IOList<labelList>>
                        (
                            new IOList<labelList>
                            (
                                IOobject
                                (
                                    "totalMagicPointsB",
                                    para->runTime.time().constant(),
                                    "../" + FolderM,
                                    para->mesh,
                                    IOobject::READ_IF_PRESENT,
                                    IOobject::NO_WRITE
                                )
                            )
                        );
    uniqueMagicPointsB = autoPtr<IOList<label>>
                         (
                             new IOList<label>
                             (
                                 IOobject
                                 (
                                     "uniqueMagicPointsB",
                                     para->runTime.time().constant(),
                                     "../" + FolderM,
                                     para->mesh,
                                     IOobject::READ_IF_PRESENT,
                                     IOobject::NO_WRITE
                                 )
                             )
                         );
    List<label> indices;
    volScalarField Indici
    (
        IOobject
        (
            MatrixName + "_B_indices",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(MatrixName + "_B_indices", dimensionSet(0, 0, 0, 0, 0, 0, 0),
                          Foam::scalar(0))
    );
    submeshB = autoPtr<fvMeshSubset>(new fvMeshSubset(mesh));

    for (label i = 0; i < magicPointsB().size(); i++)
    {
        indices = ITHACAutilities::getIndices(mesh, magicPointsB()[i], layers);
        totalMagicPointsB().append(indices);

        if (!secondTime)
        {
            ITHACAutilities::assignONE(Indici, indices);
        }
    }

    uniqueMagicPointsB() = ITHACAutilities::combineList(totalMagicPointsB());
    std::cout.setstate(std::ios_base::failbit);
#if OPENFOAM >= 1812
    submeshB->setCellSubset(uniqueMagicPointsB());
#else
    submeshB->setLargeCellSubset(uniqueMagicPointsB());
#endif
    submeshB->subMesh().fvSchemes::readOpt() = mesh.fvSchemes::readOpt();
    submeshB->subMesh().fvSolution::readOpt() = mesh.fvSolution::readOpt();
    submeshB->subMesh().fvSchemes::read();
    submeshB->subMesh().fvSolution::read();
    std::cout.clear();
    S f = submeshB->interpolate(field);
    scalar zerodot25 = 0.25;
    ITHACAutilities::assignIF(Indici, zerodot25,
                              uniqueMagicPointsB().List<label>::clone()());
    ITHACAutilities::assignONE(Indici, magicPointsB());

    if (!secondTime)
    {
        localMagicPointsB = global2local(magicPointsB(), submeshB());
        ITHACAstream::exportSolution(Indici, "1", "./ITHACAoutput/DEIM/" + MatrixName
                                    );
        totalMagicPointsB().write();
        uniqueMagicPointsB().write();
    }

    runSubMeshB = true;
    return f;
}


template<typename T>
List<label> DEIM<T>::global2local(List<label>& points,
                                  fvMeshSubset& submesh)
{
    List<label> localPoints;

    for (label i = 0; i < points.size(); i++)
    {
        for (label j = 0; j < submesh.cellMap().size(); j++)
        {
            if (submesh.cellMap()[j] == points[i])
            {
                localPoints.append(j);
                break;
            }
        }
    }

    return localPoints;
}

template<typename T>
void DEIM<T>::check3DIndices(label& ind_rowA, label&  ind_colA, label& xyz_rowA,
                             label& xyz_colA)
{
    if (ind_rowA < Ncells)
    {
        xyz_rowA = 0;
    }
    else if (ind_rowA < Ncells * 2)
    {
        xyz_rowA = 1;
        ind_rowA = ind_rowA - Ncells;
    }
    else
    {
        xyz_rowA = 2;
        ind_rowA = ind_rowA - 2 * Ncells;
    }

    if (ind_colA < Ncells )
    {
        xyz_colA = 0;
    }
    else if (ind_colA < Ncells * 2)
    {
        xyz_colA = 1;
        ind_colA = ind_colA - 2 * Ncells;
    }
    else
    {
        xyz_colA = 2;
        ind_colA = ind_colA - 2 * Ncells;
    }
};

template<typename T>
void DEIM<T>::check3DIndices(label& ind_rowA, label& xyz_rowA)
{
    if (ind_rowA < Ncells)
    {
        xyz_rowA = 0;
    }
    else if (ind_rowA < Ncells * 2)
    {
        xyz_rowA = 1;
        ind_rowA = ind_rowA - Ncells;
    }
    else
    {
        xyz_rowA = 2;
        ind_rowA = ind_rowA - 2 * Ncells;
    }
};

template<class T>
template <class F>
F DEIM<T>::generateSubFieldMatrix(F& field)
{
    M_Assert(runSubMeshA == true,
             "You have to compute the magicPointsA before calling this function, try to rerun generateSubmeshMatrix");
    F f = submeshA().interpolate(field);
    return f;
}

template<typename T>
template <class F>
F DEIM<T>::generateSubFieldVector(F& field)
{
    M_Assert(runSubMeshB == true,
             "You have to compute the magicPointsB before calling this function, try to rerun generateSubmeshVector");
    F f = submeshB().interpolate(field);
    return f;
}


template<> label DEIM<fvScalarMatrix>::getNcells(label sizeM)
{
    label Ncells = sizeM;
    return Ncells;
}

template<> label DEIM<fvVectorMatrix>::getNcells(label sizeM)
{
    label Ncells = sizeM / 3;
    return Ncells;
}

// Specialization of the constructor
template DEIM<fvScalarMatrix>::DEIM (PtrList<fvScalarMatrix>& s,
                                     label MaxModesA,
                                     label MaxModesB, word MatrixName);
template DEIM<fvVectorMatrix>::DEIM (PtrList<fvVectorMatrix>& s,
                                     label MaxModesA,
                                     label MaxModesB, word MatrixName);
template DEIM<volScalarField>::DEIM(PtrList<volScalarField>& s, label MaxModes,
                                    word FunctionName, word FieldName);
template DEIM<volVectorField>::DEIM(PtrList<volVectorField>& s, label MaxModes,
                                    word FunctionName, word FieldName);

// Specialization for generateSubFieldMatrix
template volScalarField DEIM<fvScalarMatrix>::generateSubFieldMatrix(
    volScalarField& field);
template volVectorField DEIM<fvScalarMatrix>::generateSubFieldMatrix(
    volVectorField& field);
template surfaceScalarField
DEIM<fvScalarMatrix>::generateSubFieldMatrix(surfaceScalarField& field);
template surfaceVectorField
DEIM<fvScalarMatrix>::generateSubFieldMatrix(surfaceVectorField& field);
template volScalarField DEIM<fvVectorMatrix>::generateSubFieldMatrix(
    volScalarField& field);
template volVectorField DEIM<fvVectorMatrix>::generateSubFieldMatrix(
    volVectorField& field);
template surfaceScalarField
DEIM<fvVectorMatrix>::generateSubFieldMatrix(surfaceScalarField& field);
template surfaceVectorField
DEIM<fvVectorMatrix>::generateSubFieldMatrix(surfaceVectorField& field);

// Specialization for generateSubFieldVector
template volScalarField DEIM<fvScalarMatrix>::generateSubFieldVector(
    volScalarField& field);
template volVectorField DEIM<fvScalarMatrix>::generateSubFieldVector(
    volVectorField& field);
template surfaceScalarField
DEIM<fvScalarMatrix>::generateSubFieldVector(surfaceScalarField& field);
template surfaceVectorField
DEIM<fvScalarMatrix>::generateSubFieldVector(surfaceVectorField& field);
template volScalarField DEIM<fvVectorMatrix>::generateSubFieldVector(
    volScalarField& field);
template volVectorField DEIM<fvVectorMatrix>::generateSubFieldVector(
    volVectorField& field);
template surfaceScalarField
DEIM<fvVectorMatrix>::generateSubFieldVector(surfaceScalarField& field);
template surfaceVectorField
DEIM<fvVectorMatrix>::generateSubFieldVector(surfaceVectorField& field);

// Specialization for generateSubmesh
template volScalarField DEIM<volScalarField>::generateSubmesh(
    label layers, fvMesh& mesh, volScalarField field,
    label secondTime);
template volVectorField DEIM<volVectorField>::generateSubmesh(
    label layers, fvMesh& mesh, volVectorField field,
    label secondTime);
template volScalarField DEIM<volVectorField>::generateSubmesh(
    label layers, fvMesh& mesh, volScalarField field,
    label secondTime);
template volVectorField DEIM<volScalarField>::generateSubmesh(
    label layers, fvMesh& mesh, volVectorField field,
    label secondTime);

// Specialization for generateSubmeshVector
template volScalarField DEIM<fvScalarMatrix>::generateSubmeshVector(
    label layers, fvMesh& mesh, volScalarField field, label secondTime);
template volVectorField DEIM<fvScalarMatrix>::generateSubmeshVector(
    label layers, fvMesh& mesh, volVectorField field, label secondTime);
template surfaceScalarField
DEIM<fvScalarMatrix>::generateSubmeshVector(label layers, fvMesh& mesh,
        surfaceScalarField field, label secondTime);
template surfaceVectorField
DEIM<fvScalarMatrix>::generateSubmeshVector(label layers, fvMesh& mesh,
        surfaceVectorField field, label secondTime);
template volScalarField DEIM<fvVectorMatrix>::generateSubmeshVector(
    label layers, fvMesh& mesh, volScalarField field, label secondTime);
template volVectorField DEIM<fvVectorMatrix>::generateSubmeshVector(
    label layers, fvMesh& mesh, volVectorField field, label secondTime);
template surfaceScalarField
DEIM<fvVectorMatrix>::generateSubmeshVector(label layers, fvMesh& mesh,
        surfaceScalarField field, label secondTime);
template surfaceVectorField
DEIM<fvVectorMatrix>::generateSubmeshVector(label layers, fvMesh& mesh,
        surfaceVectorField field, label secondTime);

// Specialization for generateSubmeshMatrix
template volScalarField DEIM<fvScalarMatrix>::generateSubmeshMatrix(
    label layers, fvMesh& mesh, volScalarField field, label secondTime);
template volVectorField DEIM<fvScalarMatrix>::generateSubmeshMatrix(
    label layers, fvMesh& mesh, volVectorField field, label secondTime);
template surfaceScalarField
DEIM<fvScalarMatrix>::generateSubmeshMatrix(label layers, fvMesh& mesh,
        surfaceScalarField field, label secondTime);
template surfaceVectorField
DEIM<fvScalarMatrix>::generateSubmeshMatrix(label layers, fvMesh& mesh,
        surfaceVectorField field, label secondTime);
template volScalarField DEIM<fvVectorMatrix>::generateSubmeshMatrix(
    label layers, fvMesh& mesh, volScalarField field, label secondTime);
template volVectorField DEIM<fvVectorMatrix>::generateSubmeshMatrix(
    label layers, fvMesh& mesh, volVectorField field, label secondTime);
template surfaceScalarField
DEIM<fvVectorMatrix>::generateSubmeshMatrix(label layers, fvMesh& mesh,
        surfaceScalarField field, label secondTime);
template surfaceVectorField
DEIM<fvVectorMatrix>::generateSubmeshMatrix(label layers, fvMesh& mesh,
        surfaceVectorField field, label secondTime);

