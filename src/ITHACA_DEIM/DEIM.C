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
    Pair <label> indA(ind_rowAOF, ind_colAOF);
    Pair <label> xyzA(xyz_rowA, xyz_colA);
    xyz_A.append(xyzA);
    rhoA(0) = maxA;
    magicPointsA.append(indA);
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
        Pair <label> indA(ind_rowAOF, ind_colAOF);
        Pair <label> xyzA(xyz_rowA, xyz_colA);
        xyz_A.append(xyzA);
        magicPointsA.append(indA);
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
    xyz_B.append(xyz_rowB);
    magicPointsB.append(ind_rowBOF);
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
        xyz_B.append(xyz_rowB);
        PB.conservativeResize((std::get<1>(Matrix_Modes)[i]).size(), i + 1);
        PB.insert(ind_rowB, i) = 1;
        UB.conservativeResize((std::get<1>(Matrix_Modes)[i]).size(), i + 1);
        UB.col(i) =  std::get<1>(Matrix_Modes)[i];
        rhoB.conservativeResize(i + 1);
        rhoB(i) = maxB;
        magicPointsB.append(ind_rowBOF);
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
}


template<typename T>
template<typename S>
S DEIM<T>::generateSubmeshes(label layers, fvMesh& mesh, S field,
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

        labelList a = ListListOps::combine<labelList>(totalMagicPoints(),
                      accessOp<labelList>());
#if OPENFOAM >= 1812
        inplaceUniqueSort(a);
#else
        labelList order;
        uniqueOrder(a, order);
        labelList b(order.size());

        for (label i = 0; i < order.size(); ++i)
        {
            b[i] = a[order[i]];
        }
        a.resize(order.size());
        a = b;
#endif
        uniqueMagicPoints() = a;
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
PtrList<S> DEIM<T>::generateSubmeshesMatrix(label layers, fvMesh& mesh, S field,
        label secondTime)
{
    fvMeshSubset* submesh;
    PtrList<S> fieldsA;
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

    if (!secondTime)
    {
        Indici = Indici * 0;
    }

    fieldsA.resize(0);

    for (label i = 0; i < magicPointsA.size(); i++)
    {
        submesh = new fvMeshSubset(mesh);
        indices = ITHACAutilities::getIndices(mesh, magicPointsA[i].first(), layers);
        indices.append(ITHACAutilities::getIndices(mesh, magicPointsA[i].second(),
                       layers));

        if (!secondTime)
        {
            ITHACAutilities::assignONE(Indici, indices);
        }

        std::cout.setstate(std::ios_base::failbit);
#if OPENFOAM >= 1812
        submesh->setCellSubset(indices);
#else
        submesh->setLargeCellSubset(indices);
#endif
        submesh->subMesh().fvSchemes::readOpt() = mesh.fvSchemes::readOpt();
        submesh->subMesh().fvSolution::readOpt() = mesh.fvSolution::readOpt();
        submesh->subMesh().fvSchemes::read();
        submesh->subMesh().fvSolution::read();
        std::cout.clear();
        S f = submesh->interpolate(field);
        fieldsA.append(tmp<S>(f));

        if (!secondTime)
        {
            submeshListA.append(submesh);
        }
    }

    if (!secondTime)
    {
        for (label i = 0; i < magicPointsA.size(); i++)
        {
            Indici.ref()[magicPointsA[i].first()] = 10;
            Indici.ref()[magicPointsA[i].second()] = 10;
        }
    }

    if (!secondTime)
    {
        localMagicPointsA = global2local(magicPointsA, submeshListA);
        ITHACAstream::exportSolution(Indici, "1", "./ITHACAoutput/DEIM/" + MatrixName
                                    );
    }

    runSubMeshA = true;
    return fieldsA;
}

template<typename T>
template<typename S>
PtrList<S> DEIM<T>::generateSubmeshesVector(label layers, fvMesh& mesh, S field,
        label secondTime)
{
    fvMeshSubset* submesh;
    List<label> indices;
    PtrList<S> fieldsB;
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
        dimensionedScalar(MatrixName + "_A_indices", dimensionSet(0, 0, 0, 0, 0, 0, 0),
                          Foam::scalar(0))
    );

    if (!secondTime)
    {
        Indici = Indici * 0;
    }

    fieldsB.resize(0);

    for (label i = 0; i < magicPointsB.size(); i++)
    {
        submesh = new fvMeshSubset(mesh);
        indices = ITHACAutilities::getIndices(mesh, magicPointsB[i], layers);

        if (!secondTime)
        {
            ITHACAutilities::assignONE(Indici, indices);
        }

        std::cout.setstate(std::ios_base::failbit);
#if OPENFOAM >= 1812
        submesh->setCellSubset(indices);
#else
        submesh->setLargeCellSubset(indices);
#endif
        submesh->subMesh().fvSchemes::readOpt() = mesh.fvSchemes::readOpt();
        submesh->subMesh().fvSolution::readOpt() = mesh.fvSolution::readOpt();
        submesh->subMesh().fvSchemes::read();
        submesh->subMesh().fvSolution::read();
        std::cout.clear();
        S f = submesh->interpolate(field);
        fieldsB.append(tmp<S>(f));

        if (!secondTime)
        {
            submeshListB.append(submesh);
        }
    }

    if (!secondTime)
    {
        for (label i = 0; i < magicPointsB.size(); i++)
        {
            Indici.ref()[magicPointsB[i]] = 10;
        }
    }

    if (!secondTime)
    {
        localMagicPointsB = global2local(magicPointsB, submeshListB);
        ITHACAstream::exportSolution(Indici, "1", "./ITHACAoutput/DEIM/" + MatrixName
                                    );
    }

    runSubMeshB = true;
    return fieldsB;
}


template<typename T>
List<label> DEIM<T>::global2local(List<label>& points,
                                  PtrList<fvMeshSubset>& submeshList)
{
    List<label> localPoints;

    for (label i = 0; i < points.size(); i++)
    {
        for (label j = 0; j < submeshList[i].cellMap().size(); j++)
        {
            if (submeshList[i].cellMap()[j] == points[i])
            {
                localPoints.append(j);
                break;
            }
        }
    }

    return localPoints;
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
List<Pair <label>> DEIM<T>::global2local(List<Pair <label>>& points,
                PtrList<fvMeshSubset>& submeshList)
{
    List<Pair <label>> localPoints(points.size());

    for (label i = 0; i < points.size(); i++)
    {
        for (label j = 0; j < submeshList[i].cellMap().size(); j++)
        {
            if (submeshList[i].cellMap()[j] == points[i].first())
            {
                localPoints[i].first() = j;
                break;
            }
        }

        for (label j = 0; j < submeshList[i].cellMap().size(); j++)
        {
            if (submeshList[i].cellMap()[j] == points[i].second())
            {
                localPoints[i].second() = j;
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
PtrList<F> DEIM<T>::generateSubFieldsMatrix(F& field)
{
    PtrList<F> fields;
    M_Assert(runSubMeshA == true,
             "You have to compute the magicPointsA before calling this function, try to rerun generateSubmeshesMatrix");

    for (label i = 0; i < submeshListA.size(); i++)
    {
        F f = submeshListA[i].interpolate(field);
        fields.append(tmp<F>(f));
    }

    return fields;
}

template<typename T>
template <class F>
PtrList<F> DEIM<T>::generateSubFieldsVector(F& field)
{
    M_Assert(runSubMeshB == true,
             "You have to compute the magicPointsB before calling this function, try to rerun generateSubmeshesVector");
    PtrList<F> fields;

    for (label i = 0; i < submeshListB.size(); i++)
    {
        F f = submeshListB[i].interpolate(field);
        fields.append(tmp<F>(f));
    }

    return fields;
}


template<> label DEIM<fvScalarMatrix>::getNcells(label sizeM)
{
    label Ncells = sizeM;
    return sizeM;
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

// Specialization for generateSubFieldsMatrix
template PtrList<volScalarField> DEIM<fvScalarMatrix>::generateSubFieldsMatrix(
    volScalarField& field);
template PtrList<volVectorField> DEIM<fvScalarMatrix>::generateSubFieldsMatrix(
    volVectorField& field);
template PtrList<surfaceScalarField>
DEIM<fvScalarMatrix>::generateSubFieldsMatrix(surfaceScalarField& field);
template PtrList<surfaceVectorField>
DEIM<fvScalarMatrix>::generateSubFieldsMatrix(surfaceVectorField& field);
template PtrList<volScalarField> DEIM<fvVectorMatrix>::generateSubFieldsMatrix(
    volScalarField& field);
template PtrList<volVectorField> DEIM<fvVectorMatrix>::generateSubFieldsMatrix(
    volVectorField& field);
template PtrList<surfaceScalarField>
DEIM<fvVectorMatrix>::generateSubFieldsMatrix(surfaceScalarField& field);
template PtrList<surfaceVectorField>
DEIM<fvVectorMatrix>::generateSubFieldsMatrix(surfaceVectorField& field);

// Specialization for generateSubFieldsVector
template PtrList<volScalarField> DEIM<fvScalarMatrix>::generateSubFieldsVector(
    volScalarField& field);
template PtrList<volVectorField> DEIM<fvScalarMatrix>::generateSubFieldsVector(
    volVectorField& field);
template PtrList<surfaceScalarField>
DEIM<fvScalarMatrix>::generateSubFieldsVector(surfaceScalarField& field);
template PtrList<surfaceVectorField>
DEIM<fvScalarMatrix>::generateSubFieldsVector(surfaceVectorField& field);
template PtrList<volScalarField> DEIM<fvVectorMatrix>::generateSubFieldsVector(
    volScalarField& field);
template PtrList<volVectorField> DEIM<fvVectorMatrix>::generateSubFieldsVector(
    volVectorField& field);
template PtrList<surfaceScalarField>
DEIM<fvVectorMatrix>::generateSubFieldsVector(surfaceScalarField& field);
template PtrList<surfaceVectorField>
DEIM<fvVectorMatrix>::generateSubFieldsVector(surfaceVectorField& field);

// Specialization for generateSubmeshes
template volScalarField DEIM<volScalarField>::generateSubmeshes(
    label layers, fvMesh& mesh, volScalarField field,
    label secondTime);
template volVectorField DEIM<volVectorField>::generateSubmeshes(
    label layers, fvMesh& mesh, volVectorField field,
    label secondTime);
template volScalarField DEIM<volVectorField>::generateSubmeshes(
    label layers, fvMesh& mesh, volScalarField field,
    label secondTime);
template volVectorField DEIM<volScalarField>::generateSubmeshes(
    label layers, fvMesh& mesh, volVectorField field,
    label secondTime);

// Specialization for generateSubmeshesVector
template PtrList<volScalarField> DEIM<fvScalarMatrix>::generateSubmeshesVector(
    label layers, fvMesh& mesh, volScalarField field, label secondTime);
template PtrList<volVectorField> DEIM<fvScalarMatrix>::generateSubmeshesVector(
    label layers, fvMesh& mesh, volVectorField field, label secondTime);
template PtrList<surfaceScalarField>
DEIM<fvScalarMatrix>::generateSubmeshesVector(label layers, fvMesh& mesh,
        surfaceScalarField field, label secondTime);
template PtrList<surfaceVectorField>
DEIM<fvScalarMatrix>::generateSubmeshesVector(label layers, fvMesh& mesh,
        surfaceVectorField field, label secondTime);
template PtrList<volScalarField> DEIM<fvVectorMatrix>::generateSubmeshesVector(
    label layers, fvMesh& mesh, volScalarField field, label secondTime);
template PtrList<volVectorField> DEIM<fvVectorMatrix>::generateSubmeshesVector(
    label layers, fvMesh& mesh, volVectorField field, label secondTime);
template PtrList<surfaceScalarField>
DEIM<fvVectorMatrix>::generateSubmeshesVector(label layers, fvMesh& mesh,
        surfaceScalarField field, label secondTime);
template PtrList<surfaceVectorField>
DEIM<fvVectorMatrix>::generateSubmeshesVector(label layers, fvMesh& mesh,
        surfaceVectorField field, label secondTime);
template PtrList<volScalarField> DEIM<fvScalarMatrix>::generateSubmeshesMatrix(
    label layers, fvMesh& mesh, volScalarField field, label secondTime);
template PtrList<volVectorField> DEIM<fvScalarMatrix>::generateSubmeshesMatrix(
    label layers, fvMesh& mesh, volVectorField field, label secondTime);
template PtrList<surfaceScalarField>

// Specialization for generateSubmeshesMatrix
DEIM<fvScalarMatrix>::generateSubmeshesMatrix(label layers, fvMesh& mesh,
        surfaceScalarField field, label secondTime);
template PtrList<surfaceVectorField>
DEIM<fvScalarMatrix>::generateSubmeshesMatrix(label layers, fvMesh& mesh,
        surfaceVectorField field, label secondTime);
template PtrList<volScalarField> DEIM<fvVectorMatrix>::generateSubmeshesMatrix(
    label layers, fvMesh& mesh, volScalarField field, label secondTime);
template PtrList<volVectorField> DEIM<fvVectorMatrix>::generateSubmeshesMatrix(
    label layers, fvMesh& mesh, volVectorField field, label secondTime);
template PtrList<surfaceScalarField>
DEIM<fvVectorMatrix>::generateSubmeshesMatrix(label layers, fvMesh& mesh,
        surfaceScalarField field, label secondTime);
template PtrList<surfaceVectorField>
DEIM<fvVectorMatrix>::generateSubmeshesMatrix(label layers, fvMesh& mesh,
        surfaceVectorField field, label secondTime);

