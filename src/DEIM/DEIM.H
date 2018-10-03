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
    DEIM

Description
    Implementation of the discrete empirical interpolation method

SourceFiles
    DEIM.C

\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "ITHACAPOD.H"
#include "Foam2Eigen.H"
#include "EigenFunctions.H"
#include "ITHACAutilities.H"
#include "fvMeshSubset.H"


template<typename T, typename S>
class DEIM
{
public:


    ///
    /// @brief      Construct DEIM for non-linear function
    ///
    /// @param[in]  SnapShotsMatrix  The snapshots matrix
    /// @param[in]  MaxModes         The maximum number of modes
    /// @param[in]  FunctionName     The function name
    ///
    DEIM (T& SnapShotsMatrix, int MaxModes, word FunctionName);

    ///
    /// @brief      Construct DEIM for matrix with non-linear dependency
    ///
    /// @param      SnapShotsMatrix  The snapshots matrix
    /// @param[in]  MaxModesA        The maximum number of modes for the Matrix A
    /// @param[in]  MaxModesB        The maximum number of modes for the source term b
    /// @param[in]  MatrixName       The matrix name
    ///
    DEIM (T& SnapShotsMatrix, int MaxModesA, int MaxModesB, word MatrixName);

    /// The snapshots matrix containing the nonlinear function or operator
    T SnapShotsMatrix;

    /// The POD modes of the DEIM procedure that can be 
    T modes;

    /// The name of the matrix
    word MatrixName;

    /// The name of the non-linear function
    word FunctionName;

    /// Matrix Modes
    std::tuple<List<Eigen::SparseMatrix<double> >, List<Eigen::VectorXd> > Matrix_Modes;

    /// Field used to evaluate the online non linear function
    PtrList<S> fields;
    PtrList<S> fieldsA;
    PtrList<S> fieldsB;

    /// The maximum number of modes to be considered
    int MaxModes;
    int MaxModesA;
    int MaxModesB;

    /// A given tolerance
    double tolerance;

    /// Online vector
    Eigen::VectorXd theta;
    Eigen::VectorXd thetaA;
    Eigen::VectorXd thetaB;

    /// Int size M;
    int sizeM;

    /// source
    bool source;

    /// The matrix containing the modes
    Eigen::MatrixXd MatrixModes;

    /// Online Matrix
    Eigen::MatrixXd MatrixOnline;
    List<Eigen::SparseMatrix<double> > MatrixOnlineA;
    Eigen::MatrixXd MatrixOnlineB;

    /// The U matrix of the DEIM method
    Eigen::MatrixXd U;
    List<Eigen::SparseMatrix<double> > UA;
    Eigen::MatrixXd UB;

    /// Indices of the magic points
    List<int> magicPoints;
    List<Pair <int > > magicPointsA;
    List<int> magicPointsB;
    List<Pair<int> > xyz_A;
    List<int> xyz_B;

    /// Indices of the local magic points in the subMesh
    List<int> localMagicPoints;
    List<Pair <int> > localMagicPointsA;
    List<int> localMagicPointsB;

    /// The P matrix of the DEIM method
    Eigen::SparseMatrix<double> P;
    List<Eigen::SparseMatrix<double> > PA;
    Eigen::SparseMatrix<double> PB;

    /// List of submeshes
    PtrList<fvMeshSubset> submeshList;
    PtrList<fvMeshSubset> submeshListA;
    PtrList<fvMeshSubset> submeshListB;

    /// @brief      Function to generate the submeshes
    ///
    /// @param[in]  layers  Number of layers used to generate each submesh
    /// @param      mesh    The mesh of the problem
    ///
    void generateSubmeshes(int layers, fvMesh& mesh, S field, int secondTime = 0);

    ///
    /// @brief      { function_description }
    ///
    /// @param[in]  layers  The layers
    /// @param      mesh    The mesh
    /// @param[in]  field   The field
    ///
    void generateSubmeshesMatrix(int layers, fvMesh& mesh, S field, int secondTime = 0);

    ///
    /// @brief      { function_description }
    ///
    /// @param[in]  layers  The layers
    /// @param      mesh    The mesh
    /// @param[in]  field   The field
    ///
    void generateSubmeshesVector(int layers, fvMesh& mesh, S field, int secondTime = 0);

    ///
    /// @brief      Function to get the onlineCoeffs of the DEIM approx. It is problem dependent so it must be overridden.
    ///
    void onlineCoeffs();

    ///
    /// @brief      { function_description }
    ///
    /// @param      points       The points
    /// @param      submeshList  The submesh list
    ///
    /// @return     { description_of_the_return_value }
    ///
    List<int> global2local(List<int>& points, PtrList<fvMeshSubset>& submeshList);

    ///
    /// @brief      { function_description }
    ///
    /// @param      points       The points
    /// @param      submeshList  The submesh list
    ///
    /// @return     { description_of_the_return_value }
    ///
    List<Pair <int > > global2local(List< Pair <int> >& points, PtrList<fvMeshSubset>& submeshList);

    ///
    /// @brief      { function_description }
    ///
    /// @param      ind_rowA  The ind row a
    /// @param      ind_colA  The ind col a
    /// @param      xyz_rowA  The xyz row a
    /// @param      xyz_colA  The xyz col a
    ///
    void check3D_indices(int& ind_rowA, int&  ind_colA, int& xyz_rowA, int& xyz_colA);

    ///
    /// @brief      { function_description }
    ///
    /// @param      ind_rowA  The ind row a
    /// @param      xyz_rowA  The xyz row a
    ///
    void check3D_indices(int& ind_rowA, int& xyz_rowA);

};

// Template function constructor
template<typename T, typename S>
DEIM<T, S>::DEIM (T& s, int MaxModes, word FunctionName)
    :
    SnapShotsMatrix(s),
    MaxModes(MaxModes),
    FunctionName(FunctionName)
{
    Eigen::MatrixXd A;
    Eigen::VectorXd b;
    Eigen::VectorXd c;
    Eigen::VectorXd r;
    Eigen::VectorXd rho(1);
    modes = ITHACAPOD::DEIMmodes(SnapShotsMatrix, MaxModes, FunctionName);
    MatrixModes = Foam2Eigen::PtrList2Eigen(modes);
    int ind_max, c1;
    double max = MatrixModes.cwiseAbs().col(0).maxCoeff(&ind_max, &c1);
    rho(0) = max;
    magicPoints.append(ind_max);
    U = MatrixModes.col(0);
    P.resize(MatrixModes.rows(), 1);
    P.insert(ind_max, 0) = 1;
    for (int i = 1; i < MaxModes; i++)
    {
        A = P.transpose() * U;
        b = P.transpose() * MatrixModes.col(i);
        c = A.lu().solve(b);
        r = MatrixModes.col(i) - U * c;
        max = r.cwiseAbs().maxCoeff(&ind_max, &c1);
        P.conservativeResize(MatrixModes.rows(), i + 1);
        P.insert(ind_max, i) = 1;
        U.conservativeResize(MatrixModes.rows(), i + 1);
        U.col(i) =  MatrixModes.col(i);
        rho.conservativeResize(i + 1);
        rho(i) = max;
        magicPoints.append(ind_max);
    }
    MatrixOnline = U * ((P.transpose() * U).inverse());
}

template<typename T, typename S>
DEIM< T, S>::DEIM (T& s, int MaxModesA, int MaxModesB, word MatrixName)
    :
    SnapShotsMatrix(s),
    MaxModesA(MaxModesA),
    MaxModesB(MaxModesB),
    source(source),
    MatrixName(MatrixName)

{
    Eigen::MatrixXd AA;
    Eigen::VectorXd bA;
    Eigen::MatrixXd cA;
    Eigen::SparseMatrix<double> rA;
    Eigen::VectorXd rhoA(1);
    Matrix_Modes = ITHACAPOD::DEIMmodes(SnapShotsMatrix, MaxModesA, MaxModesB, MatrixName);
    sizeM = std::get<1>(Matrix_Modes)[0].rows();
    int ind_rowA, ind_colA, xyz_rowA, xyz_colA;
    double maxA = EigenFunctions::max(std::get<0>(Matrix_Modes)[0], ind_rowA, ind_colA);
    int ind_rowAOF = ind_rowA;
    int ind_colAOF = ind_colA;
    check3D_indices(ind_rowAOF, ind_colAOF, xyz_rowA, xyz_colA);
    Pair <int> indA(ind_rowAOF, ind_colAOF);
    Pair <int> xyzA(xyz_rowA, xyz_colA);
    xyz_A.append(xyzA);
    rhoA(0) = maxA;
    magicPointsA.append(indA);
    UA.append(std::get<0>(Matrix_Modes)[0]);
    Eigen::SparseMatrix<double> Pnow(std::get<0>(Matrix_Modes)[0].rows(), std::get<0>(Matrix_Modes)[0].cols());
    Pnow.insert(ind_rowA, ind_colA) = 1;
    PA.append(Pnow);
    for (int i = 1; i < MaxModesA; i++)
    {
        AA = EigenFunctions::innerProduct(PA, UA);
        bA = EigenFunctions::innerProduct(PA, std::get<0>(Matrix_Modes)[i]);
        cA = AA.colPivHouseholderQr().solve(bA);
        rA = std::get<0>(Matrix_Modes)[i] - EigenFunctions::MVproduct(UA, cA);
        double maxA = EigenFunctions::max(rA, ind_rowA, ind_colA);
        rhoA.conservativeResize(i + 1);
        rhoA(i) = maxA;
        int ind_rowAOF = ind_rowA;
        int ind_colAOF = ind_colA;
        check3D_indices(ind_rowAOF, ind_colAOF, xyz_rowA, xyz_colA);
        Pair <int> indA(ind_rowAOF, ind_colAOF);
        Pair <int> xyzA(xyz_rowA, xyz_colA);
        xyz_A.append(xyzA);
        magicPointsA.append(indA);
        UA.append(std::get<0>(Matrix_Modes)[i]);
        Eigen::SparseMatrix<double> Pnow(std::get<0>(Matrix_Modes)[0].rows(), std::get<0>(Matrix_Modes)[0].cols());
        Pnow.insert(ind_rowA, ind_colA) = 1;
        PA.append(Pnow);
    }
    Eigen::MatrixXd Aaux = EigenFunctions::innerProduct(PA, UA).inverse();
    MatrixOnlineA = EigenFunctions::MMproduct(UA , Aaux);

    Eigen::MatrixXd AB;
    Eigen::VectorXd bB;
    Eigen::VectorXd cB;
    Eigen::VectorXd rB;
    Eigen::VectorXd rhoB(1);
    int ind_rowB, xyz_rowB, c1;
    double maxB = std::get<1>(Matrix_Modes)[0].maxCoeff(&ind_rowB, &c1);

    int ind_rowBOF = ind_rowB;
    check3D_indices(ind_rowBOF, xyz_rowB);
    rhoB(0) = maxB;
    xyz_B.append(xyz_rowB);
    magicPointsB.append(ind_rowBOF);
    UB = std::get<1>(Matrix_Modes)[0];
    PB.resize(UB.rows(), 1);
    PB.insert(ind_rowB, 0) = 1;
    for (int i = 1; i < MaxModesB; i++)
    {
        AB = PB.transpose() * UB;
        bB = PB.transpose() * std::get<1>(Matrix_Modes)[i];
        cB = AB.colPivHouseholderQr().solve(bB);
        rB = std::get<1>(Matrix_Modes)[i] - UB * cB;
        maxB = rB.cwiseAbs().maxCoeff(&ind_rowB, &c1);
        ind_rowBOF = ind_rowB;
        check3D_indices(ind_rowBOF, xyz_rowB);
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
        MatrixOnlineB = UB * ((PB.transpose() * UB).inverse());
    }
    else
    {
        MatrixOnlineB = UB;
    }
}

template<typename T, typename S>
void DEIM<T, S>::generateSubmeshes(int layers, fvMesh& mesh, S field, int secondTime)
{

    fvMeshSubset *submesh;
    List<int> indices;

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
        dimensionSet(0, 0, 0, 0, 0)
    );
    if (!secondTime)
    {
        Indici = Indici * 0;
    }
    for (int i = 0; i < magicPoints.size(); i++)
    {
        submesh = new fvMeshSubset(mesh);
        indices = ITHACAutilities::getIndices(mesh, magicPoints[i], layers);
        if (!secondTime)
        {
            ITHACAutilities::assignONE(Indici, indices);
        }
        std::cout.setstate(std::ios_base::failbit);
        submesh->setLargeCellSubset(indices);
        submesh->subMesh().fvSchemes::readOpt() = mesh.fvSchemes::readOpt();
        submesh->subMesh().fvSolution::readOpt() = mesh.fvSolution::readOpt();
        submesh->subMesh().fvSchemes::read();
        submesh->subMesh().fvSolution::read();
        std::cout.clear();
        S f = submesh->interpolate(field);
        fields.append(f);
        if (!secondTime)
        {
            submeshList.append(submesh);
        }
    }
    if (!secondTime)
    {
        localMagicPoints = global2local(magicPoints, submeshList);
        ITHACAutilities::exportSolution(Indici, "./ITHACAoutput/DEIM/" + FunctionName, "1");
    }
}

template<typename T, typename S>
void DEIM<T, S>::generateSubmeshesMatrix(int layers, fvMesh& mesh, S field, int secondTime)
{
    fvMeshSubset *submesh;
    List<int> indices;
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
        dimensionSet(0, 0, 0, 0, 0)
    );
    if (!secondTime)
    {
        Indici = Indici * 0;
    }
    fieldsA.resize(0);
    for (int i = 0; i < magicPointsA.size(); i++)
    {
        submesh = new fvMeshSubset(mesh);
        indices = ITHACAutilities::getIndices(mesh, magicPointsA[i].first(), layers);
        indices.append(ITHACAutilities::getIndices(mesh, magicPointsA[i].second(), layers));
        if (!secondTime)
        {
            ITHACAutilities::assignONE(Indici, indices);
        }        
        std::cout.setstate(std::ios_base::failbit);
        submesh->setLargeCellSubset(indices);
        submesh->subMesh().fvSchemes::readOpt() = mesh.fvSchemes::readOpt();
        submesh->subMesh().fvSolution::readOpt() = mesh.fvSolution::readOpt();
        submesh->subMesh().fvSchemes::read();
        submesh->subMesh().fvSolution::read();
        std::cout.clear();
        S f = submesh->interpolate(field);
        fieldsA.append(f);
        if (!secondTime)
        {
            submeshListA.append(submesh);
        }
    }
    if (!secondTime)
    {
        localMagicPointsA = global2local(magicPointsA, submeshListA);
        ITHACAutilities::exportSolution(Indici, "./ITHACAoutput/DEIM/" + MatrixName , "1");
    }
}

template<typename T, typename S>
void DEIM<T, S>::generateSubmeshesVector(int layers, fvMesh& mesh, S field, int secondTime)
{
    fvMeshSubset *submesh;
    List<int> indices;
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
        dimensionSet(0, 0, 0, 0, 0)
    );
    if (!secondTime)
    {
        Indici = Indici * 0;
    }
    fieldsB.resize(0);

    for (int i = 0; i < magicPointsB.size(); i++)
    {
        submesh = new fvMeshSubset(mesh);
        indices = ITHACAutilities::getIndices(mesh, magicPointsB[i], layers);
        if (!secondTime)
        {

            ITHACAutilities::assignONE(Indici, indices);
        }
        submesh->setLargeCellSubset(indices);
        submesh->subMesh().fvSchemes::readOpt() = mesh.fvSchemes::readOpt();
        submesh->subMesh().fvSolution::readOpt() = mesh.fvSolution::readOpt();
        submesh->subMesh().fvSchemes::read();
        submesh->subMesh().fvSolution::read();
        std::cout.clear();
        S f = submesh->interpolate(field);
        fieldsB.append(f);
        if (!secondTime)
        {
            submeshListB.append(submesh);
        }
    }
    if (!secondTime)
    {
        localMagicPointsB = global2local(magicPointsB, submeshListB);
        ITHACAutilities::exportSolution(Indici, "./ITHACAoutput/DEIM/" + MatrixName , "1");
    }
}

template<typename T, typename S>
List<int> DEIM<T, S>::global2local(List<int>& points, PtrList<fvMeshSubset>& submeshList)
{
    List<int> localPoints;
    for (int i = 0; i < points.size(); i++)
    {
        for (int j = 0; j < submeshList[i].cellMap().size(); j++)
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

template<typename T, typename S>
List<Pair <int > > DEIM<T, S>::global2local(List<Pair <int > >& points, PtrList<fvMeshSubset>& submeshList)
{
    List< Pair <int> > localPoints(points.size());
    for (int i = 0; i < points.size(); i++)
    {
        for (int j = 0; j < submeshList[i].cellMap().size(); j++)
        {
            if (submeshList[i].cellMap()[j] == points[i].first())
            {
                localPoints[i].first() = j;
                break;
            }
        }
        for (int j = 0; j < submeshList[i].cellMap().size(); j++)
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

template<typename T, typename S>
void DEIM<T, S>::check3D_indices(int& ind_rowA, int&  ind_colA, int& xyz_rowA, int& xyz_colA)
{
    if (ind_rowA < sizeM)
    {
        xyz_rowA = 0;
    }
    else if (ind_rowA < sizeM * 2)
    {
        xyz_rowA = 1;
        ind_rowA = ind_rowA - sizeM;
    }
    else
    {
        xyz_rowA = 2;
        ind_rowA = ind_rowA - 2 * sizeM;
    }
    if (ind_colA < sizeM )
    {
        xyz_colA = 0;
    }
    else if (ind_colA < sizeM * 2)
    {
        xyz_colA = 1;
        ind_colA = ind_colA - 2 * sizeM;
    }
    else
    {
        xyz_colA = 2;
        ind_colA = ind_colA - 2 * sizeM;
    }
};

template<typename T, typename S>
void DEIM<T, S>::check3D_indices(int& ind_rowA, int& xyz_rowA)
{
    if (ind_rowA < sizeM)
    {
        xyz_rowA = 0;
    }
    else if (ind_rowA < sizeM * 2)
    {
        xyz_rowA = 1;
        ind_rowA = ind_rowA - sizeM;
    }
    else
    {
        xyz_rowA = 2;
        ind_rowA = ind_rowA - 2 * sizeM;
    }
};