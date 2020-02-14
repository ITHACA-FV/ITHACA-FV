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

#include "ITHACAutilities.H"
#include "ITHACAstream.H"

/// \file
/// Source file of the ITHACAutilities class.

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

List<int> ITHACAutilities::getIndices(fvMesh& mesh, int index, int layers)
{
    List<int> out;
    out.resize(1);
    out[0] = index;

    for (int i = 0; i < layers; i++)
    {
        int size = out.size();

        for (int j = 0; j < size; j++)
        {
            out.append(mesh.cellCells()[out[j]]);
        }
    }

    labelList uniqueIndex;
    uniqueOrder(out, uniqueIndex);
    List<int> out2;
    forAll(uniqueIndex, i)
    {
        out2.append(out[uniqueIndex[i]]);
    }
    return out2;
}

List<int> ITHACAutilities::getIndices(fvMesh& mesh, int index_row,
                                      int index_col, int layers)
{
    List<int> out;
    out.resize(2);
    out[0] = index_row;
    out[1] = index_col;

    for (int i = 0; i < layers; i++)
    {
        int size = out.size();

        for (int j = 0; j < size; j++)
        {
            out.append(mesh.cellCells()[out[j]]);
        }
    }

    labelList uniqueIndex;
    uniqueOrder(out, uniqueIndex);
    List<int> out2;
    forAll(uniqueIndex, i)
    {
        out2.append(out[uniqueIndex[i]]);
    }
    return out2;
}


void ITHACAutilities::createSymLink(word folder)
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

Eigen::MatrixXd ITHACAutilities::rand(int rows, int cols, double min,
                                      double max)
{
    std::srand(static_cast<long unsigned int>
               (std::chrono::high_resolution_clock::now().time_since_epoch().count()));
    Eigen::MatrixXd matr = Eigen::MatrixXd::Random(rows, cols);
    matr = (matr.array() + 1) / 2;
    matr = matr.array() * (max - min);
    matr = matr.array() + min;
    return matr;
}

Eigen::MatrixXd ITHACAutilities::rand(int rows, Eigen::MatrixXd minMax)
{
    std::srand(static_cast<long unsigned int>
               (std::chrono::high_resolution_clock::now().time_since_epoch().count()));
    int cols = minMax.rows();
    Eigen::MatrixXd matr = Eigen::MatrixXd::Random(rows, cols);
    matr = (matr.array() + 1) / 2;

    for (int i = 0; i < cols; i++)
    {
        matr.col(i) = matr.col(i).array() * (minMax(i, 1) - minMax(i, 0));
        matr.col(i) = matr.col(i).array() + (minMax(i, 0));
    }

    return matr;
}


bool ITHACAutilities::isInteger(double ratio)
{
    bool checkResult = 0;

    if (abs(round(ratio) - ratio) < std::sqrt(SMALL))
    {
        checkResult = true;
    }
    else
    {
        checkResult = false;
    }

    return checkResult;
}
// Check if the modes exists
bool ITHACAutilities::check_pod()
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

    reduce(pod_exist, sumOp<int>());

    if (!pod_exist)
    {
        createSymLink("./ITHACAoutput/POD");
    }

    return pod_exist;
}

// Check if the offline data exist
bool ITHACAutilities::check_off()
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

    reduce(off_exist, sumOp<int>());

    if (!off_exist)
    {
        createSymLink("./ITHACAoutput/Offline");
    }

    return off_exist;
}

// Check if the supremizer data exist
bool ITHACAutilities::check_sup()
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

    reduce(sup_exist, sumOp<int>());

    if (!sup_exist)
    {
        createSymLink("./ITHACAoutput/supremizer");
    }

    return sup_exist;
}

bool ITHACAutilities::check_folder(word folder)
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

bool ITHACAutilities::check_file(std::string fileName)
{
    std::ifstream infile(fileName);
    return infile.good();
}

template<class TypeField>
PtrList<TypeField> ITHACAutilities::reconstruct_from_coeff(
    PtrList<TypeField>& modes, Eigen::MatrixXd& coeff_matrix, label Nmodes)
{
    PtrList<TypeField> rec_field;
    rec_field.resize(0);

    for (label k = 0; k < coeff_matrix.cols(); k++)
    {
        for (label i = 0; i < Nmodes; i++)
        {
            if ( i == 0)
            {
                rec_field.append(modes[i]*coeff_matrix(i, k));
            }
            else
            {
                rec_field[k] +=  modes[i] * coeff_matrix(i, k);
            }
        }
    }

    return rec_field;
}

template<class TypeField>
double ITHACAutilities::error_fields(TypeField& field1,
                                     TypeField& field2)
{
    double err;

    if (L2norm(field1) <= 1e-6)
    {
        err = 0;
    }
    else
    {
        err = L2norm(field1 - field2) / L2norm(field1);
    }

    return err;
}

template<>
double ITHACAutilities::error_fields(
    GeometricField<vector, fvPatchField, volMesh>& field1,
    GeometricField<vector, fvPatchField, volMesh>& field2, volScalarField& Volumes)

{
    volScalarField diffFields2 = ((field1 - field2) & (field1 - field2)) * Volumes;
    double err = Foam::sqrt(gSum(diffFields2));
    return err;
}

template<>
double ITHACAutilities::error_fields(
    GeometricField<scalar, fvPatchField, volMesh>& field1,
    GeometricField<scalar, fvPatchField, volMesh>& field2, volScalarField& Volumes)

{
    volScalarField diffFields2 = ((field1 - field2) * (field1 - field2)) * Volumes;
    double err = Foam::sqrt(gSum(diffFields2));
    return err;
}

template<class TypeField>
double ITHACAutilities::error_fields_abs(TypeField& field1,
        TypeField& field2)
{
    double err = L2norm(field1 - field2);
    return err;
}

template<class TypeField>
Eigen::MatrixXd ITHACAutilities::error_listfields(PtrList<TypeField>&
        fields1, PtrList<TypeField>& fields2)
{
    Eigen::VectorXd err;

    if (fields1.size() != fields2.size())
    {
        Info << "The two fields do not have the same size, code will abort" << endl;
        exit(0);
    }

    err.resize(fields1.size(), 1);

    for (label k = 0; k < fields1.size(); k++)
    {
        err(k, 0) = error_fields(fields1[k], fields2[k]);
        Info << " Error is " << err[k] << endl;
    }

    return err;
}

template<class TypeField>
Eigen::MatrixXd ITHACAutilities::error_listfields(
    PtrList<GeometricField<TypeField, fvPatchField, volMesh>>& fields1,
    PtrList<GeometricField<TypeField, fvPatchField, volMesh>>& fields2,
    PtrList<volScalarField>& Volumes)
{
    M_Assert(fields1.size() == fields2.size(),
             "The two fields do not have the same size, code will abort");
    M_Assert(fields1.size() == Volumes.size(),
             "The volumes field and the two solution fields do not have the same size, code will abort");
    Eigen::VectorXd err;
    err.resize(fields1.size(), 1);

    for (label k = 0; k < fields1.size(); k++)
    {
        err(k, 0) = error_fields(fields1[k], fields2[k], Volumes[k]);
        Info << " Error is " << err[k] << endl;
    }

    return err;
}

template<class TypeField>
Eigen::MatrixXd ITHACAutilities::error_listfields_abs(PtrList<TypeField>&
        fields1, PtrList<TypeField>& fields2)
{
    Eigen::VectorXd err;

    if (fields1.size() != fields2.size())
    {
        Info << "The two fields do not have the same size, code will abort" << endl;
        exit(0);
    }

    err.resize(fields1.size(), 1);

    for (label k = 0; k < fields1.size(); k++)
    {
        err(k, 0) = error_fields_abs(fields1[k], fields2[k]);
        Info << " Error is " << err[k] << endl;
    }

    return err;
}

Eigen::MatrixXd ITHACAutilities::get_mass_matrix(PtrList<volVectorField> modes,
        int Nmodes)
{
    label Msize;

    if (Nmodes == 0)
    {
        Msize =  modes.size();
    }
    else
    {
        Msize = Nmodes;
    }

    M_Assert(modes.size() >= Msize,
             "The Number of requested modes is larger then the available quantity.");
    Eigen::MatrixXd M_matrix(Msize, Msize);

    // Project everything
    for (label i = 0; i < Msize; i++)
    {
        for (label j = 0; j < Msize; j++)
        {
            M_matrix(i, j) = fvc::domainIntegrate(modes[i] & modes[j]).value();
        }
    }

    return M_matrix;
}

Eigen::MatrixXd ITHACAutilities::get_mass_matrix(PtrList<volScalarField> modes,
        int Nmodes)
{
    label Msize;

    if (Nmodes == 0)
    {
        Msize =  modes.size();
    }
    else
    {
        Msize = Nmodes;
    }

    M_Assert(modes.size() >= Msize,
             "The Number of requested modes is larger then the available quantity.");
    Eigen::MatrixXd M_matrix(Msize, Msize);

    // Project everything
    for (label i = 0; i < Msize; i++)
    {
        for (label j = 0; j < Msize; j++)
        {
            M_matrix(i, j) = fvc::domainIntegrate(modes[i] * modes[j]).value();
        }
    }

    return M_matrix;
}

template<class TypeField>
Eigen::VectorXd ITHACAutilities::get_mass_matrix_FV(
    GeometricField<TypeField, fvPatchField, volMesh>& snapshot)
{
    Eigen::MatrixXd snapEigen = Foam2Eigen::field2Eigen(snapshot);
    int dim = std::nearbyint(snapEigen.rows() / (snapshot.mesh().V()).size());
    Eigen::VectorXd volumes = Foam2Eigen::field2Eigen(snapshot.mesh().V());
    Eigen::VectorXd vol3 = volumes.replicate(dim, 1);
    return vol3;
}

Eigen::VectorXd ITHACAutilities::get_coeffs(volVectorField snapshot,
        PtrList<volVectorField>& modes, int Nmodes)
{
    label Msize;

    if (Nmodes == 0)
    {
        Msize =  modes.size();
    }
    else
    {
        Msize = Nmodes;
    }

    M_Assert(modes.size() >= Msize,
             "The Number of requested modes is larger then the available quantity.");
    Eigen::MatrixXd M_matrix = get_mass_matrix(modes, Nmodes);
    Eigen::VectorXd a(Msize);
    Eigen::VectorXd b(Msize);

    // Project everything
    for (label i = 0; i < Msize; i++)
    {
        b(i) = fvc::domainIntegrate(snapshot & modes[i]).value();
    }

    a = M_matrix.colPivHouseholderQr().solve(b);
    return a;
}

Eigen::VectorXd ITHACAutilities::get_coeffs(volScalarField snapshot,
        PtrList<volScalarField>& modes, int Nmodes)
{
    label Msize;

    if (Nmodes == 0)
    {
        Msize =  modes.size();
    }
    else
    {
        Msize = Nmodes;
    }

    M_Assert(modes.size() >= Msize,
             "The Number of requested modes is larger then the available quantity.");
    Eigen::MatrixXd M_matrix = get_mass_matrix(modes, Msize);
    Eigen::VectorXd a(Msize);
    Eigen::VectorXd b(Msize);

    // Project everything
    for (label i = 0; i < Msize; i++)
    {
        b(i) = fvc::domainIntegrate(snapshot * modes[i]).value();
    }

    a = M_matrix.colPivHouseholderQr().solve(b);
    return a;
}

template<class TypeField>
Eigen::MatrixXd ITHACAutilities::get_coeffs(PtrList<TypeField> snapshots,
        PtrList<TypeField> modes, int Nmodes)
{
    label Msize;

    if (Nmodes == 0)
    {
        Msize =  modes.size();
    }
    else
    {
        Msize = Nmodes;
    }

    Eigen::MatrixXd coeff(Msize, snapshots.size());

    for (auto i = 0; i < snapshots.size(); i++)
    {
        coeff.col(i) = get_coeffs(snapshots[i], modes, Nmodes);
    }

    return coeff;
}

template<>
Eigen::VectorXd ITHACAutilities::get_coeffs_ortho(volScalarField snapshot,
        PtrList<volScalarField>& modes, int Nmodes)
{
    label Msize;

    if (Nmodes == 0)
    {
        Msize =  modes.size();
    }
    else
    {
        Msize = Nmodes;
    }

    M_Assert(modes.size() >= Msize,
             "The Number of requested modes is larger then the available quantity.");
    Eigen::VectorXd b(Msize);

    // Project everything
    for (label i = 0; i < Msize; i++)
    {
        b(i) = fvc::domainIntegrate(snapshot * modes[i]).value();
    }

    return b;
}

template<>
Eigen::VectorXd ITHACAutilities::get_coeffs_ortho(volVectorField
        snapshot, PtrList<volVectorField>& modes, int Nmodes)
{
    label Msize;

    if (Nmodes == 0)
    {
        Msize =  modes.size();
    }
    else
    {
        Msize = Nmodes;
    }

    M_Assert(modes.size() >= Msize,
             "The Number of requested modes is larger then the available quantity.");
    Eigen::VectorXd b(Msize);

    // Project everything
    for (label i = 0; i < Msize; i++)
    {
        b(i) = fvc::domainIntegrate(snapshot & modes[i]).value();
    }

    return b;
}

template<class TypeField>
Eigen::MatrixXd ITHACAutilities::get_coeffs_ortho(PtrList<TypeField>
        snapshots, PtrList<TypeField>& modes, int Nmodes)
{
    label Msize;

    if (Nmodes == 0)
    {
        Msize =  modes.size();
    }
    else
    {
        Msize = Nmodes;
    }

    Eigen::MatrixXd coeff(Msize, snapshots.size());

    for (auto i = 0; i < snapshots.size(); i++)
    {
        coeff.col(i) = get_coeffs_ortho(snapshots[i], modes, Nmodes);
    }

    return coeff;
}

template<class TypeField>
Eigen::MatrixXd ITHACAutilities::getCoeffsFrobenius(PtrList<TypeField>
        snapshots,
        PtrList<TypeField>& modes, int nModes)
{
    label Msize;

    if (nModes == 0)
    {
        Msize =  modes.size();
    }
    else
    {
        Msize = nModes;
    }

    Eigen::MatrixXd ModesE = (Foam2Eigen::PtrList2Eigen(modes)).leftCols(nModes);
    Eigen::MatrixXd SnapsE = Foam2Eigen::PtrList2Eigen(snapshots);
    Eigen::MatrixXd Mass = ModesE.transpose() * ModesE;
    Eigen::MatrixXd rhs = ModesE.transpose() * SnapsE;
    Eigen::MatrixXd coeffs;
    coeffs.resize(ModesE.cols(), SnapsE.cols());

    for (int j = 0; j < SnapsE.cols(); j++)
    {
        coeffs.col(j) = Mass.fullPivLu().solve(rhs.col(j));
    }

    return coeffs;
}

double ITHACAutilities::L2norm(volScalarField field)
{
    double a;
    a = Foam::sqrt(fvc::domainIntegrate(field * field).value());
    return a;
}

double ITHACAutilities::L2norm(volVectorField field)
{
    double a;
    a = Foam::sqrt(fvc::domainIntegrate(field & field).value());
    return a;
}

double ITHACAutilities::H1seminorm(volScalarField field)
{
    double a;
    a = Foam::sqrt(fvc::domainIntegrate(fvc::grad(field) & fvc::grad(
                                            field)).value());
    return a;
}

double ITHACAutilities::H1seminorm(volVectorField field)
{
    double a;
    a = Foam::sqrt(fvc::domainIntegrate(fvc::grad(field)
                                        && fvc::grad(field)).value());
    return a;
}

template<class TypeField>
double frobNorm(TypeField field)
{
    double norm(0);
    Eigen::VectorXd vF = Foam2Eigen::field2Eigen(field);
    norm = vF.norm();
    return norm;
}

void ITHACAutilities::setBoxToValue(volScalarField& field, Eigen::MatrixXd Box,
                                    double value)
{
    M_Assert(Box.rows() == 2
             && Box.cols() == 3,
             "The box must be a 2*3 matrix shaped in this way: \nBox = \t|x0, y0, z0|\n\t|x1, yi, z1|\n");

    for (label i = 0; i < field.internalField().size(); i++)
    {
        auto cx = field.mesh().C()[i].component(vector::X);
        auto cy = field.mesh().C()[i].component(vector::Y);
        auto cz = field.mesh().C()[i].component(vector::Z);

        if (cx >= Box(0, 0) && cy >= Box(0, 1) && cz >= Box(0, 2) && cx <= Box(1, 0)
                && cy <= Box(1, 1) && cz <= Box(1, 2) )
        {
            field.ref()[i] = value;
        }
    }

    for (label i = 0; i < field.boundaryField().size(); i++)
    {
        for (label j = 0; j < field.boundaryField()[i].size(); j++)
        {
            if (field.boundaryField()[i].type() == "fixedValue"
                    || field.boundaryField()[i].type() == "calculated")
            {
                auto cx = field.mesh().C().boundaryField()[i][j][0];
                auto cy = field.mesh().C().boundaryField()[i][j][1];
                auto cz = field.mesh().C().boundaryField()[i][j][2];

                if (cx >= Box(0, 0) && cy >= Box(0, 1) && cz >= Box(0, 2) && cx <= Box(1, 0)
                        && cy <= Box(1, 1) && cz <= Box(1, 2) )
                {
                    field.boundaryFieldRef()[i][j] = value;
                }
            }
        }
    }
}

void ITHACAutilities::assignONE(volScalarField& s, List<int>& L)
{
    for (label i = 0; i < L.size(); i++)
    {
        s.ref()[L[i]] = 1;
    }
}

// Assign a BC for a vector field
void ITHACAutilities::assignBC(volScalarField& s, label BC_ind, double& value)
{
    if (s.boundaryField()[BC_ind].type() == "fixedValue"
            || s.boundaryField()[BC_ind].type() == "calculated")
    {
        for (label i = 0; i < s.boundaryField()[BC_ind].size(); i++)
        {
            s.boundaryFieldRef()[BC_ind][i] = value;
        }
    }
    else if (s.boundaryField()[BC_ind].type() == "fixedGradient")
    {
        fixedGradientFvPatchScalarField& Tpatch =
            refCast<fixedGradientFvPatchScalarField>(s.boundaryFieldRef()[BC_ind]);
        scalarField& gradTpatch = Tpatch.gradient();
        forAll(gradTpatch, faceI)
        {
            gradTpatch[faceI] = value;
        }
    }
    else if (s.boundaryField()[BC_ind].type() == "fixedFluxPressure")
    {
        for (label i = 0; i < s.boundaryField()[BC_ind].size(); i++)
        {
            s.boundaryFieldRef()[BC_ind][i] = value;
        }
    }
    else if (s.boundaryField()[BC_ind].type() == "freestream")
    {
        for (label i = 0; i < s.boundaryField()[BC_ind].size(); i++)
        {
            s.boundaryFieldRef()[BC_ind][i] = value;
        }

        freestreamFvPatchField<scalar>& Tpatch =
            refCast<freestreamFvPatchField<scalar>>(s.boundaryFieldRef()[BC_ind]);
        scalarField& gradTpatch = Tpatch.freestreamValue();
        forAll(gradTpatch, faceI)
        {
            gradTpatch[faceI] = value;
        }
    }
    else if (s.boundaryField()[BC_ind].type() == "calculated")
    {
        for (label i = 0; i < s.boundaryField()[BC_ind].size(); i++)
        {
            s.boundaryFieldRef()[BC_ind][i] = value;
        }
    }
    else if (s.boundaryField()[BC_ind].type() == "processor")
    {
        for (label i = 0; i < s.boundaryField()[BC_ind].size(); i++)
        {
            s.boundaryFieldRef()[BC_ind][i] = value;
        }
    }
    else if (s.boundaryField()[BC_ind].type() == "empty")
    {
    }
}

// Assign a BC for a scalar field
void ITHACAutilities::assignBC(volScalarField& s, label BC_ind,
                               List<double> value)
{
    if (s.boundaryField()[BC_ind].type() == "fixedValue")
    {
        s.boundaryFieldRef()[BC_ind] = value;
    }
    else if (s.boundaryField()[BC_ind].type() == "fixedGradient")
    {
        fixedGradientFvPatchScalarField& Tpatch =
            refCast<fixedGradientFvPatchScalarField>(s.boundaryFieldRef()[BC_ind]);
        scalarField& gradTpatch = Tpatch.gradient();
        gradTpatch = value;
    }
    else if (s.boundaryField()[BC_ind].type() == "empty")
    {
    }
    else
    {
        Info << "This type of boundary condition is not yet implemented, code will abort"
             << endl;
        exit(0);
    }
}

// Assign a BC for a scalar field
void ITHACAutilities::assignBC(volVectorField& s, label BC_ind,
                               Vector<double>& value)
{
    if (s.boundaryField()[BC_ind].type() == "fixedValue"
            || s.boundaryField()[BC_ind].type() == "processor"
            || s.boundaryField()[BC_ind].type() == "calculated")
    {
        for (label i = 0; i < s.boundaryField()[BC_ind].size(); i++)
        {
            s.boundaryFieldRef()[BC_ind][i] = value;
        }
    }
    else if (s.boundaryField()[BC_ind].type() == "fixedGradient")
    {
        Info << "This Feature is not implemented for this boundary condition" << endl;
        exit(0);
    }
    else if (s.boundaryField()[BC_ind].type() == "freestream")
    {
        for (label i = 0; i < s.boundaryField()[BC_ind].size(); i++)
        {
            s.boundaryFieldRef()[BC_ind][i] = value;
        }

        freestreamFvPatchField<vector>& Tpatch =
            refCast<freestreamFvPatchField<vector>>(s.boundaryFieldRef()[BC_ind]);
        vectorField& gradTpatch = Tpatch.freestreamValue();
        forAll(gradTpatch, faceI)
        {
            gradTpatch[faceI] = value;
        }
    }
}

void ITHACAutilities::assignBC(volScalarField& s, label BC_ind,
                               Eigen::MatrixXd valueVec)
{
    word typeBC = s.boundaryField()[BC_ind].type();

    if (typeBC == "fixedValue" || typeBC == "calculated"
            || typeBC == "fixedFluxPressure" ||  typeBC == "processor")
    {
        for (label i = 0; i < s.boundaryField()[BC_ind].size(); i++)
        {
            double value = valueVec(i);
            s.boundaryFieldRef()[BC_ind][i] = value;
        }
    }
    else if (s.boundaryField()[BC_ind].type() == "fixedGradient")
    {
        fixedGradientFvPatchScalarField& Tpatch =
            refCast<fixedGradientFvPatchScalarField>(s.boundaryFieldRef()[BC_ind]);
        scalarField& gradTpatch = Tpatch.gradient();
        forAll(gradTpatch, faceI)
        {
            double value = valueVec(faceI);
            gradTpatch[faceI] = value;
        }
    }
    else if (s.boundaryField()[BC_ind].type() == "freestream")
    {
        for (label i = 0; i < s.boundaryField()[BC_ind].size(); i++)
        {
            double value = valueVec(i);
            s.boundaryFieldRef()[BC_ind][i] = value;
        }

        freestreamFvPatchField<scalar>& Tpatch =
            refCast<freestreamFvPatchField<scalar>>(s.boundaryFieldRef()[BC_ind]);
        scalarField& gradTpatch = Tpatch.freestreamValue();
        forAll(gradTpatch, faceI)
        {
            double value = valueVec(faceI);
            gradTpatch[faceI] = value;
        }
    }
    else if (s.boundaryField()[BC_ind].type() == "empty")
    {
    }
}

// Assign a BC for a scalar field
void ITHACAutilities::assignBC(volVectorField& s, label BC_ind,
                               Eigen::MatrixXd valueVec)
{
    word typeBC = s.boundaryField()[BC_ind].type();
    int sizeBC = s.boundaryField()[BC_ind].size();

    if (typeBC == "fixedValue" || typeBC == "calculated" || typeBC == "processor")
    {
        for (label i = 0; i < sizeBC; i++)
        {
            Vector<double> value(valueVec(i), valueVec(i + sizeBC),
                                 valueVec(i + sizeBC * 2));
            s.boundaryFieldRef()[BC_ind][i] = value;
        }
    }
    else if (s.boundaryField()[BC_ind].type() == "fixedGradient")
    {
        Info << "This Feature is not implemented for this boundary condition" << endl;
        exit(0);
    }
    else if (s.boundaryField()[BC_ind].type() == "freestream")
    {
        for (label i = 0; i < s.boundaryField()[BC_ind].size(); i++)
        {
            Vector<double> value(valueVec(i), valueVec(i + sizeBC),
                                 valueVec(i + sizeBC * 2));
            s.boundaryFieldRef()[BC_ind][i] = value;
        }

        freestreamFvPatchField<vector>& Tpatch =
            refCast<freestreamFvPatchField<vector>>(s.boundaryFieldRef()[BC_ind]);
        vectorField& gradTpatch = Tpatch.freestreamValue();
        forAll(gradTpatch, faceI)
        {
            Vector<double> value(valueVec(faceI), valueVec(faceI + sizeBC),
                                 valueVec(faceI + sizeBC * 2));
            gradTpatch[faceI] = value;
        }
    }
}
template<class TypeField>
PtrList<TypeField> ITHACAutilities::averageSubtract(PtrList<TypeField>
        fields, Eigen::MatrixXd ind, PtrList<TypeField>& ave)
{
    PtrList<TypeField> aveSubtracted;
    Eigen::VectorXd newInd;
    newInd.resize(ind.size() + 1);
    newInd.head(ind.size()) = ind;
    newInd(ind.size()) = fields.size();

    for (label i = 0; i < ind.size(); i++)
    {
        TypeField aveTemp("nut", fields[0] * 0);

        for (label j = newInd(i); j < newInd(i + 1); j++)
        {
            aveTemp += fields[j];
        }

        aveTemp /= newInd(i + 1) - newInd(i);
        ave.append(aveTemp);
    }

    for (label i = 0; i < ind.size(); i++)
    {
        for (label j = newInd(i); j < newInd(i + 1); j++)
        {
            TypeField newfield("nut", fields[0] * 0);
            newfield = fields[j] - ave[i];
            aveSubtracted.append(newfield);
        }
    }

    return aveSubtracted;
}

template<class TypeField>
TypeField ITHACAutilities::computeAverage(PtrList<TypeField>& fields)
{
    TypeField av(fields[0]);

    for (int i = 1; i < fields.size(); i++)
    {
        av += fields[i];
    }

    av = av / fields.size();
    return av;
}



Eigen::MatrixXd ITHACAutilities::parTimeCombMat(List<Eigen::VectorXd>
        acquiredSnapshotsTimes,
        Eigen::MatrixXd parameters)
{
    int parsNum = parameters.cols();
    int parsSamplesNum = parameters.rows();
    M_Assert(parsSamplesNum == acquiredSnapshotsTimes.size(),
             "The list of time instants does not have the same number of vectors as the number of parameters samples");
    Eigen::MatrixXd comb;
    int totalSnapshotsNum = 0;

    for (label k = 0; k < acquiredSnapshotsTimes.size(); k++)
    {
        totalSnapshotsNum += acquiredSnapshotsTimes[k].size();
    }

    comb.resize(totalSnapshotsNum, parsNum + 1);
    label i = 0;

    for (label j = 0; j < acquiredSnapshotsTimes.size(); j++)
    {
        for (label k = 0; k < acquiredSnapshotsTimes[j].size(); k++)
        {
            comb(i, parsNum) = (acquiredSnapshotsTimes[j])(k, 0);
            comb.block(i, 0, 1, parsNum) = parameters.row(j);
            i = i + 1;
        }
    }

    return comb;
}

template<class TypeField>
void ITHACAutilities::changeBCtype(
    GeometricField<TypeField, fvPatchField, volMesh>& field, word BCtype,
    label BC_ind)
{
    field.boundaryFieldRef().set(BC_ind, fvPatchField<TypeField>::New(BCtype,
                                 field.mesh().boundary()[BC_ind], field));
}

List<vector> ITHACAutilities::displacedSegment(List<vector> x0, double mux1,
        double mux2, double muy1, double muy2)
{
    vector minimum = min(x0);
    vector maximum = max(x0);
    double l = Foam::sqrt(pow(minimum[0] - maximum[0],
                              2) + pow(minimum[1] - maximum[1], 2));
    double limin;
    double limax;
    List<vector> xdef(x0);
    List<vector> displ(x0);

    for (int i = 0; i < x0.size(); i++)
    {
        limin = Foam::sqrt(pow(x0[i][0] - minimum[0], 2) + pow(x0[i][1] - minimum[1],
                           2));
        limax = Foam::sqrt(pow(x0[i][0] - maximum[0], 2) + pow(x0[i][1] - maximum[1],
                           2));
        xdef[i][0] = x0[i][0] + mux1 * (1 - limin / l) + mux2 * (1 - limax / l);
        xdef[i][1] = x0[i][1] + muy1 * (1 - limin / l) + muy2 * (1 - limax / l);
        xdef[i][2] = x0[i][2];
    }

    displ = xdef - x0;
    return displ;
}

template<>
void ITHACAutilities::setIndices2Value(labelList& ind2set,
                                       List<vector>& value2set, labelList& movingIDS, List<vector>& originalList)
{
    M_Assert(ind2set.size() == value2set.size(),
             "The size of the indices must be equal to the size of the values list");
    M_Assert(originalList.size() >= value2set.size(),
             "The size of the original list of values must be bigger than the size of the list of values you want to set");
    labelList ind_ok(ind2set);

    for (int i = 0; i < ind2set.size(); i++)
    {
        for (int k = 0; k < movingIDS.size(); k++)
        {
            if (ind2set[i] == movingIDS[k])
            {
                ind_ok[i] = k;
                break;
            }
        }
    }

    for (int i = 0; i < ind2set.size(); i++)
    {
        originalList[ind_ok[i]] = value2set[i];
    }
}

volScalarField ITHACAutilities::meshNonOrtho(fvMesh& mesh,
        volScalarField& NonOrtho)
{
    scalarField sno = (polyMeshTools::faceOrthogonality(mesh, mesh.Sf(), mesh.C()));

    for (int i = 0; i < sno.size(); i++)
    {
        sno[i] = Foam::acos(min(1, sno[i])) * 180 / constant::mathematical::pi;
    }

    surfaceScalarField pippo = mesh.magSf();
    const fvPatchList& patches = mesh.boundary();

    for (int i = 0; i < pippo.internalField().size(); i++)
    {
        pippo.ref()[i] = sno[i];
    }

    for (int i = 0; i < patches.size(); i++)
    {
        if ( patches[i].type() != "empty" )
        {
            label start = patches[i].patch().start();
            label n = patches[i].patch().size();

            for (int k = 0; k < n; k++)
            {
                pippo.boundaryFieldRef()[i][k] = sno[start + k];
            }
        }
    }

    NonOrtho = fvc::average(pippo);
    return NonOrtho;
}

void ITHACAutilities::getPointsFromPatch(fvMesh& mesh, label ind,
        List<vector>& points, labelList& indices)
{
    pointField meshPoints(mesh.points());
    const polyPatch& patchFound = mesh.boundaryMesh()[ind];
    labelList labelPatchFound(patchFound.meshPoints());
    points.resize(labelPatchFound.size());

    for (int i = 0; i < labelPatchFound.size(); i++)
    {
        points[i] = meshPoints[labelPatchFound[i]];
    }

    indices = labelPatchFound;
}

vector ITHACAutilities::displacePoint(vector x0, vector x_low, vector x_up,
                                      double mux_low, double mux_up, double muy_low, double muy_up, double muz_low,
                                      double muz_up)
{
    vector direction = x_up - x_low;
    double t0;
    double t1;
    double t2;
    vector x_low_def(x_low[0] + mux_low, x_low[1] + muy_low, x_low[2] + muz_low);
    vector x_up_def(x_up[0] + mux_up, x_up[1] + muy_up, x_up[2] + muz_up);
    vector direction_def = x_up_def - x_low_def;

    if (abs(direction[0]) > 1e-16)
    {
        t0 = (x0[0] - x_low[0]) / direction[0];
    }
    else
    {
        t0 = 0;
    }

    if (abs(direction[1]) > 1e-16)
    {
        t1 = (x0[1] - x_low[1]) / direction[1];

        if (t0 == 0)
        {
            t0 = t1;
        }
    }
    else
    {
        t1 = t0;
    }

    if (abs(direction[2]) > 1e-16)
    {
        t2 = (x0[2] - x_low[2]) / direction[2];

        if (t1 == 0)
        {
            t1 = t2;
        }

        if (t0 == 0)
        {
            t0 = t2;
        }
    }
    else
    {
        t2 = t1;
    }

    Info << abs((x_low + t0 * direction - x0)[0]) << endl;
    Info << abs((x_low + t1 * direction - x0)[1]) << endl;
    Info << abs((x_low + t2 * direction - x0)[2]) << endl;
    M_Assert(abs(abs((x_low + t0 * direction - x0)[0]) + abs((
                     x_low + t0 * direction - x0)[1]) + abs((x_low + t0 * direction - x0)[2])) <
             1e-6, "The givent point is not on the segment");
    vector def_point = x_low_def + t1 * direction_def;
    Info << def_point << endl;
    return def_point;
}

List<vector> ITHACAutilities::rotatePoints(const List<vector>& originalPoints,
        vector AxisOfRotation, double AngleOfRotation)
{
    double theta = AngleOfRotation / 180 *  constant::mathematical::pi;
    quaternion q(AxisOfRotation, theta);
    List<vector> rotatedPoints(originalPoints);

    for (int i = 0; i < rotatedPoints.size(); i++)
    {
        rotatedPoints[i] = q.transform(rotatedPoints[i]);
    }

    return rotatedPoints;
}

Field<vector> ITHACAutilities::rotateMesh(fvMesh& mesh, double r1, double r2,
        vector axis,
        double alpha, labelList movingPointsIDs,
        List<double> radii, word angleVariationMethod, double v)
{
    M_Assert(angleVariationMethod == "Linear"
             || angleVariationMethod == "Sinusoidal" || angleVariationMethod == "Sigmoid",
             "The variation function of the angle from the inner radius to the outer radius must be either Linear or Sinusoidal or Sigmoid");
    Field<vector> pointRot(mesh.points());

    for (int i = 0; i < movingPointsIDs.size(); i++)
    {
        double l = radii[i];
        vector pointNow = pointRot[movingPointsIDs[i]];

        if (l <= r1)
        {
            double theta = alpha / 180 *  constant::mathematical::pi;
            quaternion q(axis, theta);
            pointRot[movingPointsIDs[i]] = q.transform(pointNow);
        }
        else if (l > r1 && l <= r2)
        {
            if (angleVariationMethod == "Linear")
            {
                double theta = alpha / 180 *  constant::mathematical::pi * (r2 - l) / (r2 - r1);
                quaternion q(axis, theta);
                pointRot[movingPointsIDs[i]] = q.transform(pointNow);
            }
            else if (angleVariationMethod == "Sinusoidal")
            {
                double theta = alpha / 180 *  constant::mathematical::pi * std::sin(
                                   constant::mathematical::pi / 2 * (r2 - l) / (r2 - r1));
                quaternion q(axis, theta);
                pointRot[movingPointsIDs[i]] = q.transform(pointNow);
            }
            else if (angleVariationMethod == "Sigmoid")
            {
                double theta = alpha / 180 *  constant::mathematical::pi * (1 - 1 /
                               (1 + std::exp(
                                    -v * (l - (r1 + r2) / 2))));
                quaternion q(axis, theta);
                pointRot[movingPointsIDs[i]] = q.transform(pointNow);
            }
        }
    }

    return pointRot;
}

Eigen::MatrixXd ITHACAutilities::rotationMatrix(vector AxisOfRotation,
        double AngleOfRotation)
{
    Eigen::MatrixXd R(3, 3);
    double theta = AngleOfRotation / 180 *  constant::mathematical::pi;
    scalar di = mag(AxisOfRotation);
    scalar ux = AxisOfRotation[0] / di;
    scalar uy = AxisOfRotation[1] / di;
    scalar uz = AxisOfRotation[2] / di;
    R(0, 0) = Foam::cos(theta) + ux * ux * (1 - Foam::cos(theta));
    R(1, 0) = uy * ux * (1 - Foam::cos(theta)) + uz * Foam::sin(theta);
    R(2, 0) = uz * ux * (1 - Foam::cos(theta)) - uy * Foam::sin(theta);
    R(0, 1) = ux * uz * (1 - Foam::cos(theta)) - uz * Foam::sin(theta);
    R(1, 1) = Foam::cos(theta) + uy * uy * (1 - Foam::cos(theta));
    R(2, 1) = ux * uy * (1 - Foam::cos(theta)) + ux * Foam::sin(theta);
    R(0, 2) = ux * uz * (1 - Foam::cos(theta)) + uy * Foam::sin(theta);
    R(1, 2) = uy * uz * (1 - Foam::cos(theta)) - ux * Foam::sin(theta);
    R(2, 2) = Foam::cos(theta) + uz * uz * (1 - Foam::cos(theta));
    return R;
}

Eigen::VectorXd ITHACAutilities::boudaryFaceToCellDistance(
    fvMesh& mesh, label BC_ind)
{
    Eigen::VectorXd cellFaceDistance;
    const polyPatch& cPatch = mesh.boundaryMesh()[BC_ind];
    //Starting index of the face in a patch
    label faceId_start = cPatch.start() ;
    //List of cells close to a boundary
    const labelUList& faceCells = cPatch.faceCells();
    cellFaceDistance.conservativeResize(cPatch.size());
    forAll(cPatch, faceI)
    {
        // index of each face
        label faceID = faceId_start + faceI;
        //id of the owner cell having the face
        label faceOwner = faceCells[faceI] ;
        cellFaceDistance(faceI) = mag(mesh.C()[faceOwner] - mesh.Cf()[faceID]);
    }
    return (cellFaceDistance);
}


template<>
void ITHACAutilities::assignMixedBC(
    GeometricField<scalar, fvPatchField, volMesh>& field, label BC_ind,
    Eigen::MatrixXd& value, Eigen::MatrixXd& grad, Eigen::MatrixXd& valueFrac)
{
    std::string message = "Patch is NOT mixed. It is of type: " +
                          field.boundaryField()[BC_ind].type();
    M_Assert(field.boundaryField()[BC_ind].type() == "mixed", message.c_str());

    if (field.boundaryField()[BC_ind].type() == "mixed")
    {
        mixedFvPatchScalarField& Tpatch =
            refCast<mixedFvPatchScalarField>(field.boundaryFieldRef()[BC_ind]);
        scalarField& valueTpatch = Tpatch.refValue();
        scalarField& gradTpatch = Tpatch.refGrad();
        scalarField& valueFracTpatch = Tpatch.valueFraction();
        Foam2Eigen::Eigen2field(valueTpatch, value);
        Foam2Eigen::Eigen2field(gradTpatch, grad);
        Foam2Eigen::Eigen2field(valueFracTpatch, valueFrac);
    }
}

template<>
void ITHACAutilities::assignMixedBC(
    GeometricField<vector, fvPatchField, volMesh>& field, label BC_ind,
    Eigen::MatrixXd& value, Eigen::MatrixXd& grad, Eigen::MatrixXd& valueFrac)
{
    std::string message = "Patch is NOT mixed. It is of type: " +
                          field.boundaryField()[BC_ind].type();
    M_Assert(field.boundaryField()[BC_ind].type() == "mixed", message.c_str());

    if (field.boundaryField()[BC_ind].type() == "mixed")
    {
        mixedFvPatchVectorField& Tpatch =
            refCast<mixedFvPatchVectorField>(field.boundaryFieldRef()[BC_ind]);
        vectorField& valueTpatch = Tpatch.refValue();
        vectorField& gradTpatch = Tpatch.refGrad();
        scalarField& valueFracTpatch = Tpatch.valueFraction();
        Foam2Eigen::Eigen2field(valueTpatch, value);
        Foam2Eigen::Eigen2field(gradTpatch, grad);
        Foam2Eigen::Eigen2field(valueFracTpatch, valueFrac);
    }
}

template<>
void ITHACAutilities::assignMixedBC(
    GeometricField<scalar, fvPatchField, volMesh>& field, label BC_ind,
    List<scalar>& value, List<scalar>& grad, List<scalar>& valueFrac)
{
    std::string message = "Patch is NOT mixed. It is of type: " +
                          field.boundaryField()[BC_ind].type();
    M_Assert(field.boundaryField()[BC_ind].type() == "mixed", message.c_str());

    if (field.boundaryField()[BC_ind].type() == "mixed")
    {
        mixedFvPatchScalarField& Tpatch =
            refCast<mixedFvPatchScalarField>(field.boundaryFieldRef()[BC_ind]);
        scalarField& valueTpatch = Tpatch.refValue();
        scalarField& gradTpatch = Tpatch.refGrad();
        scalarField& valueFracTpatch = Tpatch.valueFraction();
        valueTpatch = value;
        gradTpatch = grad;
        valueFracTpatch = valueFrac;
    }
}

template<>
void ITHACAutilities::assignMixedBC(
    GeometricField<vector, fvPatchField, volMesh>& field, label BC_ind,
    List<vector>& value, List<vector>& grad, List<scalar>& valueFrac)
{
    std::string message = "Patch is NOT mixed. It is of type: " +
                          field.boundaryField()[BC_ind].type();
    M_Assert(field.boundaryField()[BC_ind].type() == "mixed", message.c_str());

    if (field.boundaryField()[BC_ind].type() == "mixed")
    {
        mixedFvPatchVectorField& Tpatch =
            refCast<mixedFvPatchVectorField>(field.boundaryFieldRef()[BC_ind]);
        vectorField& valueTpatch = Tpatch.refValue();
        vectorField& gradTpatch = Tpatch.refGrad();
        scalarField& valueFracTpatch = Tpatch.valueFraction();
        valueTpatch = value;
        gradTpatch = grad;
        valueFracTpatch = valueFrac;
    }
}

List<int> ITHACAutilities::getIndicesFromDisc(fvMesh& mesh, double radius,
        vector origin, vector axis, List<double>& radii)
{
    pointField meshPoints(mesh.points());
    List<int> indices(meshPoints.size());
    radii.resize(meshPoints.size());
    int k = 0;

    for (int i = 0; i < meshPoints.size(); i++)
    {
        vector pointNow = meshPoints[i];
        vector r = pointNow - origin;
        vector projComponent = (r & axis) / (pow(mag(axis), 2)) * axis;
        vector d = r - projComponent;
        double l = Foam::sqrt(pow(d[0],
                                  2) + pow(d[1], 2) + pow(d[2], 2));

        if (l < radius)
        {
            indices[k] = i;
            radii[k] = l;
            k++;
        }
    }

    indices.resize(k);
    radii.resize(k);
    return indices;
}

template<typename type_f>
List<int> ITHACAutilities::getIndicesFromBox(
    GeometricField<type_f, fvPatchField, volMesh>& field, Eigen::MatrixXd Box)
{
    M_Assert(Box.rows() == 2
             && Box.cols() == 3,
             "The box must be a 2*3 matrix shaped in this way: \nBox = \t|x0, y0, z0|\n\t|x1, yi, z1|\n");
    List<int> indices(field.internalField().size());
    int k = 0;

    for (label i = 0; i < field.internalField().size(); i++)
    {
        auto cx = field.mesh().C()[i][0];
        auto cy = field.mesh().C()[i][1];
        auto cz = field.mesh().C()[i][2];

        if (cx >= Box(0, 0) && cy >= Box(0, 1) && cz >= Box(0, 2) && cx <= Box(1, 0)
                && cy <= Box(1, 1) && cz <= Box(1, 2) )
        {
            indices[k] = i;
            k++;
        }
    }

    indices.resize(k);
    return indices;
}

template<typename type_f>
fvMeshSubset* ITHACAutilities::getSubMeshFromBox(
    GeometricField<type_f, fvPatchField, volMesh>& field, Eigen::MatrixXd Box)
{
    List<int> indices = getIndicesFromBox(field, Box);
    fvMeshSubset* sub;
    sub = new fvMeshSubset(field.mesh());
    (field.mesh());
#if OPENFOAM >= 1812
    sub->setCellSubset(indices);
#else
    sub->setLargeCellSubset(indices);
#endif
    return sub;
}

template<typename type_f>
void ITHACAutilities::normalizeFields(
    PtrList<GeometricField<type_f, fvPatchField, volMesh>>& fields)
{
    Eigen::MatrixXd eigenFields = Foam2Eigen::PtrList2Eigen(fields);
    List<Eigen::MatrixXd> eigenFieldsBC = Foam2Eigen::PtrList2EigenBC(fields);

    for (label i = 0; i < fields.size(); i++)
    {
        double norm = L2norm(fields[i]);
        GeometricField<type_f, fvPatchField, volMesh> tmp(fields[0].name(),
                fields[0] * 0);
        Eigen::VectorXd vec = eigenFields.col(i) / norm;
        tmp = Foam2Eigen::Eigen2field(tmp, vec);

        // Adjusting boundary conditions
        for (int k = 0; k < tmp.boundaryField().size(); k++)
        {
            Eigen::MatrixXd vec = eigenFieldsBC[k].col(i) / norm;
            ITHACAutilities::assignBC(tmp, k, vec);
        }

        fields.set(i, tmp);
    }
}

template void ITHACAutilities::normalizeFields(
    PtrList<GeometricField<scalar, fvPatchField, volMesh>>& fields);
template void ITHACAutilities::normalizeFields(
    PtrList<GeometricField<vector, fvPatchField, volMesh>>& fields);

template Eigen::MatrixXd ITHACAutilities::getCoeffsFrobenius(
    PtrList<volScalarField> snapshots,
    PtrList<volScalarField>& modes, int nModes);

template Eigen::MatrixXd ITHACAutilities::getCoeffsFrobenius(
    PtrList<volVectorField> snapshots,
    PtrList<volVectorField>& modes, int nModes);

template fvMeshSubset* ITHACAutilities::getSubMeshFromBox(
    GeometricField<scalar, fvPatchField, volMesh>& field, Eigen::MatrixXd Box);
template fvMeshSubset* ITHACAutilities::getSubMeshFromBox(
    GeometricField<vector, fvPatchField, volMesh>& field, Eigen::MatrixXd Box);

template List<int> ITHACAutilities::getIndicesFromBox(
    GeometricField<scalar, fvPatchField, volMesh>& field, Eigen::MatrixXd Box);
template List<int> ITHACAutilities::getIndicesFromBox(
    GeometricField<vector, fvPatchField, volMesh>& field, Eigen::MatrixXd Box);

template PtrList<volScalarField> ITHACAutilities::reconstruct_from_coeff(
    PtrList<volScalarField>& modes, Eigen::MatrixXd& coeff_matrix, label Nmodes);
template PtrList<volVectorField> ITHACAutilities::reconstruct_from_coeff(
    PtrList<volVectorField>& modes, Eigen::MatrixXd& coeff_matrix, label Nmodes);

template PtrList<volScalarField> ITHACAutilities::averageSubtract(
    PtrList<volScalarField>
    fields, Eigen::MatrixXd ind, PtrList<volScalarField>& ave);
template PtrList<volVectorField> ITHACAutilities::averageSubtract(
    PtrList<volVectorField>
    fields, Eigen::MatrixXd ind, PtrList<volVectorField>& ave);
template volVectorField ITHACAutilities::computeAverage(
    PtrList<volVectorField>& fields);
template volScalarField ITHACAutilities::computeAverage(
    PtrList<volScalarField>& fields);

template double ITHACAutilities::error_fields(volScalarField& field1,
        volScalarField& field2);
template double ITHACAutilities::error_fields(volVectorField& field1,
        volVectorField& field2);

template double ITHACAutilities::error_fields_abs(volScalarField& field1,
        volScalarField& field2);
template double ITHACAutilities::error_fields_abs(volVectorField& field1,
        volVectorField& field2);

template Eigen::MatrixXd ITHACAutilities::error_listfields(
    PtrList<GeometricField<scalar, fvPatchField, volMesh>>& fields1,
    PtrList<GeometricField<scalar, fvPatchField, volMesh>>& fields2,
    PtrList<volScalarField>& Volumes);
template Eigen::MatrixXd ITHACAutilities::error_listfields(
    PtrList<GeometricField<vector, fvPatchField, volMesh>>& fields1,
    PtrList<GeometricField<vector, fvPatchField, volMesh>>& fields2,
    PtrList<volScalarField>& Volumes);


template Eigen::MatrixXd ITHACAutilities::error_listfields(
    PtrList<volScalarField>& fields1,
    PtrList<volScalarField>& fields2);
template Eigen::MatrixXd ITHACAutilities::error_listfields(
    PtrList<volVectorField>& fields1,
    PtrList<volVectorField>& fields2);

template Eigen::MatrixXd ITHACAutilities::error_listfields_abs(
    PtrList<volScalarField>& fields1,
    PtrList<volScalarField>& fields2);
template Eigen::MatrixXd ITHACAutilities::error_listfields_abs(
    PtrList<volVectorField>& fields1,
    PtrList<volVectorField>& fields2);

template Eigen::VectorXd ITHACAutilities::get_mass_matrix_FV(
    GeometricField<scalar, fvPatchField, volMesh>& snapshot);
template Eigen::VectorXd ITHACAutilities::get_mass_matrix_FV(
    GeometricField<vector, fvPatchField, volMesh>& snapshot);

template Eigen::MatrixXd ITHACAutilities::get_coeffs(PtrList<volScalarField>
        snapshots, PtrList<volScalarField> modes, int Nmodes);
template Eigen::MatrixXd ITHACAutilities::get_coeffs(PtrList<volVectorField>
        snapshots, PtrList<volVectorField> modes, int Nmodes);

template Eigen::MatrixXd ITHACAutilities::get_coeffs_ortho(
    PtrList<volScalarField> snapshots, PtrList<volScalarField>& modes, int Nmodes);
template Eigen::MatrixXd ITHACAutilities::get_coeffs_ortho(
    PtrList<volVectorField> snapshots, PtrList<volVectorField>& modes, int Nmodes);

template void ITHACAutilities::changeBCtype<scalar>
(GeometricField<scalar, fvPatchField, volMesh>& field, word BCtype,
 label BC_ind);
template void ITHACAutilities::changeBCtype<vector>
(GeometricField<vector, fvPatchField, volMesh>& field, word BCtype,
 label BC_ind);

