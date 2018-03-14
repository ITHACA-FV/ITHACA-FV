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

List<int> ITHACAutilities::getIndices(fvMesh& mesh, int index_row, int index_col, int layers)
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
    mkDir(folder);
    word command1("ln -rs  $(readlink -f constant/) " + folder + "/" + " >/dev/null 2>&1");
    word command2("ln -rs  $(readlink -f system/) " + folder + "/" + " >/dev/null 2>&1");
    word command3("ln -rs  $(readlink -f 0/) " + folder + "/" + " >/dev/null 2>&1");
    std::cout.setstate(std::ios_base::failbit);
    system(command1);
    system(command2);
    system(command3);
    std::cout.clear();
    
}

Eigen::MatrixXd ITHACAutilities::rand(int rows, int cols, double min, double max)
{
    Eigen::MatrixXd matr = Eigen::MatrixXd::Random(rows, cols);
    matr = (matr.array() + 1) / 2;
    matr = matr.array() * (max - min);
    matr = matr.array() + min;
    return matr;
}

Eigen::MatrixXd ITHACAutilities::rand(int rows, Eigen::MatrixXd minMax)
{
    int cols = minMax.rows();
    Eigen::MatrixXd matr = Eigen::MatrixXd::Random(rows, cols);
    matr = (matr.array() + 1) / 2;
    for (int i = 0; i<cols; i++)
    {
        matr.col(i) = matr.col(i).array()*(minMax(i,1)-minMax(i,0));
        matr.col(i) = matr.col(i).array()+(minMax(i,0));
    }
    return matr; 

}

// Check if the modes exists
bool ITHACAutilities::check_pod()
{
    struct stat sb;
    bool pod_exist;
    int c = 0;
    if (stat("./ITHACAoutput/POD", &sb) == 0 && S_ISDIR(sb.st_mode))
    {
        pod_exist = true;
        Info << "POD data already exist, reading existing modes" << endl;
    }
    else
    {
        pod_exist = false;
        Info << "POD don't exist, performing a POD decomposition" << endl;
        mkDir("./ITHACAoutput/POD");
        c += abs(system("ln -s ../../constant ./ITHACAoutput/POD/constant"));
        c += abs(system("ln -s ../../0 ./ITHACAoutput/POD/0"));
        c += abs(system("ln -s ../../system ./ITHACAoutput/POD/system"));
    }
    if (c > 0)
    {
        Info << "Error in system operation" << endl;
        exit(0);
    }
    return pod_exist;
}

// Check if the offline data exist
bool ITHACAutilities::check_off()
{
    struct stat sb;
    bool off_exist;
    int c = 0;
    if (stat("./ITHACAoutput/Offline", &sb) == 0 && S_ISDIR(sb.st_mode))
    {
        off_exist = true;
        Info << "Offline data already exist, reading existing data" << endl;
    }
    else
    {
        off_exist = false;
        Info << "Offline don't exist, performing the Offline Solve" << endl;
        mkDir("./ITHACAoutput/Offline");
        c += abs(system("ln -s ../../constant ./ITHACAoutput/Offline/constant"));
        c += abs(system("ln -s ../../0 ./ITHACAoutput/Offline/0"));
        c += abs(system("ln -s ../../system ./ITHACAoutput/Offline/system"));
    }
    if (c > 0)
    {
        Info << "Error in system operation" << endl;
        exit(0);
    }
    return off_exist;
}

// Check if the supremizer data exist
bool ITHACAutilities::check_sup()
{
    struct stat sb;
    bool sup_exist;
    int c = 0;
    if (stat("./ITHACAoutput/supremizer", &sb) == 0 && S_ISDIR(sb.st_mode))
    {
        sup_exist = true;
        Info << "Supremizer data already exist, reading existing data" << endl;
    }
    else
    {
        sup_exist = false;
        Info << "Supremizers don't exist, performing a POD decomposition" << endl;
        mkDir("./ITHACAoutput/supremizer");
        c += abs(system("ln -s ../../constant ./ITHACAoutput/supremizer/constant"));
        c += abs(system("ln -s ../../0 ./ITHACAoutput/supremizer/0"));
        c += abs(system("ln -s ../../system ./ITHACAoutput/supremizer/system"));
    }
    if (c > 0)
    {
        Info << "Error in system operation" << endl;
        exit(0);
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


double ITHACAutilities::error_fields(volVectorField & field1, volVectorField & field2)
{
    double err = L2norm(field1 - field2) / L2norm(field1);
    return err;
}

Eigen::MatrixXd ITHACAutilities::error_listfields(PtrList<volVectorField>& fields1, PtrList<volVectorField>& fields2)
{
    Eigen::VectorXd err;
    if (fields1.size() != fields2.size())
    {
        Info << "The two fields does not have the same size, code will abort" << endl;
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

double ITHACAutilities::error_fields(volScalarField & field1, volScalarField & field2)
{
    double err = L2norm(field1 - field2) / L2norm(field1);
    return err;
}


Eigen::MatrixXd ITHACAutilities::error_listfields(PtrList<volScalarField>& fields1, PtrList<volScalarField>& fields2)
{
    Eigen::VectorXd err;
    if (fields1.size() != fields2.size())
    {
        Info << "The two fields does not have the same size, code will abort" << endl;
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

Eigen::MatrixXd ITHACAutilities::get_mass_matrix(PtrList<volVectorField> modes)
{
    label Msize = modes.size();
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

Eigen::MatrixXd ITHACAutilities::get_mass_matrix(PtrList<volScalarField> modes)
{
    label Msize = modes.size();
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

Eigen::VectorXd ITHACAutilities::get_coeffs(volVectorField snapshot, PtrList<volVectorField>& modes)
{
    label Msize = modes.size();
    Eigen::MatrixXd M_matrix = get_mass_matrix(modes);
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

Eigen::VectorXd ITHACAutilities::get_coeffs(volScalarField snapshot, PtrList<volScalarField>& modes)
{
    label Msize = modes.size();
    Eigen::MatrixXd M_matrix = get_mass_matrix(modes);
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

Eigen::MatrixXd ITHACAutilities::get_coeffs(PtrList<volScalarField>  snapshots, PtrList<volScalarField> modes)
{
    label Msize = modes.size();
    Eigen::MatrixXd M_matrix = get_mass_matrix(modes);
    Eigen::VectorXd a(Msize);
    Eigen::VectorXd b(Msize);

    Eigen::MatrixXd coeff(modes.size(), snapshots.size());
    for (auto i = 0; i < snapshots.size(); i++)
    {
        for (label k = 0; k < Msize; k++)
        {
            b(k) = fvc::domainIntegrate(snapshots[i] * modes[k]).value();
        }
        coeff.col(i) = M_matrix.ldlt().solve(b);;
    }
    return coeff;
}

Eigen::MatrixXd ITHACAutilities::get_coeffs( PtrList<volVectorField>  snapshots, PtrList<volVectorField> modes)
{
    label Msize = modes.size();
    Eigen::MatrixXd M_matrix = get_mass_matrix(modes);
    Eigen::VectorXd a(Msize);
    Eigen::VectorXd b(Msize);

    Eigen::MatrixXd coeff(modes.size(), snapshots.size());
    for (auto i = 0; i < snapshots.size(); i++)
    {
        for (label k = 0; k < Msize; k++)
        {
            b(k) = fvc::domainIntegrate(snapshots[i] & modes[k]).value();
        }
        coeff.col(i) = M_matrix.ldlt().solve(b);;
    }
    return coeff;
}


Eigen::MatrixXd ITHACAutilities::get_coeffs_ortho(PtrList<volScalarField> snapshots, PtrList<volScalarField>& modes)
{
    Eigen::MatrixXd coeff(modes.size(), snapshots.size());
    for (auto i = 0; i < modes.size(); i++)
    {
        for (auto j = 0; j < snapshots.size(); j++)
        {
            coeff(i, j) = fvc::domainIntegrate(snapshots[j] * modes[i]).value();
        }
    }
    return coeff;
}


Eigen::MatrixXd ITHACAutilities::get_coeffs_ortho(PtrList<volVectorField> snapshots, PtrList<volVectorField>& modes)
{
    Eigen::MatrixXd coeff(modes.size(), snapshots.size());
    for (auto i = 0; i < modes.size(); i++)
    {
        for (auto j = 0; j < snapshots.size(); j++)
        {
            coeff(i, j) = fvc::domainIntegrate(snapshots[j] & modes[i]).value();
        }
    }
    return coeff;
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
    a = Foam::sqrt(fvc::domainIntegrate(fvc::grad(field) & fvc::grad(field)).value());
    return a;
}

double ITHACAutilities::H1seminorm(volVectorField field)
{
    double a;
    a = Foam::sqrt(fvc::domainIntegrate(fvc::grad(field) && fvc::grad(field)).value());
    return a;
}

void ITHACAutilities::setBoxToValue(volScalarField& field, Eigen::MatrixXd Box, double value)
{
    for (label i = 0; i < field.internalField().size(); i++)
    {
        auto cx = field.mesh().C()[i].component(vector::X);
        auto cy = field.mesh().C()[i].component(vector::Y);
        auto cz = field.mesh().C()[i].component(vector::Z);
        if (cx >= Box(0, 0) && cy >= Box(0, 1) && cz >= Box(0, 2) && cx <= Box(1, 0) && cy <= Box(1, 1) && cz <= Box(1, 2) )
        {
            field.ref()[i] = value;
        }
    }
    for (label i = 0; i < field.boundaryField().size(); i++)
    {
        for (label j = 0; j < field.boundaryField()[i].size(); j++)
        {
            if (field.boundaryField()[i].type() == "fixedValue" || field.boundaryField()[i].type() == "calculated")
            {
                auto cx = field.mesh().C().boundaryField()[i][j][0];
                auto cy = field.mesh().C().boundaryField()[i][j][1];
                auto cz = field.mesh().C().boundaryField()[i][j][2];
                if (cx >= Box(0, 0) && cy >= Box(0, 1) && cz >= Box(0, 2) && cx <= Box(1, 0) && cy <= Box(1, 1) && cz <= Box(1, 2) )
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
    if (s.boundaryField()[BC_ind].type() == "fixedValue")
    {
        for (label i = 0; i < s.boundaryField()[BC_ind].size(); i++)
        {
            s.boundaryFieldRef()[BC_ind][i] = value;
        }
    }
    else if (s.boundaryField()[BC_ind].type() == "fixedGradient")
    {
        fixedGradientFvPatchScalarField& Tpatch = refCast<fixedGradientFvPatchScalarField>(s.boundaryFieldRef()[BC_ind]);
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
        freestreamFvPatchField<scalar>& Tpatch = refCast<freestreamFvPatchField<scalar> >(s.boundaryFieldRef()[BC_ind]);
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
    else if (s.boundaryField()[BC_ind].type() == "empty")
    {

    }
}

// Assign a BC for a scalar field
void ITHACAutilities::assignBC(volVectorField& s, label BC_ind, Vector<double>& value)
{
    if (s.boundaryField()[BC_ind].type() == "fixedValue")
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
        freestreamFvPatchField<vector>& Tpatch = refCast<freestreamFvPatchField<vector> >(s.boundaryFieldRef()[BC_ind]);
        vectorField& gradTpatch = Tpatch.freestreamValue();
        forAll(gradTpatch, faceI)
        {
            gradTpatch[faceI] = value;
        }
    }
}

void ITHACAutilities::assignBC(volScalarField& s, label BC_ind, Eigen::MatrixXd valueVec)
{
    word typeBC = s.boundaryField()[BC_ind].type();
    if (typeBC == "fixedValue" || typeBC == "calculated" || typeBC == "fixedFluxPressure" )
    {
        for (label i = 0; i < s.boundaryField()[BC_ind].size(); i++)
        {
            double value = valueVec(i);
            s.boundaryFieldRef()[BC_ind][i] = value;
        }
    }
    else if (s.boundaryField()[BC_ind].type() == "fixedGradient")
    {
        fixedGradientFvPatchScalarField& Tpatch = refCast<fixedGradientFvPatchScalarField>(s.boundaryFieldRef()[BC_ind]);
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
        freestreamFvPatchField<scalar>& Tpatch = refCast<freestreamFvPatchField<scalar> >(s.boundaryFieldRef()[BC_ind]);
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
void ITHACAutilities::assignBC(volVectorField& s, label BC_ind, Eigen::MatrixXd valueVec)
{
    word typeBC = s.boundaryField()[BC_ind].type();
    int sizeBC = s.boundaryField()[BC_ind].size();
    if (typeBC == "fixedValue" || typeBC == "calculated")
    {
        for (label i = 0; i < sizeBC; i++)
        {
            Vector<double> value(valueVec(i),valueVec(i+sizeBC),valueVec(i+sizeBC*2)); 
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
            Vector<double> value(valueVec(i),valueVec(i+sizeBC),valueVec(i+sizeBC*2)); 
            s.boundaryFieldRef()[BC_ind][i] = value;
        }
        freestreamFvPatchField<vector>& Tpatch = refCast<freestreamFvPatchField<vector> >(s.boundaryFieldRef()[BC_ind]);
        vectorField& gradTpatch = Tpatch.freestreamValue();
        forAll(gradTpatch, faceI)
        {
            Vector<double> value(valueVec(faceI),valueVec(faceI+sizeBC),valueVec(faceI+sizeBC*2)); 
            gradTpatch[faceI] = value;
        }
    }
}





// * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //

