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
#include "ITHACAparameters.H"
#include "turbulentTransportModel.H"

/// \file
/// Source file of the ITHACAutilities namespace.

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

namespace ITHACAutilities
{

double L2norm(volScalarField field)
{
    double a;
    a = Foam::sqrt(fvc::domainIntegrate(field * field).value());
    return a;
}

double L2norm(volVectorField field)
{
    double a;
    a = Foam::sqrt(fvc::domainIntegrate(field & field).value());
    return a;
}



template<class TypeField>
double frobNorm(TypeField& field)
{
    double norm(0);
    Eigen::VectorXd vF = Foam2Eigen::field2Eigen(field);
    norm = vF.norm();
    return norm;
}

template double frobNorm(volScalarField& field);
template double frobNorm(volVectorField& field);








Eigen::MatrixXd rand(int rows, int cols, double min,
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

Eigen::MatrixXd rand(int rows, Eigen::MatrixXd minMax)
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


bool isInteger(double ratio)
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

bool isTurbulent()
{
    bool checkTurb;
    ITHACAparameters* para = ITHACAparameters::getInstance();
    auto& tur =
        para->mesh.lookupObject<incompressible::turbulenceModel>("turbulenceProperties");

    if (tur.type() == "Stokes" || tur.type() == "Maxwell"
            || tur.type() == "laminarModel")
    {
        checkTurb = false;
    }
    else
    {
        checkTurb = true;
    }

    return checkTurb;
}

template<class TypeField>
PtrList<TypeField> reconstruct_from_coeff(
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
double errorFieldsFrob(TypeField& field1,
                       TypeField& field2)
{
    double err;
    TypeField errField = field1 - field2;

    if (frobNorm(field1) <= 1e-6)
    {
        err = 0;
    }
    else
    {
        err = frobNorm(errField) / frobNorm(field1);
    }

    return err;
}

template<class TypeField>
double error_fields(TypeField& field1,
                    TypeField& field2)
{
    double err;

    if (L2norm(field1) <= 1e-6)
    {
        err = 0;
    }
    else
    {
        err = L2norm(field1 - field2) / L2norm(
                  field1);
    }

    return err;
}

template double error_fields(volScalarField& field1, volScalarField& field2);
template double error_fields(volVectorField& field1, volVectorField& field2);

template<>
double error_fields(
    GeometricField<vector, fvPatchField, volMesh>& field1,
    GeometricField<vector, fvPatchField, volMesh>& field2, volScalarField& Volumes)

{
    volScalarField diffFields2 = ((field1 - field2) & (field1 - field2)) * Volumes;
    double err = Foam::sqrt(gSum(diffFields2));
    return err;
}

template<>
double error_fields(
    GeometricField<scalar, fvPatchField, volMesh>& field1,
    GeometricField<scalar, fvPatchField, volMesh>& field2, volScalarField& Volumes)

{
    volScalarField diffFields2 = ((field1 - field2) * (field1 - field2)) * Volumes;
    double err = Foam::sqrt(gSum(diffFields2));
    return err;
}

template<class TypeField>
double error_fields_abs(TypeField& field1,
                        TypeField& field2)
{
    double err = L2norm(field1 - field2);
    return err;
}

template<class TypeField>
Eigen::MatrixXd error_listfields(PtrList<TypeField>&
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
Eigen::MatrixXd errorListFieldsFrob(PtrList<TypeField>&
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
        err(k, 0) = errorFieldsFrob(fields1[k], fields2[k]);
        Info << " Error is " << err[k] << endl;
    }

    return err;
}

template<class TypeField>
Eigen::MatrixXd error_listfields(
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
Eigen::MatrixXd error_listfields_abs(PtrList<TypeField>&
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

Eigen::MatrixXd get_mass_matrix(PtrList<volVectorField> modes,
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

Eigen::MatrixXd get_mass_matrix(PtrList<volScalarField> modes,
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
Eigen::VectorXd get_mass_matrix_FV(
    GeometricField<TypeField, fvPatchField, volMesh>& snapshot)
{
    Eigen::MatrixXd snapEigen = Foam2Eigen::field2Eigen(snapshot);
    int dim = std::nearbyint(snapEigen.rows() / (snapshot.mesh().V()).size());
    Eigen::VectorXd volumes = Foam2Eigen::field2Eigen(snapshot.mesh().V());
    Eigen::VectorXd vol3 = volumes.replicate(dim, 1);
    return vol3;
}

Eigen::VectorXd get_coeffs(volVectorField snapshot,
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

Eigen::VectorXd get_coeffs(volScalarField snapshot,
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
Eigen::MatrixXd get_coeffs(PtrList<TypeField> snapshots,
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
Eigen::VectorXd get_coeffs_ortho(volScalarField snapshot,
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
Eigen::VectorXd get_coeffs_ortho(volVectorField
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
Eigen::MatrixXd get_coeffs_ortho(PtrList<TypeField>
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
Eigen::MatrixXd getCoeffsFrobenius(PtrList<TypeField>
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

    Eigen::MatrixXd ModesE = (Foam2Eigen::PtrList2Eigen(modes)).leftCols(Msize);
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

double H1seminorm(volScalarField field)
{
    double a;
    a = Foam::sqrt(fvc::domainIntegrate(fvc::grad(field) & fvc::grad(
                                            field)).value());
    return a;
}

double H1seminorm(volVectorField field)
{
    double a;
    a = Foam::sqrt(fvc::domainIntegrate(fvc::grad(field)
                                        && fvc::grad(field)).value());
    return a;
}

void setBoxToValue(volScalarField& field, Eigen::MatrixXd Box,
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

void assignONE(volScalarField& s, List<int>& L)
{
    for (label i = 0; i < L.size(); i++)
    {
        s.ref()[L[i]] = 1;
    }
}

// Assign a BC for a vector field
void assignBC(volScalarField& s, label BC_ind, double& value)
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
void assignBC(volScalarField& s, label BC_ind,
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
void assignBC(volVectorField& s, label BC_ind,
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

void assignBC(volScalarField& s, label BC_ind,
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
void assignBC(volVectorField& s, label BC_ind,
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
PtrList<TypeField> averageSubtract(PtrList<TypeField>
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
TypeField computeAverage(PtrList<TypeField>& fields)
{
    TypeField av(fields[0]);

    for (int i = 1; i < fields.size(); i++)
    {
        av += fields[i];
    }

    av = av / fields.size();
    return av;
}



Eigen::MatrixXd parTimeCombMat(List<Eigen::VectorXd>
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
void changeBCtype(
    GeometricField<TypeField, fvPatchField, volMesh>& field, word BCtype,
    label BC_ind)
{
    field.boundaryFieldRef().set(BC_ind, fvPatchField<TypeField>::New(BCtype,
                                 field.mesh().boundary()[BC_ind], field));
}



void setIndices2Value(labelList& ind2set,
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

volScalarField meshNonOrtho(fvMesh& mesh,
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


List<vector> rotatePoints(const List<vector>& originalPoints,
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




template<typename T>
void assignIF(GeometricField<T, fvPatchField, volMesh>& s,
              T& value)
{
    for (label i = 0; i < s.internalField().size(); i++)
    {
        s.ref()[i] = value;
    }
}

template void assignIF(
    GeometricField<scalar, fvPatchField, volMesh>& field, scalar& value);
template void assignIF(
    GeometricField<vector, fvPatchField, volMesh>& field, vector& value);

template<typename T>
void assignIF(GeometricField<T, fvPatchField, volMesh>& s,
              T& value, List<int>& indices)
{
    for (label i = 0; i < indices.size(); i++)
    {
        s.ref()[indices[i]] = value;
    }
}

template<typename T>
void assignIF(GeometricField<T, fvPatchField, volMesh>& s, T& value, int index)
{
    s.ref()[index] = value;
}

template void assignIF(GeometricField<scalar, fvPatchField, volMesh>& field,
                       scalar& value, int index);
template void assignIF(GeometricField<vector, fvPatchField, volMesh>& field,
                       vector& value, int index);


template<typename T>
Eigen::MatrixXd get_mass_matrix_Eigen(PtrList<T>& fields,
                                      bool consider_volumes)
{
    Eigen::MatrixXd F = Foam2Eigen::PtrList2Eigen(fields);
    Eigen::VectorXd V = Foam2Eigen::field2Eigen(fields[0].mesh().V());
    Eigen::MatrixXd VM = V.asDiagonal();
    Eigen::MatrixXd M;

    if (consider_volumes)
    {
        M = F.transpose() * VM * F;
    }
    else
    {
        M = F.transpose() * F;
    }

    return M;
}


void assignMixedBC(
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

void assignMixedBC(
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

void assignMixedBC(
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

void assignMixedBC(
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


template<typename type_f>
void normalizeFields(
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
            assignBC(tmp, k, vec);
        }

        fields.set(i, tmp);
    }
}

template PtrList<volScalarField> reconstruct_from_coeff(
    PtrList<volScalarField>& modes, Eigen::MatrixXd& coeff_matrix, label Nmodes);
template PtrList<volVectorField> reconstruct_from_coeff(
    PtrList<volVectorField>& modes, Eigen::MatrixXd& coeff_matrix, label Nmodes);

template void normalizeFields(
    PtrList<GeometricField<scalar, fvPatchField, volMesh>>& fields);
template void normalizeFields(
    PtrList<GeometricField<vector, fvPatchField, volMesh>>& fields);

template Eigen::MatrixXd getCoeffsFrobenius(
    PtrList<volScalarField> snapshots,
    PtrList<volScalarField>& modes, int nModes);

template Eigen::MatrixXd getCoeffsFrobenius(
    PtrList<volVectorField> snapshots,
    PtrList<volVectorField>& modes, int nModes);

template PtrList<volScalarField> averageSubtract(
    PtrList<volScalarField>
    fields, Eigen::MatrixXd ind, PtrList<volScalarField>& ave);
template PtrList<volVectorField> averageSubtract(
    PtrList<volVectorField>
    fields, Eigen::MatrixXd ind, PtrList<volVectorField>& ave);
template volVectorField computeAverage(
    PtrList<volVectorField>& fields);
template volScalarField computeAverage(
    PtrList<volScalarField>& fields);

template double errorFieldsFrob(volScalarField& field1,
                                volScalarField& field2);
template double errorFieldsFrob(volVectorField& field1,
                                volVectorField& field2);

template double error_fields_abs(volScalarField& field1,
                                 volScalarField& field2);
template double error_fields_abs(volVectorField& field1,
                                 volVectorField& field2);

template Eigen::MatrixXd error_listfields(
    PtrList<GeometricField<scalar, fvPatchField, volMesh>>& fields1,
    PtrList<GeometricField<scalar, fvPatchField, volMesh>>& fields2,
    PtrList<volScalarField>& Volumes);
template Eigen::MatrixXd error_listfields(
    PtrList<GeometricField<vector, fvPatchField, volMesh>>& fields1,
    PtrList<GeometricField<vector, fvPatchField, volMesh>>& fields2,
    PtrList<volScalarField>& Volumes);


template Eigen::MatrixXd error_listfields(
    PtrList<volScalarField>& fields1,
    PtrList<volScalarField>& fields2);
template Eigen::MatrixXd error_listfields(
    PtrList<volVectorField>& fields1,
    PtrList<volVectorField>& fields2);

template Eigen::MatrixXd errorListFieldsFrob(
    PtrList<volScalarField>& fields1,
    PtrList<volScalarField>& fields2);
template Eigen::MatrixXd errorListFieldsFrob(
    PtrList<volVectorField>& fields1,
    PtrList<volVectorField>& fields2);

template Eigen::MatrixXd error_listfields_abs(
    PtrList<volScalarField>& fields1,
    PtrList<volScalarField>& fields2);
template Eigen::MatrixXd error_listfields_abs(
    PtrList<volVectorField>& fields1,
    PtrList<volVectorField>& fields2);

template Eigen::VectorXd get_mass_matrix_FV(
    GeometricField<scalar, fvPatchField, volMesh>& snapshot);
template Eigen::VectorXd get_mass_matrix_FV(
    GeometricField<vector, fvPatchField, volMesh>& snapshot);

template Eigen::MatrixXd get_coeffs(PtrList<volScalarField>
                                    snapshots, PtrList<volScalarField> modes, int Nmodes);
template Eigen::MatrixXd get_coeffs(PtrList<volVectorField>
                                    snapshots, PtrList<volVectorField> modes, int Nmodes);

template Eigen::MatrixXd get_coeffs_ortho(
    PtrList<volScalarField> snapshots, PtrList<volScalarField>& modes, int Nmodes);
template Eigen::MatrixXd get_coeffs_ortho(
    PtrList<volVectorField> snapshots, PtrList<volVectorField>& modes, int Nmodes);

template void changeBCtype<scalar>
(GeometricField<scalar, fvPatchField, volMesh>& field, word BCtype,
 label BC_ind);
template void changeBCtype<vector>
(GeometricField<vector, fvPatchField, volMesh>& field, word BCtype,
 label BC_ind);

}

