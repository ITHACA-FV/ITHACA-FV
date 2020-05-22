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







template PtrList<volScalarField> reconstruct_from_coeff(
    PtrList<volScalarField>& modes, Eigen::MatrixXd& coeff_matrix, label Nmodes);
template PtrList<volVectorField> reconstruct_from_coeff(
    PtrList<volVectorField>& modes, Eigen::MatrixXd& coeff_matrix, label Nmodes);



template Eigen::MatrixXd getCoeffsFrobenius(
    PtrList<volScalarField> snapshots,
    PtrList<volScalarField>& modes, int nModes);

template Eigen::MatrixXd getCoeffsFrobenius(
    PtrList<volVectorField> snapshots,
    PtrList<volVectorField>& modes, int nModes);




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



}

