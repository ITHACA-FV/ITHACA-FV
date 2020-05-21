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

#include "ITHACAerror.H"

/// \file
/// Source file of the ITHACAerror file.

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

namespace ITHACAutilities
{

template<typename T>
double errorFrobRel(GeometricField<T, fvPatchField, volMesh>& field1,
                    GeometricField<T, fvPatchField, volMesh>& field2)
{
    double err;
    GeometricField<T, fvPatchField, volMesh> errField = field1 - field2;

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

template double errorFrobRel(GeometricField<scalar, fvPatchField, volMesh>&
                             field1,
                             GeometricField<scalar, fvPatchField, volMesh>& field2);
template double errorFrobRel(GeometricField<vector, fvPatchField, volMesh>&
                             field1,
                             GeometricField<vector, fvPatchField, volMesh>& field2);

template<typename T>
double errorL2Rel(GeometricField<T, fvPatchField, volMesh>& field1,
                  GeometricField<T, fvPatchField, volMesh>& field2)
{
    double err;
    GeometricField<T, fvPatchField, volMesh> errField = field1 - field2;

    if (L2Norm(field1) <= 1e-6)
    {
        err = 0;
    }
    else
    {
        err = L2Norm(errField) / L2Norm(
                  field1);
    }

    return err;
}

template double errorL2Rel(GeometricField<scalar, fvPatchField, volMesh>&
                           field1,
                           GeometricField<scalar, fvPatchField, volMesh>& field2);
template double errorL2Rel(GeometricField<vector, fvPatchField, volMesh>&
                           field1,
                           GeometricField<vector, fvPatchField, volMesh>& field2);

template<>
double errorL2Abs(GeometricField<vector, fvPatchField, volMesh>& field1,
                  GeometricField<vector, fvPatchField, volMesh>& field2, volScalarField& Volumes)

{
    volScalarField diffFields2 = ((field1 - field2) & (field1 - field2)) * Volumes;
    double err = Foam::sqrt(gSum(diffFields2));
    return err;
}

template<>
double errorL2Abs(GeometricField<scalar, fvPatchField, volMesh>& field1,
                  GeometricField<scalar, fvPatchField, volMesh>& field2, volScalarField& Volumes)

{
    volScalarField diffFields2 = ((field1 - field2) * (field1 - field2)) * Volumes;
    double err = Foam::sqrt(gSum(diffFields2));
    return err;
}

template<typename T>
double errorL2Abs(GeometricField<T, fvPatchField, volMesh>& field1,
                  GeometricField<T, fvPatchField, volMesh>& field2)
{
    GeometricField<T, fvPatchField, volMesh> errField = field1 - field2;
    double err = L2Norm(errField);
    return err;
}

template double errorL2Abs(GeometricField<scalar, fvPatchField, volMesh>&
                           field1,
                           GeometricField<scalar, fvPatchField, volMesh>& field2);
template double errorL2Abs(GeometricField<vector, fvPatchField, volMesh>&
                           field1,
                           GeometricField<vector, fvPatchField, volMesh>& field2);

template<typename T>
Eigen::MatrixXd errorL2Rel(PtrList<GeometricField<T, fvPatchField, volMesh>>&
                           fields1, PtrList<GeometricField<T, fvPatchField, volMesh>>& fields2)
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
        err(k, 0) = errorL2Rel(fields1[k], fields2[k]);
        Info << " Error is " << err[k] << endl;
    }

    return err;
}

template Eigen::MatrixXd errorL2Rel(
    PtrList<GeometricField<scalar, fvPatchField, volMesh>>& fields1,
    PtrList<GeometricField<scalar, fvPatchField, volMesh>>& fields2);
template Eigen::MatrixXd errorL2Rel(
    PtrList<GeometricField<vector, fvPatchField, volMesh>>& fields1,
    PtrList<GeometricField<vector, fvPatchField, volMesh>>& fields2);

template<typename T>
Eigen::MatrixXd errorFrobRel(PtrList<GeometricField<T, fvPatchField, volMesh>>&
                             fields1,
                             PtrList<GeometricField<T, fvPatchField, volMesh>>& fields2)
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
        err(k, 0) = errorFrobRel(fields1[k], fields2[k]);
        Info << " Error is " << err[k] << endl;
    }

    return err;
}

template Eigen::MatrixXd errorFrobRel(
    PtrList<GeometricField<scalar, fvPatchField, volMesh>>& fields1,
    PtrList<GeometricField<scalar, fvPatchField, volMesh>>& fields2);
template Eigen::MatrixXd errorFrobRel(
    PtrList<GeometricField<vector, fvPatchField, volMesh>>& fields1,
    PtrList<GeometricField<vector, fvPatchField, volMesh>>& fields2);


template<typename T>
Eigen::MatrixXd errorL2Abs(
    PtrList<GeometricField<T, fvPatchField, volMesh>>& fields1,
    PtrList<GeometricField<T, fvPatchField, volMesh>>& fields2,
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
        err(k, 0) = errorL2Abs(fields1[k], fields2[k], Volumes[k]);
        Info << " Error is " << err[k] << endl;
    }

    return err;
}

template Eigen::MatrixXd errorL2Abs(
    PtrList<GeometricField<scalar, fvPatchField, volMesh>>& fields1,
    PtrList<GeometricField<scalar, fvPatchField, volMesh>>& fields2,
    PtrList<volScalarField>& Volumes);
template Eigen::MatrixXd errorL2Abs(
    PtrList<GeometricField<vector, fvPatchField, volMesh>>& fields1,
    PtrList<GeometricField<vector, fvPatchField, volMesh>>& fields2,
    PtrList<volScalarField>& Volumes);

template<typename T>
Eigen::MatrixXd errorL2Abs(PtrList<GeometricField<T, fvPatchField, volMesh>>&
                           fields1,
                           PtrList<GeometricField<T, fvPatchField, volMesh>>& fields2)
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
        err(k, 0) = errorL2Abs(fields1[k], fields2[k]);
        Info << " Error is " << err[k] << endl;
    }

    return err;
}

template Eigen::MatrixXd errorL2Abs(
    PtrList<GeometricField<scalar, fvPatchField, volMesh>>& fields1,
    PtrList<GeometricField<scalar, fvPatchField, volMesh>>& fields2);
template Eigen::MatrixXd errorL2Abs(
    PtrList<GeometricField<vector, fvPatchField, volMesh>>& fields1,
    PtrList<GeometricField<vector, fvPatchField, volMesh>>& fields2);

template<>
double L2Norm(GeometricField<scalar, fvPatchField, volMesh>& field)
{
    double a;
    a = Foam::sqrt(fvc::domainIntegrate(field * field).value());
    return a;
}

template<>
double L2Norm(GeometricField<vector, fvPatchField, volMesh>& field)
{
    double a;
    a = Foam::sqrt(fvc::domainIntegrate(field & field).value());
    return a;
}

template<>
double H1Seminorm(GeometricField<scalar, fvPatchField, volMesh>& field)
{
    double a;
    a = Foam::sqrt(fvc::domainIntegrate(fvc::grad(field) & fvc::grad(
                                            field)).value());
    return a;
}

template<>
double H1Seminorm(GeometricField<vector, fvPatchField, volMesh>& field)
{
    double a;
    a = Foam::sqrt(fvc::domainIntegrate(fvc::grad(field)
                                        && fvc::grad(field)).value());
    return a;
}

template<class T>
double frobNorm(GeometricField<T, fvPatchField, volMesh>& field)
{
    double norm(0);
    Eigen::VectorXd vF = Foam2Eigen::field2Eigen(field);
    norm = vF.norm();
    return norm;
}

template double frobNorm(GeometricField<scalar, fvPatchField, volMesh>& field);
template double frobNorm(GeometricField<vector, fvPatchField, volMesh>& field);

}