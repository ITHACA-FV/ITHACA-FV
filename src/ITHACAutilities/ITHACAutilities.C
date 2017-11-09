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

\*---------------------------------------------------------------------------*/

#include "ITHACAutilities.H"

/// \file
/// Source file of the ITHACAutilities class.

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Check if the modes exists
bool ITHACAutilities::check_pod()
{
    struct stat sb;
    bool pod_exist;
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
        system("ln -s ../../constant ./ITHACAoutput/POD/constant");
        system("ln -s ../../0 ./ITHACAoutput/POD/0");
        system("ln -s ../../system ./ITHACAoutput/POD/system");
    }
    return pod_exist;
}

// Check if the offline data exist
bool ITHACAutilities::check_off()
{
    struct stat sb;
    bool off_exist;
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
        system("ln -s ../../constant ./ITHACAoutput/Offline/constant");
        system("ln -s ../../0 ./ITHACAoutput/Offline/0");
        system("ln -s ../../system ./ITHACAoutput/Offline/system");
    }
    return off_exist;
}

// Check if the supremizer data exist
bool ITHACAutilities::check_sup()
{
    struct stat sb;
    bool sup_exist;
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
        system("ln -s ../../constant ./ITHACAoutput/supremizer/constant");
        system("ln -s ../../0 ./ITHACAoutput/supremizer/0");
        system("ln -s ../../system ./ITHACAoutput/supremizer/system");
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
    double err = L2norm(field1 - field2)/L2norm(field1);
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
    err.resize(fields1.size(),1);
    for (label k = 0; k < fields1.size(); k++)
    {
        err(k,0) = error_fields(fields1[k], fields2[k]);
        Info << " Error is " << err[k] << endl;
    }
    return err;
}

double ITHACAutilities::error_fields(volScalarField & field1, volScalarField & field2)
{
    double err = L2norm(field1 - field2)/L2norm(field1);
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
    err.resize(fields1.size(),1);
    for (label k = 0; k < fields1.size(); k++)
    {
        err(k,0) = error_fields(fields1[k], fields2[k]);
        Info << " Error is " << err[k] << endl;
    }
    return err;
}

Eigen::MatrixXd ITHACAutilities::get_mass_matrix(PtrList<volVectorField>& modes)
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

Eigen::MatrixXd ITHACAutilities::get_mass_matrix(PtrList<volScalarField>& modes)
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

Eigen::MatrixXd ITHACAutilities::foam2eigen(PtrList<volVectorField>& fields1)
{
    Eigen::MatrixXd out;
    out.resize(int(fields1[0].size()*3),fields1.size());

    for(int k=0; k<fields1.size(); k++)
    {
        for(int l=0; l<fields1[k].size(); l++)
        {
            out(3*l,k) = fields1[k][l][0];
            out(3*l+1,k) = fields1[k][l][1];
            out(3*l+2,k) = fields1[k][l][2];                 
        }
    }
    return out;
}

Eigen::MatrixXd ITHACAutilities::foam2eigen(PtrList<volScalarField>& fields1)
{
    Eigen::MatrixXd out;
    out.resize(int(fields1[0].size()),fields1.size());

    for(int k=0; k<fields1.size(); k++)
    {
        for(int l=0; l<fields1[k].size(); l++)
        {
            out(l,k) = fields1[k][l];                
        }
    }
    return out;
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

// * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //

