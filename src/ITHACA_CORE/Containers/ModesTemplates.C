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
    EigenFunctions
Description
    Container that contains a list of modes with additional operations
SourceFiles
    Modes.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file containing the implementation of the template functions.

template<class T>
List<Eigen::MatrixXd> Modes<T>::toEigen()
{
    NBC = 0;

    for (int i = 0; i < (this->first()).boundaryFieldRef().size(); i++)
    {
        if ((this->first()).boundaryFieldRef()[i].type() != "processor")
        {
            NBC++;
        }
    }

    EigenModes.resize(NBC + 1);
    EigenModes[0] = Foam2Eigen::PtrList2Eigen(this->toPtrList());
    List<Eigen::MatrixXd> BC = Foam2Eigen::PtrList2EigenBC(this->toPtrList());

    for (int i = 0; i < NBC; i++)
    {
        EigenModes[i + 1] = BC[i];
    }

    return EigenModes;
}

template<class T>
List<Eigen::MatrixXd> Modes<T>::project(fvMatrix<T>& Af)
{
    List<Eigen::MatrixXd> LinSys;
    LinSys.resize(2);

    if (EigenModes.size() == 0)
    {
        toEigen();
    }

    Eigen::SparseMatrix<double> Ae;
    Eigen::VectorXd be;
    Foam2Eigen::fvMatrix2Eigen(Af, Ae, be);
    LinSys[0] = EigenModes[0].transpose() * Ae * EigenModes[0];
    LinSys[1] = EigenModes[0].transpose() * be;
    return LinSys;
}


template<class T>
GeometricField<T, fvPatchField, volMesh> Modes<T>::reconstruct(
    Eigen::MatrixXd& Coeff,
    word Name)
{
    if (EigenModes.size() == 0)
    {
        toEigen();
    }

    GeometricField<T, fvPatchField, volMesh> field(Name,
            (this->toPtrList())[0] * 0);
    int Nmodes = Coeff.rows();
    Eigen::VectorXd InField = EigenModes[0].leftCols(Nmodes) * Coeff;
    field = Foam2Eigen::Eigen2field(field, InField);

    for (int i = 0; i < 1; i++)
    {
        Eigen::VectorXd BF = EigenModes[i + 1].leftCols(Nmodes) * Coeff;
        ITHACAutilities::assignBC(field, i, BF);
    }

    return field;
}

template<class T>
PtrList<GeometricField<T, fvPatchField, volMesh>>
        Modes<T>::projectSnapshots(
            PtrList<GeometricField<T, fvPatchField, volMesh>> snapshots,
            int numberOfModes,
            word innerProduct)
{
    if (numberOfModes == 0)
    {
        numberOfModes == this->size();
    }

    if (EigenModes.size() == 0)
    {
        toEigen();
    }

    M_Assert(numberOfModes <= this->size(),
             "The number of Modes used for the projection cannot be bigger than the number of available modes");
    Eigen::MatrixXd M_vol;
    Eigen::MatrixXd M;
    PtrList<GeometricField<T, fvPatchField, volMesh>> projSnap;
    Eigen::MatrixXd projSnapI;
    Eigen::MatrixXd projSnapCoeff;

    for (int i = 0; i < snapshots.size(); i++)
    {
        Eigen::MatrixXd F_eigen = Foam2Eigen::field2Eigen(snapshots[i]);

        if (innerProduct == "L2")
        {
            M_vol = ITHACAutilities::get_mass_matrix_FV(snapshots[i]);
        }
        else if (innerProduct == "Frobenius")
        {
            M_vol =  Eigen::MatrixXd::Identity(F_eigen.rows(), F_eigen.rows());
        }
        else
        {
            std::cout << "Inner product not defined" << endl;
            exit(0);
        }

        M = EigenModes[0].transpose() * M_vol * EigenModes[0];
        projSnapI = EigenModes[0].transpose() * M_vol * F_eigen;
        projSnapCoeff = M.fullPivLu().solve(projSnapI);
        projSnap.append(reconstruct(projSnapCoeff, "projSnap"));
    }

    return projSnap;
}
