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

/// \file
/// source file for the ITHACADMD class

#include "ITHACADMD.H"

template<class Type, template<class> class PatchField, class GeoMesh>
ITHACADMD<Type, PatchField, GeoMesh>::ITHACADMD(
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& snapshots, double dt)
    :
    snapshotsDMD(snapshots),
    NSnaps(snapshots.size()),
    originalDT(dt)
{
    ITHACAparameters* para(ITHACAparameters::getInstance());
    redSVD = para->ITHACAdict->lookupOrDefault<bool>("redSVD", false);
}


template<class Type, template<class> class PatchField, class GeoMesh>
void ITHACADMD<Type, PatchField, GeoMesh>::getModes(label SVD_rank, bool exact,
        bool exportDMDmodes)
{
    // Check the rank if rank < 0, full rank.
    if (SVD_rank < 0)
    {
        SVD_rank = NSnaps - 1;
    }

    std::string assertMessage = "The SVD_rank is equal to: " + name(
                                    SVD_rank) +
                                ", it cannot be bigger than the number of snapshots - 1. NSnapshots is equal to: "
                                + name(NSnaps);
    SVD_rank_public = SVD_rank;
    M_Assert(SVD_rank < NSnaps, assertMessage.c_str());
    // Convert the OpenFoam Snapshots to Matrix
    Eigen::MatrixXd SnapEigen = Foam2Eigen::PtrList2Eigen(snapshotsDMD);
    List<Eigen::MatrixXd> SnapEigenBC = Foam2Eigen::PtrList2EigenBC(snapshotsDMD);
    Eigen::MatrixXd Xm = SnapEigen.leftCols(NSnaps - 1);
    Eigen::MatrixXd Ym = SnapEigen.rightCols(NSnaps - 1);
    List<Eigen::MatrixXd> XmBC(SnapEigenBC.size());
    List<Eigen::MatrixXd> YmBC(SnapEigenBC.size());

    for (label i = 0 ; i < SnapEigenBC.size() ; i++)
    {
        XmBC[i] = SnapEigenBC[i].leftCols(NSnaps - 1);
        YmBC[i] = SnapEigenBC[i].rightCols(NSnaps - 1);
    }

    Eigen::MatrixXcd U;
    Eigen::MatrixXcd V;
    Eigen::VectorXd S;

    if (redSVD)
    {
        Info << "SVD using Randomized method" << endl;
        RedSVD::RedSVD<Eigen::MatrixXd> svd(Xm, SVD_rank);
        U = svd.matrixU();
        V = svd.matrixV();
        S = svd.singularValues().array().cwiseInverse();
    }
    else
    {
        Info << "SVD using Divide and Conquer method" << endl;
        Eigen::BDCSVD<Eigen::MatrixXd> svd(Xm,
                                           Eigen::ComputeThinU | Eigen::ComputeThinV);
        U = svd.matrixU().leftCols(SVD_rank);
        V = svd.matrixV().leftCols(SVD_rank);
        S = svd.singularValues().head(SVD_rank).array().cwiseInverse();
    }

    Eigen::MatrixXcd A_tilde = U.transpose().conjugate() * Ym *
                               V *
                               S.asDiagonal();
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> esEg(A_tilde);
    eigenValues = esEg.eigenvalues();
    Eigen::VectorXd ln = eigenValues.array().log().imag().abs();
    // Sort based on The Frequencies
    typedef std::pair<double, label> mypair;
    std::vector<mypair> sortedList(ln.size());

    for (label i = 0; i < ln.size(); i++)
    {
        sortedList[i].first = ln(i);
        sortedList[i].second = i;
    }

    std::sort(sortedList.begin(), sortedList.end(),
              std::less<std::pair<double, label>>());

    if (exact)
    {
        DMDEigenModes = Ym * V * S.asDiagonal() * esEg.eigenvectors();
        DMDEigenModesBC.resize(YmBC.size());

        for (label i = 0; i < YmBC.size(); i++)
        {
            DMDEigenModesBC[i] = YmBC[i] * V * S.asDiagonal() * esEg.eigenvectors();
        }
    }
    else
    {
        Eigen::VectorXd eigenValueseigLam =
            S.array().sqrt();
        PODm = (Xm * V) * eigenValueseigLam.asDiagonal();
        PODmBC.resize(XmBC.size());

        for (label i = 0; i < XmBC.size(); i++)
        {
            PODmBC[i] = (XmBC[i] * V) * eigenValueseigLam.asDiagonal();
        }

        DMDEigenModes = PODm * esEg.eigenvectors();
        DMDEigenModesBC.resize(PODmBC.size());

        for (label i = 0; i < PODmBC.size(); i++)
        {
            DMDEigenModesBC[i] = PODmBC[i] * esEg.eigenvectors();
        }
    }

    Amplitudes = DMDEigenModes.real().bdcSvd(Eigen::ComputeThinU |
                 Eigen::ComputeThinV).solve(Xm.col(0));

    if (exportDMDmodes)
    {
        convert2Foam();
        Info << "exporting the DMDmodes for " << snapshotsDMD[0].name() << endl;
        ITHACAstream::exportFields(DMDmodesReal.toPtrList(), "ITHACAoutput/DMD/",
                                   snapshotsDMD[0].name() + "_Modes_" + name(SVD_rank) + "_Real");
        ITHACAstream::exportFields(DMDmodesImag.toPtrList(), "ITHACAoutput/DMD/",
                                   snapshotsDMD[0].name() + "_Modes_" + name(SVD_rank) + "_Imag");
    }
}

template<class Type, template<class> class PatchField, class GeoMesh>
void ITHACADMD<Type, PatchField, GeoMesh>::getDynamics(double tStart,
        double tFinal, double dt)
{
    Eigen::VectorXcd omega = eigenValues.array().log() / originalDT;
    label ncols = static_cast<label>((tFinal - tStart) / dt ) + 1;
    dynamics.resize(SVD_rank_public, ncols);
    label i = 0;

    for (double t = tStart; t <= tFinal; t += dt)
    {
        Eigen::VectorXcd coli = (omega * t).array().exp() * Amplitudes.array();
        dynamics.col(i) = coli;
        i++;
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void ITHACADMD<Type, PatchField, GeoMesh>::exportEigs(word exportFolder)
{
    Eigen::MatrixXcd eigs = eigenValues;
    mkDir(exportFolder);
    std::string path = exportFolder + "/eigs.npy";
    cnpy::save(eigs, path);
    return;
}


template<class Type, template<class> class PatchField, class GeoMesh>
void ITHACADMD<Type, PatchField, GeoMesh>::reconstruct(word exportFolder,
        word fieldName)
{
    PtrList<GeometricField<Type, PatchField, GeoMesh>> snapshotsrec;

    for (label i = 0; i < dynamics.cols(); i++)
    {
        Eigen::MatrixXcd col = DMDEigenModes * dynamics.col(i);
        Eigen::VectorXd vec = col.real();
        GeometricField<Type, PatchField, GeoMesh> tmp2("TMP",
                snapshotsDMD[0] * 0);
        tmp2 = Foam2Eigen::Eigen2field(tmp2, vec);

        for (label k = 0; k < tmp2.boundaryField().size(); k++)
        {
            Eigen::VectorXd vecBC = (DMDEigenModesBC[k] * dynamics.col(i)).real();
            ITHACAutilities::assignBC(tmp2, k, vecBC);
        }

        snapshotsrec.append(tmp<GeometricField<Type, PatchField, GeoMesh>>(tmp2));
    }

    ITHACAstream::exportFields(snapshotsrec, exportFolder, fieldName);
}
template<class Type, template<class> class PatchField, class GeoMesh>
void ITHACADMD<Type, PatchField, GeoMesh>::convert2Foam()
{
    DMDmodesReal.resize(DMDEigenModes.cols());
    DMDmodesImag.resize(DMDEigenModes.cols());
    GeometricField<Type, PatchField, GeoMesh> tmp2Real(
        snapshotsDMD[0].name(), snapshotsDMD[0] * 0);
    GeometricField<Type, PatchField, GeoMesh> tmp2Imag(
        snapshotsDMD[0].name(), snapshotsDMD[0] * 0);

    for (label i = 0; i < DMDmodesReal.size(); i++)
    {
        Eigen::VectorXd vecReal = DMDEigenModes.col(i).real();
        Eigen::VectorXd vecImag = DMDEigenModes.col(i).imag();
        tmp2Real = Foam2Eigen::Eigen2field(tmp2Real, vecReal);
        tmp2Imag = Foam2Eigen::Eigen2field(tmp2Imag, vecImag);

        // Adjusting boundary conditions
        for (label k = 0; k < tmp2Real.boundaryField().size(); k++)
        {
            ITHACAutilities::assignBC(tmp2Real, k, DMDEigenModesBC[k].col(i).real());
            ITHACAutilities::assignBC(tmp2Imag, k, DMDEigenModesBC[k].col(i).imag());
        }

        DMDmodesReal.set(i, tmp<GeometricField<Type, PatchField, GeoMesh>>(tmp2Real));
        DMDmodesImag.set(i, tmp<GeometricField<Type, PatchField, GeoMesh>>(tmp2Imag));
    }
}
template class ITHACADMD<scalar, fvPatchField, volMesh>;
template class ITHACADMD<vector, fvPatchField, volMesh>;
