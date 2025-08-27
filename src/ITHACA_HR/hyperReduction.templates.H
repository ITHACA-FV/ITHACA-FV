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
    it under the terms of the GN3U Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    ITHACA-FV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with ITHACA-FV. If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#ifndef hyperReduction_templates_H
#define hyperReduction_templates_H

#include "hyperReduction.H"

// TODO increase modularity to avoid heavy compiling costs

// Template function constructor
template <typename... SnapshotsLists>
HyperReduction<SnapshotsLists...>::HyperReduction(label n_modes,
        label n_nodes,
        Eigen::VectorXi initialSeeds,
        word problemName,
        SnapshotsLists &&...snapshotsLists)
    : vectorial_dim{compute_vectorial_dim(snapshotsLists...)},
      n_modes{n_modes},
      n_nodes{n_nodes},
      initialSeeds{initialSeeds},
      problemName{problemName},
      snapshotsListTuple{std::forward_as_tuple(snapshotsLists...)}
{
    Info << "Init HyperReduction class with vectorial dim: " << vectorial_dim <<
         endl;
    para = ITHACAparameters::getInstance();
    folderProblem = "ITHACAoutput/" + problemName;
    mkDir(folderProblem);
    // TODO check snapshotsLists is not empty
    // TODO check that the first snapshotsList is not empty (includes above)
    // TODO initialize vectorial_dim only from SnapshotsLists type
    // vectorial_dim = std::apply(compute_vectorial_dim<SnapshotsLists&&...>, snapshotsListTuple);
    // Get fields' names
    std::apply([this](auto & ...snapList)
    {
        (..., stackNames(snapList));
    }, snapshotsListTuple);

    // Get fields' dimensions
    std::apply([this](auto & ...snapList)
    {
        (..., stackDimensions(snapList));
    }, snapshotsListTuple);

    sumFieldsDim = 0;
    std::apply([this](auto & ...snapList)
    {
        (..., sumDimensions(this->sumFieldsDim, snapList));
    }, snapshotsListTuple);

    Info << "Fields names: ";

    for (unsigned int ith_field = 0; ith_field < fieldNames.size(); ith_field++)
    {
        Info << fieldNames[ith_field] << " (dim=" << fieldDims[ith_field] << "); ";
    }

    Info << endl;
    n_snapshots = std::get<0>(snapshotsListTuple).size();
    Info << "The number of snapshots is: " << n_snapshots << endl;
    n_cells = std::get<0>(snapshotsListTuple)[0].size();
    Info << "The number of cells is: " << n_cells << endl;
    Info << "Initial seeds length: " << initialSeeds.rows() << endl;
    // get boundaries info
    n_boundary_patches = std::get<0>(snapshotsListTuple)[0].boundaryField().size();
    n_boundary_cells_list.resize(n_boundary_patches);
    n_boundary_cells = 0;

    for (unsigned int ith_boundary_patch = 0;
            ith_boundary_patch < n_boundary_patches; ith_boundary_patch++)
    {
        n_boundary_cells_list[ith_boundary_patch] = std::get<0>
            (snapshotsListTuple)[0].boundaryField()[ith_boundary_patch].size();
        n_boundary_cells = n_boundary_cells + n_boundary_cells_list[ith_boundary_patch];
    }

    Info << "Number of boundary cells: " << n_boundary_cells << endl;
}


// Template function constructor
template <typename... SnapshotsLists>
HyperReduction<SnapshotsLists...>::HyperReduction(label n_modes,
        label n_nodes,
        unsigned int vectorialDim,
        label n_cells,
        Eigen::VectorXi initialSeeds,
        word problemName)
    : n_modes{n_modes},
      n_nodes{n_nodes},
      vectorial_dim{vectorialDim},
      n_cells{n_cells},
      initialSeeds{initialSeeds},
      problemName{problemName}
{
    Info << "Init HyperReduction class with vectorial dim: " << vectorial_dim <<
         endl;
    para = ITHACAparameters::getInstance();
    folderProblem = problemName;
    mkDir(folderProblem);
    Info << "Initial seeds length: " << initialSeeds.rows() << endl;
}

template <typename... SnapshotsLists>
void HyperReduction<SnapshotsLists...>::getModesSVD(
    SnapshotsListTuple& snapshotsListTuple, Eigen::MatrixXd& modesSVD,
    Eigen::VectorXd& fieldWeights, bool saveModesFlag)
{
    // TODO move inside ITHACAPOD
    ITHACAparameters* para(ITHACAparameters::getInstance());
    word folderSVD = "ITHACAoutput/" + problemName + "/ModesSVD/";

    if (!ITHACAutilities::check_folder(folderSVD))
    {
        mkDir(folderSVD);
        // Initialize stacked snapshotsMatrix and normalizing weights
        Eigen::MatrixXd snapMatrix;
        getSnapMatrix(snapMatrix, fieldWeights);
        word dimensionReductionMethod =
            para->ITHACAdict->lookupOrDefault<word>("DimensionReductionMethod", "SVD");
        Eigen::MatrixXd eigenVectoreig;

        if (dimensionReductionMethod == "RSVD")
        {
            Info << "####### Performing RSVD for " << problemName << " #######" << endl;
            unsigned int r_modeSVD =
                para->ITHACAdict->lookupOrDefault<unsigned int>("RSVDdim", snapMatrix.cols());
            RedSVD::RedSVD svd(fieldWeights.array().cwiseInverse().matrix().asDiagonal() *
                               snapMatrix, r_modeSVD);
            Info << "####### End of RSVD for " << problemName << " #######" <<
                 endl;
            eigenValueseig = svd.singularValues().real();
            eigenVectoreig = svd.matrixU().real();
        }
        else
        {
            Info << "####### Performing SVD for " << problemName << " #######" << endl;
            Eigen::JacobiSVD<Eigen::MatrixXd> svd(
                fieldWeights.array().cwiseInverse().matrix().asDiagonal() * snapMatrix,
                Eigen::ComputeThinU | Eigen::ComputeThinV);
            Info << "####### End of SVD for " << problemName << " #######" <<
                 endl;
            eigenValueseig = svd.singularValues().real();
            eigenVectoreig = svd.matrixU().real();
        }

        modesSVD = eigenVectoreig;
        // TODO correct boundary conditions
        cnpy::save(eigenValueseig,
                   "ITHACAoutput/" + problemName + "/ModesSVD/eigenvalues.npy");
        eigenValueseig = eigenValueseig / eigenValueseig.sum();
        Eigen::VectorXd cumEigenValues(eigenValueseig);

        for (label j = 1; j < cumEigenValues.size(); ++j)
        {
            cumEigenValues(j) += cumEigenValues(j - 1);
        }

        Eigen::saveMarketVector(eigenValueseig,
                                "ITHACAoutput/" + problemName + "/ModesSVD/Eigenvalues_" + problemName,
                                para->precision,
                                para->outytpe);
        Eigen::saveMarketVector(cumEigenValues,
                                "ITHACAoutput/" + problemName + "/ModesSVD/CumEigenvalues_" + problemName,
                                para->precision,
                                para->outytpe);
        cnpy::save(modesSVD, "ITHACAoutput/" + problemName + "/ModesSVD/modes.npy");
        cnpy::save(fieldWeights,
                   "ITHACAoutput/" + problemName + "/ModesSVD/normalizingWeights.npy");

        if (saveModesFlag)
        {
            unsigned int rowIndex{0};

            for (unsigned int modeIndex = 0; modeIndex < modesSVD.cols(); modeIndex++)
            {
                std::apply([this, & modesSVD, & rowIndex, & modeIndex,
                            & folderSVD](auto & ...snapList)
                {
                    (..., saveModes(snapList, modesSVD, rowIndex, modeIndex, folderSVD));
                }, snapshotsListTuple);

                rowIndex = 0;
            }
        }
    }
    else
    {
        Info << "Reading the existing modes" << endl;
        cnpy::load(modesSVD, folderSVD + "/modes.npy");
        cnpy::load(fieldWeights, folderSVD + "/normalizingWeights.npy");
    }
}

template <typename... SnapshotsLists>
void HyperReduction<SnapshotsLists...>::getModesSVD(
    SnapshotsListTuple& snapshotsListTuple, Eigen::MatrixXd& modesSVD,
    Eigen::VectorXd& fieldWeights, Eigen::MatrixXd& modesSVDBoundary,
    Eigen::VectorXd& fieldWeightsBoundary, bool saveModesFlag)
{
    // TODO move inside ITHACAPOD
    ITHACAparameters* para(ITHACAparameters::getInstance());
    word folderSVD = "ITHACAoutput/" + problemName + "/ModesSVD/";

    if (!ITHACAutilities::check_folder(folderSVD))
    {
        mkDir(folderSVD);
        // Initialize stacked snapshotsMatrix and normalizing weights
        Eigen::MatrixXd snapMatrix;
        List<Eigen::MatrixXd> snapMatrixBoundary;
        List<Eigen::VectorXd> fieldWeightsBoundaryList;
        getSnapMatrix(snapMatrix, fieldWeights, snapMatrixBoundary,
                      fieldWeightsBoundaryList);

        for (unsigned int id = 0; id < n_boundary_patches; id++)
        {
            snapMatrix.conservativeResize(snapMatrix.rows() + snapMatrixBoundary[id].rows(),
                                          snapMatrix.cols());
            snapMatrix.bottomRows(snapMatrixBoundary[id].rows()) = snapMatrixBoundary[id];
            fieldWeightsBoundary.conservativeResize(fieldWeightsBoundary.rows() +
                                                    snapMatrixBoundary[id].rows());
            fieldWeightsBoundary.tail(fieldWeightsBoundaryList[id].rows()) =
                                    fieldWeightsBoundaryList[id];
        }

        cnpy::save(fieldWeights,
                   "ITHACAoutput/" + problemName + "/ModesSVD/normalizingWeights.npy");
        cnpy::save(fieldWeightsBoundary,
                   "ITHACAoutput/" + problemName + "/ModesSVD/normalizingWeightsBoundary.npy");
        fieldWeights.conservativeResize(fieldWeights.rows() +
                                        fieldWeightsBoundary.rows());
        fieldWeights.tail(fieldWeightsBoundary.rows()) = fieldWeightsBoundary;
        word dimensionReductionMethod =
            para->ITHACAdict->lookupOrDefault<word>("DimensionReductionMethod", "SVD");
        Eigen::MatrixXd eigenVectoreig;

        if (dimensionReductionMethod == "RSVD")
        {
            Info << "####### Performing RSVD for " << problemName << " #######" << endl;
            unsigned int r_modeSVD =
                para->ITHACAdict->lookupOrDefault<unsigned int>("RSVDdim", snapMatrix.cols());
            RedSVD::RedSVD svd(fieldWeights.array().cwiseInverse().matrix().asDiagonal() *
                               snapMatrix, r_modeSVD);
            Info << "####### End of RSVD for " << problemName << " #######" << endl;
            eigenValueseig = svd.singularValues().real();
            eigenVectoreig = svd.matrixU().real();
        }
        else
        {
            Info << "####### Performing SVD for " << problemName << " #######" << endl;
            Eigen::JacobiSVD<Eigen::MatrixXd> svd(
                fieldWeights.array().cwiseInverse().matrix().asDiagonal() * snapMatrix,
                Eigen::ComputeThinU | Eigen::ComputeThinV);
            Info << "####### End of SVD for " << problemName << " #######" <<
                 endl;
            eigenValueseig = svd.singularValues().real();
            eigenVectoreig = svd.matrixU().real();
        }

        modesSVD = eigenVectoreig.bottomRows(vectorial_dim * n_cells);
        modesSVDBoundary = eigenVectoreig.topRows(n_boundary_cells);
        fieldWeights.conservativeResize(fieldWeights.rows() -
                                        fieldWeightsBoundary.rows());
        cnpy::save(eigenValueseig,
                   "ITHACAoutput/" + problemName + "/ModesSVD/eigenvalues.npy");
        eigenValueseig = eigenValueseig / eigenValueseig.sum();
        Eigen::VectorXd cumEigenValues(eigenValueseig);

        for (label j = 1; j < cumEigenValues.size(); ++j)
        {
            cumEigenValues(j) += cumEigenValues(j - 1);
        }

        Eigen::saveMarketVector(eigenValueseig,
                                "ITHACAoutput/" + problemName + "/ModesSVD/Eigenvalues_" + problemName,
                                para->precision,
                                para->outytpe);
        Eigen::saveMarketVector(cumEigenValues,
                                "ITHACAoutput/" + problemName + "/ModesSVD/CumEigenvalues_" + problemName,
                                para->precision,
                                para->outytpe);
        cnpy::save(modesSVD, "ITHACAoutput/" + problemName + "/ModesSVD/modes.npy");
        cnpy::save(modesSVDBoundary,
                   "ITHACAoutput/" + problemName + "/ModesSVD/modesBoundary.npy");

        if (saveModesFlag)
        {
            unsigned int rowIndex{0};
            unsigned int rowIndexBoundary{0};

            for (unsigned int modeIndex = 0; modeIndex < modesSVD.cols(); modeIndex++)
            {
                std::apply([this, & modesSVD, & modesSVDBoundary, & rowIndex, & modeIndex,
                            & rowIndexBoundary, & folderSVD](auto & ...snapList)
                {
                    (..., saveModes(snapList, modesSVD, modesSVDBoundary, rowIndex,
                                    rowIndexBoundary, modeIndex, folderSVD));
                }, snapshotsListTuple);

                rowIndex = 0;
            }
        }
    }
    else
    {
        Info << "Reading the existing modes" << endl;
        cnpy::load(modesSVD, folderSVD + "/modes.npy");
        cnpy::load(modesSVDBoundary,
                   "ITHACAoutput/" + problemName + "/ModesSVD/modesBoundary.npy");
        cnpy::load(fieldWeights, folderSVD + "/normalizingWeights.npy");
        cnpy::load(fieldWeightsBoundary,
                   "ITHACAoutput/" + problemName + "/ModesSVD/normalizingWeightsBoundary.npy");
    }
}

template <typename... SnapshotsLists>
void HyperReduction<SnapshotsLists...>::getSnapMatrix(Eigen::MatrixXd&
        snapMatrix, Eigen::VectorXd& fieldWeights)
{
    std::apply([this, & fieldWeights, & snapMatrix](auto & ...snapList)
    {
        (..., stackSnapshots(snapList, snapMatrix, fieldWeights));
    }, snapshotsListTuple);
}

template <typename... SnapshotsLists>
void HyperReduction<SnapshotsLists...>::getSnapMatrix(Eigen::MatrixXd&
        snapMatrix, Eigen::VectorXd& fieldWeights,
        List<Eigen::MatrixXd>& snapMatrixBoundary,
        List<Eigen::VectorXd>& fieldWeightsBoundary)
{
    fieldWeightsBoundary.resize(n_boundary_patches);
    snapMatrixBoundary.resize(n_boundary_patches);
    std::apply([this, & fieldWeights, & snapMatrix, & fieldWeightsBoundary,
                & snapMatrixBoundary](auto & ...snapList)
    {
        (..., stackSnapshots(snapList, snapMatrix, fieldWeights));
        (..., stackSnapshotsBoundary(snapList, snapMatrixBoundary,
                                     fieldWeightsBoundary));
    },

    snapshotsListTuple);
}

template <typename... SnapshotsLists>
template <typename SnapshotsList>
void HyperReduction<SnapshotsLists...>::saveModes(SnapshotsList sList,
        Eigen::MatrixXd& snapshotsMatrix, unsigned int& rowIndex,
        unsigned int& modeIndex, word folder)
{
    unsigned int field_dim = get_field_dim<typename SnapshotsList::value_type>();
    auto fieldName = sList[0].name();
    unsigned int fieldSize = field_dim * sList[0].size();
    Eigen::VectorXd fieldBlock = snapshotsMatrix.block(rowIndex, modeIndex,
                                 fieldSize, 1);
    auto fieldOut = Foam2Eigen::Eigen2field(sList[0], fieldBlock, true);
    rowIndex += fieldSize;
    ITHACAstream::exportSolution(fieldOut, name(modeIndex), folder);
}

template <typename... SnapshotsLists>
template <typename SnapshotsList>
void HyperReduction<SnapshotsLists...>::saveModes(SnapshotsList sList,
        Eigen::MatrixXd& snapshotsMatrix, Eigen::MatrixXd& snapshotsMatrixBoundary,
        unsigned int& rowIndex,  unsigned int& rowIndexBoundary,
        unsigned int& modeIndex, word folder)
{
    unsigned int field_dim = get_field_dim<typename SnapshotsList::value_type>();
    auto fieldName = sList[0].name();
    unsigned int fieldSize = field_dim * sList[0].size();
    Eigen::VectorXd fieldBlock = snapshotsMatrix.block(rowIndex, modeIndex,
                                 fieldSize, 1);
    List<Eigen::VectorXd> fieldBlockBoundary;
    fieldBlockBoundary.resize(n_boundary_patches);

    for (unsigned int id = 0; id < n_boundary_patches; id++)
    {
        unsigned int bfieldDim = field_dim * int(n_boundary_cells_list[id]);
        fieldBlockBoundary[id] = snapshotsMatrixBoundary.block(rowIndexBoundary,
                                 modeIndex, bfieldDim, 1);
        rowIndexBoundary += bfieldDim;
    }

    auto fieldOut = Foam2Eigen::Eigen2field(sList[0], fieldBlock,
                                            fieldBlockBoundary);
    rowIndex += fieldSize;
    ITHACAstream::exportSolution(fieldOut, name(modeIndex), folder);
}

template <typename... SnapshotsLists>
void HyperReduction<SnapshotsLists...>::offlineGappyDEIM(
    Eigen::MatrixXd& snapshotsModes, Eigen::VectorXd& weights,
    word folderMethodName)
{
    folderMethod = folderMethodName;
    Info << "FolderMethod : " << folderMethod << endl;
    mkDir(folderMethod);
    normalizingWeights = weights;
    nodePoints = autoPtr<IOList<label >> (new IOList<label>(
            IOobject("nodePoints", para->runTime.time().constant(), "../" + folderMethod,
                     para->mesh, IOobject::READ_IF_PRESENT, IOobject::NO_WRITE)));
    auto greedyMetric = para->ITHACAdict->lookupOrDefault<word>("GreedyMetric",
                        "L2");
    bool offlineStage = !nodePoints().headerOk();

    if (offlineStage)
    {
        assert(n_modes > 0);
        assert(n_nodes >= n_modes);
        Eigen::VectorXd mp_not_mask = Eigen::VectorXd::Constant(n_cells * vectorial_dim,
                                      1);
        std::set<label> nodePointsSet;

        // set initialSeeds
        if (initialSeeds.rows() > 0)
        {
            initSeeds(mp_not_mask, nodePointsSet);
        }

        int na = n_nodes - initialSeeds.rows();

        if (na > 0)
        {
            Eigen::SparseMatrix<double> reshapeMat;
            initReshapeMat(reshapeMat);

            if (greedyMetric == "L2")
            {
                Info << "####### Begin Greedy-L2 GappyDEIM #######" << endl;
                Info << "####### Modes=" <<  n_modes
                     << ", nodePoints=" << n_nodes << " #######"
                     << endl;
                int nb = 0;
                int nit = std::min(n_modes, na);
                int ncimin = std::floor(n_modes / nit);
                int naimin = std::floor(na / n_modes);
                Eigen::MatrixXd A;
                Eigen::VectorXd b;
                Eigen::VectorXd c;
                Eigen::VectorXd r;
                label ind_max, c1;
                double max;

                for (int i = 1; i <= nit; i++)
                {
                    int nci = ncimin;

                    // add the remaining modes in the quotient n_modes / nite
                    if (i <= n_modes % nit)
                    {
                        nci = ncimin + 1;
                    }

                    int nai = naimin;

                    // add the remaining nodePoints in the quotient na / n_modes
                    if (i <= na % n_modes)
                    {
                        nai = naimin + 1;
                    }

                    Eigen::MatrixXd V;

                    // select basis
                    if (i == 1)
                    {
                        V = snapshotsModes.leftCols(nci);
                    }
                    else
                    {
                        for (int q = 1; q <= nci; q++)
                        {
                            A = P.transpose() * basisMatrix;
                            b = P.transpose() * snapshotsModes.col(nb + q - 1);
                            c = A.fullPivLu().solve(b);
                            r = snapshotsModes.col(nb + q - 1) - basisMatrix * c;
                            V.conservativeResize(snapshotsModes.rows(), q);
                            V.col(q - 1) = r;
                        }
                    }

                    for (int j = 1; j <= nai; j++)
                    {
                        if (P.cols() > 0)
                        {
                            max = (reshapeMat * (mp_not_mask.asDiagonal() * V)
                                   .rowwise()
                                   .lpNorm<2>().array().square().matrix())
                                  .maxCoeff(& ind_max, & c1);
                        }
                        else
                        {
                            max = (reshapeMat * V.rowwise().lpNorm<2>().array().square().matrix()).maxCoeff(
                                      & ind_max, & c1);
                        }

                        updateNodes(P, ind_max, mp_not_mask);
                    }

                    nb += nci;
                    basisMatrix = snapshotsModes.leftCols(nb);
                }

                Info << "####### End of greedy GappyDEIM #######\n";
            }
            else if (greedyMetric == "SOPT")
            {
                Info << "####### Begin SOPT GappyDEIM #######" << endl;
                Info << "####### Modes=" <<  n_modes
                     << ", nodePoints=" << n_nodes << " #######"
                     << endl;
                label ind_max, c1;
                double max;
                Eigen::MatrixXd V = snapshotsModes.leftCols(n_modes);

                for (unsigned int ith_node = 0; ith_node < na; ith_node++)
                {
                    if (P.cols() > 0)
                    {
                        // TODO make efficient
                        Eigen::MatrixXd tmp = P.transpose() * V;
                        Eigen::MatrixXd R = mp_not_mask.asDiagonal() * V;
                        Eigen::MatrixXd inv = (tmp.transpose() * tmp).fullPivLu().inverse();
                        Eigen::VectorXd num = 1 + (reshapeMat * ((R * inv).array() *
                                                   R.array()).rowwise().sum().matrix()).array();
                        Eigen::VectorXd norma = (tmp.array().square()).colwise().sum();
                        Eigen::MatrixXd adding = reshapeMat * R.array().square().matrix();
                        adding.rowwise() += norma.transpose();
                        Eigen::VectorXd denom = (adding.array().pow(1. /
                                                 V.cols()).matrix()).rowwise().prod();
                        max = ((mp_not_mask.head(n_cells).asDiagonal() * num).array() /
                               denom.array()).maxCoeff(& ind_max, & c1);
                    }
                    else
                    {
                        max = (reshapeMat * V.rowwise().lpNorm<2>().array().square().matrix()).maxCoeff(
                                  & ind_max, & c1);
                    }

                    updateNodes(P, ind_max, mp_not_mask);
                }

                basisMatrix = V;
                Info << "####### End of SOPT GappyDEIM #######\n";
            }
            else if (greedyMetric == "SOPTE")
            {
                Info << "####### Begin SOPT-exact GappyDEIM #######" << endl;
                Info << "####### Modes=" <<  n_modes
                     << ", nodePoints=" << n_nodes << " #######"
                     << endl;
                label ind_max, c1;
                double max;
                Eigen::MatrixXd A;
                Eigen::VectorXd b;
                Eigen::VectorXd c;
                Eigen::VectorXd r;
                Eigen::MatrixXd V = snapshotsModes.col(0);
                max = (reshapeMat * (mp_not_mask.asDiagonal() * snapshotsModes.leftCols(
                                         n_modes)).rowwise().lpNorm<2>().array().square().matrix()).maxCoeff(& ind_max,
                                             & c1);
                updateNodes(P, ind_max, mp_not_mask);

                for (unsigned int ith_node = 1; ith_node < na; ith_node++)
                {
                    if (ith_node < n_modes)
                    {
                        A = P.transpose() * V;
                        b = P.transpose() * snapshotsModes.col(ith_node);
                        c = A.fullPivLu().solve(b);
                        r = snapshotsModes.col(ith_node) - V * c;
                        V.conservativeResize(snapshotsModes.rows(), ith_node + 1);
                        V.col(ith_node) =  snapshotsModes.col(ith_node);
                        max = (reshapeMat * r.rowwise().lpNorm<2>().array().square().matrix()).maxCoeff(
                                  & ind_max, & c1);
                    }
                    else
                    {
                        // TODO optimize
                        Eigen::MatrixXd tmp = P.transpose() * snapshotsModes.leftCols(n_modes);
                        Info << "SOPT: " << s_optimality(tmp) << endl;
                        tmp.conservativeResize(tmp.rows() + vectorial_dim, tmp.cols());
                        Eigen::MatrixXd masked = mp_not_mask.asDiagonal() * snapshotsModes.leftCols(
                                                     n_modes);
                        Eigen::VectorXd results(n_cells);

                        for (unsigned int ith_cell = 0; ith_cell < n_cells; ith_cell++)
                        {
                            for (unsigned int ith_field = 0; ith_field < vectorial_dim; ith_field++)
                            {
                                tmp.row(tmp.rows() - vectorial_dim + ith_field) = masked.row(
                                        ith_cell +  ith_field * n_cells);
                            }

                            results(ith_cell) = s_optimality(tmp);
                        }

                        max = results.maxCoeff(& ind_max, & c1);
                    }

                    updateNodes(P, ind_max, mp_not_mask);
                }

                basisMatrix = snapshotsModes.leftCols(n_modes);
                Info << "####### End of SOPTE GappyDEIM #######\n";
            }
        }
        else
        {
            basisMatrix = snapshotsModes.leftCols(n_modes);
            Info << "####### InitialSeeds equal to n_nodes #######\n";
        }

        evaluatePinv(P, basisMatrix, normalizingWeights);
        renormalizedBasisMatrix = normalizingWeights.asDiagonal() * basisMatrix;
        MatrixOnline = renormalizedBasisMatrix * pinvPU;
        P.makeCompressed();
        cnpy::save(basisMatrix, folderMethod + "/basisMatrix.npy");
        cnpy::save(P, folderMethod + "/projectionMatrix.npz");
        cnpy::save(normalizingWeights, folderMethod + "/normalizingWeights.npy");
        cnpy::save(nodes, folderMethod + "/mp.npy");
        cnpy::save(pinvPU, folderMethod + "/pinvPU.npy");
        cnpy::save(MatrixOnline, folderMethod + "/MatrixOnline.npy");
        nodePoints().write();
        Info << "Projection Matrix shape: " << P.rows() << " " << P.cols() << endl;
        Info << "Basis Matrix shape: " << basisMatrix.rows() << " " <<
             basisMatrix.cols() << endl;
        Info << "Pseudo Inverse Matrix shape: " << pinvPU.rows() << " " << pinvPU.cols()
             << endl;
    }
    else
    {
        Info << "Read GappyDEIM\n";
        cnpy::load(basisMatrix, folderMethod + "/basisMatrix.npy");
        cnpy::load(MatrixOnline, folderMethod + "/MatrixOnline.npy");
        cnpy::load(normalizingWeights, folderMethod + "/normalizingWeights.npy");
        cnpy::load(P, folderMethod + "/projectionMatrix.npz");
        cnpy::load(nodes, folderMethod + "/mp.npy");
        cnpy::load(pinvPU, folderMethod + "/pinvPU.npy");
        renormalizedBasisMatrix = normalizingWeights.asDiagonal() * basisMatrix;
        Info << "Projection Matrix shape: " << P.rows() << " " << P.cols() << endl;
        Info << "Basis Matrix shape: " << basisMatrix.rows() << " " <<
             basisMatrix.cols() << endl;
        Info << "Pseudo Inverse Matrix shape: " << pinvPU.rows() << " " << pinvPU.cols()
             << endl;
    }
}

template <typename... SnapshotsLists>
void HyperReduction<SnapshotsLists...>::offlineECP(Eigen::MatrixXd&
        snapshotsModes, Eigen::VectorXd& weights, word folderMethodName)
{
    folderMethod = folderMethodName;
    Info << "FolderMethod : " << folderMethod << endl;
    mkDir(folderMethod);
    n_nodes = n_modes + 1;
    normalizingWeights = weights;
    Info << "####### Begin ECP #######" << endl;
    Info << "####### Modes=" <<  n_modes
         << ", nodePoints=" << n_nodes << " #######"
         << endl;
    nodePoints = autoPtr<IOList<label >> (new IOList<label>(
            IOobject("nodePoints", para->runTime.time().constant(), "../" + folderMethod,
                     para->mesh, IOobject::READ_IF_PRESENT, IOobject::NO_WRITE)));
    bool offlineStage = !nodePoints().headerOk();

    if (offlineStage)
    {
        assert(n_modes > 0);
        assert(n_nodes >= n_modes);
        // matrices for quadratureWeights evaluation
        Eigen::MatrixXd J;//shape: (vectorial_dim * (n_modes + 1), nodePoints->size())
        Eigen::MatrixXd Jwhole(vectorial_dim * (n_modes + 1), n_cells);
        Eigen::VectorXd q(vectorial_dim * (n_modes + 1));
        // matrices for greedy selection of the nodes
        Eigen::MatrixXd A(vectorial_dim * n_modes, n_cells);
        Eigen::VectorXd b = Eigen::VectorXd::Constant(vectorial_dim * n_modes, 1);
        Eigen::VectorXd volumes = ITHACAutilities::getMassMatrixFV(std::get<0>
                                  (snapshotsListTuple)[0]);
        double volume = volumes.array().sum();

        for (unsigned int ith_field = 0; ith_field < vectorial_dim; ith_field++)
        {
            Eigen::MatrixXd block = snapshotsModes.block(ith_field * n_cells, 0, n_cells,
                                    n_modes).transpose();
            q.segment(ith_field * (n_modes + 1), n_modes) = block.rowwise().sum();
            q(ith_field * (n_modes + 1) + n_modes) = volume;
            Jwhole.middleRows(ith_field * (n_modes + 1), n_modes) = block;
            Jwhole.row(ith_field * (n_modes + 1) + n_modes) = Eigen::VectorXd::Constant(
                      n_cells, 1);
            Eigen::VectorXd mean = block.rowwise().mean();
            block.colwise() -= mean;
            Eigen::VectorXd Anorm = block.colwise().lpNorm<2>();

            if (Anorm.minCoeff() > 0)
            {
                block.array().rowwise() /= Anorm.transpose().array();
            }

            A.middleRows(ith_field * n_modes, n_modes) = block;
        }

        Eigen::VectorXd mp_not_mask = Eigen::VectorXd::Constant(n_cells * vectorial_dim,
                                      1);
        std::set<label> nodePointsSet;

        // set initialSeeds
        if (initialSeeds.rows() > 0)
        {
            initialSeeds(mp_not_mask, nodePointsSet);
            computeLS(J, Jwhole, b, q);
        }

        int na = n_nodes - initialSeeds.rows();
        Eigen::SparseMatrix<double> reshapeMat;
        initReshapeMat(reshapeMat);

        if (na > 0)
        {
            label ind_max, r1;
            double max;

            for (unsigned int ith_node = 0; ith_node < na; ith_node++)
            {
                max = ((b.transpose() * A * mp_not_mask.head(
                            n_cells).asDiagonal()).colwise().sum()).maxCoeff(& r1, & ind_max);
                updateNodes(P, ind_max, mp_not_mask);
                computeLS(J, Jwhole, b, q);
            }
        }

        Info << "####### End ECP #######" << endl;
        basisMatrix = snapshotsModes.leftCols(n_modes);
        evaluateWPU(P, basisMatrix, normalizingWeights, quadratureWeights);
        cnpy::save(basisMatrix, folderMethod + "/basisMatrix.npy");
        cnpy::save(quadratureWeights, folderMethod + "/quadratureWeights.npy");
        cnpy::save(normalizingWeights, folderMethod + "/normalizingWeights.npy");
        cnpy::save(nodes, folderMethod + "/mp.npy");
        cnpy::save(wPU, folderMethod + "/wPU.npy");
        nodePoints().write();
    }
    else
    {
        Info << "Read ECP\n";
        cnpy::load(basisMatrix, folderMethod + "/basisMatrix.npy");
        cnpy::load(quadratureWeights, folderMethod + "/quadratureWeights.npy");
        cnpy::load(normalizingWeights, folderMethod + "/normalizingWeights.npy");
        cnpy::load(nodes, folderMethod + "/mp.npy");
        cnpy::load(wPU, folderMethod + "/wPU.npy");
    }
}

template<typename... SnapshotsLists>
void HyperReduction<SnapshotsLists...>::initSeeds(Eigen::VectorXd mp_not_mask,
        std::set<label> nodePointsSet)
{
    P.resize(n_cells * vectorial_dim, initialSeeds.rows() * vectorial_dim);
    P.reserve(Eigen::VectorXi::Constant(initialSeeds.rows() *
                                        vectorial_dim /*n_cols*/, 1 /*n_non_zero_elements*/));
    unsigned int trueSeeds{0};

    for (int i = 0; i < initialSeeds.rows(); i++)
    {
        unsigned int index = initialSeeds(i) % n_cells;

        // check that there are not repeated nodes in initialSeeds
        if (nodePointsSet.find(index) == nodePointsSet.end())
        {
            nodePointsSet.insert(index);
            nodePoints->append(index);

            for (unsigned int ith_field = 0; ith_field < vectorial_dim; ith_field++)
            {
                P.insert(index + ith_field * n_cells,
                         vectorial_dim * trueSeeds + ith_field) = 1;
                mp_not_mask(index + ith_field * n_cells) = 0;
            }

            trueSeeds++;
        }
    }

    P.conservativeResize(P.rows(), nodePoints->size() * vectorial_dim);
    P.makeCompressed();
    M_Assert(nodePoints->size() <= n_nodes,
             "Size of 'initialSeeds' is greater than 'n_nodes'");
}

template <typename... SnapshotsLists>
void HyperReduction<SnapshotsLists...>::updateNodes(Eigen::SparseMatrix<double>
        & P, label& ind, Eigen::VectorXd& mp_not_mask)
{
    nodePoints->append(ind);
    nodes.conservativeResize(nodes.rows() + 1);
    nodes(nodes.rows() - 1) = ind;
    unsigned int last_col{0};

    if (P.rows() == 0)
    {
        P.resize(n_cells * vectorial_dim, vectorial_dim);
    }
    else
    {
        last_col = P.cols();
        P.resize(P.rows(), P.cols() + vectorial_dim);
    }

    int step = P.cols() / vectorial_dim - 1;

    for (int ith_node = 0; ith_node <= step; ith_node++)
    {
        for (unsigned int ith_field = 0; ith_field < vectorial_dim; ith_field++)
        {
            P.insert(nodes[ith_node] + ith_field * n_cells,
                     ith_node + ith_field * (step + 1)) = 1;
        }
    }

    for (unsigned int ith_field = 0; ith_field < vectorial_dim; ith_field++)
    {
        mp_not_mask(ind + ith_field * n_cells) = 0;
    }
}

template<typename... SnapshotsLists>
void HyperReduction<SnapshotsLists...>::computeLS(Eigen::MatrixXd& J,
        Eigen::MatrixXd& Jwhole, Eigen::VectorXd& b, Eigen::VectorXd& q)
{
    // TODO fix with only positive weights (NNLS)
    J.resize(vectorial_dim * (n_modes + 1), nodePoints->size());
    J = Jwhole(Eigen::indexing::all, nodes);
    quadratureWeights.resize(nodePoints->size() * vectorial_dim);

    for (unsigned int ith_field = 0; ith_field < vectorial_dim; ith_field++)
    {
        quadratureWeights.segment(ith_field * nodePoints->size(),
                                  nodePoints->size()) = J.middleRows(ith_field * (n_modes + 1),
                                      n_modes + 1).colPivHouseholderQr().solve(q.segment(ith_field * (n_modes + 1),
                                          n_modes + 1));
        b.segment(ith_field * n_modes, n_modes) = q.segment(ith_field * (n_modes + 1),
            n_modes) - J.middleRows(ith_field * (n_modes + 1),
                                    n_modes) * quadratureWeights.segment(ith_field * nodePoints->size(),
                                        nodePoints->size());
    }

    b = b / b.lpNorm<2>();
}

template<typename... SnapshotsLists>
void HyperReduction<SnapshotsLists...>::initReshapeMat(
    Eigen::SparseMatrix<double>& reshapeMat)
{
    reshapeMat.resize(n_cells, n_cells * vectorial_dim);
    reshapeMat.reserve(Eigen::VectorXi::Constant(n_cells * vectorial_dim /*n_cols*/,
                       1 /*n_non_zero_elements*/));

    for (unsigned int i = 0; i < n_cells; i++)
    {
        for (unsigned int ith_field = 0; ith_field < vectorial_dim; ith_field++)
        {
            reshapeMat.insert(i, i + n_cells * ith_field) = 1;
        }
    }

    reshapeMat.makeCompressed();
}

template <typename... SnapshotsLists>
template <typename SnapshotsList>
void HyperReduction<SnapshotsLists...>::stackSnapshots(SnapshotsList sList,
        Eigen::MatrixXd& snapshotsMatrix, Eigen::VectorXd& fieldWeights)
{
    unsigned int field_dim = get_field_dim<typename SnapshotsList::value_type>();
    Eigen::MatrixXd tmpSnapshots = Foam2Eigen::PtrList2Eigen(sList);
    // get volumes
    Eigen::VectorXd V = ITHACAutilities::getMassMatrixFV(sList[0]);
    double maxVal = std::sqrt(tmpSnapshots.colwise().lpNorm<2>().maxCoeff());
    fieldWeights.conservativeResize(fieldWeights.rows() + field_dim * n_cells);
    fieldWeights.tail(n_cells * field_dim) = V.array().sqrt().cwiseInverse() *
                maxVal;
    snapshotsMatrix.conservativeResize(snapshotsMatrix.rows() + n_cells * field_dim,
                                       n_snapshots);
    snapshotsMatrix.bottomRows(n_cells * field_dim) = tmpSnapshots;
}

template <typename... SnapshotsLists>
template <typename SnapshotsList>
void HyperReduction<SnapshotsLists...>::stackSnapshotsBoundary(
    SnapshotsList sList, List<Eigen::MatrixXd>& snapshotsMatrixBoundary,
    List<Eigen::VectorXd>& fieldWeightsBoundary)
{
    unsigned int field_dim = get_field_dim<typename SnapshotsList::value_type>();
    List<Eigen::MatrixXd> tmpBoundarySnapshots = Foam2Eigen::PtrList2EigenBC(sList);
    List<double> maxVal;
    maxVal.resize(n_boundary_patches);

    for (label id = 0; id < n_boundary_patches; id++)
    {
        maxVal[id] = std::sqrt(
                         tmpBoundarySnapshots[id].colwise().lpNorm<2>().maxCoeff());
        // get surface areas
        Eigen::VectorXd S = Foam2Eigen::field2Eigen(
                                sList[0].mesh().magSf().boundaryField()[id]);
        S = S.replicate(field_dim, 1);
        unsigned int bSize = field_dim * int(n_boundary_cells_list[id]);
        fieldWeightsBoundary[id].conservativeResize(fieldWeightsBoundary[id].rows() +
                bSize);
        fieldWeightsBoundary[id].tail(bSize) = S.array().pow(3 / 4) * maxVal[id];
        snapshotsMatrixBoundary[id].conservativeResize(
            snapshotsMatrixBoundary[id].rows() + bSize, n_snapshots);
        snapshotsMatrixBoundary[id].bottomRows(bSize) = tmpBoundarySnapshots[id];
    }
}

template<typename... SnapshotsLists>
void HyperReduction<SnapshotsLists...>::evaluatePinv(Eigen::SparseMatrix<double>
        & Projector, Eigen::MatrixXd& Modes, Eigen::VectorXd& fieldWeights)
{
    assert(Projector.cols() > 0);
    Eigen::MatrixXd restricted = Projector.transpose() * Modes;
    Info << "S-Optimalty: " << s_optimality(restricted) << endl;
    pinvPU = ITHACAutilities::invertMatrix(restricted,
                                           para->ITHACAdict->lookupOrDefault<word>("InversionMethod",
                                               "completeOrthogonalDecomposition"));
    pinvPU = pinvPU * (Projector.transpose() *
                       fieldWeights.array().cwiseInverse().matrix()).asDiagonal();
}

template<typename... SnapshotsLists>
void HyperReduction<SnapshotsLists...>::evaluateWPU(Eigen::SparseMatrix<double>
        & Projector, Eigen::MatrixXd& Modes, Eigen::VectorXd& fieldWeights,
        Eigen::VectorXd& quadratureWeights)
{
    assert(Projector.cols() > 0);
    Eigen::VectorXd quadratureWeightsOrderedAsProjector(quadratureWeights.rows());
    unsigned int n_weightsPerField = quadratureWeights.rows() / vectorial_dim;
    unsigned int reorderIndex{0};

    for (unsigned int ith_field = 0; ith_field < vectorial_dim; ith_field++)
    {
        for (unsigned int ith_weight = 0; ith_weight < n_weightsPerField; ith_weight++)
        {
            quadratureWeightsOrderedAsProjector(reorderIndex) = quadratureWeights(
                    ith_weight + ith_field * n_weightsPerField);
            reorderIndex++;
        }
    }

    wPU = quadratureWeightsOrderedAsProjector.transpose() *
          (Projector.transpose() *
           fieldWeights.array().cwiseInverse().matrix()).asDiagonal();
}

template<typename... SnapshotsLists>
void HyperReduction<SnapshotsLists...>::generateSubmesh(label layers,
        const fvMesh& mesh)
{
    Info << "####### Extract submesh #######\n";
    ITHACAparameters* para(ITHACAparameters::getInstance());
    totalNodePoints = autoPtr<IOList<labelList >> (new IOList<labelList>(IOobject(
                          "totalNodePoints", para->runTime.time().constant(),  "../" + folderMethod,
                          para->mesh, IOobject::READ_IF_PRESENT, IOobject::NO_WRITE)));
    uniqueNodePoints = autoPtr<IOList<label >> (new IOList<label>(IOobject(
                           "uniqueNodePoints", para->runTime.time().constant(),  "../" + folderMethod,
                           para->mesh, IOobject::READ_IF_PRESENT, IOobject::NO_WRITE)));
    volScalarField Indici
    (
        IOobject
        (
            problemName + "_indices",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(problemName + "_indices",
                          dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0));
    submesh = autoPtr<fvMeshSubset>(new fvMeshSubset(mesh));
    bool offlineStage = !totalNodePoints().headerOk();

    if (offlineStage)
    {
        List<label> indices;

        for (label i = 0; i < nodePoints().size(); i++)
        {
            indices = ITHACAutilities::getIndices(mesh, nodePoints()[i], layers);
            totalNodePoints().append(indices);
        }

        labelList a = ListListOps::combine<labelList>(totalNodePoints(),
                      accessOp<labelList>());
        inplaceUniqueSort(a);
        uniqueNodePoints() = a;
        scalar zerodot25 = 0.25;
        ITHACAutilities::assignIF(Indici, zerodot25,
                                  uniqueNodePoints().List<label>::clone()());
        ITHACAutilities::assignONE(Indici, nodePoints());
        totalNodePoints().write();
        uniqueNodePoints().write();
        ITHACAstream::exportSolution(Indici, "1", folderMethod);
    }
    submesh->setCellSubset(uniqueNodePoints());
    submesh->subMesh().fvSchemes::readOpt() = mesh.fvSchemes::readOpt();
    submesh->subMesh().fvSolution::readOpt() = mesh.fvSolution::readOpt();
    submesh->subMesh().fvSchemes::read();
    submesh->subMesh().fvSolution::read();
    std::cout.clear();
    localNodePoints = global2local(nodePoints(), submesh());
    n_cellsSubfields = submesh().cellMap().size();
    Info << "####### End extract submesh size = " << n_cellsSubfields <<
         " #######\n";
    createMasks(offlineStage);
    Info << "####### End create masks #######\n";
}

template<typename... SnapshotsLists>
void HyperReduction<SnapshotsLists...>::createMasks(bool offlineStage)
{
    if (offlineStage)
    {
        field2submesh.resize(submesh().cellMap().size() * vectorial_dim,
                             n_cells * vectorial_dim);
        field2submesh.reserve(Eigen::VectorXi::Constant(n_cells * vectorial_dim, 1));

        for (unsigned int ith_subCell{0} ; ith_subCell < submesh().cellMap().size();
                ith_subCell++)
        {
            for (unsigned int ith_field = 0; ith_field < vectorial_dim; ith_field++)
            {
                field2submesh.insert(ith_subCell + ith_field * submesh().cellMap().size(),
                                     submesh().cellMap()[ith_subCell] + n_cells * ith_field) = 1;
            }
        }

        field2submesh.makeCompressed();
        submesh2nodes.resize(nodePoints().size() * vectorial_dim,
                             submesh().cellMap().size() * vectorial_dim);
        submesh2nodes.reserve(Eigen::VectorXi::Constant(submesh().cellMap().size() *
                              vectorial_dim, 1));
        submesh2nodesMask.resize(nodePoints().size() * vectorial_dim);
        int index_col = 0;

        for (unsigned int ith_node{0} ; ith_node < nodePoints().size(); ith_node++)
        {
            for (unsigned int ith_subCell{0} ; ith_subCell < submesh().cellMap().size();
                    ith_subCell++)
            {
                if (nodePoints()[ith_node] == submesh().cellMap()[ith_subCell])
                {
                    for (unsigned int ith_field = 0; ith_field < vectorial_dim; ith_field++)
                    {
                        submesh2nodes.insert(ith_node + nodePoints().size() * ith_field,
                                             index_col + ith_field * submesh().cellMap().size()) = 1;
                        submesh2nodesMask(ith_node + nodePoints().size() * ith_field) = index_col +
                            ith_field * submesh().cellMap().size();
                    }

                    break;
                }
                else
                {
                    index_col++;
                }
            }

            index_col = 0;
        }

        submesh2nodes.makeCompressed();
        cnpy::save(field2submesh, folderMethod + "/field2submesh.npz");
        cnpy::save(submesh2nodes, folderMethod + "/submesh2nodes.npz");
        cnpy::save(submesh2nodesMask, folderMethod + "/submesh2nodesMask.npy");
    }
    else
    {
        cnpy::load(field2submesh, folderMethod + "/field2submesh.npz");
        cnpy::load(submesh2nodes, folderMethod + "/submesh2nodes.npz");
        cnpy::load(submesh2nodesMask, folderMethod + "/submesh2nodesMask.npy");
    }
}

template<typename... SnapshotsLists>
List<label> HyperReduction<SnapshotsLists...>::global2local(
    List<label>& points, fvMeshSubset& submesh)
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

#endif
