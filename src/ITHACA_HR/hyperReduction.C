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

#include "hyperReduction.H"

// Using the Eigen library, using the SVD decomposition method to solve the
// matrix pseudo-inverse, the default error er is 0
Eigen::MatrixXd pinv_eigen_based(Eigen::MatrixXd &origin, const float er = 0) {
  // perform svd decomposition
  Eigen::JacobiSVD<Eigen::MatrixXd> svd_holder(origin, Eigen::ComputeThinU |
                                                           Eigen::ComputeThinV);
  // Build SVD decomposition results
  Eigen::MatrixXd U = svd_holder.matrixU();
  Eigen::MatrixXd V = svd_holder.matrixV();
  Eigen::MatrixXd D = svd_holder.singularValues();

  // Build the S matrix
  Eigen::MatrixXd S(V.cols(), U.cols());
  S.setZero();

  for (unsigned int i = 0; i < D.size(); ++i) {

    if (D(i, 0) > er) {
      S(i, i) = 1 / D(i, 0);
    } else {
      S(i, i) = 0;
    }
  }

  std::cout << "Condition number of inverse: " <<  D(D.size()-1, 0)/D(0, 0) << std::endl;
  // pinv_matrix = V * S * U^T
  return V * S * U.transpose();
}

// Template function constructor
template <typename... SnapshotsLists>
HyperReduction<SnapshotsLists...>::HyperReduction(HyperReductionMethod hrMethod,
                                                  label maxModes,
                                                  label maxNodes,
                                                  Eigen::VectorXi initialSeeds,
                                                  word functionName,
                                                  SnapshotsLists &&...snapshotsLists)
    :  methodName{methodNames[int(hrMethod)]},
      vectorial_dim{compute_vectorial_dim(snapshotsLists...)},
      maxModes{maxModes},
      maxNodes{maxNodes},
      initialSeeds{initialSeeds},
      functionName{functionName},
      snapshotsTuple{std::forward_as_tuple(snapshotsLists...)}
{
    // check snapshotsLists is not empty
    // check that the first snapshotsList is not empty (includes above)

    // TODO initialize vectorial_dim only from SnapshotsLists type
    // vectorial_dim = std::apply(compute_vectorial_dim<SnapshotsLists&&...>, snapshotsTuple);

    ITHACAparameters *para(ITHACAparameters::getInstance());
    Folder = "ITHACAoutput/" + methodName + "/" + functionName;

    nodePoints = autoPtr<IOList<label>>(new IOList<label>(
        IOobject("nodePoints", para->runTime.time().constant(), "../" + Folder,
                 para->mesh, IOobject::READ_IF_PRESENT, IOobject::NO_WRITE)));

    nodePointsFields = autoPtr<IOList<label>>(new IOList<label>(IOobject(
        "nodePointsFields", para->runTime.time().constant(), "../" + Folder,
        para->mesh, IOobject::READ_IF_PRESENT, IOobject::NO_WRITE)));

    n_snapshots = std::get<0>(snapshotsTuple).size();
    Info << "The number of snapshots is: " << n_snapshots << endl;
    n_cells = std::get<0>(snapshotsTuple)[0].size();
    Info << "The number of cells is: " << n_cells << endl;
    
    std::apply([this](auto& ...x){(..., compute_matrix_modes(x));}, snapshotsTuple);

    if (!nodePoints().headerOk())
    {

        assert(maxModes > 0);
        assert(maxNodes >= maxModes);

        Eigen::VectorXd mp_not_mask = Eigen::VectorXd::Constant(n_cells * vectorial_dim, 1);

        // set initialSeeds
        if (initialSeeds.rows() > 0)
        {

            P.resize(n_cells * vectorial_dim, initialSeeds.rows() * vectorial_dim);
            P.reserve(Eigen::VectorXi::Constant(initialSeeds.rows() * vectorial_dim /*n_cols*/, 1 /*n_non_zero_elements*/));

            int col_index{0};
            for (int i = 0; i < initialSeeds.rows(); i++)
            {
                for (unsigned int ith_field = 0; ith_field < vectorial_dim; ith_field++)
                {
                    P.insert(initialSeeds(i)+ith_field*n_cells, col_index++) = 1;
                    mp_not_mask(initialSeeds(i) + ith_field * n_cells) = 0;
                }
                
                nodePointsSet.insert(initialSeeds(i));
            }

            P.makeCompressed();
        }

        // std::cout << "Initial seeds: " << P << std::endl;

        for (int k = 0; k < initialSeeds.outerSize(); ++k) {
          for (Eigen::SparseMatrix<double>::InnerIterator it(P, k); it; ++it) {
            assert(it.value() == 1);
            nodePointsSet.insert(it.row() % n_cells);
          }
        }

        int na = maxNodes - initialSeeds.rows();
        if (na > 0)
        {
            int nb = 0;
            int nit = std::min(maxModes, na);
            int ncimin = std::floor(maxModes / nit);
            int naimin = std::floor(na / maxModes);

            Eigen::MatrixXd A;
            Eigen::VectorXd b;
            Eigen::VectorXd c;
            Eigen::VectorXd r;
            label ind_max, c1;

            double max;

            Info << "######## Begin HR training ########" << endl;
            Info << "######## Modes=" <<  maxModes
                 << ", nodePoints=" << maxNodes << " ########"
                 << endl;

            for (int i = 1; i <= nit; i++)
            {
                int nci = ncimin;
                // add the remaining modes in the quotien maxModes / nite
                if (i <= maxModes % nit)
                {
                    nci = ncimin + 1;
                }
                int nai = naimin;
                // add the remaining nodePoints_xyz in the quotien na / maxModes
                if (i <= na % maxModes)
                {
                    nai = naimin + 1;
                }

                Eigen::MatrixXd V;

                // select basis
                if (i == 1)
                {
                    V = matrixModes.leftCols(nci);
                }
                else
                {
                    for (int q = 1; q <= nci; q++)
                    {
                        A = P.transpose() * U;

                        b = P.transpose() * matrixModes.col(nb + q - 1);

                        c = A.fullPivLu().solve(b);

                        r = matrixModes.col(nb + q - 1) - U * c;

                        V.conservativeResize(matrixModes.rows(), q);
                        V.col(q - 1) = r;
                    }
                }

                // select nodePoints_xyz
                for (int j = 1; j <= nai; j++)
                {
                    // TODO: implement as in Carlberg paper taking into account degrees
                    // of freedom associated to each node
                    // auto Id = Eigen::MatrixXd::Identity(P.rows(), P.rows());

                    if (P.cols() > 0)
                    {
                        // std::cout << mp_not_mask.asDiagonal().rows() << " " << V.rows() << V.cols() << std::endl;
                        max = (mp_not_mask.asDiagonal() * V)
                                  .rowwise()
                                  .lpNorm<2>()
                                  .maxCoeff(&ind_max, &c1);
                    }
                    else
                    {
                        // max = V.topRows(2*n_cells).rowwise().lpNorm<2>().maxCoeff(&ind_max, &c1);
                        max = V.rowwise().lpNorm<2>().maxCoeff(&ind_max, &c1);
                    }

                    // std::cout << "Maximum found: " << max << " at " << ind_max
                    //           << std::endl;

                    int idxes = P.cols();
                    int maxFirstMesh = ind_max % n_cells;
                    P.conservativeResize(matrixModes.rows(), idxes + vectorial_dim);
                    
                    for (unsigned int ith_field = 0; ith_field < vectorial_dim; ith_field++)
                    {
                        P.insert(maxFirstMesh + ith_field * n_cells, idxes + ith_field) = 1;
                        mp_not_mask(maxFirstMesh + ith_field * n_cells) = 0;
                    }

                    nodePointsSet.insert(maxFirstMesh);
                }

                nb += nci;
                U = matrixModes.leftCols(nb);
            }
        }
        else
        {
            U = matrixModes;
        }

        // convert nodePointsSet to nodePoints i.e. list
        for (const int &x : nodePointsSet)
        {
            nodePoints().append(x);
        }

        // int i{0};
        // nodePoints->resize(nodePoints().size());
        // for (int &x : nodePoints())
        // {
        //     std::cout << x << " ";
        //     nodePoints()(i) = x;
        //     i++;
        // }
        // cnpy::save(nodePoints(), "./numpy_datasets/mp.npy");
        // Info << "\n # DEBUG GNAT.C, line 740 # " << endl;
        // std::cout << std::endl;

        // evaluatePinv(P, U);
        // std::cout << "pinv_PU: " << pinv_PU.rows() << std::endl;

        mkDir(Folder);
        cnpy::save(pinv_PU, Folder + "/MatrixOnline.npy");
        // cnpy::save(pinv_PUV, Folder + "/MatrixOnlineV.npy");
        // cnpy::save(pinv_PUS, Folder + "/MatrixOnlineS.npy");
        cnpy::save(maskAux, Folder + "/MaskAux.npy");

        // std::cout << "nodePoints size: " << nodePoints().size() << "\n" << std::endl;
        // nodePoints.resize(nodePoints().size(), 1);
        // for (int i = 0; i < nodePoints().size(); i++)
        // {
        //   nodePoints.coeffRef(i, 0) = nodePoints()[i];
        // }
        // cnpy::save(nodePoints, "./nummpy_datasets/mp.npy");
        nodePoints().write();
    }
    else
    {
        cnpy::load(pinv_PU, Folder + "/MatrixOnline.npy");
        // cnpy::load(pinv_PUV, Folder + "/MatrixOnlineV.np
        // cnpy::load(pinv_PUS, Folder + "/MatrixOnlineS.npy");y");
        // cnpy::load(pinv_PUS, Folder + "/MatrixOnlineS.npy");
        cnpy::load(maskAux, Folder + "/MaskAux.npy");
        Info << "Loaded mask and UpinvPU with size: "
                  << pinv_PU.rows() << "x" << pinv_PU.cols() << endl;
    }
}

template <typename... SnapshotsLists>
template <typename SnapshotsList>
void HyperReduction<SnapshotsLists...>::compute_matrix_modes(SnapshotsList sList)
{
    // TODO normalization
    // Eigen::MatrixXf::Index maxIndex;
    // double maxVal = modes.col(0).cwiseAbs().maxCoeff(&maxIndex);
    // std::cout << "MaxValues: " << maxVal << std::endl;

    // modes = modes / maxVal;

    // maxVal = modes.col(0).cwiseAbs().maxCoeff(&maxIndex);
    // std::cout << "MaxValues: " << maxVal << std::endl;

    word fieldName = sList[0].name();
    unsigned int field_dim = get_field_dim<typename SnapshotsList::value_type>();

    auto tmp_modes = ITHACAPOD::DEIMmodes(sList, maxModes, "GNAT/" + functionName + "V", fieldName);

    matrixModes.conservativeResize(n_cells * field_dim, n_snapshots);
    matrixModes.bottomRows(n_cells * field_dim) = Foam2Eigen::PtrList2Eigen(tmp_modes);
}

template<typename... SnapshotsLists>
void HyperReduction<SnapshotsLists...>::evaluatePinv(Eigen::SparseMatrix<double> &Projector, Eigen::MatrixXd &Modes) {
  assert(Projector.cols() > 0);
  Eigen::MatrixXd restricted = Projector.transpose() * Modes;

//   Info << "Restricted modes matrix: " << restricted.rows() << " x "
//             << restricted.cols() << endl;

  int ssize = restricted.rows()/3;
  int vsize = 2*ssize;
  
  pinv_PU = pinv_eigen_based(restricted);
//   Eigen::MatrixXd PUV = restricted.topRows(vsize);
//   Eigen::MatrixXd PUS = restricted.bottomRows(ssize);
//   cnpy::save(restricted, "./numpy_datasets/restricted.npy");
//   pinv_PUV = pinv_eigen_based(PUV, 1e-10);
//   pinv_PUS = pinv_eigen_based(PUS, 1e-10);
//   Info << "PUV: " << PUV.rows() << " " << PUV.cols() << endl;
//   Info << "PUS: " << PUS.rows() << " " << PUS.cols() << endl;

//   pinv_PU.resize(restricted.cols(), ssize*3);
//   pinv_PU.leftCols(vsize) = pinv_PUV;
//   pinv_PU.rightCols(ssize) = pinv_PUS;
  // Info << "pinv_PU: " << pinv_PU << endl;

//   Info << "Pseudo inverse of restricted modes: " << pinv_PU.rows() << " x "
//             << pinv_PU.cols() << endl;

  Eigen::MatrixXd Upinv_PU_ = restricted * pinv_PU;
  Info << "Pseudo inverse of restricted modes: " << Upinv_PU_.rows()
            << " x " << Upinv_PU_.cols() << endl;

  maskingAux = Modes * pinv_PU;
//   Info << "Masking auxiliary matrix: " << maskingAux.rows() << " x "
//             << maskingAux.cols() << endl;

  // auto ret = std::make_pair<Eigen::MatrixXd, Eigen::MatrixXd>(Upinv_PU, maskingAux);
}

template<typename... SnapshotsLists>
void HyperReduction<SnapshotsLists...>::generateSubmeshes(label layers, fvMesh &mesh, label secondTime) {

  ITHACAparameters *para(ITHACAparameters::getInstance());

  totalMagicPoints = autoPtr<IOList<labelList>>(new IOList<labelList>(IOobject(
      "totalMagicPoints", para->runTime.time().constant(),  Folder,
      para->mesh, IOobject::READ_IF_PRESENT, IOobject::NO_WRITE)));

  uniqueMagicPoints = autoPtr<IOList<label>>(new IOList<label>(IOobject(
      "uniqueMagicPoints", para->runTime.time().constant(),  Folder,
      para->mesh, IOobject::READ_IF_PRESENT, IOobject::NO_WRITE)));

  volScalarField Indici(
      IOobject(FunctionName + "_indices", mesh.time().timeName(), mesh,
               IOobject::NO_READ, IOobject::NO_WRITE),
      mesh,
      dimensionedScalar(FunctionName + "_indices",
                        dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0));

  submesh = autoPtr<fvMeshSubset>(new fvMeshSubset(mesh));

  if (!totalMagicPoints().headerOk()) {
    List<label> indices;

    for (label i = 0; i < magicPoints().size(); i++) {
      indices = ITHACAutilities::getIndices(mesh, magicPoints()[i], layers);
      totalMagicPoints().append(indices);
    }

    labelList a = ListListOps::combine<labelList>(totalMagicPoints(),
                                                  accessOp<labelList>());

    inplaceUniqueSort(a);
    uniqueMagicPoints() = a;

    scalar zerodot25 = 0.25;
    ITHACAutilities::assignIF(Indici, zerodot25,
                              uniqueMagicPoints().List<label>::clone()());
    ITHACAutilities::assignONE(Indici, magicPoints());

    totalMagicPoints().write();
    uniqueMagicPoints().write();
  }

  submesh->setCellSubset(uniqueMagicPoints());
  submesh->subMesh().fvSchemes::readOpt() = mesh.fvSchemes::readOpt();
  submesh->subMesh().fvSolution::readOpt() = mesh.fvSolution::readOpt();
  submesh->subMesh().fvSchemes::read();
  submesh->subMesh().fvSolution::read();
  std::cout.clear();

  if (!secondTime) {
      localMagicPoints = global2local(magicPoints(), submesh());
      ITHACAstream::exportSolution(Indici, "1",
                                 "./ITHACAoutput/" + methodName + "/");
  }

  NcellsSubfields = submesh().cellMap().size();
  createMasks();
}

template<typename... SnapshotsLists>
List<label> HyperReduction<SnapshotsLists...>::global2local(List<label> &points, fvMeshSubset &submesh) {
  List<label> localPoints;

  for (label i = 0; i < points.size(); i++) {
    for (label j = 0; j < submesh.cellMap().size(); j++) {
      if (submesh.cellMap()[j] == points[i]) {
        localPoints.append(j);
        break;
      }
    }
  }

  return localPoints;
}

template<typename... SnapshotsLists>
template<typename FieldType>
void HyperReduction<SnapshotsLists...>::interpolateField(FieldType& field) {
    return autoPtr<FieldType>(new FieldType(submesh->interpolate(field)));
}

template<typename... SnapshotsLists>
void HyperReduction<SnapshotsLists...>::createMasks() {
  assert(uniqueMagicPoints().size() > 0);

  // initialize nonRepeatedSubmeshes
  nonRepeatedSubmeshes.resize(uniqueMagicPoints().size());
  int vector_idx{0};
  for (int &x : uniqueMagicPoints()) {
    nonRepeatedSubmeshes(vector_idx) = x;
    vector_idx++;
  }

  // initialize nonRepeatedSubmeshes_xy
  nonRepeatedSubmeshes_xy.resize(uniqueMagicPoints().size() * 3);
  nonRepeatedSubmeshes_xy.head(uniqueMagicPoints().size()) =
      nonRepeatedSubmeshes;
  nonRepeatedSubmeshes_xy.segment(uniqueMagicPoints().size(), uniqueMagicPoints().size()) = (Ncells + nonRepeatedSubmeshes.array()).matrix();
  nonRepeatedSubmeshes_xy.tail(uniqueMagicPoints().size()) =
      (2 * Ncells + nonRepeatedSubmeshes.array()).matrix();


  // initialize P_xymp_mp
  // due to the fact that x, y coordinates of the residual in a single cell are
  // computed at the same time in OF, if two magicPoints_xy refer to the same
  // cell, the residual can be evaluated just once. In case only one coordinate
  // is a magicPoint the following mask will be used to distinguish the two
  // cases
  cnpy::load(magicpoints, "./numpy_datasets/mp.npy");
  // std::cout << magicpoints << std::endl;

  P_xymp_submesh.resize(2 * magicpoints.rows(), 2 * NcellsSubfields);
  P_xymp_submesh_scalar.resize(magicpoints.rows(), NcellsSubfields);


  // ITHACAparameters *para(ITHACAparameters::getInstance());
  // magicPoints.clear();
  // magicPoints = autoPtr<IOList<label>>(new IOList<label>(
  //     IOobject("magicPoints", para->runTime.time().constant(), "../ITHACAoutput/GNAT/" + FunctionName,
  //              para->mesh, IOobject::READ_IF_PRESENT, IOobject::NO_WRITE)));

  // for (int i{0}; i< magicPoints().size(); i++){
  //   std::cout << magicPoints()[i] << " ";
  // }

  int index_row = 0;
  int index_col = 0;
  for (const int &cell : submesh().cellMap()) {
    for (int mpxy{0}; mpxy < magicpoints.rows(); mpxy++) {
      // std::cout << "cell: "<< cell << " mp: " << magicpoints(mpxy) << std::endl;
      if (cell == int(magicpoints(mpxy))) {
        P_xymp_submesh.insert(index_row, index_col) = 1;
        P_xymp_submesh.insert(index_row + magicpoints.rows(), index_col + NcellsSubfields) = 1;

        P_xymp_submesh_scalar.insert(index_row, index_col) = 1;
      }

      index_row++;volVectorField fieldV, volScalarField fieldS
    }
    index_row = 0;
    index_col++;
  }
//   Info << " # DEBUG GNAT.C, line 429 # " << P_xymp_submesh.rows() << " " << P_xymp_submesh.cols() << endl;
// //   Info << " # DEBUG GNAT.C, line 430 # " << P_xymp_submesh_scalar.rows() << " " << P_xymp_submesh_scalar.cols() << endl;

  Eigen::MatrixXd P_xymp_submesh_ = P_xymp_submesh;
  Eigen::MatrixXd P_xymp_submesh_scalar_ = P_xymp_submesh_scalar;
  cnpy::save(P_xymp_submesh_, "./numpy_datasets/submesh_vector.npy");
  cnpy::save(P_xymp_submesh_scalar_, "./numpy_datasets/submesh_scalar.npy");
}

template HyperReduction<PtrList<volScalarField>&>::HyperReduction(HyperReductionMethod hrMethod,
                                                        label maxModes,
                                                        label maxNodes,
                                                        Eigen::VectorXi initialSeeds,
                                                        word functionName,
                                                        PtrList<volScalarField>& snapshotsLists1);

// template HyperReduction<PtrList<volVectorField>&>::HyperReduction(HyperReductionMethod hrMethod,
//                                                         label maxModes,
//                                                         label maxNodes,
//                                                         Eigen::VectorXi initialSeeds,
//                                                         word functionName,
//                                                         PtrList<volVectorField>& snapshotsLists1);

// template HyperReduction<PtrList<volVectorField>&, PtrList<volScalarField>&>::HyperReduction(HyperReductionMethod hrMethod,
//                                                         label maxModes,
//                                                         label maxNodes,
//                                                         Eigen::VectorXi initialSeeds,
//                                                         word functionName,
//                                                         PtrList<volVectorField>& snapshotsLists1,
//                                                         PtrList<volScalarField>&snapshotsLists2);