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

#ifndef hyperReduction_templates_H
#define hyperReduction_templates_H

#include "hyperReduction.H"

// Using the Eigen library, using the SVD decomposition method to solve the
// matrix pseudo-inverse, the default error er is 0
Eigen::MatrixXd pinv_eigen_based(Eigen::MatrixXd &origin, const float er = 0) {
  // perform svd decomposition
  Eigen::JacobiSVD<Eigen::MatrixXd> svd_holder(origin, Eigen::ComputeThinU | Eigen::ComputeThinV);
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
  return V * S * U.transpose();
}

template <typename... SnapshotsLists>
void HyperReduction<SnapshotsLists...>::getModesSVD(
    SnapshotsTuple& snapshotsTuple, Eigen::MatrixXd& modesSVD, Eigen::VectorXd& fieldWeights)
{
  ITHACAparameters* para(ITHACAparameters::getInstance());
  if (!ITHACAutilities::check_folder(Folder + "/Modes/"))
  {
    // Initialize stacked snapshotsMatrix and normalizing weights
    Eigen::MatrixXd snapMatrix;
    std::apply([this, &snapMatrix, &fieldWeights](auto& ...snapList){(..., stackSnapshots(snapList, snapMatrix, fieldWeights));}, snapshotsTuple);

    Info << "####### Performing SVD for " << problemName << " #######" << endl;

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(fieldWeights.array().cwiseInverse().matrix().asDiagonal()*snapMatrix,
                                          Eigen::ComputeThinU | Eigen::ComputeThinV);

    Info << "####### End of SVD for " << problemName << " #######" <<
          endl;

    Eigen::VectorXd eigenValueseig;
    Eigen::MatrixXd eigenVectoreig;
    eigenValueseig = svd.singularValues().real();
    eigenVectoreig = svd.matrixU().real();

    modesSVD = eigenVectoreig;

    // TODO correct boundary conditions

    eigenValueseig = eigenValueseig / eigenValueseig.sum();
    Eigen::VectorXd cumEigenValues(eigenValueseig);

    for (label j = 1; j < cumEigenValues.size(); ++j)
    {
        cumEigenValues(j) += cumEigenValues(j - 1);
    }

    mkDir(Folder + "/Modes/");
    Eigen::saveMarketVector(eigenValueseig,
                            Folder+"/Modes/Eigenvalues_" + problemName, para->precision,
                            para->outytpe);
    Eigen::saveMarketVector(cumEigenValues,
                            Folder+"/Modes/CumEigenvalues_" + problemName, para->precision,
                            para->outytpe);
    cnpy::save(modesSVD, Folder + "/Modes/modes.npy");
    cnpy::save(fieldWeights, Folder + "/Modes/normalizingWeights.npy");
  }
  else
  {
    Info << "Reading the existing modes" << endl;
    cnpy::load(modesSVD, Folder + "/Modes/modes.npy");
    cnpy::load(fieldWeights, Folder + "/Modes/normalizingWeights.npy");
  }
}

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
      snapshotsTuple{std::forward_as_tuple(snapshotsLists...)}
{
  Info << "Init HyperReduction class with vectorial dim: " << vectorial_dim << endl;

  // TODO check snapshotsLists is not empty
  // TODO check that the first snapshotsList is not empty (includes above)

  // TODO initialize vectorial_dim only from SnapshotsLists type
  // vectorial_dim = std::apply(compute_vectorial_dim<SnapshotsLists&&...>, snapshotsTuple);

  // Get fields' names
  std::apply([this](auto& ...snapList){(..., stackNames(snapList));}, snapshotsTuple);

  // Get fields' dimensions
  std::apply([this](auto& ...snapList){(..., stackDimensions(snapList));}, snapshotsTuple);

  Info << "Fields names: ";
  for (unsigned int ith_field = 0; ith_field < fieldNames.size(); ith_field++)
  {
    Info << fieldNames[ith_field] << " (dim=" << fieldDims[ith_field] << "); ";
  }
  Info << endl;
  
  n_snapshots = std::get<0>(snapshotsTuple).size();
  Info << "The number of snapshots is: " << n_snapshots << endl;
  n_cells = std::get<0>(snapshotsTuple)[0].size();
  Info << "The number of cells is: " << n_cells << endl;
  Info << "Initial seeds length: " << initialSeeds.rows() << endl;
}

template <typename... SnapshotsLists>
void HyperReduction<SnapshotsLists...>::offlineStage()
{
  ITHACAparameters *para(ITHACAparameters::getInstance());
  word methodName = para->ITHACAdict->lookupOrDefault<word>("HyperReduction", "GappyDEIM");
  
  if (methodName=="GappyDEIM")
  {
    Folder = "ITHACAoutput/" + problemName + "/" + methodName;
    mkDir(Folder);

    nodePoints = autoPtr<IOList<label>>(new IOList<label>(
        IOobject("nodePoints", para->runTime.time().constant(), "../" + Folder,
                  para->mesh, IOobject::READ_IF_PRESENT, IOobject::NO_WRITE)));

    Eigen::MatrixXd snapshotsModes;
    getModesSVD(snapshotsTuple, snapshotsModes, normalizingWeights);  

    if (!nodePoints().headerOk())
    {
        assert(n_modes > 0);
        assert(n_nodes >= n_modes);

        Eigen::VectorXd mp_not_mask = Eigen::VectorXd::Constant(n_cells * vectorial_dim, 1);
        std::set<label> nodePointsSet;
        std::deque<label> nodePointsDeque;

        // set initialSeeds
        if (initialSeeds.rows() > 0)
        {
          Info << "Add initialSeeds\n";
          P.resize(n_cells * vectorial_dim, initialSeeds.rows() * vectorial_dim);
          P.reserve(Eigen::VectorXi::Constant(initialSeeds.rows() * vectorial_dim /*n_cols*/, 1 /*n_non_zero_elements*/));

          for (int i = 0; i < initialSeeds.rows(); i++)
          {
            unsigned int index = initialSeeds(i) % n_cells;
            for (unsigned int ith_field = 0; ith_field < vectorial_dim; ith_field++)
            {
              P.insert(index+ith_field*n_cells, vectorial_dim*i+ith_field) = 1;
              mp_not_mask(index + ith_field * n_cells) = 0;
            }
            
            // check that are not repeated nodes in initialSeeds
            if (nodePointsSet.find(index) == nodePointsSet.end()){
              nodePointsSet.insert(index);
              nodePointsDeque.push_back(index);
            }
          }

          P.makeCompressed();
        }

        int na = n_nodes - initialSeeds.rows();
        if (na > 0)
        {
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

            Info << "####### Begin Greedy GappyDEIM #######" << endl;
            Info << "####### Modes=" <<  n_modes
                  << ", nodePoints=" << n_nodes << " #######"
                  << endl;
            


            Eigen::SparseMatrix<double> reshapeMat;
            reshapeMat.resize(n_cells, n_cells*vectorial_dim);
            reshapeMat.reserve(Eigen::VectorXi::Constant(n_cells * vectorial_dim /*n_cols*/, 1 /*n_non_zero_elements*/));

            for (unsigned int i = 0; i < n_cells; i++)
            {
              for (unsigned int ith_field = 0; ith_field < vectorial_dim; ith_field++)
              {
                reshapeMat.insert(i, i+n_cells*ith_field) = 1;
              }
            }
              
            reshapeMat.makeCompressed();

            for (int i = 1; i <= nit; i++)
            {
                int nci = ncimin;
                // add the remaining modes in the quotient n_modes / nite
                if (i <= n_modes % nit)
                {
                    nci = ncimin + 1;
                }
                int nai = naimin;
                // add the remaining nodePoints_xyz in the quotient na / n_modes
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
                    max = (reshapeMat*(mp_not_mask.asDiagonal() * V)
                              .rowwise()
                              .lpNorm<2>())
                              .maxCoeff(&ind_max, &c1);
                  }
                  else
                  {
                    max = V.rowwise().lpNorm<2>().maxCoeff(&ind_max, &c1);
                  }
                  int idxes = P.cols();
                  int maxFirstMesh = ind_max % n_cells;
                  P.conservativeResize(snapshotsModes.rows(), idxes + vectorial_dim);
                  
                  for (unsigned int ith_field = 0; ith_field < vectorial_dim; ith_field++)
                  {
                    P.insert(maxFirstMesh + ith_field * n_cells, idxes + ith_field) = 1;
                    mp_not_mask(maxFirstMesh + ith_field * n_cells) = 0;
                  }

                  nodePointsDeque.push_back(maxFirstMesh);
                }

                nb += nci;
                basisMatrix = snapshotsModes.leftCols(nb);
            }
        }
        else
        {
            basisMatrix = snapshotsModes;
        }

        Info << "####### End of greedy GappyDEIM #######\n";
        cnpy::save(basisMatrix, Folder +"/basisMatrix.npy");
        cnpy::save(P, Folder +"/projectionMatrix.npy");
        Info << "Projection Matrix shape: " << P.rows() << " " << P.cols() << endl;
        Info << "Basis Matrix shape: " << basisMatrix.rows() << " " << basisMatrix.cols() << endl;

        evaluatePinv(P, basisMatrix, normalizingWeights);
        renormalizedBasisMatrix = normalizingWeights.asDiagonal()*basisMatrix;
        Info << "Pseudo Inverse Matrix shape: " << pinvPU.rows() << " " << pinvPU.cols() << endl;

        // convert nodePointsSet to nodePoints i.e. list
        Eigen::VectorXi nodes(nodePointsDeque.size());
        unsigned int i{0};
        for (auto &x : nodePointsDeque)
        {
          nodePoints().append(x);
          nodes(i) = x;
          i++;
        }
        cnpy::save(nodes, Folder + "/mp.npy");
        cnpy::save(pinvPU, Folder + "/pinvPU.npy");
        nodePoints().write();
    }
    else
    {
        cnpy::load(pinvPU, Folder + "/pinvPU.npy");
        cnpy::load(basisMatrix, Folder +"/basisMatrix.npy");
    }
  }
  else
  {
    Info << "Offline stage not performed. Invalid HR method.\n";
  }
}

template <typename... SnapshotsLists>
template <typename SnapshotsList>
void HyperReduction<SnapshotsLists...>::stackSnapshots(SnapshotsList sList, Eigen::MatrixXd& snapshotsMatrix, Eigen::VectorXd& fieldWeights)
{
  Info << "####### Compute hyper-reduction basis #######\n";

  unsigned int field_dim = get_field_dim<typename SnapshotsList::value_type>();

  Eigen::MatrixXd tmpSnapshots = Foam2Eigen::PtrList2Eigen(sList);
  auto fieldName = sList[0].name();

  double maxVal = std::sqrt(tmpSnapshots.colwise().lpNorm<2>().maxCoeff());

  // get volumes
  Eigen::VectorXd V = ITHACAutilities::getMassMatrixFV(sList[0]);

  fieldWeights.conservativeResize(fieldWeights.rows()+field_dim*n_cells);
  fieldWeights.tail(n_cells * field_dim) = V.array().sqrt().cwiseInverse()*maxVal;

  snapshotsMatrix.conservativeResize(snapshotsMatrix.rows() + n_cells * field_dim, n_snapshots);
  snapshotsMatrix.bottomRows(n_cells * field_dim) = tmpSnapshots;

  Info << "####### End of computing hyper-reduction basis #######\n";
}

template <typename... SnapshotsLists>
template<typename SnapshotsList>
void HyperReduction<SnapshotsLists...>::stackNames(SnapshotsList sList)
{
  fieldNames.append(sList[0].name());
}

template <typename... SnapshotsLists>
template<typename SnapshotsList>
void HyperReduction<SnapshotsLists...>::stackDimensions(SnapshotsList sList)
{
  fieldDims.append(get_field_dim<typename SnapshotsList::value_type>());
}

template<typename... SnapshotsLists>
void HyperReduction<SnapshotsLists...>::evaluatePinv(Eigen::SparseMatrix<double> &Projector, Eigen::MatrixXd &Modes, Eigen::VectorXd &fieldWeights) {
  assert(Projector.cols() > 0);
  Eigen::MatrixXd restricted = Projector.transpose() * Modes;
  pinvPU = pinv_eigen_based(restricted);
  pinvPU = pinvPU*(Projector.transpose() * fieldWeights.array().cwiseInverse().matrix()).asDiagonal();
}

template<typename... SnapshotsLists>
void HyperReduction<SnapshotsLists...>::generateSubmesh(label layers, fvMesh &mesh) {

  Info << "####### Extract submesh #######\n";
  ITHACAparameters *para(ITHACAparameters::getInstance());

  totalNodePoints = autoPtr<IOList<labelList>>(new IOList<labelList>(IOobject(
      "totalNodePoints", para->runTime.time().constant(),  "../" + Folder,
      para->mesh, IOobject::READ_IF_PRESENT, IOobject::NO_WRITE)));

  uniqueNodePoints = autoPtr<IOList<label>>(new IOList<label>(IOobject(
      "uniqueNodePoints", para->runTime.time().constant(),  "../" + Folder,
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

  if (!totalNodePoints().headerOk()) {
    List<label> indices;
    for (label i = 0; i < nodePoints().size(); i++) {
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
  }

  submesh->setCellSubset(uniqueNodePoints());
  submesh->subMesh().fvSchemes::readOpt() = mesh.fvSchemes::readOpt();
  submesh->subMesh().fvSolution::readOpt() = mesh.fvSolution::readOpt();
  submesh->subMesh().fvSchemes::read();
  submesh->subMesh().fvSolution::read();
  std::cout.clear();

  localNodePoints = global2local(nodePoints(), submesh());
  ITHACAstream::exportSolution(Indici, "1", Folder);

  n_cellsSubfields = submesh().cellMap().size();
  Info << "####### End extract submesh #######\n";
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
autoPtr<FieldType> HyperReduction<SnapshotsLists...>::interpolateField(FieldType& field) {
    return autoPtr<FieldType>(new FieldType(submesh->interpolate(field)));
}

#endif
