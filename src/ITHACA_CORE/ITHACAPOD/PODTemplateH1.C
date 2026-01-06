#include "PODTemplateH1.H"

namespace ITHACAPOD {

template <typename T, typename G>
PODTemplateH1<T,G>::PODTemplateH1(Parameters* myParameters, const word& myfield_name) :
  PODTemplate<T>::PODTemplate(myParameters, myfield_name),
  gradfield_name("grad"+field_name)
{
}

template <typename T, typename G>
Eigen::MatrixXd PODTemplateH1<T,G>::buildCovMatrix()
{
  precomputeGradients();

  l_hilbertSp = "L2";
  PODTemplate<T>::define_paths();
  Eigen::MatrixXd covMatrix = PODTemplate<T>::buildCovMatrix();
  l_hilbertSp = m_parameters->get_hilbertSpacePOD()[field_name];
  PODTemplate<T>::define_paths();

  m_parameters->set_hilbertSpacePOD(gradfield_name, "L2");
  m_parameters->set_nModes(gradfield_name, l_nmodes);
  PODTemplate<G> ithacaFVPODgrad(m_parameters,gradfield_name);
  ithacaFVPODgrad.computeMeanField();
  Eigen::MatrixXd covMatrixGrad = ithacaFVPODgrad.buildCovMatrix();

  if (l_hilbertSp == "wH1")
  {
    weightH1 = covMatrix.trace()/covMatrixGrad.trace();
    covMatrixGrad *= weightH1;
  }
  else
  {
    weightH1 = 1.0;
  }
  m_parameters->set_weightH1(weightH1);

  covMatrix += covMatrixGrad;

  if (Pstream::master())
  {
    mkDir(folder_covMatrix);
    cnpy::save(covMatrix, folder_covMatrix + name_covMatrix+ ".npy");
  }

  double varyingEnergy = covMatrix.trace()/l_nSnapshot;
  m_parameters->set_varyingEnergy( field_name + "_" + l_hilbertSp, varyingEnergy);
  Info << "Total varying energy " << l_hilbertSp << " : " << varyingEnergy << endl;
  Info << endl;

  return covMatrix;
}

template <typename T, typename G>
void PODTemplateH1<T,G>::precomputeGradients()
{
  T snapshotj = *f_field;
  G gradSnapshotsj = fvc::grad(*f_field);
  word local_file;
  bool exist_precomputed_fields = true;

  word pathProcessor("");
  if (Pstream::parRun())
  {
      pathProcessor = "processor" + name(Pstream::myProcNo()) + "/";
  }

  local_file = runTime2.caseName() + pathProcessor + runTime2.times()[1].name()
    +  "/" + gradfield_name;
  exist_precomputed_fields = exist_precomputed_fields && ITHACAutilities::check_file(local_file);
  for (label j = 0; j < l_nSnapshot + l_nSnapshotSimulation - 1 ; j++)
  {
    // Read the j-th field
    local_file = runTime2.caseName() + pathProcessor + runTime2.times()[l_startTime + j].name()
      +  "/" + gradfield_name;
    exist_precomputed_fields = exist_precomputed_fields && ITHACAutilities::check_file(local_file);
  }

  if (!exist_precomputed_fields)
  {
    Info << "Evaluating gradient fields" << endl;
    local_file = runTime2.caseName() + ".";
    ITHACAstream::read_snapshot(snapshotj, "0", runTime2.caseName());
    gradSnapshotsj = fvc::grad(snapshotj);
    ITHACAstream::exportSolution(gradSnapshotsj,
                                 "0",
                                 local_file,
                                 gradfield_name);
    for (label j = 0; j < l_nSnapshot + l_nSnapshotSimulation - 1; j++)
      {
        label index = l_startTime + j;
        ITHACAstream::read_snapshot(snapshotj, timeFolders[index].name(), runTime2.caseName());
        gradSnapshotsj = fvc::grad(snapshotj);
        // Write the j-th field
        string idxTimeStr( runTime2.times()[ index ].name() );
        fileName subfolder = idxTimeStr;
        ITHACAstream::exportSolution(gradSnapshotsj,
                                     subfolder,
                                     local_file,
                                     gradfield_name);
      }
  }
} 

template <typename T, typename G>
void PODTemplateH1<T,G>::changeEigenFolderUnitTest(const word& name) {
        
    l_nmodes = l_nSnapshot -1;
    word pathCentered("");
    if (this->b_centeredOrNot && !(this->field_name=="U"))
    {
        pathCentered = "_centered";
    }
    this->folder_eigen = "./ITHACAoutput/AllEigenValues_" + name
    + pathCentered + "_" + std::to_string(this->l_nmodes) + "snapshots/";

    this->exist_eigenDecomposition = ITHACAutilities::check_file(this->folder_eigen + this->name_eigenValues+".npy")
        && ITHACAutilities::check_file(this->folder_eigen + this->name_eigenVector+".npy");

}

}
