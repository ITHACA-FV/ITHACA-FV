#ifndef IthacaFVPODVectorH1_H
#define IthacaFVPODVectorH1_H

#include "PODTemplate.H"
#include "fvCFD.H"


namespace ITHACAPOD {

// Class that compute POD (spatial modes and temporal modes)
template<typename T, typename G>
class PODTemplateH1 : public PODTemplate<T>
{
  using PODTemplate<T>::ithacaFVParameters;
  using PODTemplate<T>::field_name;
  using PODTemplate<T>::casenameData;
  using PODTemplate<T>::l_nSnapshot;
  using PODTemplate<T>::l_nBlocks;
  using PODTemplate<T>::l_nmodes;
  using PODTemplate<T>::l_hilbertSp;
  using PODTemplate<T>::weightH1;
  using PODTemplate<T>::l_startTime;
  using PODTemplate<T>::l_endTime;
  using PODTemplate<T>::l_nSnapshotSimulation;
  using PODTemplate<T>::l_endTimeSimulation;
  using PODTemplate<T>::timeFolders;
  using PODTemplate<T>::name_covMatrix;
  using PODTemplate<T>::folder_covMatrix;
  using PODTemplate<T>::exist_covMatrix;
  using PODTemplate<T>::name_eigenValues;
  using PODTemplate<T>::name_eigenValuesNormalized;
  using PODTemplate<T>::name_cumEigenValues;
  using PODTemplate<T>::name_eigenVector;
  using PODTemplate<T>::folder_eigen;
  using PODTemplate<T>::exist_eigenDecomposition;
  using PODTemplate<T>::folder_spatialModes;
  using PODTemplate<T>::exist_spatialModes;
  using PODTemplate<T>::folder_temporalModes;
  using PODTemplate<T>::exist_temporalModes;
  using PODTemplate<T>::folder_temporalModesSimulation;
  using PODTemplate<T>::exist_temporalModesSimulation;
  using PODTemplate<T>::f_field;
  using PODTemplate<T>::b_centeredOrNot;
  using PODTemplate<T>::f_meanField;
  using PODTemplate<T>::folder_mean;
  using PODTemplate<T>::exist_noMean;
  using PODTemplate<T>::eigenValueseigLam;
  using PODTemplate<T>::lambda;
  using PODTemplate<T>::w_eigensolver;
  using PODTemplate<T>::i_precision;
  using PODTemplate<T>::ios_outytpe;
  using PODTemplate<T>::runTime2;

  public :
    word gradfield_name;

    // Constructor
    PODTemplateH1<T,G>(Parameters* myParameters, const word& myfield_name);

    ~PODTemplateH1<T,G>() override {};

    // First step : computation of covariance matrix
    // -------------------------------------------------------------------------------------------
    // method to compute the covariance matrix
    // add a L2 covariances of U and of grad(U)
    // param[out] : the covariance matrix
    Eigen::MatrixXd buildCovMatrix() override;

    // Precompute grad(U) and save it
    void precomputeGradients();


    // Change the result folder name and the number of snapshots 
    // for the eigenValues UnitTest
    void changeEigenFolderUnitTest(const word& name = "");

};

template class PODTemplateH1<volVectorField,volTensorField>;

template class PODTemplateH1<volScalarField,volVectorField>;
}

#endif
