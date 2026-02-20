#ifndef PODTemplate_H
#define PODTemplate_H

#include <Time.H>
#include <volFieldsFwd.H>
#include <FieldField.H>
#include "Foam2Eigen.H"

#include "ITHACAstream.H"
#include "ITHACAutilities.H"
#include "StoredParameters.H"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#include <Eigen/Eigen>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/SymEigsSolver.h>
#pragma GCC diagnostic pop

namespace ITHACAPOD
{

// Structure to store the indexes of the added block to the cov matrix
// Triangular block
struct indexTri
{
    Foam::label index_start;
    Foam::label index_end;
};

// Structure to store the indexes of the added block to the cov matrix
// Square block
struct indexSquare
{
    Foam::label index1_start;
    Foam::label index1_end;
    Foam::label index2_start;
    Foam::label index2_end;
};


// Class that compute POD (spatial modes and temporal modes)
template<typename T>
class PODTemplate
{
    protected :
        StoredParameters* m_parameters;
        // Path where OpenFOAM simulation outputs are
        Foam::word field_name;
        // Type of data used in the POD
        const Foam::fileName& casenameData;
        // Number of snapshot, also number of time step
        Foam::label l_nSnapshot;
        // number of blocks to load snapshots
        const Foam::label l_nBlocks;
        // number of modes
        Foam::label l_nmodes;
        // Hilbert space used for the POD scalar product
        // (and the associated norm to sort the modes)
        Foam::word l_hilbertSp;
        //Weight of the weighted H1 scalar product
        double weightH1;
        //Weight of the boundary conditions
        const double weightBC;
        //Patch for the boundary conditions
        const Foam::word patchBC;
        // start time of the POD decomposition
        Foam::label l_startTime;
        // end time of the POD decomposition
        Foam::label l_endTime;

        // Number of Simulation snapshot, also number of time step
        Foam::label l_nSnapshotSimulation;
        // end time of the Simulation
        const Foam::label l_endTimeSimulation;

        // Path to the directory where the snapshots are stored 
        // (usually casenameData, but can be different through 3rd arg in constructor)
        fileName snapshotsPath;

        // List linking indices and time for the snapshots data directories
        Foam::instantList timeFolders;

        // name and folder to save the covMatrix
        Foam::fileName name_covMatrix;
        Foam::word folder_covMatrix;
        bool exist_covMatrix;

        // name and folder to save eigenValues et eigenVectors
        Foam::fileName name_eigenValues;
        Foam::fileName name_eigenValuesNormalized;
        Foam::fileName name_cumEigenValues;
        Foam::fileName name_eigenVector;
        Foam::word folder_eigen;
        bool exist_eigenDecomposition;

        // folder of the spatialModes
        Foam::word folder_spatialModes;
        bool exist_spatialModes;

        // folder of the temporal Modes
        Foam::word folder_temporalModes;
        bool exist_temporalModes;

        // folder of the temporal Modes
        Foam::word folder_temporalModesSimulation;
        bool exist_temporalModesSimulation;

        // template of field readed, useful for the ITHACAstream::read_fields mlethod
        T* f_field;

        // if true, the temporal mean of snapshot is computed and subtract to snapshot fields
        // else, temporal mean is not computed and snapshot fields are not modified
        bool b_centeredOrNot;
        T* f_meanField;
        Foam::word folder_mean;
        bool exist_noMean;

        // Temporal modes energies
        Eigen::VectorXd eigenValueseigLam;
        Eigen::VectorXd lambda;

        const Foam::word w_eigensolver;
        const int i_precision;
        const std::_Ios_Fmtflags ios_outytpe;

        const Foam::Time runTime2;

        Eigen::MatrixXi inletIndex;
        bool lifting;
        Foam::PtrList<T> liftfield;

    public :
        // Constructor
        PODTemplate(Parameters* myParameters, const Foam::word& myfield_name, const word& mySnapshots_path = "default_path");

        // Destructor
        virtual ~PODTemplate();

        void define_paths();

        // Optionnal preliminary computation
        // -------------------------------------------------------------------------------------------
        // method to compute the temporal mean of the snapshots
        // it modifies the private attribute f_meanField
        // if centeredOrNot is true, the mean is computed
        // else the mean is a zero field
        void computeMeanField();

        // in case the mean has been computed
        // the mean field is appended to spatialModes
        // param[inout] : spatialModes without mean field, spatialModes with mean field
        void appendMeanfieldtoSpatialModes(Foam::PtrList<T>& spatialModes);
        // -------------------------------------------------------------------------------------------

        // if a covMatrix temporary file is found in the covMat folder
        // it is loaded and the process continues where it stopped
        void findTempFile(Eigen::MatrixXd* covMat, int* index1, int* index2);

        // temporary covMatrix filename
        Foam::word nameTempCovMatrix( int i, int j);

        // save the temporary covMatrix file
        void saveTempCovMatrix(Eigen::MatrixXd& covMatrix, int i, int j);

        // delete the temporary covMatrix file
        void deleteTempCovMatrix(int i, int j);

        // delete previous [step-N] covMatrix temp file -> keep last N covMatrix temp file
        void deletePreviousTempCovMatrix_N(int* valI, int* valJ, int i, int j, int N);

        // First step : computation of covariance matrix
        // -------------------------------------------------------------------------------------------
        // method to compute the covariance matrix
        // param[out] : the covariance matrix
        virtual Eigen::MatrixXd buildCovMatrix();

        // method to add a new triangular block to the covariance matrix (near the diagonal)
        // used in buildCovMatrix method
        // param[inout] : covMatrix, covariance matrix
        // param[in] : snapshots sample of the field
        // param[in] : indexes of the added block
        void addCovMatrixTriCoeff(Eigen::MatrixXd& covMatrix,
                                  Foam::PtrList<T>& snapshots,
                                  indexTri& indTri);

        // method to add a new square block to the covariance matrix
        // used in buildCovMatrix method
        // param[inout] : covMatrix, covariance matrix
        // param[in] : snapshots1 and snapshots2, sample of the field
        // param[in] : indexes of the added block
        void addCovMatrixSquareCoeff(Eigen::MatrixXd& covMatrix,
                                     Foam::PtrList<T>& snapshots1,
                                     Foam::PtrList<T>& snapshots2,
                                     indexSquare& indSquare);
        // -------------------------------------------------------------------------------------------


        // Second step : diagonalisation of the covariance matrix
        // -------------------------------------------------------------------------------------------
        // method to diagonalize the covariance matrix
        // param[in] covMatrix : covariance matrix (compute in the buildCovMatrix method)
        // param[out] eigenValueseig : eigen values of the diagonalisation
        // param[out] eigenVectoreig : eigen vector of the diagonalisation
        void diagonalisation(Eigen::MatrixXd& covMatrix,
                             Eigen::VectorXd& eigenValueseig, Eigen::MatrixXd& eigenVectoreig);
        // -------------------------------------------------------------------------------------------


        // Third step : computation of spatial modes (alse named topos)
        // -------------------------------------------------------------------------------------------
        // compute the spatial from the eigen values and vectors of the diagonalisation
        // param[in] eigenValueseig and eigenVectoreig : eigen values and vectors (compute in the diagonalisation method)
        // param[out] : spatialModes
        Foam::PtrList<T> computeSpatialModes(Eigen::VectorXd& eigenValueseig,
                                             Eigen::MatrixXd& eigenVectoreig);

        // // get the modes on Eigen format from the eigen values and vectors and snapshots
        // // param[in] : eigenValueseigLam and eigenVectoreig
        // // param[out] : modesEig, the Eigen formatting modes
        // // used in computeSpatialModes method
        // List<Eigen::MatrixXd> computeModesEig(Eigen::VectorXd& eigenValueseigLam, Eigen::MatrixXd& eigenVectoreig);
        // -------------------------------------------------------------------------------------------


        // Forth step : computation of temporal modes (also named chronos)
        // -------------------------------------------------------------------------------------------
        // compute the temporals from the spatial modes and snapshots
        // param[in] : spatial modes
        // param[out] : temporals modes
        Eigen::MatrixXd computeTemporalModes(Eigen::VectorXd& eigenValueseig,
                                             Eigen::MatrixXd& eigenVectoreig);
        // -------------------------------------------------------------------------------------------


        // Function which execute all the steps implemented in the previous functions
        // -------------------------------------------------------------------------------------------
        // get the temporal and spatial modes
        // param[out] : spatialModes
        // param[out] : temporalModes
        void getModes(Foam::PtrList<T>& spatialModes, Eigen::MatrixXd& temporalModes,
                      Eigen::MatrixXd& temporalModesSimulation, Eigen::MatrixXd& covMatrix);
        // -------------------------------------------------------------------------------------------

        // computation of test temporal modes
        // -------------------------------------------------------------------------------------------
        // compute the temporals from the spatial modes and snapshots
        // param[in] : spatial modes
        // param[out] : temporals modes
        Eigen::MatrixXd computeSimulationTemporalModes(Foam::PtrList<T>& spatialModes);

        // computation of test temporal modes energies
        void compute_lambda(Eigen::MatrixXd& temporalModes);

        const T& get_mean() const
        {
            return (*f_meanField);
        }

        void set_b_centeredOrNot(const Foam::label& input_b_centeredOrNot)
        {
            b_centeredOrNot = input_b_centeredOrNot;
            define_paths();
        }

        //--------------------------------------------------------------------------
        // Homogenize the snapshot matrix, it works with PtrList of volVectorField and volScalarField
        //
        // @param[in]      Lfield     The list of snapshots to be homogenized.
        void lift(Foam::PtrList<T>& snapshots);
        void lift(T& snapshot);


        // Modify parameters for the snapshots. Useful if they are different from the OpenFOAM simulation
        void set_snapFolderParams(label nSnapshot=-1, label nSnapshotSimulation=-1, label startTime=-1, label endTime=-1)
        {
            if (nSnapshot != -1) {l_nSnapshot = nSnapshot;}
            if (nSnapshotSimulation != -1) {l_nSnapshotSimulation = nSnapshotSimulation;}
            if (startTime != -1) {l_startTime = startTime;}
            if (endTime != -1) {l_endTime = endTime;}
        };

};

void computeLift(Foam::PtrList<Foam::volScalarField>& Lfield,
                 Foam::PtrList<Foam::volScalarField>& liftfield,
                 Foam::PtrList<Foam::volScalarField>& omfield, Eigen::MatrixXi inletIndex);
void computeLift(Foam::PtrList<Foam::volVectorField>& Lfield,
                 Foam::PtrList<Foam::volVectorField>& liftfield,
                 Foam::PtrList<Foam::volVectorField>& omfield, Eigen::MatrixXi inletIndex);
void computeLift(Foam::PtrList<Foam::volTensorField>& Lfield,
                 Foam::PtrList<Foam::volTensorField>& liftfield,
                 Foam::PtrList<Foam::volTensorField>& omfield, Eigen::MatrixXi inletIndex);

}

#endif
