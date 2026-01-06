
#include "PODTemplate.H"

namespace ITHACAPOD
{

template<typename T>
PODTemplate<T>::PODTemplate(Parameters* myParameters,
                            const word& myfield_name, const word& mySnapshots_path) :
    m_parameters(static_cast<StoredParameters*>(myParameters)),
    field_name(myfield_name),
    casenameData(m_parameters->get_casenameData()),
    l_nSnapshot(m_parameters->get_nSnapshots()),
    snapshotsPath(mySnapshots_path),
    l_nBlocks(m_parameters->get_nBlocks()),
    l_nmodes(m_parameters->get_nModes()[field_name]),
    l_hilbertSp(m_parameters->get_hilbertSpacePOD()[field_name]),
    weightH1(m_parameters->get_weightH1()),
    weightBC(m_parameters->get_weightPOD()),
    patchBC(m_parameters->get_patchBC()),
    l_startTime(m_parameters->get_startTime()),
    l_endTime(m_parameters->get_endTime()),
    l_nSnapshotSimulation(m_parameters->get_nSnapshotsSimulation()),
    l_endTimeSimulation(m_parameters->get_endTimeSimulation()),
    b_centeredOrNot(m_parameters->get_centeredOrNot()),
    lambda(Eigen::VectorXd::Zero(l_nmodes)),
    w_eigensolver(m_parameters->get_eigensolver()),
    i_precision(m_parameters->get_precision()),
    ios_outytpe(m_parameters->get_outytpe()),
    runTime2(Foam::Time::controlDictName, ".",
             m_parameters->get_casenameData())
{
    if (snapshotsPath == "default_path")
    {
        snapshotsPath = casenameData;
    }

    word pathProcessor("");
    if (Pstream::parRun())
    {
        pathProcessor = "processor" + name(Pstream::myProcNo()) + "/";
    }
    timeFolders = runTime2.findTimes(snapshotsPath + pathProcessor);

    l_startTime = Time::findClosestTimeIndex(timeFolders,std::stoi(runTime2.times()[l_startTime].name()));
    l_endTime = l_startTime + l_nSnapshot - 1;

    f_field = new T(
        IOobject
        (
            field_name,
            snapshotsPath + timeFolders[1].name(),
            m_parameters->get_mesh(),
            IOobject::MUST_READ
        ),
        m_parameters->get_mesh()
    );
    f_meanField = new T(f_field->name(), *f_field);
    // // Set the inlet boundaries where we have non homogeneous boundary conditions
    // inletIndex.resize(1, 2);
    // inletIndex(0, 0) = 3;
    // inletIndex(0, 1) = 2;
    // if (field_name=="U")
    // {
    //   lifting = true;
    // }
    // else
    // {
    //   lifting = false;
    // }
    // b_centeredOrNot = false;
    lifting = false;
    define_paths();
}

template<typename T>
PODTemplate<T>::~PODTemplate()
{
    delete f_field;
    delete f_meanField;
}


template<typename T>
void PODTemplate<T>::define_paths()
{
    word pathCentered("");

    if (b_centeredOrNot && !(field_name == "U"))
    {
        pathCentered = "_centered";
    }

    if (lifting)
    {
        pathCentered += "Lifted";
    }

    word pathProcessor("");
    if (Pstream::parRun())
    {
        pathProcessor = "processor" + name(Pstream::myProcNo()) + "/";
    }

    word pathHilbertSpace(m_parameters->get_pathHilbertSpace_fromHS(
                              l_hilbertSp));
    // name and folder of the covariance
    // check if the matrix was already computed
    name_covMatrix = "covMatrix" + f_field->name();
    folder_covMatrix = "./ITHACAoutput/CovMatrices" + pathHilbertSpace
                       + pathCentered + "/";
    exist_covMatrix = ITHACAutilities::check_file(folder_covMatrix + name_covMatrix
                      + ".npy");
    // name and folder of the eigen decomposition
    name_eigenValues = "Eigenvalues_" + f_field->name();
    name_eigenValuesNormalized = "EigenvaluesNormalized_" + f_field->name();
    name_cumEigenValues = "/CumEigenvalues_" + f_field->name();
    name_eigenVector = "/Eigenvector_" + f_field->name() ;
    folder_eigen = "./ITHACAoutput/EigenValuesandVector"
                   + pathCentered + "_" + std::to_string(l_nmodes) + "modes/";
    exist_eigenDecomposition = ITHACAutilities::check_file(folder_eigen +
                               name_eigenValues + ".npy")
                               && ITHACAutilities::check_file(folder_eigen + name_eigenVector + ".npy");
    // folder of the spatial modes and check if every modes were already computed
    folder_spatialModes = "./ITHACAoutput/spatialModes"
                          + pathCentered + "_" + std::to_string(l_nmodes) + "modes/";
    exist_spatialModes = ITHACAutilities::check_file(folder_spatialModes + pathProcessor + "1/" +
                         f_field->name());
    // folder of the temporal modes and check if modes were already computed
    folder_temporalModes = "./ITHACAoutput/temporalModes"
                           + pathCentered + "_" + std::to_string(l_nmodes) + "modes/";
    exist_temporalModes = ITHACAutilities::check_file(folder_temporalModes +
                          f_field->name() + ".npy");
    // folder of the temporal modes and check if modes were already computed
    folder_temporalModesSimulation = "./ITHACAoutput/temporalModesSimulation"
                                     + pathCentered + "_" + std::to_string(l_nmodes) + "modes/";
    exist_temporalModesSimulation = ITHACAutilities::check_file(
                                        folder_temporalModesSimulation + f_field->name() + ".npy");
    // folder of mean field and check if mean was already computed
    folder_mean = "./ITHACAoutput/mean/";
    exist_noMean = !ITHACAutilities::check_file(folder_mean + "/" + pathProcessor +
                    std::to_string(1) + "/" + f_field->name());
}


template<typename T>
void PODTemplate<T>::computeMeanField()
{
    if (lifting)
    {
        Info << "Loading lifting functions" << endl;
        liftfield.resize(inletIndex.rows());

        for (label k = 0; k < inletIndex.rows(); k++)
        {
            T* f_lift = new T
            (
                IOobject
                (
                    f_field->name() + "lift" + std::to_string(k),
                    snapshotsPath + timeFolders[1].name(),
                    m_parameters->get_mesh(),
                    IOobject::MUST_READ
                ),
                m_parameters->get_mesh()
            );
            liftfield.set(k, f_lift );
        }
    }

    if (b_centeredOrNot)
    {
        if (exist_noMean)
        {
            ITHACAutilities::setToZero(*f_meanField);
            Info << "Computing the mean of " << f_field->name() << " field" << endl;
            T snapshotj = *f_field;
            b_centeredOrNot = false;

            for (label j = 0; j < l_nSnapshot; j++)
            {
                // Read the j-th field
                ITHACAstream::read_snapshot(snapshotj, timeFolders[l_startTime+j].name(), snapshotsPath);
                lift(snapshotj);
                // add j-th field to meanfield
                ITHACAutilities::addFields(*f_meanField, snapshotj);
            }

            b_centeredOrNot = true;
            ITHACAutilities::multField(*f_meanField, 1 / double(l_nSnapshot));
            PtrList<T> meanExport(1);
            meanExport.set(0, new T(f_meanField->name(), *f_meanField));
            ITHACAstream::exportFields(meanExport, folder_mean, f_field->name());
        }
        else
        {
            PtrList<T> meanRead;
            Info << "Reading the mean of " << f_field->name() << " field" << endl;
            ITHACAstream::read_fields(meanRead, (*f_field), folder_mean);
            *f_meanField = meanRead[0];
        }

        double energyMean = ITHACAutilities::dot_product_L2(*f_meanField, *f_meanField);
        double energyHilbertMean = ITHACAutilities::dot_product_POD(*f_meanField,
                                   *f_meanField, l_hilbertSp);
        m_parameters->set_meanEnergy( f_meanField->name(), energyMean);
        m_parameters->set_meanEnergy( f_meanField->name() + "_" + l_hilbertSp,
                                            energyHilbertMean);
    }
    else
    {
        ITHACAutilities::setToZero(*f_meanField);
    }

    Info << endl;
}


template<typename T>
void PODTemplate<T>::appendMeanfieldtoSpatialModes(PtrList<T>& spatialModes)
{
    if (b_centeredOrNot)
    {
        spatialModes.set(l_nmodes, f_meanField);
    }
}

template<typename T>
void PODTemplate<T>::findTempFile(Eigen::MatrixXd* covMat, int* index1,
                                  int* index2)
{
    word pathTemp = name_covMatrix + "_temp_";
    DIR* dir;
    struct dirent* entry;
    dir = opendir(folder_covMatrix.c_str());
    word extTemp = ".npy";
    Info << "looking for " << pathTemp.c_str() << "*" << extTemp.c_str() << " in "
         << folder_covMatrix.c_str() << endl;
    // INDEX : index1 and index2 of oldest and most recent *_temp_*.npy file
    // INDEX[O][*] : index min (for example with *_temp_0_1.npy and *_temp_3_2.npy => INDEX[0][0]=0 and INDEX[0][1]=1 )
    // INDEX[1][*] : index max (for example with *_temp_0_1.npy and *_temp_3_2.npy => INDEX[1][0]=3 and INDEX[1][1]=2 )
    int INDEX[2][2], N_INDEX = 0;

    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            INDEX[i][j] = -1;
        }
    }

    if (dir != NULL)
    {
        while ((entry = readdir(dir)) != NULL)
        {
            char ext_name[4];
            strncpy(ext_name, "    ", 4);

            if ( strlen(entry->d_name) >= 4 )
            {
                for (int i = 0; i < 4; i++)
                {
                    ext_name[i] = entry->d_name[strlen(entry->d_name) - (4 - i)];
                }
            }

            if ((strncmp(entry->d_name, pathTemp.c_str(), pathTemp.size()) == 0)
                    && (strncmp(ext_name, extTemp.c_str(), 4) == 0))
            {
                int num1, num2;
                char endFileName[12];
                strncpy(endFileName, entry->d_name+strlen(entry->d_name)-12, 12);
                sscanf(endFileName, "%*[^0-9]%d_%d", &num1, &num2);
                *index1 = num1;
                *index2 = num2;

                if ( N_INDEX == 0 )
                {
                    INDEX[0][0] = *index1;
                    INDEX[0][1] = *index2;
                    INDEX[1][0] = *index1;
                    INDEX[1][1] = *index2;
                }
                else
                {
                    // updating INDEX from oldest *_temp_*.npy file ( -> INDEX[0][*])
                    if ( *index1 < INDEX[0][0] )
                    {
                        INDEX[0][0] = *index1;
                        INDEX[0][1] = *index2;
                    }
                    else
                    {
                        if (( *index1 == INDEX[0][0] ) && ( *index2 < INDEX[0][1] ))
                        {
                            INDEX[0][1] = *index2;
                        }
                    }

                    // updating INDEX from most recent *_temp_*.npy file ( -> INDEX[1][*])
                    if ( *index1 > INDEX[1][0] )
                    {
                        INDEX[1][0] = *index1;
                        INDEX[1][1] = *index2;
                    }
                    else
                    {
                        if (( *index1 == INDEX[1][0] ) && ( *index2 > INDEX[1][1] ))
                        {
                            INDEX[1][1] = *index2;
                        }
                    }
                }

                N_INDEX = N_INDEX + 1;
                Info << "  " << entry->d_name << " FOUND => covMat  can be updated" << endl;
            }
        }

        closedir(dir);
    }

    // restart from oldest *_temp_*.npy file or most recent *_temp_*.npy file
    if ( N_INDEX > 0 )
    {
        // updating matrix from oldest *_temp_*.npy file ( -> INDEX[0][*])
        *index1 = INDEX[0][0];
        *index2 = INDEX[0][1];
        Info << "    -> RESTART from INDEX : " << *index1 << " " << *index2 <<
             " (oldest *_temp_*.npy file)" << endl;
        char str1[10];
        sprintf(str1, "%d", *index1);
        char str2[10];
        sprintf(str2, "%d", *index2);
        word suffix = static_cast<word>(str1) + "_" + static_cast<word>(str2);
        cnpy::load(*covMat, folder_covMatrix + pathTemp + suffix + extTemp);
        Info << "       with covMat size=(" << covMat->rows() << ", " << covMat->cols()
             << ")" << endl;
    }
}

template<typename T>
word PODTemplate<T>::nameTempCovMatrix( int i, int j)
{
    char str1[10];
    sprintf(str1, "%d", i);
    char str2[10];
    sprintf(str2, "%d", j);
    word suffix = "_temp_" + static_cast<word>(str1) + "_" + static_cast<word>
                  (str2);
    word filename = folder_covMatrix + name_covMatrix + suffix + ".npy";
    return filename;
}

template<typename T>
void PODTemplate<T>::saveTempCovMatrix(Eigen::MatrixXd& covMatrix, int i, int j)
{
    word filename = nameTempCovMatrix(i, j);
    if (Pstream::master())
    {
        cnpy::save(covMatrix, filename);
    }
}

template<typename T>
void PODTemplate<T>::deleteTempCovMatrix(int i, int j)
{
    word filename = nameTempCovMatrix(i, j);
    std::ifstream fileStr(filename.c_str());

    if (fileStr)
    {
        remove(filename.c_str());
    }
}

template<typename T>
void PODTemplate<T>::deletePreviousTempCovMatrix_N(int* valI, int* valJ, int i,
        int j, int N)
{
    // NUM = (sum_{n=1}^{n=i} n) = ( (m+(m+i))+(m+1+(m+i)-1)+...+((m+i)+m))/2 = (2*m*i+i+i*i)/2
    //  => NUM(m=0) = (i+i*i)/2
    //     NUM(m=0)+j+1 = (sum_{n=1}^{n=i} n) +j+1 = (1+i)*i/2 +j
    int NUM = (1 + i) * i / 2 + j;
    int dNUM = NUM - N;
    *valI = 0;
    int S = 0;

    while ( dNUM >= S )
    {
        *valI = *valI + 1;
        S = (1 + *valI)** valI / 2;
    }

    *valI = *valI - 1;
    S = (1 + *valI)** valI / 2;
    *valJ = dNUM - S;

    if ( *valI != i )
    {
        *valJ = *valJ - 1;
    }

    if ( *valJ == -1 )
    {
        *valI = *valI - 1;
        *valJ = *valJ + 1;
    }

    if ((*valI >= 0) && (*valJ >= 0))
    {
        deleteTempCovMatrix(*valI, *valJ);
    }

    *valJ = j;
}

template<typename T>
Eigen::MatrixXd PODTemplate<T>::buildCovMatrix()
{
    // the covariance matrix is initialized to a unreal value
    double CovUnrealValue = 999999.999999;
    Eigen::MatrixXd covMatrix;
    covMatrix.setConstant(l_nSnapshot, l_nSnapshot, CovUnrealValue);
    bool exist_covMatrix_bin = ITHACAutilities::check_file(folder_covMatrix +
                               name_covMatrix);
    int valI, valJ;
    //initializing value of etapeI and etapeJ : no temp file found
    int etapeI, etapeJ;
    etapeI = -1;
    etapeJ = -1;
    // create folder
    mkDir(folder_covMatrix);

    // in case the cov matrix was not computed
    if (!exist_covMatrix && !exist_covMatrix_bin)
    {
        // updating covMatrix from previous step (when a *_temp_* file exist)
        findTempFile(&covMatrix, &etapeI, &etapeJ);
        Info << "Computing the covariance matrix of the " << f_field->name() << " field"
             << endl;
        // the size is cov matrix l_nSnapshot, it has been subdivised in l_nBlocks
        label q = l_nSnapshot / l_nBlocks;
        label r = l_nSnapshot % l_nBlocks;
        Info << "q = " << q << endl;
        Info << "r = " << r << endl;
        Info << "nSnapshot = " << l_nSnapshot << endl;
        Info << "nBlocks = " << l_nBlocks << endl;
        Info << endl;
        PtrList<T> snapshots;

        // If a temp file is found load it and go to the next step
        if (etapeI != -1)
        {
            if (etapeJ == etapeI - 1)
            {
                etapeJ = 0;
                etapeI += 1;
            }
            else
            {
                etapeJ += 1;
            }
        }
        else
        {
            etapeI = 0;
            etapeJ = 0;
        }

        // Number of previous tempory CovMatrix to keep in directory
        // (parameter used by [deletePreviousTempCovMatrix_N])
        int N_previous_temp_Mat = 3;

        for (label i = etapeI; i < l_nBlocks; i++)
        {
            ITHACAstream::read_fields(snapshots, (*f_field), snapshotsPath,
                                      l_startTime - 2 + i * q, q);
            lift(snapshots);
            indexTri indTri;
            indTri.index_start = i * q;
            indTri.index_end = (i + 1) * q;
            addCovMatrixTriCoeff(covMatrix, snapshots, indTri);
            PtrList<T> snapshots2;
            int valI = etapeI, valJ = etapeJ;

            for (label j = etapeJ; j < i; j++)
            {
                ITHACAstream::read_fields(snapshots2, (*f_field), snapshotsPath,
                                          l_startTime - 2 + j * q, q);
                lift(snapshots2);
                indexSquare indSquare;
                indSquare.index1_start = i * q;
                indSquare.index1_end = (i + 1) * q;
                indSquare.index2_start = j * q;
                indSquare.index2_end = (j + 1) * q;
                addCovMatrixSquareCoeff(covMatrix, snapshots, snapshots2, indSquare);
                // Clear the pointer to snapshots2
                snapshots2.clear();
                // Save cove matrix temp file
                saveTempCovMatrix(covMatrix, i, j);
                // Delete previous covMatrix temp file after saving the current one
                deletePreviousTempCovMatrix_N(&valI, &valJ, i, j, N_previous_temp_Mat);
            }

            // Clear the pointer to snapshots
            snapshots.clear();
            valI = i;
            etapeJ = 0;
        }

        if (r != 0)
        {
            PtrList<T> snapshotsEnd;
            ITHACAstream::read_fields(snapshotsEnd, (*f_field), snapshotsPath,
                                      l_startTime - 2 + l_nBlocks * q, r);
            lift(snapshotsEnd);
            indexTri indTri;
            indTri.index_start = l_nBlocks * q;
            indTri.index_end = l_nSnapshot;
            addCovMatrixTriCoeff(covMatrix, snapshotsEnd, indTri);
            PtrList<T> snapshotsEnd2;

            for (label j = 0; j < l_nBlocks; j++)
            {
                ITHACAstream::read_fields(snapshotsEnd2, (*f_field), snapshotsPath,
                                          l_startTime - 2 + j * q, q);
                lift(snapshotsEnd2);
                indexSquare indSquare;
                indSquare.index1_start = l_nBlocks * q;
                indSquare.index1_end =  l_nSnapshot;
                indSquare.index2_start = j * q;
                indSquare.index2_end = (j + 1) * q;
                addCovMatrixSquareCoeff(covMatrix, snapshotsEnd, snapshotsEnd2, indSquare);
                snapshotsEnd2.clear();
                // Save cove matrix temp file
                saveTempCovMatrix(covMatrix, l_nBlocks, j);
                // Delete previous covMatrix temp file after saving the current one
                deletePreviousTempCovMatrix_N(&valI, &valJ, l_nBlocks, j, N_previous_temp_Mat);
            }

            snapshotsEnd.clear();
        }
        // covMatrix is symetric, the lower part is used to build the upper part
        covMatrix = covMatrix.selfadjointView<Eigen::Lower>();
        if (Pstream::master())
        {
            cnpy::save(covMatrix, folder_covMatrix + name_covMatrix + ".npy");

        }
        // Delete previous covMatrix temp file after saving the current one
        if (r == 0)
        {
            deletePreviousTempCovMatrix_N(&valI, &valJ, l_nBlocks - 1, l_nBlocks - 1, 1);
        }

        if (r != 0)
        {
            deletePreviousTempCovMatrix_N(&valI, &valJ, l_nBlocks, l_nBlocks, 1);
        }
    }
    // in case the cov matrix was already computed, it reads in the hard disk
    else if (exist_covMatrix)
    {
        Info << "Reading the covariance matrix of the " << f_field->name() << " field"
             << endl;
        cnpy::load(covMatrix, folder_covMatrix + name_covMatrix + ".npy");
    }
    // in case the cov matrix was already computed, it reads in the hard disk
    else
    {
        Info << "Reading (binary) the covariance matrix of the " << f_field->name() <<
             " field" << endl;
        ITHACAstream::ReadDenseMatrix(covMatrix, folder_covMatrix, name_covMatrix);
        if (Pstream::master())
        {
            cnpy::save(covMatrix, folder_covMatrix + name_covMatrix + ".npy");
        }
    }

    // looking for CovUnrealValue
    int covMatrixOK = 1;
    int NbCovUnrealValue = 0;

    for (int i = 0; i < l_nSnapshot; i++)
    {
        for (int j = i; j < l_nSnapshot; j++)
        {
            if (covMatrix(i, j) == CovUnrealValue)
            {
                covMatrixOK = 0;
                NbCovUnrealValue += 1;
            }
        }
    }

    if (covMatrixOK == 0)
    {
        Info << "\n!!! OUPS !!! Unreal value [" << CovUnrealValue << "] found " <<
             NbCovUnrealValue << " times in triangular up part of " << name_covMatrix <<
             " !!!\n" << endl;
        abort();
    }

    double varyingEnergy = covMatrix.trace() / l_nSnapshot;

    if (l_hilbertSp == "L2" || l_hilbertSp == "dL2")
    {
        m_parameters->set_varyingEnergy( f_field->name(), varyingEnergy);
    }

    Info << "Total varying " << l_hilbertSp << " energy for " << field_name << " : "
         << varyingEnergy << endl;
    Info << endl;
    return covMatrix;
}


template<typename T>
void PODTemplate<T>::addCovMatrixTriCoeff(Eigen::MatrixXd& covMatrix,
        PtrList<T>& snapshots, indexTri& indTri)
{
    Info << "Adding the triangular block [" << indTri.index_start << ":" <<
         indTri.index_end - 1 << "]x["
         << indTri.index_start << ":" << indTri.index_end - 1 <<
         "] to the covariance matrix" << endl;
    Eigen::MatrixXd covMatrixTemp(ITHACAutilities::dot_product_POD(snapshots,
                                  snapshots, l_hilbertSp, weightBC, patchBC));

    for (label i = 0; i < indTri.index_end - indTri.index_start; i++)
    {
        for (label j = 0; j <= i; j++)
        {
            covMatrix(i + indTri.index_start, j + indTri.index_start) =
                covMatrixTemp(i, j);
        }
    }

    Info << endl;
}


template<typename T>
void PODTemplate<T>::addCovMatrixSquareCoeff(Eigen::MatrixXd& covMatrix,
        PtrList<T>& snapshots1,
        PtrList<T>& snapshots2,
        indexSquare& indSquare)
{
    Info << "Adding the square block [" << indSquare.index1_start << ":" <<
         indSquare.index1_end - 1 << "]x["
         << indSquare.index2_start << ":" << indSquare.index2_end - 1 <<
         "] to the covariance matrix" << endl;
    Eigen::MatrixXd covMatrixTemp(ITHACAutilities::dot_product_POD(snapshots1,
                                  snapshots2, l_hilbertSp, weightBC, patchBC));

    for (label i = 0; i < indSquare.index1_end - indSquare.index1_start; i++)
    {
        for (label j = 0; j < indSquare.index2_end - indSquare.index2_start; j++)
        {
            covMatrix(i + indSquare.index1_start, j + indSquare.index2_start) =
                covMatrixTemp(i, j);
        }
    }

    Info << endl;
}


template<typename T>
void PODTemplate<T>::diagonalisation(Eigen::MatrixXd& covMatrix,
                                     Eigen::VectorXd& eigenValueseig, Eigen::MatrixXd& eigenVectoreig)
{
    if (!exist_eigenDecomposition)
    {
        mkDir(folder_eigen);
        Info << "Performing the eigen decomposition" << endl;

        if (w_eigensolver == "spectra")
        {
            Spectra::DenseSymMatProd<double> op(covMatrix);
            Spectra::SymEigsSolver<Spectra::DenseSymMatProd<double>> es(op, l_nmodes, l_nSnapshot);
            std::cout << "Using Spectra EigenSolver " << std::endl;
            es.init();
            es.compute(Spectra::SortRule::LargestAlge);
            M_Assert(es.info() == Spectra::CompInfo::Successful,
                     "The Eigenvalue Decomposition did not succeed");
            eigenVectoreig = es.eigenvectors().real();
            eigenValueseig = es.eigenvalues().real();
        }
        else if (w_eigensolver == "eigen")
        {
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esEg;
            std::cout << "Using Eigen EigenSolver " << std::endl;
            esEg.compute(covMatrix);
            M_Assert(esEg.info() == Eigen::Success,
                     "The Eigenvalue Decomposition did not succeed");
            eigenVectoreig = esEg.eigenvectors().real().rowwise().reverse().leftCols(
                                 l_nmodes);
            eigenValueseig = esEg.eigenvalues().real().reverse().head(l_nmodes);
        }
        else
        {
        }

        // Compute the norm of each Modes
        eigenValueseigLam = eigenValueseig.real().array().abs().sqrt();
        // save the eigen values

        if (Pstream::master())
        {
            cnpy::save(eigenValueseig, folder_eigen + name_eigenValues + ".npy");
            // save the eigen vectors
            cnpy::save(eigenVectoreig,
                    folder_eigen + "/Eigenvector_" + f_field->name() + ".npy");
            // save the norm of each Modes
            cnpy::save(eigenValueseigLam,
                    folder_eigen + "/EigenvectorLambda_" + f_field->name() + ".npy");
            Eigen::VectorXd eigenValueseigNormalized = eigenValueseig /
                    eigenValueseig.sum();
            Eigen::VectorXd cumEigenValues(eigenValueseigNormalized);

            for (int j = 1; j < cumEigenValues.size(); ++j)
            {
                cumEigenValues(j) += cumEigenValues(j - 1);
            }

            // save eigen values normalized
            cnpy::save(eigenValueseigNormalized,
                    folder_eigen + name_eigenValuesNormalized + ".npy");
            // save the cumulated eigen values
            cnpy::save(cumEigenValues, folder_eigen + name_cumEigenValues + ".npy");
        }
    }
    else
    {
        // in case the eigen decomposition was already performed
        // load eigen values
        cnpy::load(eigenValueseig, folder_eigen + name_eigenValues + ".npy");
        // load eigen vectors
        cnpy::load(eigenVectoreig,
                   folder_eigen + "/Eigenvector_" + f_field->name() + ".npy");
        // load the norm of each Modes
        cnpy::load(eigenValueseigLam,
                   folder_eigen + "/EigenvectorLambda_" + f_field->name() + ".npy");
    }

    double resolvedVaryingEnergy = eigenValueseig.sum() / l_nSnapshot;

    if (l_hilbertSp == "L2")
    {
        m_parameters->set_resolvedVaryingEnergy( f_field->name(),
                resolvedVaryingEnergy);
    }
    else
    {
        m_parameters->set_resolvedVaryingEnergy( f_field->name() + "_" +
                l_hilbertSp, resolvedVaryingEnergy);
    }

    Info << "% of varying " << l_hilbertSp << " energy captures by the modes for "
         << field_name
         << " : " << 100.0 * eigenValueseig.sum() / covMatrix.trace() << "%" << endl;
    Info << endl;
}


template<typename T>
PtrList<T> PODTemplate<T>::computeSpatialModes(Eigen::VectorXd& eigenValueseig,
        Eigen::MatrixXd& eigenVectoreig)
{
    PtrList<T> spatialModes;

    // In case the modes were not computed
    if (!exist_spatialModes)
    {
        Info << "Computing the spatial modes of the " << f_field->name() << " field" <<
             endl;
        spatialModes.resize(l_nmodes);

        for (int k = 0; k < l_nmodes; k++)
        {
            spatialModes.set(k, new T(f_field->name(), *f_field) );
            ITHACAutilities::setToZero(spatialModes[k]);
        }

        for (label j = 0; j < l_nSnapshot; j++)
        {
            T snapshotj = *f_field;
            ITHACAstream::read_snapshot(snapshotj, timeFolders[l_startTime+j].name(), snapshotsPath);

            if ((m_parameters->get_DEIMInterpolatedField() == "nut" 
                || ITHACAutilities::containsSubstring(m_parameters->get_DEIMInterpolatedField(), "reducedNut")) 
                && l_hilbertSp == "dL2")
            {
                ITHACAutilities::multField(snapshotj, m_parameters->get_deltaWeight());
            }

            lift(snapshotj);

            for (label k = 0; k < l_nmodes; k++)
            {
                ITHACAutilities::addFields(spatialModes[k], snapshotj, eigenVectoreig(j, k));
            }
        }

        for (int k = 0; k < l_nmodes; k++)
        {
            ITHACAutilities::multField(spatialModes[k], 1 / eigenValueseigLam(k));
        }

        mkDir(folder_spatialModes);
        ITHACAstream::exportFields(spatialModes, folder_spatialModes, f_field->name());
    }
    // in the case the spatial modes was already computed, it reads in the hard disk
    else
    {
        Info << "Reading the spatial modes of the " << f_field->name() << " field" <<
             endl;
        ITHACAstream::read_fields(spatialModes, (*f_field), folder_spatialModes);
    }

    Info << endl;
    return spatialModes;
}


template<typename T>
Eigen::MatrixXd PODTemplate<T>::computeTemporalModes(Eigen::VectorXd&
        eigenValueseig, Eigen::MatrixXd& eigenVectoreig)
{
    Eigen::MatrixXd temporalModes(l_nSnapshot, l_nmodes);

    if (!exist_temporalModes)
    {
        Info << "Computing the temporal modes of the " << f_field->name() << " field" <<
             endl;
        mkDir( folder_temporalModes );

        for (int i = 0; i < l_nmodes; i++)
        {
            temporalModes.col(i) = eigenVectoreig.col(i) * eigenValueseigLam(i);
        }

        if (Pstream::master())
        {
            cnpy::save(temporalModes, folder_temporalModes + f_field->name() + ".npy");
        }
    }
    else
    {
        Info << "Reading the temporal modes of the " << f_field->name() << " field" <<
             endl;
        cnpy::load(temporalModes, folder_temporalModes + f_field->name() + ".npy");
    }

    Info << endl;
    return temporalModes;
}

template<typename T>
void PODTemplate<T>::getModes(PtrList<T>& spatialModes,
                              Eigen::MatrixXd& temporalModes, Eigen::MatrixXd& temporalModesSimulation,
                              Eigen::MatrixXd& covMatrix)
{
    Info << "-------------------------------------------------------------------------------------------"
         << endl;
    Info << "The POD is performing with the following parameters :" << endl;
    Info << "Field : " << f_field->name() << endl;
    Info << "Number of modes : " << l_nmodes << endl;
    Info << "POD Hilbert space : " << l_hilbertSp << endl;
    Info << "Weight for the boundary conditions (0 for L2 POD Hilbert space) : " <<
         weightBC << endl;
    Info << "Patch for the boundary conditions(inlet by default) used in L2wBC : "
         << patchBC << endl;
    Info << "start time : " << l_startTime << endl;
    Info << "End time : " << l_endTime << endl;
    Info << "Number of snapshots : " << l_nSnapshot << endl;
    Info << "Number of test snapshots : " << l_nSnapshotSimulation << endl;
    Info << "Number of blocks : " << l_nBlocks << endl;
    Info << "Centered datas or not : " << b_centeredOrNot <<
         " (1 centered, 0 not centered)" << endl;
    Info << "Name of eigensolver used : " << w_eigensolver << endl;
    Info << "Results folder : " << "ITHACAoutput" << endl;
    computeMeanField();
    covMatrix = this->buildCovMatrix();
    // Reduction (parallelization)
    Eigen::VectorXd eigenValueseig = Eigen::VectorXd::Zero(l_nmodes);
    Eigen::MatrixXd eigenVectoreig = Eigen::MatrixXd::Zero(l_nSnapshot, l_nmodes);
    diagonalisation(covMatrix, eigenValueseig, eigenVectoreig);
    spatialModes = computeSpatialModes(eigenValueseig, eigenVectoreig);
    temporalModes = computeTemporalModes(eigenValueseig, eigenVectoreig);
    compute_lambda(temporalModes);
    temporalModesSimulation = computeSimulationTemporalModes(spatialModes);
    // Reduction (parallelization)

    if (field_name == "U")
    {
        m_parameters->set_eigenValues_U(eigenValueseig);
        m_parameters->set_lambda(lambda);
    }

    Info << "-------------------------------------------------------------------------------------------"
         << endl;
}

template<typename T>
Eigen::MatrixXd PODTemplate<T>::computeSimulationTemporalModes(
    PtrList<T>& f_spatialModes)
{
    Eigen::MatrixXd temporalModesSimulation(l_nSnapshotSimulation, l_nmodes);

    if (!exist_temporalModesSimulation)
    {
        Info << "Computing the Simulation temporal modes of the " << f_field->name() <<
             " field" << endl;
        mkDir( folder_temporalModesSimulation );
        label l_startTimeSimulation(l_endTime);

        for (label j = 0; j < l_nSnapshotSimulation; j++)
        {
            T snapshotj = *f_field;
            ITHACAstream::read_snapshot(snapshotj, timeFolders[l_startTimeSimulation+j].name(), snapshotsPath);

            if ((m_parameters->get_DEIMInterpolatedField() == "nut" 
                 || ITHACAutilities::containsSubstring(m_parameters->get_DEIMInterpolatedField(), "reducedNut")) 
                 && l_hilbertSp == "dL2")
            {
                ITHACAutilities::multField(snapshotj, m_parameters->get_deltaWeight());
            }

            lift(snapshotj);

            for (label i = 0; i < l_nmodes; i++)
            {
                temporalModesSimulation(j, i) = ITHACAutilities::dot_product_POD(snapshotj,
                                                f_spatialModes[i], l_hilbertSp, weightBC, patchBC, weightH1);
            }
        }

        if (Pstream::master())
        {
            cnpy::save(temporalModesSimulation,
                   folder_temporalModesSimulation + f_field->name() + ".npy");
        }
    }
    else
    {
        Info << "Reading the Simulation temporal modes of the " << f_field->name() <<
             " field" << endl;
        cnpy::load(temporalModesSimulation,
                   folder_temporalModesSimulation + f_field->name() + ".npy");
    }

    Info << endl;
    return temporalModesSimulation;
}

template<typename T>
void PODTemplate<T>::compute_lambda(Eigen::MatrixXd& temporalModes)
{
    for (int p = 0; p < l_nmodes; p++)
    {
        for (int i = 0; i < l_nSnapshot; i++)
        {
            lambda(p) += temporalModes(i, p) * temporalModes(i, p);
        }

        lambda(p) /= l_nSnapshot;
    }
}


template<typename T>
void PODTemplate<T>::lift(PtrList<T>& snapshots)
{
    if (lifting)
    {
        PtrList<T> omfield;
        computeLift(snapshots, liftfield, omfield, inletIndex);
        snapshots = omfield;
    }

    if (b_centeredOrNot)
    {
        for (label k = 0; k < snapshots.size(); k++)
        {
            ITHACAutilities::subtractFields(snapshots[k], *f_meanField);
        }
    }
}

template<typename T>
void PODTemplate<T>::lift(T& snapshot)
{
    PtrList<T> snapshots(1);
    snapshots.set(0, new T(snapshot.name(), snapshot));
    lift(snapshots);
    snapshot = snapshots[0];
    List<Eigen::VectorXd> snapshotsBC = Foam2Eigen::field2EigenBC(snapshots[0]);

    for (int k = 0; k < snapshotsBC.size(); k++)
    {
        ITHACAutilities::assignBC(snapshot, k, snapshotsBC[k]);
    }
}

// Specialisation
template PODTemplate<volTensorField>::PODTemplate(Parameters*
        myParameters, const word& myfield_name, const word& mySnapshots_path);
// template PODTemplate<volTensorField>::PODTemplate(m_parameters* myParameters, const word& myfield_name);
// template PODTemplate<volTensorField>::PODTemplate(Parameters* myParameters, const word& myfield_name, bool b_centeredOrNot);
template PODTemplate<volTensorField>::~PODTemplate();
template void PODTemplate<volTensorField>::define_paths();
template void PODTemplate<volTensorField>::computeMeanField();
template void PODTemplate<volTensorField>::appendMeanfieldtoSpatialModes(
    PtrList<volTensorField>& spatialModes);
template void PODTemplate<volTensorField>::findTempFile(Eigen::MatrixXd* covMat,
        int* index1, int* index2);
template void PODTemplate<volTensorField>::deleteTempCovMatrix(int i, int j);

template word PODTemplate<volTensorField>::nameTempCovMatrix(int i, int j);
template void PODTemplate<volTensorField>::saveTempCovMatrix(
    Eigen::MatrixXd& covMatrix, int i, int j);
template void PODTemplate<volTensorField>::deletePreviousTempCovMatrix_N(
    int* valI, int* valJ, int i, int j, int N);

template Eigen::MatrixXd PODTemplate<volTensorField>::buildCovMatrix();
template void PODTemplate<volTensorField>::addCovMatrixTriCoeff(
    Eigen::MatrixXd& covMatrix,
    PtrList<volTensorField>& snapshots,
    indexTri& indTri);
template void PODTemplate<volTensorField>::addCovMatrixSquareCoeff(
    Eigen::MatrixXd& covMatrix,
    PtrList<volTensorField>& snapshots1,
    PtrList<volTensorField>& snapshots2,
    indexSquare& indSquare);
template void PODTemplate<volTensorField>::diagonalisation(
    Eigen::MatrixXd& covMatrix, Eigen::VectorXd& eigenValueseig,
    Eigen::MatrixXd& eigenVectoreig);
template PtrList<volTensorField>
PODTemplate<volTensorField>::computeSpatialModes(Eigen::VectorXd&
        eigenValueseig, Eigen::MatrixXd& eigenVectoreig);
// template List<Eigen::MatrixXd> PODTemplate<volTensorField>::computeModesEig(Eigen::VectorXd& eigenValueseigLam, Eigen::MatrixXd& eigenVectoreig);
template Eigen::MatrixXd PODTemplate<volTensorField>::computeTemporalModes(
    Eigen::VectorXd& eigenValueseig, Eigen::MatrixXd& eigenVectoreig);
template void PODTemplate<volTensorField>::getModes(PtrList<volTensorField>&
        spatialModes, Eigen::MatrixXd& temporalModes,
        Eigen::MatrixXd& temporalModesSimulation, Eigen::MatrixXd& covMatrix);
template Eigen::MatrixXd
PODTemplate<volTensorField>::computeSimulationTemporalModes(
    PtrList<volTensorField>& spatialModes);
template void PODTemplate<volTensorField>::compute_lambda(
    Eigen::MatrixXd& temporalModes);
template void PODTemplate<volTensorField>::lift(PtrList<volTensorField>&
        snapshots);
template void PODTemplate<volTensorField>::lift(volTensorField& snapshot);

// Specialisation
template PODTemplate<volVectorField>::PODTemplate(Parameters*
        myParameters, const word& myfield_name, const word& mySnapshots_path);
// template PODTemplate<volVectorField>::PODTemplate(m_parameters* myParameters, const word& myfield_name);
// template PODTemplate<volVectorField>::PODTemplate(Parameters* myParameters, const word& myfield_name, bool b_centeredOrNot);
template PODTemplate<volVectorField>::~PODTemplate();
template void PODTemplate<volVectorField>::define_paths();
template void PODTemplate<volVectorField>::computeMeanField();
template void PODTemplate<volVectorField>::appendMeanfieldtoSpatialModes(
    PtrList<volVectorField>& spatialModes);
template Eigen::MatrixXd PODTemplate<volVectorField>::buildCovMatrix();
template void PODTemplate<volVectorField>::addCovMatrixTriCoeff(
    Eigen::MatrixXd& covMatrix,
    PtrList<volVectorField>& snapshots,
    indexTri& indTri);
template void PODTemplate<volVectorField>::addCovMatrixSquareCoeff(
    Eigen::MatrixXd& covMatrix,
    PtrList<volVectorField>& snapshots1,
    PtrList<volVectorField>& snapshots2,
    indexSquare& indSquare);
template void PODTemplate<volVectorField>::diagonalisation(
    Eigen::MatrixXd& covMatrix, Eigen::VectorXd& eigenValueseig,
    Eigen::MatrixXd& eigenVectoreig);
template PtrList<volVectorField>
PODTemplate<volVectorField>::computeSpatialModes(Eigen::VectorXd&
        eigenValueseig, Eigen::MatrixXd& eigenVectoreig);
// template List<Eigen::MatrixXd> PODTemplate<volVectorField>::computeModesEig(Eigen::VectorXd& eigenValueseigLam, Eigen::MatrixXd& eigenVectoreig);
template Eigen::MatrixXd PODTemplate<volVectorField>::computeTemporalModes(
    Eigen::VectorXd& eigenValueseig, Eigen::MatrixXd& eigenVectoreig);
template void PODTemplate<volVectorField>::getModes(PtrList<volVectorField>&
        spatialModes, Eigen::MatrixXd& temporalModes,
        Eigen::MatrixXd& temporalModesSimulation, Eigen::MatrixXd& covMatrix);
template Eigen::MatrixXd
PODTemplate<volVectorField>::computeSimulationTemporalModes(
    PtrList<volVectorField>& spatialModes);
template void PODTemplate<volVectorField>::compute_lambda(
    Eigen::MatrixXd& temporalModes);
template void PODTemplate<volVectorField>::lift(PtrList<volVectorField>&
        snapshots);
template void PODTemplate<volVectorField>::lift(volVectorField& snapshot);

// Specialisation
template PODTemplate<volScalarField>::PODTemplate(Parameters*
        myParameters, const word& myfield_name, const word& mySnapshots_path);
// template PODTemplate<volScalarField>::PODTemplate(m_parameters* myParameters, const word& myfield_name);
// template PODTemplate<volScalarField>::PODTemplate(Parameters* myParameters, const word& myfield_name, bool b_centeredOrNot);
template PODTemplate<volScalarField>::~PODTemplate();
template void PODTemplate<volScalarField>::define_paths();
template void PODTemplate<volScalarField>::computeMeanField();
template void PODTemplate<volScalarField>::appendMeanfieldtoSpatialModes(
    PtrList<volScalarField>& spatialModes);
template Eigen::MatrixXd PODTemplate<volScalarField>::buildCovMatrix();
template void PODTemplate<volScalarField>::addCovMatrixTriCoeff(
    Eigen::MatrixXd& covMatrix,
    PtrList<volScalarField>& snapshots,
    indexTri& indTri);
template void PODTemplate<volScalarField>::addCovMatrixSquareCoeff(
    Eigen::MatrixXd& covMatrix,
    PtrList<volScalarField>& snapshots1,
    PtrList<volScalarField>& snapshots2,
    indexSquare& indSquare);
template void PODTemplate<volScalarField>::diagonalisation(
    Eigen::MatrixXd& covMatrix, Eigen::VectorXd& eigenValueseig,
    Eigen::MatrixXd& eigenVectoreig);
template PtrList<volScalarField>
PODTemplate<volScalarField>::computeSpatialModes(Eigen::VectorXd&
        eigenValueseig, Eigen::MatrixXd& eigenVectoreig);
// template List<Eigen::MatrixXd> PODTemplate<volScalarField>::computeModesEig(Eigen::VectorXd& eigenValueseigLam, Eigen::MatrixXd& eigenVectoreig);
template Eigen::MatrixXd PODTemplate<volScalarField>::computeTemporalModes(
    Eigen::VectorXd& eigenValueseig, Eigen::MatrixXd& eigenVectoreig);
template void PODTemplate<volScalarField>::getModes(PtrList<volScalarField>&
        spatialModes, Eigen::MatrixXd& temporalModes,
        Eigen::MatrixXd& temporalModesSimulation, Eigen::MatrixXd& covMatrix);
template Eigen::MatrixXd
PODTemplate<volScalarField>::computeSimulationTemporalModes(
    PtrList<volScalarField>& spatialModes);
template void PODTemplate<volScalarField>::compute_lambda(
    Eigen::MatrixXd& temporalModes);
template void PODTemplate<volScalarField>::lift(PtrList<volScalarField>&
        snapshots);
template void PODTemplate<volScalarField>::lift(volScalarField& snapshot);



void computeLift(PtrList<volTensorField>& Lfield,
                 PtrList<volTensorField>& liftfield, PtrList<volTensorField>& omfield,
                 Eigen::MatrixXi inletIndex)
{
    scalar u_bc;
    scalar u_lf;
    scalar area;

    for (label k = 0; k < inletIndex.rows(); k++)
    {
        label p = inletIndex(k, 0);
        label l = inletIndex(k, 1);
        area = gSum(Lfield[0].mesh().magSf().boundaryField()[p]);
        u_lf = gSum(liftfield[k].mesh().magSf().boundaryField()[p] *
                    liftfield[k].boundaryField()[p]).component(l) / area;
        M_Assert(std::abs(u_lf) > 1e-5,
                 "The lift cannot be computed. Please, check your inletIndex definition");

        for (label j = 0; j < Lfield.size(); j++)
        {
            if (k == 0)
            {
                u_bc = gSum(Lfield[j].mesh().magSf().boundaryField()[p] *
                            Lfield[j].boundaryField()[p]).component(l) / area;
                volTensorField C(Lfield[0].name(), Lfield[j] - liftfield[k]*u_bc / u_lf);
                omfield.append(C.clone());
            }
            else
            {
                u_bc = gSum(omfield[j].mesh().magSf().boundaryField()[p] *
                            omfield[j].boundaryField()[p]).component(l) / area;
                volTensorField C(Lfield[0].name(), omfield[j] - liftfield[k]*u_bc / u_lf);
                omfield.set(j, C.clone());
            }
        }
    }
}

void computeLift(PtrList<volVectorField>& Lfield,
                 PtrList<volVectorField>& liftfield, PtrList<volVectorField>& omfield,
                 Eigen::MatrixXi inletIndex)
{
    scalar u_bc;
    scalar u_lf;
    scalar area;

    for (label k = 0; k < inletIndex.rows(); k++)
    {
        label p = inletIndex(k, 0);
        label l = inletIndex(k, 1);
        area = gSum(Lfield[0].mesh().magSf().boundaryField()[p]);
        u_lf = gSum(liftfield[k].mesh().magSf().boundaryField()[p] *
                    liftfield[k].boundaryField()[p]).component(l) / area;
        M_Assert(std::abs(u_lf) > 1e-5,
                 "The lift cannot be computed. Please, check your inletIndex definition");

        for (label j = 0; j < Lfield.size(); j++)
        {
            if (k == 0)
            {
                u_bc = gSum(Lfield[j].mesh().magSf().boundaryField()[p] *
                            Lfield[j].boundaryField()[p]).component(l) / area;
                volVectorField C(Lfield[0].name(), Lfield[j] - liftfield[k]*u_bc / u_lf);
                omfield.append(C.clone());
            }
            else
            {
                u_bc = gSum(omfield[j].mesh().magSf().boundaryField()[p] *
                            omfield[j].boundaryField()[p]).component(l) / area;
                volVectorField C(Lfield[0].name(), omfield[j] - liftfield[k]*u_bc / u_lf);
                omfield.set(j, C.clone());
            }
        }
    }
}

void computeLift(PtrList<volScalarField>& Lfield,
                 PtrList<volScalarField>& liftfield, PtrList<volScalarField>& omfield,
                 Eigen::MatrixXi inletIndex)
{
    scalar t_bc;
    scalar t_lf;
    scalar area;

    for (label k = 0; k < inletIndex.rows(); k++)
    {
        label p = inletIndex(k, 0);
        area = gSum(Lfield[0].mesh().magSf().boundaryField()[p]);
        t_lf = gSum(liftfield[k].mesh().magSf().boundaryField()[p] *
                    liftfield[k].boundaryField()[p]) / area;

        for (label j = 0; j < Lfield.size(); j++)
        {
            if (k == 0)
            {
                t_bc = gSum(Lfield[j].mesh().magSf().boundaryField()[p] *
                            Lfield[j].boundaryField()[p]) / area;
                volScalarField C(Lfield[0].name(), Lfield[j] - liftfield[k]*t_bc / t_lf);
                omfield.append(C.clone());
            }
            else
            {
                t_bc = gSum(omfield[j].mesh().magSf().boundaryField()[p] *
                            omfield[j].boundaryField()[p]) / area;
                volScalarField C(Lfield[0].name(), omfield[j] - liftfield[k]*t_bc / t_lf);
                omfield.set(j, C.clone());
            }
        }
    }
}

}
