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

#include "TurbDiffusionHyperreduction.H"

/// \file
/// Source file of the TurbDiffusionHyperreduction class.



template<typename volF, typename T, typename S>
TurbDiffusionHyperreduction<volF,T,S>::TurbDiffusionHyperreduction
        (Parameters* myParameters, bool mgPoints_0_Neighborhoods_1, T& template_HRInterpField, 
         S& template_HRSnapshotsField, volF*& meanU, PtrList<volF>& spatialModesU):
    m_parameters(static_cast<StoredParameters*>(myParameters)),
    ithacaUnsteadyNSTurb(autoPtr<UnsteadyNSTurb>(new UnsteadyNSTurb(m_parameters))),
    l_nSnapshot(m_parameters->get_nSnapshots()),  
    l_startTime(m_parameters->get_startTime()),
    l_nSnapshotSimulation(m_parameters->get_nSnapshotsSimulation()),
    template_HRInterpField(template_HRInterpField),
    template_HRSnapshotsField(template_HRSnapshotsField),
    f_spatialModesU(spatialModesU),
    f_meanU(meanU),
    folder_HR(m_parameters->get_folder_DEIM()),
    runTime2(Foam::Time::controlDictName, ".", m_parameters->get_casenameData()),
    snapshotsFolder(m_parameters->get_casenameData()),
    mgPoints_0_Neighborhoods_1(mgPoints_0_Neighborhoods_1),
    HRMethod(m_parameters->get_HRMethod()),
    HRInterpField(m_parameters->get_DEIMInterpolatedField()),
    HRSnapshotsField(m_parameters->get_HRSnapshotsField()),
    ECPAlgo(m_parameters->get_ECPAlgo()),
    interpFieldCenteredOrNot(m_parameters->get_interpFieldCenteredOrNot())
{
    //
}



template<typename volF, typename T, typename S>
TurbDiffusionHyperreduction<volF,T,S>::~TurbDiffusionHyperreduction()
{
    //
}


template<typename volF, typename T, typename S>
template<typename volGradF>
void TurbDiffusionHyperreduction<volF,T,S>::computeDefTensor(volGradF tensorMeanU)
{
    int nModesU = f_spatialModesU.size();
    PtrList<volGradF> defTensorOfModesOnMgPointsORMgNeighborhoods;
    defTensorOfModesOnMgPointsORMgNeighborhoods.setSize(nModesU + 1);

    volGradF temp = ithacaHyperreduction->interpolateField(tensorMeanU);
    defTensorOfModesOnMgPointsORMgNeighborhoods.set(0, new volGradF(temp));

    for (int c=0 ; c<nModesU ; c++)
    {
        volGradF defTensorMode(UnsteadyNSTurb::computeS_fromU(f_spatialModesU[c]));
        volGradF temp = ithacaHyperreduction->interpolateField(defTensorMode);
        defTensorOfModesOnMgPointsORMgNeighborhoods.set(c+1, new volGradF(temp));
    }
    if(!mgPoints_0_Neighborhoods_1)
    {
        m_parameters->set_deformationTensorOfModesOnMagicPoints(defTensorOfModesOnMgPointsORMgNeighborhoods);
    }
    else
    {
        m_parameters->set_deformationTensorOfModesOnMagicNeighborhoods(defTensorOfModesOnMgPointsORMgNeighborhoods);
    }
}

template<typename volF, typename T, typename S>
void TurbDiffusionHyperreduction<volF,T,S>::common_MeanHRInterpField()
{
      T meanHR(template_HRInterpField.name(), template_HRInterpField);
      if (interpFieldCenteredOrNot)
      {
        PtrList<T> meanRead;  
        ITHACAstream::read_fields(meanRead, HRInterpField, "./ITHACAoutput/mean/");
        meanHR = meanRead[0];
      }
      else
      {
        ITHACAutilities::setToZero(meanHR);
      }
      T meanHRMagic = ithacaHyperreduction->interpolateField(meanHR);
      m_parameters->set_meanDEIMMagic(meanHRMagic);  
}

template<typename volF, typename T, typename S>
void TurbDiffusionHyperreduction<volF,T,S>::precomputeTurbDiffusionFunctions(word& fieldToCompute)
{
    label nTotalSnapshots = l_nSnapshot;
    if (fieldToCompute.substr(0,HRInterpField.size()) == HRInterpField)
    {
        nTotalSnapshots += (l_nSnapshotSimulation-1);
    }

    double startingTime = std::stod(runTime2.times()[l_startTime].name());
    double saveTime = m_parameters->get_saveTime();

    if (HRMethod == "ECP" && ITHACAutilities::containsSubstring(fieldToCompute, HRSnapshotsField))
    {
        if (ECPAlgo == "Global")
        {
            string pathCentered = "";
            if (interpFieldCenteredOrNot) {pathCentered = "_centered";}
            snapshotsFolder = "./ITHACAoutput/Hyperreduction/" + HRSnapshotsField.substr(0, HRSnapshotsField.find("_")) 
                                + "/Global" + pathCentered + "/";
            nTotalSnapshots = f_spatialModesU.size()*(l_nSnapshot);
            if (HRInterpField == "nut" || ITHACAutilities::containsSubstring(HRInterpField, "reducedNut"))
            {
                nTotalSnapshots = (f_spatialModesU.size()*(f_spatialModesU.size() + 1)/2 + f_spatialModesU.size()) * l_nSnapshot;
            }
            startingTime = 1;
            saveTime = 1;
            POD_nSnapshot = nTotalSnapshots;
            POD_startTime = 2;   // Subfolders indices are 0:constant, 1:0, 2:1
            POD_endTime = POD_startTime + POD_nSnapshot - 1;
        }
        else
        {
            snapshotsFolder = "./ITHACAoutput/Hyperreduction/" + HRSnapshotsField.substr(0, HRSnapshotsField.find("_")) + "/EachMode/";
        }
        l_nSnapshotSimulation = 0;
    }

    word pathProcessor("");
    if (Pstream::parRun())
    {
        pathProcessor = "processor" + name(Pstream::myProcNo()) + "/";
    }

    string local_file = snapshotsFolder + pathProcessor + runTime2.times()[1].name() +  "/" + fieldToCompute;
    bool exist_precomputed_fields = ITHACAutilities::check_file(local_file);

    for (label j = 0; j < nTotalSnapshots ; j++)
    {
        // Read the j-th field
        local_file = snapshotsFolder + pathProcessor + ITHACAutilities::double2ConciseString(startingTime + j*saveTime)
                     + "/" + fieldToCompute;
        exist_precomputed_fields = exist_precomputed_fields && ITHACAutilities::check_file(local_file);
    }

    string pathCentered = "";
    if (interpFieldCenteredOrNot) {pathCentered = "_centered";}
    bool exist_covMatrix = ITHACAutilities::check_file("./ITHACAoutput/CovMatrices" 
                           + pathCentered + "/covMatrix" + fieldToCompute + ".npy")
                           && ITHACAutilities::check_file("./ITHACAoutput/mean/1/" + fieldToCompute);
    exist_precomputed_fields = exist_precomputed_fields || exist_covMatrix;

    if (!exist_precomputed_fields)
    {
        label offset_1 = 0;
        if (fieldToCompute == HRInterpField || fieldToCompute == "fullStressFunction" || fieldToCompute == "nut")
        {
            Info << "Evaluating " << fieldToCompute << " term" << endl;
            T nonLinearSnapshotsj = template_HRInterpField;

            if (fieldToCompute == "fullStressFunction" 
                || ITHACAutilities::containsSubstring(fieldToCompute, "reducedFullStressFunction")
                || ITHACAutilities::containsSubstring(fieldToCompute, "reducedNut"))
            {
              ITHACAstream::exportSolution(nonLinearSnapshotsj, "0", snapshotsFolder, fieldToCompute);
              ITHACAstream::exportSolution(nonLinearSnapshotsj, "0", m_parameters->get_casenameData(), fieldToCompute);
            }
            else if(fieldToCompute == "nut")
            {
                // /!\ Nut requires Snapshot at time 0 to be initialized, not fullStressFunction field
                if (l_startTime == 1)
                {
                    offset_1 = 1;
                }
            }
            else
            {
                Info << "Error : HRInterpolatedField not valid : " << HRInterpField << endl;
                Info << "HR is available for fullStressFunction and nut only."  << endl;
                abort();
            }

            for (label j = offset_1; j < nTotalSnapshots; j++)
            {
                word snap_time = runTime2.times()[l_startTime + j].name();
                string idxTimeStr = ITHACAutilities::double2ConciseString(startingTime + j*saveTime);
                if (!ITHACAutilities::containsSubstring(fieldToCompute, "reduced"))
                {
                  ithacaUnsteadyNSTurb->computeNonLinearSnapshot_at_time(snap_time, nonLinearSnapshotsj) ;
                }
                else
                {
                  ithacaUnsteadyNSTurb->computeNonLinearSnapshot_at_time(snap_time, nonLinearSnapshotsj, f_spatialModesU) ;
                }
                ITHACAstream::exportSolution(nonLinearSnapshotsj, idxTimeStr, snapshotsFolder, fieldToCompute);
            }

        } 

        // Compute snapshots field if ECP, ie projected field
        else if (ITHACAutilities::containsSubstring(fieldToCompute, HRSnapshotsField))
        {
            if(ITHACAutilities::containsSubstring(HRSnapshotsField, "projFullStressFunction")
               || ITHACAutilities::containsSubstring(HRSnapshotsField, "projReducedFullStressFunction"))
            {
                int nUsedModesU = 0;
                int currentModeU = 0;
                T meanHRInterpField = template_HRInterpField;
                S snapshotECPcj(fieldToCompute, template_HRSnapshotsField);
                S meanECPModec = template_HRSnapshotsField;
                ITHACAutilities::setToZero(meanECPModec);

                if (ECPAlgo == "Global")
                {
                    nUsedModesU = f_spatialModesU.size();
                    PtrList<T> meanRead;  
                    ITHACAstream::read_fields(meanRead, HRInterpField, "./ITHACAoutput/mean/");
                    meanHRInterpField = meanRead[0];
                }
                else
                {
                    nUsedModesU = 1;
                    std::vector<int> intInFieldName = ITHACAutilities::extractIntFromString(fieldToCompute);
                    currentModeU = intInFieldName[intInFieldName.size()-1]-1;
                }
            
                Info << "Evaluating " << fieldToCompute << " term" << endl;
                ITHACAstream::exportSolution(snapshotECPcj, "0", snapshotsFolder, fieldToCompute);
                ITHACAstream::exportSolution(snapshotECPcj, "0", m_parameters->get_casenameData(), fieldToCompute);

                for (label c = 0; c < nUsedModesU; c++)
                {
                    if (ECPAlgo == "Global")
                    {
                        currentModeU = c;
                        if (interpFieldCenteredOrNot)
                        {
                            // Computing the mean of projFullStressFunction for mode c
                            ithacaUnsteadyNSTurb->computeProjSmagTerm_fromSmag_fromMode
                                (meanECPModec, meanHRInterpField, f_spatialModesU[currentModeU]);
                        }
                    }

                    for (label j = 0; j < nTotalSnapshots/nUsedModesU; j++)       
                    {
                        word snap_time = runTime2.times()[l_startTime + j].name();
                        string idxStr = ITHACAutilities::double2ConciseString(startingTime + (c*l_nSnapshot + j)*saveTime);
                        ithacaUnsteadyNSTurb->computeNonLinearSnapshot_at_time(snap_time, snapshotECPcj, f_spatialModesU[currentModeU]);
                        ITHACAutilities::subtractFields(snapshotECPcj, meanECPModec);
                        ITHACAstream::exportSolution(snapshotECPcj, idxStr, snapshotsFolder, fieldToCompute);
                        ITHACAstream::printProgress(double(c*l_nSnapshot + j) / (nTotalSnapshots-1));
                    }
                }
                Info << endl;
            }

            else if(ITHACAutilities::containsSubstring(HRSnapshotsField, "projSmagFromNut")
                    || ITHACAutilities::containsSubstring(HRSnapshotsField, "projSmagFromReducedNut"))
            {
                std::vector<int> currentModesU;
                T meanHRInterpField = template_HRInterpField;
                S snapshotECPckj(fieldToCompute, template_HRSnapshotsField);
                S meanECPModeck = template_HRSnapshotsField;
                ITHACAutilities::setToZero(meanECPModeck);

                PtrList<volVectorField> f_mean_and_spatialModesU(f_spatialModesU.size() + 1);
                for (label c = 0; c < f_spatialModesU.size(); c++)      
                {
                    f_mean_and_spatialModesU.set(c+1, new volF(f_spatialModesU[c]));
                }
                f_mean_and_spatialModesU.set(0, new volF(*f_meanU));

                if (ECPAlgo == "Global")
                {
                    PtrList<T> meanRead;  
                    ITHACAstream::read_fields(meanRead, HRInterpField, "./ITHACAoutput/mean/");
                    meanHRInterpField = meanRead[0];
                }
                else
                {
                    currentModesU = ITHACAutilities::extractIntFromString(fieldToCompute);
                }
          
                Info << "Evaluating " << fieldToCompute << " term" << endl;
                ITHACAstream::exportSolution(snapshotECPckj, "0", snapshotsFolder, fieldToCompute);
                ITHACAstream::exportSolution(snapshotECPckj, "0", m_parameters->get_casenameData(), fieldToCompute);

                for (label c = 0; c < f_spatialModesU.size(); c++)
                {
                    for (label k = 0; k <= c+1; k++)
                    {
                        if (ECPAlgo == "EachMode")
                        {
                            c = currentModesU[currentModesU.size()-2] - 1;
                            k = currentModesU[currentModesU.size()-1];
                        }
                        else if (ECPAlgo == "Global" && interpFieldCenteredOrNot)
                        {
                            // Computing the mean of projSmagFromNut for modes c,k
                            ithacaUnsteadyNSTurb->computeProjSmagFromNut_fromNut_fromModes(meanECPModeck, meanHRInterpField, 
                                                                            f_spatialModesU[c], f_mean_and_spatialModesU[k]);
                        }

                        for (label j = 0; j < l_nSnapshot; j++)       
                        {
                            word snap_time = runTime2.times()[l_startTime + j].name();
                            string idxStr = ITHACAutilities::double2ConciseString(startingTime + (((c*(c+1)/2 + c)*l_nSnapshot
                                                        + k*l_nSnapshot)*(ECPAlgo=="Global") + j)*saveTime);
                            ithacaUnsteadyNSTurb->computeNonLinearSnapshot_at_time(snap_time, snapshotECPckj,
                                                             f_spatialModesU[c], f_mean_and_spatialModesU[k]);
                            ITHACAutilities::subtractFields(snapshotECPckj, meanECPModeck);
                            ITHACAstream::exportSolution(snapshotECPckj, idxStr, snapshotsFolder, fieldToCompute);
                            ITHACAstream::printProgress(double(((c*(c+1)/2 + c)*l_nSnapshot + k*l_nSnapshot)
                                                            *(ECPAlgo=="Global") + j) / (nTotalSnapshots-1));
                        }

                        if (ECPAlgo == "EachMode")
                        {
                            goto endloops;
                        }
                    }
                }
                endloops:
                Info << endl;
            }
        }
        else
        {
            Info << "Error : Unexpected field for precomputeNonPolynomialFunction : " << fieldToCompute << endl;
            abort();
        }
    }


    // Compute/Read and store in m_parameters the mean of the non polynomial field from high-fidelity velocity
    if (fieldToCompute == HRInterpField || fieldToCompute == "fullStressFunction" || fieldToCompute == "nut")  
    {
        if (interpFieldCenteredOrNot)
        {
            T meanField = template_HRInterpField;
            if (ITHACAutilities::check_file("./ITHACAoutput/mean/" + pathProcessor + "1/" + fieldToCompute))
            {
                Info << "Reading the mean of " << fieldToCompute << " field" << endl;
                PtrList<T> meanRead;  
                ITHACAstream::read_fields(meanRead, fieldToCompute, "./ITHACAoutput/mean/");
                meanField = meanRead[0];
            }
            else  // Computing the mean
            {
                ITHACAPOD::PODTemplate<T> ithacaFVPOD_interp_field(m_parameters, fieldToCompute, snapshotsFolder);
                ithacaFVPOD_interp_field.set_b_centeredOrNot(1);
                ithacaFVPOD_interp_field.computeMeanField();
                meanField = ithacaFVPOD_interp_field.get_mean();
            }

            if (fieldToCompute == "fullStressFunction" || fieldToCompute == "nut")
            {
                m_parameters->set_meanDEIM(meanField);
            }
        }  
        else
        {
            T meanZero(fieldToCompute, template_HRInterpField);
            ITHACAutilities::setToZero(meanZero);
            m_parameters->set_meanDEIM(meanZero);
        }
    }
}

template<typename volF, typename T, typename S>
void TurbDiffusionHyperreduction<volF,T,S>::computeTurbDiffusionHyperreduction()
{
    // Compute the mean of the nonpolynomial field from high-fidelity velocity if HR learnt from reduced velocity
    if (ITHACAutilities::containsSubstring(HRInterpField, "reduced"))
    {
        word nameNotReduced = HRInterpField.substr(7, HRInterpField.find("_")-7);
        nameNotReduced[0] = tolower(nameNotReduced[0]);
        precomputeTurbDiffusionFunctions(nameNotReduced);
        snapshotsFolder = "./ITHACAoutput/Hyperreduction/" + HRInterpField.substr(0, HRInterpField.find("_")) + "/";
    }

    // Evaluate and save the nonpolynomial field for HR learning
    precomputeTurbDiffusionFunctions(HRInterpField);

    std::vector<word> HRPODFields;
    // Evaluate and save the projected nonpolynomial field if ECP
    if (HRMethod == "DEIM")
    {
        HRPODFields.push_back(HRInterpField);
    }
    else if (HRMethod == "ECP")
    {
        word HRSnapshotsField_copy = HRSnapshotsField;
        std::vector<word> HRFields_temp{std::begin(m_parameters->get_field_name()), 
                                           std::end(m_parameters->get_field_name())};
        HRFields_temp.erase(std::remove_if(HRFields_temp.begin(), HRFields_temp.end(),
            [HRSnapshotsField_copy](const word x) {return !ITHACAutilities::containsSubstring(x,HRSnapshotsField_copy);}),
             HRFields_temp.end());
        HRPODFields = HRFields_temp;

        for (auto fieldIterator: HRPODFields)
        {
            precomputeTurbDiffusionFunctions(fieldIterator);
        }                  
    }

    for (auto fieldIterator: HRPODFields)
    {
        // POD of nonpolynomial fields
        PtrList<S> f_subsetSpatialModesHR;
        ITHACAPOD::PODTemplate<S> ithacaFVPOD_template_field(m_parameters, fieldIterator, snapshotsFolder);
        ithacaFVPOD_template_field.set_snapFolderParams(POD_nSnapshot,l_nSnapshotSimulation,POD_startTime,POD_endTime);
        ithacaFVPOD_template_field.set_b_centeredOrNot(interpFieldCenteredOrNot);
        ithacaFVPOD_template_field.getModes(f_subsetSpatialModesHR, m_temporalModesHR,
                                            m_temporalModesHRSimulation, covMatrixHR);
        forAll(f_subsetSpatialModesHR, i){f_spatialModesHR.append(tmp<S>(f_subsetSpatialModesHR[i]));}
    }

    if (Pstream::parRun())
    {
        bool waitMaster = false;
        if (Pstream::master()){ waitMaster = true; }
        Pstream::scatter(waitMaster);

        Info << "Hyperreduction POD done in parallel. Magic points selection must run sequentially. Aborting." << endl << endl;
        std::exit(111); 
    }

    // Convertion to Eigen and weight of the HR modes
    Eigen::VectorXd weights = ITHACAutilities::getMassMatrixFV(f_spatialModesHR[0]).array().sqrt(); 
    Eigen::MatrixXd spatialModesHR = weights.asDiagonal() * Foam2Eigen::PtrList2Eigen(f_spatialModesHR);
    Eigen::VectorXi initSeeds(0);

    // Constructor of the hyperReduction object
    ithacaHyperreduction = autoPtr<HyperReduction<PtrList<S>>>(
            new HyperReduction<PtrList<S>>(f_spatialModesHR.size(), 
                                           m_parameters->get_nMagicPoints(), 
                                           std::is_same<S, volScalarField>::value ? 1 : 3,
                                           f_spatialModesHR[0].size(),
                                           initSeeds, 
                                           folder_HR,
                                           weights.array().square().matrix()));

    Eigen::VectorXd fieldWeightsOnes = Eigen::VectorXd::Ones(ithacaHyperreduction->n_cells * ithacaHyperreduction->vectorial_dim);
    List<Eigen::VectorXd> quadWeights(HRPODFields.size()); 

   
    if (HRMethod == "DEIM")
    {
        ithacaHyperreduction->offlineGappyDEIM(spatialModesHR, fieldWeightsOnes, folder_HR);
    }

    else if (HRMethod == "ECP")
    {
        // ECP snapshots are weighted by cells' volume not sqrt
        spatialModesHR = weights.asDiagonal() * spatialModesHR;

        bool weightModesWithEigenval = 0;

        if (weightModesWithEigenval)
        {
            for (unsigned int c = 0; c < HRPODFields.size(); c++)
            {
                label nModesSubspace = m_parameters->get_nModes()[HRInterpField];
                Eigen::VectorXd eigenValuesModesHR;
                std::string pathCentered = "";
                if (interpFieldCenteredOrNot){pathCentered = "_centered";}
                std::string folder_eigen = "./ITHACAoutput/EigenValuesandVector" + pathCentered + "_" + 
                                            std::to_string(nModesSubspace) + "modes/";
                std::string eigen_name = "EigenvectorLambda_" + HRPODFields[c];
                cnpy::load(eigenValuesModesHR, folder_eigen + eigen_name + ".npy");
                spatialModesHR.middleCols(nModesSubspace*c, nModesSubspace).array().rowwise() *= 
                            (eigenValuesModesHR.transpose()/eigenValuesModesHR(0)).array();
            }
        }

        if (ECPAlgo == "Global")
        {
            quadWeights.resize(f_spatialModesU.size());
            ithacaHyperreduction->offlineECP(spatialModesHR, fieldWeightsOnes, folder_HR);
            for (label c = 0; c < f_spatialModesU.size(); c++)
            {
                quadWeights[c] = ithacaHyperreduction->quadratureWeights; 
            }
        }
        else if (ECPAlgo == "EachMode")
        {
            List<Eigen::MatrixXd> listspatialModesHR(HRPODFields.size());
            label nModesSubspace = m_parameters->get_nModes()[HRInterpField];
            for (unsigned int c = 0; c < HRPODFields.size(); c++)
            {
                listspatialModesHR[c] = spatialModesHR.middleCols(nModesSubspace*c, nModesSubspace);
            }
            ithacaHyperreduction->offlineECP(listspatialModesHR, fieldWeightsOnes, folder_HR);
            quadWeights = ithacaHyperreduction->quadratureWeightsSubspaces;
        }
    }


    if (HRMethod == "DEIM")
    {
        // Condition number 
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(ithacaHyperreduction->pinvPU); 
        double cond = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size() - 1);
        Info << "Condition number of P^T U = " << cond << endl;

        // Saves
        word folder_K_DEIM = folder_HR + "K_DEIM/";
        mkDir(folder_K_DEIM);
        cnpy::save(ithacaHyperreduction->pinvPU, folder_K_DEIM + "/P_U_inv.npy");
        cnpy::save(ithacaHyperreduction->MatrixOnline, folder_K_DEIM + "/MatrixOnline.npy");

        m_parameters->set_K_DEIM(ithacaHyperreduction->MatrixOnline);
    }

    else if (HRMethod == "ECP" && ITHACAutilities::containsSubstring(HRInterpField, "fullStressFunction"))
    {
        Eigen::VectorXd vectorWeights = ITHACAutilities::getMassMatrixFV(f_spatialModesU[0]).array().sqrt(); 
        Eigen::MatrixXd spatialModesU = vectorWeights.asDiagonal() * Foam2Eigen::PtrList2Eigen(f_spatialModesU);
        Eigen::SparseMatrix<double> P = ithacaHyperreduction->maskToOtherDim(3);
        Eigen::MatrixXd K_DEIM = spatialModesU.transpose() * P;
        for (label m = 0; m < ithacaHyperreduction->n_nodes; m++)
        {
            for (label c = 0; c < f_spatialModesU.size(); c++)
            {
                K_DEIM(c, 3*m) *= quadWeights[c][m];
                K_DEIM(c, 3*m + 1) *= quadWeights[c][m];
                K_DEIM(c, 3*m + 2) *= quadWeights[c][m];
            }
        }
        cnpy::save(K_DEIM, folder_HR + "/K_DEIM.npy");
        m_parameters->set_K_DEIM(K_DEIM);
    }

    else if (HRMethod == "ECP" && ITHACAutilities::containsSubstring(HRInterpField, "nut"))
    {
        Eigen::MatrixXd K_DEIM(HRPODFields.size(), ithacaHyperreduction->n_nodes);
        Eigen::VectorXd sparseWeights = weights.transpose() * ithacaHyperreduction->P;
        for (unsigned int c = 0; c < HRPODFields.size(); c++)
        {
            K_DEIM.row(c) = (sparseWeights.array() * quadWeights[c].array()).matrix();
        }
        cnpy::save(K_DEIM, folder_HR + "/K_DEIM.npy");
        m_parameters->set_K_DEIM(K_DEIM);
    }

        
    // Generate submeshes
    m_parameters->set_magicPoints(ithacaHyperreduction->nodePoints);
    ithacaHyperreduction->problemName = HRInterpField;
    ithacaHyperreduction->generateSubmesh(2, m_parameters->get_mesh());

    m_parameters->set_localMagicPoints(ithacaHyperreduction->localNodePoints);
    m_parameters->set_submesh(ithacaHyperreduction->submesh->subMesh());

    // Set delta on submesh
    m_parameters->set_magicDelta(ithacaHyperreduction->interpolateField(m_parameters->get_delta()));
         
    // Get deformation tensor on submesh
    computeDefTensor(UnsteadyNSTurb::computeS_fromU(*f_meanU));

    // Get temporal mean of interpolated nonpolynomial fields on mesh and submesh
    common_MeanHRInterpField();
}


template TurbDiffusionHyperreduction<volVectorField,volVectorField,volVectorField>::TurbDiffusionHyperreduction(
        Parameters* myParameters,
        bool mgPoints_0_Neighborhoods_1,
        volVectorField& template_HRInterpField,
        volVectorField& template_HRSnapshotsField,
        volVectorField*& meanU, PtrList<volVectorField>& spatialModesU);
template TurbDiffusionHyperreduction<volVectorField,volVectorField,volScalarField>::TurbDiffusionHyperreduction(
        Parameters* myParameters,
        bool mgPoints_0_Neighborhoods_1,
        volVectorField& template_HRInterpField,
        volScalarField& template_HRSnapshotsField,
        volVectorField*& meanU, PtrList<volVectorField>& spatialModesU);
template TurbDiffusionHyperreduction<volVectorField,volScalarField,volScalarField>::TurbDiffusionHyperreduction(
        Parameters* myParameters,
        bool mgPoints_0_Neighborhoods_1,
        volScalarField& template_HRInterpField,
        volScalarField& template_HRSnapshotsField,
        volVectorField*& meanU, PtrList<volVectorField>& spatialModesU);
        
template TurbDiffusionHyperreduction<volVectorField,volVectorField,volVectorField>::~TurbDiffusionHyperreduction();
template TurbDiffusionHyperreduction<volVectorField,volVectorField,volScalarField>::~TurbDiffusionHyperreduction();
template TurbDiffusionHyperreduction<volVectorField,volScalarField,volScalarField>::~TurbDiffusionHyperreduction();

template void TurbDiffusionHyperreduction<volVectorField,volVectorField,volVectorField>::computeDefTensor(volTensorField defTensorMean);
template void TurbDiffusionHyperreduction<volVectorField,volVectorField,volScalarField>::computeDefTensor(volTensorField defTensorMean);
template void TurbDiffusionHyperreduction<volVectorField,volScalarField,volScalarField>::computeDefTensor(volTensorField defTensorMean);

template void TurbDiffusionHyperreduction<volVectorField,volVectorField,volVectorField>::computeTurbDiffusionHyperreduction();
template void TurbDiffusionHyperreduction<volVectorField,volVectorField,volScalarField>::computeTurbDiffusionHyperreduction();
template void TurbDiffusionHyperreduction<volVectorField,volScalarField,volScalarField>::computeTurbDiffusionHyperreduction();

template void TurbDiffusionHyperreduction<volVectorField,volVectorField,volVectorField>::precomputeTurbDiffusionFunctions(word& fieldToCompute);
template void TurbDiffusionHyperreduction<volVectorField,volVectorField,volScalarField>::precomputeTurbDiffusionFunctions(word& fieldToCompute);
template void TurbDiffusionHyperreduction<volVectorField,volScalarField,volScalarField>::precomputeTurbDiffusionFunctions(word& fieldToCompute);

template void TurbDiffusionHyperreduction<volVectorField,volVectorField,volVectorField>::common_MeanHRInterpField();
template void TurbDiffusionHyperreduction<volVectorField,volVectorField,volScalarField>::common_MeanHRInterpField();
template void TurbDiffusionHyperreduction<volVectorField,volScalarField,volScalarField>::common_MeanHRInterpField();


