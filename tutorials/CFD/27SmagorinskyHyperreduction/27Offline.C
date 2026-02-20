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
Description
    Example of the hyperreduction of the Smagorinsky term in a ROM
SourceFiles
    27SmagorinskyHyperreduction.C
\*---------------------------------------------------------------------------*/

#include "27Offline.H"


tutorial27_offline::tutorial27_offline(int argc, char* argv[]):
    m_parameters(new StoredParameters(argc, argv)),
    m_UnsteadyNSTurb(new UnsteadyNSTurb(m_parameters)),
    nModesU(m_parameters->get_nModes()["U"]),
    interpolatedField(m_parameters->get_DEIMInterpolatedField()),
    nModesHR(m_parameters->get_nModes()[interpolatedField]),
    nMagicPoints(m_parameters->get_nMagicPoints()),
    magicPoints(m_parameters->get_magicPoints())
{
    if (!ITHACAutilities::containsSubstring(interpolatedField, "reduced"))
    {
        Info << "Error: Hyperreduction InterpolatedField not valid: " << interpolatedField << endl;
        Info << "Expects reducedFullStressFunction or reducedNut" << endl;
        abort();
    }
}


void tutorial27_offline::decompose()
{
    ITHACAPOD::PODTemplate<volVectorField> podU(m_parameters, "U");
    
    // compute spatial modes, temporal modes and cov matrix for U
    podU.getModes(spatialModesU, temporalModesU, temporalModesUSimulation, covMatrixU);
    
    // get mean of U field computed when getModes was used
    meanU = new volVectorField(podU.get_mean());
    m_parameters->set_meanU(*meanU);

    volVectorField template_HRInterpField(interpolatedField, m_UnsteadyNSTurb->initSmagFunction());
    m_parameters->set_template_field_fullStressFunction(template_HRInterpField);

    if (ITHACAutilities::containsSubstring(interpolatedField, "fullStressFunction"))
    {
        bool mgPoints_0_Neighborhoods_1 = true;

        if (m_parameters->get_HRMethod() == "DEIM")
        {
            TurbDiffusionHyperreduction<volVectorField,volVectorField,volVectorField>*
                turbDiffHR_fullStressFunction = new TurbDiffusionHyperreduction<volVectorField,volVectorField,volVectorField>
                (m_parameters, mgPoints_0_Neighborhoods_1, template_HRInterpField, template_HRInterpField, meanU, spatialModesU);

            turbDiffHR_fullStressFunction->computeTurbDiffusionHyperreduction();
        }
        else if (m_parameters->get_HRMethod() == "ECP")
        {
            volScalarField template_HRSnapshotsField(m_parameters->get_HRSnapshotsField(),
                                                        m_UnsteadyNSTurb->initProjSmagFunction());
            TurbDiffusionHyperreduction<volVectorField,volVectorField,volScalarField>*
                turbDiffHR_projFullStressFunction = new TurbDiffusionHyperreduction<volVectorField,volVectorField,volScalarField>
                (m_parameters, mgPoints_0_Neighborhoods_1, template_HRInterpField, template_HRSnapshotsField, meanU, spatialModesU);

            turbDiffHR_projFullStressFunction->computeTurbDiffusionHyperreduction();
        }
    }
    else if (ITHACAutilities::containsSubstring(interpolatedField, "nut"))
    {
        bool mgPoints_0_Neighborhoods_1 = false;
        volScalarField template_HRInterpField(interpolatedField, m_UnsteadyNSTurb->initNutFunction());
        volScalarField template_HRSnapshotsField(m_UnsteadyNSTurb->initProjSmagFromNutFunction());      
        TurbDiffusionHyperreduction<volVectorField,volScalarField,volScalarField>*
            turbDiffHR_nut = new TurbDiffusionHyperreduction<volVectorField,volScalarField,volScalarField>
            (m_parameters, mgPoints_0_Neighborhoods_1, template_HRInterpField, template_HRSnapshotsField, meanU, spatialModesU);
        
        turbDiffHR_nut->computeTurbDiffusionHyperreduction();
    }
    else
    {
        Info << "Error: Hyperreduction InterpolatedField not valid: " << interpolatedField << endl;
        Info << "Expects reducedFullStressFunction or reducedNut" << endl;
        abort();
    }
}


void tutorial27_offline::project()
{
    word folder_DEIM = m_parameters->get_folder_DEIM();

    int d = 3;
    Eigen::MatrixXd K_proj(nModesU, d*(nMagicPoints+1));
    Eigen::MatrixXd K_DEIM = m_parameters->get_K_DEIM();

    if (ITHACAutilities::containsSubstring(interpolatedField, "fullStressFunction"))
    {
        
        word file_projK_DEIM = folder_DEIM + "projected_K_DEIM_" + std::to_string(nModesU) + ".npy";
        bool exist_projK_DEIM = ITHACAutilities::check_file(file_projK_DEIM);

        if (!exist_projK_DEIM)
        {
            
            const volScalarField weight = sqrt(m_parameters->get_volume());
            volVectorField meanSmag(m_parameters->get_meanVectorDEIM());

            for (int i = 0 ; i < nModesU ; i++)
            {
                Eigen::VectorXd K_proj_line_i(d*nMagicPoints);
                if (m_parameters->get_HRMethod() == "DEIM")
                {
                    volVectorField weightedSpatialModesU_OF = weight * spatialModesU[i];
                    Eigen::VectorXd weightedSpatialModesU = Foam2Eigen::field2Eigen(weightedSpatialModesU_OF);
                    K_proj_line_i = K_DEIM.transpose() * weightedSpatialModesU;
                }
                else if (m_parameters->get_HRMethod() == "ECP")
                {
                    K_proj_line_i = K_DEIM.row(i);
                }

                K_proj.row(i).head(d*nMagicPoints) = K_proj_line_i;
                K_proj(i, d*nMagicPoints) = ITHACAutilities::dot_product_L2(spatialModesU[i], meanSmag);
                K_proj(i, d*nMagicPoints + 1) = 0.0;
                K_proj(i, d*nMagicPoints + 2) = 0.0;
            }

            cnpy::save(K_proj, file_projK_DEIM);
            m_parameters->set_projected_K_DEIM(K_proj);
        }
        else
        {
            cnpy::load(K_proj, file_projK_DEIM);
            m_parameters->set_projected_K_DEIM(K_proj);
        }
    }
    else if (ITHACAutilities::containsSubstring(interpolatedField, "nut"))
    {
        Eigen::Tensor<double, 3> xi_onMagicPts(nModesU + 1, nMagicPoints + 1, nModesU);
        word xi_onMagicPts_file("xi_onMagicPts.npy");
        bool exist_xi_onMagicPts = ITHACAutilities::check_file(folder_DEIM + xi_onMagicPts_file);

        if (!exist_xi_onMagicPts)
        {
            volScalarField meanNut = m_parameters->get_meanScalarDEIM();

            if (m_parameters->get_HRMethod() == "DEIM")
            {
                K_DEIM.array().colwise() /= Foam2Eigen::field2Eigen(m_parameters->get_volume()).array().sqrt();
                for (label i = 0; i < nModesU; i++)
                {
                    for (label q = 0; q < nModesU+1; q++)
                    {
                        for (label j = 0; j < nMagicPoints; j++)
                        {
                            Eigen::VectorXd K_DEIM_col_j = K_DEIM.col(j);
                            if (q == 0)
                            {
                                xi_onMagicPts(q,j,i) = ITHACAutilities::dot_product_L2(spatialModesU[i],
                                    m_UnsteadyNSTurb->diffusion(Foam2Eigen::Eigen2field(meanNut,K_DEIM_col_j),*meanU));
                            }
                            else
                            {
                                xi_onMagicPts(q,j,i) = ITHACAutilities::dot_product_L2(spatialModesU[i],
                                    m_UnsteadyNSTurb->diffusion(Foam2Eigen::Eigen2field(meanNut,K_DEIM_col_j),spatialModesU[q-1]));
                            }
                        }
                        if (q == 0)
                        {
                            xi_onMagicPts(q,nMagicPoints,i) = ITHACAutilities::dot_product_L2(spatialModesU[i],
                                m_UnsteadyNSTurb->diffusion(meanNut,*meanU));
                        }
                        else
                        {
                            xi_onMagicPts(q,nMagicPoints,i) = ITHACAutilities::dot_product_L2(spatialModesU[i],
                                m_UnsteadyNSTurb->diffusion(meanNut,spatialModesU[q-1]));
                        }
                    }
                }
            }
            else if (m_parameters->get_HRMethod() == "ECP")
            {
                PtrList<volTensorField> defTensorModesOnMgPts = m_parameters->get_deformationTensorOfModesOnMagicPoints();
                labelList localMagicPoints(m_parameters->get_localMagicPoints());
                
                for (label i = 0; i < nModesU; i++)
                {
                    for (label q = 0; q < nModesU+1; q++)
                    {
                        for (label j = 0; j < nMagicPoints; j++)
                        {
                            label row_id = ((pow(max(i,q-1),2)+3*max(i,q-1))/2 + min(i+1,q)) * (m_parameters->get_ECPAlgo()!="Global");
                            xi_onMagicPts(q,j,i) = -2*K_DEIM(row_id,j)*(dev(defTensorModesOnMgPts[i+1][localMagicPoints[j]]) 
                                                    && dev(defTensorModesOnMgPts[q][localMagicPoints[j]]));
                        }

                        if (q == 0)
                        {
                            xi_onMagicPts(q,nMagicPoints,i) = -fvc::domainIntegrate(m_UnsteadyNSTurb->
                                    projDiffusionIBP(meanNut, *meanU, spatialModesU[i])).value();
                        }
                        else
                        {
                            xi_onMagicPts(q,nMagicPoints,i) = -fvc::domainIntegrate(m_UnsteadyNSTurb->
                                    projDiffusionIBP(meanNut, spatialModesU[q-1], spatialModesU[i])).value();
                        }
                    }
                }
            }
            cnpy::save(xi_onMagicPts, folder_DEIM + xi_onMagicPts_file);
            m_parameters->set_xi_onMagicPts(xi_onMagicPts);
        }
        else
        {
            cnpy::load(xi_onMagicPts, folder_DEIM + xi_onMagicPts_file);
            m_parameters->set_xi_onMagicPts(xi_onMagicPts);
        }
    }
}

