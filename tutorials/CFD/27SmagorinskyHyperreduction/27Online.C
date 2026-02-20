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

#include "27Online.H"


tutorial27_online::tutorial27_online(StoredParameters* parameters):
    m_parameters(parameters),
    m_UnsteadyNSTurb(new UnsteadyNSTurb(m_parameters)),
    runTime2(Foam::Time::controlDictName, ".", m_parameters->get_casenameData()),
    meanU(new volVectorField(m_parameters->get_meanU())),
    nModesU(m_parameters->get_nModes()["U"]),
    interpolatedField(m_parameters->get_DEIMInterpolatedField()),
    nModesHR(m_parameters->get_nModes()[interpolatedField]),
    nMagicPoints(m_parameters->get_nMagicPoints()),
    magicPoints(m_parameters->get_magicPoints()),
    localMagicPoints(m_parameters->get_localMagicPoints()),
    Ck(m_parameters->get_Ck()),
    Ce(m_parameters->get_Ce())
{
    ITHACAstream::read_fields(spatialModesU, "U", "./ITHACAoutput/spatialModes_" + std::to_string(nModesU) + "modes/");
    massMatrixInv = ITHACAutilities::getMassMatrix(spatialModesU).completeOrthogonalDecomposition().pseudoInverse();
    cnpy::load(temporalModesUSimulation, "./ITHACAoutput/temporalModesSimulation_" + std::to_string(nModesU) + "modes/U.npy");

    weightAtMg = Eigen::VectorXd(nMagicPoints);
    for (int i = 0; i < nMagicPoints; i++)
    {
        weightAtMg(i) = std::pow(m_parameters->get_submesh().V()[localMagicPoints[i]], 0.5);
    }

    if (ITHACAutilities::containsSubstring(interpolatedField, "fullStressFunction"))
    {
        d = 3;
        meanSmagOnMagicNeighborhoods = new volVectorField(m_parameters->get_meanVectorDEIMMagic());
        for (auto& defTensorOfMode:m_parameters->get_deformationTensorOfModesOnMagicNeighborhoods())
        {
            defTensorOfModesAtMg.push_back(defTensorOfMode);
        }
    }
    else if (ITHACAutilities::containsSubstring(interpolatedField, "nut"))
    {
        d = 1;
        meanNutOnMagicPoints = new volScalarField(m_parameters->get_meanScalarDEIMMagic());
        for (auto& defTensorOfMode:m_parameters->get_deformationTensorOfModesOnMagicPoints())
        {
            defTensorOfModesAtMg.push_back(defTensorOfMode);
        }
    }

    nut0 = new volScalarField("nut", volScalarField(
        IOobject("nut0", m_parameters->get_casenameData() + "0", m_parameters->get_submesh(), IOobject::MUST_READ),
        m_parameters->get_submesh()));

    delta = new volScalarField(m_parameters->get_magicDelta());
    aaa = new volScalarField(Ce/(*delta));
    inv2aaa = new volScalarField(1/(2*(*aaa)));

    fullDefField = new volTensorField(defTensorOfModesAtMg[0]);
    bbb = new volScalarField((2.0/3.0)*tr(*fullDefField));
    ccc = new volScalarField(2*Ck*(*delta)*(dev(*fullDefField) && (*fullDefField)));
    sqrtk = new volScalarField((-(*bbb) + sqrt(sqr(*bbb) + 4*(*aaa)*(*ccc)))*(*inv2aaa));
    nut = new volScalarField(Ck*(*delta)*(*sqrtk));
}


void tutorial27_online::evaluateApproxNut(const Eigen::VectorXd& reducedCoeffs)
{
    *fullDefField = defTensorOfModesAtMg[0];
    for (int c = 0 ; c < reducedCoeffs.size() ; c++)
    {
        *fullDefField += reducedCoeffs(c)*defTensorOfModesAtMg[c+1];
    }
    *bbb = (2.0/3.0) * tr(*fullDefField);
    *ccc = 2*Ck * (*delta) * ((dev(*fullDefField) && (*fullDefField)));
    *sqrtk = mag(-(*bbb) + sqrt(sqr(*bbb) + 4*(*aaa) * (*ccc))) * (*inv2aaa);
    *nut = Ck * (*delta) * (*sqrtk);
    for (int i = 0; i < nut0->size(); i++)
    {
        (*nut0)[i] = (*nut)[i];
    }
    nut0->correctBoundaryConditions();
}

Eigen::VectorXd tutorial27_online::computeApproxSmagMg(const Eigen::VectorXd& reducedCoeffs)
{
    evaluateApproxNut(reducedCoeffs);
    volVectorField stressField = fvc::div(2*(*nut0) * dev(*fullDefField)) - *meanSmagOnMagicNeighborhoods;

    Eigen::VectorXd output = Eigen::VectorXd::Zero(d * (nMagicPoints + 1));
    for (int i = 0; i < nMagicPoints; i++)
    {
        vector stressVector = weightAtMg(i) * (stressField)[localMagicPoints[i]];

        output(i * d) = stressVector.x();
        output(i * d + 1) = stressVector.y();
        output(i * d + 2) = stressVector.z();
    }
    output(d * nMagicPoints) = 1.0;
    return output;
}

Eigen::VectorXd tutorial27_online::computeApproxNutMg(const Eigen::VectorXd& reducedCoeffs)
{
    evaluateApproxNut(reducedCoeffs);
    *nut0 -= *meanNutOnMagicPoints;

    Eigen::VectorXd output(nMagicPoints + 1);
    for (int i = 0; i < nMagicPoints; i++)
    {
        output(i) = weightAtMg(i) * (*nut0)[localMagicPoints[i]];
    }
    output(nMagicPoints) = 1.0;
    return output;
}

Eigen::VectorXd tutorial27_online::predictSmagROMCoeffs(const Eigen::VectorXd& reducedCoeffs)
{
    if (ITHACAutilities::containsSubstring(interpolatedField, "fullStressFunction"))
    {
        Eigen::MatrixXd projected_K_DEIM = m_parameters->get_projected_K_DEIM();
        Eigen::VectorXd smagVect(computeApproxSmagMg(reducedCoeffs));

        return massMatrixInv*projected_K_DEIM*smagVect;
    }

    else if (ITHACAutilities::containsSubstring(interpolatedField, "nut"))
    {
        Eigen::VectorXd nutVect(computeApproxNutMg(reducedCoeffs));
        Eigen::Tensor<double, 3 > xi_complete_onMagicPts(m_parameters->get_xi_onMagicPts());
        Eigen::VectorXd b_extended(reducedCoeffs.size() + 1);
        b_extended << 1.0,reducedCoeffs;

        Eigen::VectorXd output(nModesU);
        for (int i = 0; i < nModesU; i++)
        {
            output(i) = 0;
            for (int j = 0 ; j < nMagicPoints +1; j++)
            {
                for (int q = 0 ; q < nModesU + 1 ; q++)
                {
                    output(i) += xi_complete_onMagicPts(q,j,i) * b_extended(q) * nutVect(j);
                }
            }
        }
        return massMatrixInv*output;
    }
}

// Compute a vector field onto velocity modes from its reduced temporal modes. Homogenize to a stress if specified.
volVectorField tutorial27_online::computeROMproj_fromCoeffs(const Eigen::VectorXd& reducedCoeffs, bool stressUnit)
{
    Eigen::MatrixXd reducedCoeffsMatrix(nModesU, 1);
    for (label k = 0; k < nModesU; k++)
    {
        reducedCoeffsMatrix(k,0) = reducedCoeffs[k];
    }
    PtrList<volVectorField> reconstruction = ITHACAutilities::reconstructFromCoeff(spatialModesU, reducedCoeffsMatrix, nModesU);

    if (stressUnit)
    {
        // Dummy scalar to homogenized the reconstruction
        dimensionedScalar time_inv("time_inv", dimensionSet(0,0,-1,0,0), scalar(1.0));
        return reconstruction[0]*time_inv;
    }
    else
    {
        return reconstruction[0];
    }
}

// Compute the coefficients of a vector field in the velocity basis (V^T f)
Eigen::VectorXd tutorial27_online::computeROMcoeffs_fromFullDim(volVectorField& f_full)
{
    bool consider_volumes = true;
    Eigen::MatrixXd coeffsVelocityPOD = ITHACAutilities::getCoeffs(f_full, spatialModesU, nModesU, consider_volumes);

    return coeffsVelocityPOD.col(0);
}

// Compute f
volVectorField tutorial27_online::computeSmagTerm_fromChronos(const Eigen::VectorXd& reducedCoeffs)
{
    return m_UnsteadyNSTurb->computeSmagTerm_fromU(*meanU + computeROMproj_fromCoeffs(reducedCoeffs));
}


void tutorial27_online::prediction()
{
    word folder_results = "./ITHACAoutput/Online/" + m_parameters->get_HRMethod() + "_" + interpolatedField.substr(7)
                          + "U_" + std::to_string(nMagicPoints) + "MgPts/";
    mkDir(folder_results);
    Eigen::MatrixXd predictedCoeffs(m_parameters->get_nSnapshotsSimulation(), nModesU);
    Eigen::MatrixXd exactCoeffs(m_parameters->get_nSnapshotsSimulation(), nModesU);
    Eigen::VectorXd relProjError(m_parameters->get_nSnapshotsSimulation());

    for (label t = 0; t < m_parameters->get_nSnapshotsSimulation(); t++)
    {
        label index_time = t + m_parameters->get_startTime() + m_parameters->get_nSnapshots() - 1;

        // Computing prediction
        predictedCoeffs.row(t) = predictSmagROMCoeffs(temporalModesUSimulation.row(t));
        volVectorField predictedSmagProj = computeROMproj_fromCoeffs(predictedCoeffs.row(t), true);

        // Computing references from reduced velocity
        volVectorField exactSmag = computeSmagTerm_fromChronos(temporalModesUSimulation.row(t));
        exactCoeffs.row(t) = computeROMcoeffs_fromFullDim(exactSmag);
        volVectorField exactSmagProj = computeROMproj_fromCoeffs(exactCoeffs.row(t), true);

        // Evaluating errors with reference
        relProjError(t) = ITHACAutilities::errorL2Rel(exactSmagProj, predictedSmagProj);

        std::string idxTimeStr(runTime2.times()[index_time].name());
        ITHACAstream::exportSolution(exactSmagProj, idxTimeStr, folder_results + "/exactSmagProj","fullStressFunction");
        ITHACAstream::exportSolution(predictedSmagProj, idxTimeStr, folder_results + "/predictedSmagProj","fullStressFunction");
    }
    cnpy::save(predictedCoeffs, folder_results + "predictedCoeffs.npy");
    cnpy::save(exactCoeffs, folder_results + "exactCoeffs.npy");
    cnpy::save(relProjError, folder_results + "relativeSmagProjError.npy");
}

