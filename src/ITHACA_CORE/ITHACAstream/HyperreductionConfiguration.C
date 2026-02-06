#include "HyperreductionConfiguration.H"
#include "GeometricField.H"
#include "fvPatchField.H"
#include "volFieldsFwd.H"

void HyperreductionConfiguration::setTracerGradOnMagicNeighborhoods(Foam::PtrList<Foam::volVectorField>& defT, Foam::label nModes)
{
    m_tracerGradOnMagicNeighborhoods.resize(nModes);
    for (Foam::label k = 0; k < nModes; k++)
    {
        m_tracerGradOnMagicNeighborhoods.set(k, new Foam::volVectorField(defT[k]));
    }
}

void HyperreductionConfiguration::setTracerGradOnMagicPoints(Foam::PtrList<Foam::volVectorField>& defT, Foam::label nModes)
{
    m_tracerGradOnMagicPoints.resize(nModes);
    for (Foam::label k = 0; k < nModes; k++)
    {
        m_tracerGradOnMagicPoints.set(k, new Foam::volVectorField(defT[k]));
    }
}

void HyperreductionConfiguration::setDeformationTensorOnMagicNeighborhoods(Foam::PtrList<Foam::volTensorField>& defT, Foam::label nModes)
{
    m_deformationTensorOnMagicNeighborhoods.resize(nModes);
    for (Foam::label k = 0; k < nModes; k++)
    {
        m_deformationTensorOnMagicNeighborhoods.set(k, new Foam::volTensorField(defT[k]));
    }
}

void HyperreductionConfiguration::setDeformationTensorOnMagicPoints(Foam::PtrList<Foam::volTensorField>& defT, Foam::label nModes)
{
    m_deformationTensorOnMagicPoints.resize(nModes);
    for (Foam::label k = 0; k < nModes; k++)
    {
        m_deformationTensorOnMagicPoints.set(k, new Foam::volTensorField(defT[k]));
    }
}

void HyperreductionConfiguration::setMeanVectorDEIM(const Foam::volVectorField& mean)
{
    m_meanVectorDEIM = new Foam::volVectorField(mean);
}

void HyperreductionConfiguration::setMeanScalarDEIM(const Foam::volScalarField& mean)
{
    m_meanScalarDEIM = new Foam::volScalarField(mean);
}

void HyperreductionConfiguration::setMeanVectorDEIMMagic(const Foam::volVectorField& mean)
{
    m_meanVectorDEIMMagic = new Foam::volVectorField(mean);
}

void HyperreductionConfiguration::setMeanScalarDEIMMagic(const Foam::volScalarField& mean)
{
    m_meanScalarDEIMMagic = new Foam::volScalarField(mean);
}

void HyperreductionConfiguration::setDeltaWeight(const Eigen::VectorXd& dw)
{
    m_deltaWeight = dw;
}

void HyperreductionConfiguration::set_magicDelta(const Foam::volScalarField& mD) 
{
    m_magicDelta = new Foam::GeometricField<Foam::scalar, Foam::fvPatchField, Foam::volMesh>(mD); 
}

void HyperreductionConfiguration::initializeHyperreduction(std::unique_ptr<PODConfiguration>& PODConfig)
{
    if (m_HRMethod == "GappyDEIM")
    {
        m_HRMethod = "DEIM";
    }

    if (m_interpFieldCentered)
    {
        m_folderDEIM = "./ITHACAoutput/Hyperreduction/" + m_HRMethod + "_centered/";
    }
    else
    {
        m_folderDEIM = "./ITHACAoutput/Hyperreduction/" + m_HRMethod + "/";
    }

    std::string nbModesUInFolderDEIM = "";
    std::string ECPAlgoInFolderDEIM = "";

    bool isReducedField = ITHACAutilities::containsSubstring(m_HRInterpolatedField, "reduced");
    if (isReducedField)
    {
        initializeReducedField(m_HRInterpolatedField, PODConfig);
        nbModesUInFolderDEIM = std::to_string(PODConfig->nModes()["U"]) + "modesU/";
    }
    m_nMagicPoints = PODConfig->nModes()[m_HRInterpolatedField];
    m_HRSnapshotsField = m_HRInterpolatedField;

    if (m_HRMethod == "ECP")
    {
        if (!(m_ECPAlgo == "Global" || m_ECPAlgo == "EachMode"))
        {
            Foam::Info << "Error: ECPAlgo must be Global or EachMode" << Foam::endl;
            abort();
        }

        if (m_HRInterpolatedField == "fullStressFunction")
        {
            m_HRSnapshotsField = "projFullStressFunction_" + std::to_string(PODConfig->nModes()["U"]) + "modes";
        }
        else if (ITHACAutilities::containsSubstring(m_HRInterpolatedField, "reducedFullStressFunction"))
        {
            m_HRSnapshotsField = "projReducedFullStressFunction_" + std::to_string(PODConfig->nModes()["U"]) + "modes";
        }
        else if (m_HRInterpolatedField == "nut")
        {
            m_HRSnapshotsField = "projSmagFromNut_" + std::to_string(PODConfig->nModes()["U"]) + "modes";
        }
        else if (ITHACAutilities::containsSubstring(m_HRInterpolatedField, "reducedNut"))
        {
            m_HRSnapshotsField = "projSmagFromReducedNut_" + std::to_string(PODConfig->nModes()["U"]) + "modes";
        }
        else
        {
            Foam::Info << "Error: ECP method is coded only for fullStressFunction and nut" << Foam::endl;
            abort();
        }

        ECPAlgoInFolderDEIM = m_ECPAlgo + "/";
        nbModesUInFolderDEIM = std::to_string(PODConfig->nModes()["U"]) + "modesU/";

        int nECPFields = 1;
        if (m_ECPAlgo == "EachMode")
        {
            nECPFields = PODConfig->nModes()["U"];
        }

        for (int c = 0; c < nECPFields; c++)
        {
            for (int k = 0; k <= (c + 1) * (m_ECPAlgo == "EachMode") * (ITHACAutilities::containsSubstring(m_HRInterpolatedField, "nut")); k++)
            {
                Foam::word fieldNameModec = m_HRSnapshotsField;
                if (m_ECPAlgo == "EachMode")
                {
                    fieldNameModec += "_" + std::to_string(c + 1);
                    if (ITHACAutilities::containsSubstring(m_HRInterpolatedField, "nut"))
                    {
                        fieldNameModec += "_" + std::to_string(k);
                    }
                }
                PODConfig->appendField(fieldNameModec, "scalar");
                PODConfig->insert_nModes(fieldNameModec, PODConfig->nModes()[m_HRInterpolatedField]);
                PODConfig->insert_hilbertSpacePOD(PODConfig->field_name().last(), "L2");
                PODConfig->set_varyingEnergy(PODConfig->field_name().last(), 0);
                PODConfig->set_resolvedVaryingEnergy(PODConfig->field_name().last(), 0);
            }
        }
    }

    m_folderDEIM += m_HRInterpolatedField.substr(0, m_HRInterpolatedField.find("_")) + "/" + nbModesUInFolderDEIM;
    m_folderDEIM += std::to_string(m_nMagicPoints) + "magicPoints/" + ECPAlgoInFolderDEIM;
}

void HyperreductionConfiguration::initializeReducedField(const Foam::word origField,
                                       std::unique_ptr<PODConfiguration>& PODConfig)
{
    Foam::word nameToReplace = origField.substr(7);
    nameToReplace[0] = tolower(nameToReplace[0]);
    Foam::word modifiedFieldName = createReducedFieldName(origField, PODConfig);
    m_HRInterpolatedField = modifiedFieldName;
    addReducedFieldToPOD(nameToReplace, modifiedFieldName, PODConfig);
}

void HyperreductionConfiguration::addReducedFieldToPOD(const Foam::word& nameToReplace,
     const Foam::word& modifiedFieldName, std::unique_ptr<PODConfiguration>& PODConfig)
{
    const auto& field_type = PODConfig->field_type();
    const auto& field_name = PODConfig->field_name();

    PODConfig->appendField(modifiedFieldName , field_type[find(field_name.begin(), field_name.end(), nameToReplace) - field_name.begin()]);
    PODConfig->insert_nModes(modifiedFieldName, PODConfig->nModes()[nameToReplace]);
    PODConfig->insert_hilbertSpacePOD(modifiedFieldName, PODConfig->hilbertSpacePOD()[nameToReplace]);
    PODConfig->set_varyingEnergy(modifiedFieldName, 0);
    PODConfig->set_resolvedVaryingEnergy(modifiedFieldName, 0);
}

Foam::word HyperreductionConfiguration::createReducedFieldName(const Foam::word& fieldName,
                                              std::unique_ptr<PODConfiguration>& PODConfig)
{
    return fieldName + "_" + std::to_string(PODConfig->nModes()["U"]) + "modes";
}
