#include "DEIMConfiguration.H"
#include "GeometricField.H"
#include "fvPatchField.H"
#include "volFieldsFwd.H"

void ITHACAPOD::DEIMConfiguration::setTracerGradOnMagicNeighborhoods(Foam::PtrList<Foam::volVectorField>& defT, Foam::label nModes)
{
    m_tracerGradOnMagicNeighborhoods.resize(nModes);
    for (Foam::label k = 0; k < nModes; k++)
    {
        m_tracerGradOnMagicNeighborhoods.set(k, new Foam::volVectorField(defT[k]));
    }
}

void ITHACAPOD::DEIMConfiguration::setTracerGradOnMagicPoints(Foam::PtrList<Foam::volVectorField>& defT, Foam::label nModes)
{
    m_tracerGradOnMagicPoints.resize(nModes);
    for (Foam::label k = 0; k < nModes; k++)
    {
        m_tracerGradOnMagicPoints.set(k, new Foam::volVectorField(defT[k]));
    }
}

void ITHACAPOD::DEIMConfiguration::setDeformationTensorOnMagicNeighborhoods(Foam::PtrList<Foam::volTensorField>& defT, Foam::label nModes)
{
    m_deformationTensorOnMagicNeighborhoods.resize(nModes);
    for (Foam::label k = 0; k < nModes; k++)
    {
        m_deformationTensorOnMagicNeighborhoods.set(k, new Foam::volTensorField(defT[k]));
    }
}

void ITHACAPOD::DEIMConfiguration::setDeformationTensorOnMagicPoints(Foam::PtrList<Foam::volTensorField>& defT, Foam::label nModes)
{
    m_deformationTensorOnMagicPoints.resize(nModes);
    for (Foam::label k = 0; k < nModes; k++)
    {
        m_deformationTensorOnMagicPoints.set(k, new Foam::volTensorField(defT[k]));
    }
}

void ITHACAPOD::DEIMConfiguration::setMeanVectorDEIM(const Foam::volVectorField& mean)
{
    m_meanVectorDEIM = new Foam::volVectorField(mean);
}

void ITHACAPOD::DEIMConfiguration::setMeanScalarDEIM(const Foam::volScalarField& mean)
{
    m_meanScalarDEIM = new Foam::volScalarField(mean);
}

void ITHACAPOD::DEIMConfiguration::setMeanVectorDEIMMagic(const Foam::volVectorField& mean)
{
    m_meanVectorDEIMMagic = new Foam::volVectorField(mean);
}

void ITHACAPOD::DEIMConfiguration::setMeanScalarDEIMMagic(const Foam::volScalarField& mean)
{
    m_meanScalarDEIMMagic = new Foam::volScalarField(mean);
}
