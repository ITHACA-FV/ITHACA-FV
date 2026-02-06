#include "ROMExecutionConfig.H"
#include "GeometricField.H"
#include "fvPatchField.H"

void ROMExecutionConfig::setROMTemporalScheme(const Foam::word& scheme)
{
    m_ROMTemporalScheme = scheme;
}

void ROMExecutionConfig::setPressureResolutionKind(PressureResolutionKind kind)
{
    m_pressureResolutionKind = kind;
}

void ROMExecutionConfig::setMeanU(const Foam::volVectorField& mean)
{
    m_meanU = new Foam::volVectorField(mean);
}
