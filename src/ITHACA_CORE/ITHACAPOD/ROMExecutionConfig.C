#include "ROMExecutionConfig.H"
#include "GeometricField.H"
#include "fvPatchField.H"

void ITHACAPOD::ROMExecutionConfig::setROMTemporalScheme(const Foam::word& scheme)
{
    m_ROMTemporalScheme = scheme;
}

void ITHACAPOD::ROMExecutionConfig::setPressureResolutionKind(PressureResolutionKind kind)
{
    m_pressureResolutionKind = kind;
}

void ITHACAPOD::ROMExecutionConfig::setMeanU(const Foam::volVectorField& mean)
{
    m_meanU = new Foam::volVectorField(mean);
}
