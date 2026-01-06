#pragma once
#ifndef ROMExecutionConfig_H
#define ROMExecutionConfig_H

#include <Eigen/Eigen>
#include <volFieldsFwd.H>
#include "dimensionedScalarFwd.H"


/**
 * @brief Runtime settings for ROM execution
 *
 * Controls temporal schemes, pressure resolution, export options,
 * and various simulation flags.
 */
enum class PressureResolutionKind {
    FullOrder,
    ReducedOrder,
    Neglected,
    Undefined
};

class ROMExecutionConfig
{
public:
    ROMExecutionConfig():
        m_nu(nullptr), m_meanU(nullptr) { }

    // Temporal scheme
    const Foam::word& ROMTemporalScheme() const { return m_ROMTemporalScheme; }
    void setROMTemporalScheme(const Foam::word& scheme);

    // Pressure handling
    const PressureResolutionKind& pressureResolutionKind() const { return m_pressureResolutionKind; }
    void setPressureResolutionKind(PressureResolutionKind kind);

    // Physical properties
    const Foam::dimensionedScalar& nu() const { return *m_nu; }
    void setNu(Foam::dimensionedScalar* nu) { m_nu = nu; }

    // Mean field
    const Foam::volVectorField* meanU() const { return m_meanU; }
    void setMeanU(const Foam::volVectorField& mean);

    // SOTA method
    const Foam::word& useSOTA() const { return m_useSOTA; }
    void setUseSOTA(const Foam::word& sota) { m_useSOTA = sota; }


private:
    Foam::word m_ROMTemporalScheme;
    PressureResolutionKind m_pressureResolutionKind;

    Foam::dimensionedScalar* m_nu;
    Foam::volVectorField* m_meanU;

    Foam::word m_useSOTA;

};

#endif // ROMExecutionConfig_H
