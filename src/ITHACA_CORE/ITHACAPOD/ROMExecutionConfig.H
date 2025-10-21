#pragma once
#ifndef ROMExecutionConfig_H
#define ROMExecutionConfig_H

#include <Eigen/Core>
#include <volFieldsFwd.H>


#include "dimensionedScalarFwd.H"


namespace ITHACAPOD
{
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
        m_nu(nullptr), m_meanU(nullptr), m_deltaWeight(nullptr) { }

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

    // Hyper-reduction

    const Foam::word& HRSnapshotsField() const { return m_HRSnapshotsField; }
    const Foam::word& HRInterpolatedField() const { return m_HRInterpolatedField; }
    const Foam::word& ECPAlgo() const { return m_ECPAlgo; }
    const Foam::word& useSOTA() const { return m_useSOTA; }

    void setHRInterpolatedField(const Foam::word& field) { m_HRInterpolatedField = field; }
    void setHRSnapshotsField(const Foam::word& field) { m_HRSnapshotsField = field; }
    void setECPAlgo(const Foam::word& algo) { m_ECPAlgo = algo; }
    void setUseSOTA(const Foam::word& sota) { m_useSOTA = sota; }

    // Delta weight for hyper-reduction
    const Eigen::VectorXd* deltaWeight() const { return m_deltaWeight; }
    void setDeltaWeight(const Eigen::VectorXd& dw);

private:
    Foam::word m_ROMTemporalScheme;
    PressureResolutionKind m_pressureResolutionKind;

    Foam::dimensionedScalar* m_nu;
    Foam::volVectorField* m_meanU;

    // Hyper-reduction settings
    Foam::word m_HRInterpolatedField;
    Foam::word m_HRSnapshotsField;
    Foam::word m_ECPAlgo;
    Foam::word m_useSOTA;

    Eigen::VectorXd* m_deltaWeight;
};


} // namespace ITHACAPOD

#endif // ROMExecutionConfig_H
