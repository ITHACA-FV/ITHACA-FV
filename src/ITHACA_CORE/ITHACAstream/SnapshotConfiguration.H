#ifndef SNAPSHOTCONFIGURATION_H
#define SNAPSHOTCONFIGURATION_H
#include "fileName.H"
#include "labelFwd.H"
#pragma once



// ============================================================================
// Snapshot Configuration
// ============================================================================
/**
 * @brief Manages snapshot and time-related parameters
 */
class SnapshotConfiguration
{
public:
    SnapshotConfiguration():
        m_startTime(0), m_endTime(0), m_nSnapshots(0),
        m_endTimeSimulation(0), m_nSnapshotsSimulation(0),
        m_InitialTime(0.0), m_FinalTime(0.0),
        m_FinalTimeSimulation(0.0), m_writeInterval(0.0), m_nSimu(1) { }

    const Foam::fileName& get_casenameData() const { return m_casenameData; }

    const Foam::label& get_startTime() const { return m_startTime; }
    const Foam::label& get_endTime() const { return m_endTime; }
    const Foam::label& get_nSnapshots() const { return m_nSnapshots; }

    const Foam::label& get_endTimeSimulation() const { return m_endTimeSimulation; }
    const Foam::label& get_nSnapshotsSimulation() const { return m_nSnapshotsSimulation; }

    const double& get_initialTime() const { return m_InitialTime; }
    const double& get_finalTime() const { return m_FinalTime; }
    const double& get_finalTimeSimulation() const { return m_FinalTimeSimulation; }
    const double& get_writeInterval() const { return m_writeInterval; }
    const Foam::label& get_nSimu() const { return m_nSimu; }

    void set_casenameData(const Foam::fileName& cn) { m_casenameData = cn; }
    void set_startTime(const Foam::label& st) { m_startTime = st; }
    void set_endTime(const Foam::label& et) { m_endTime = et; }
    void set_nSnapshots(const Foam::label& ns) { m_nSnapshots = ns; }
    void set_endTimeSimulation(const Foam::label& ets) { m_endTimeSimulation = ets; }
    void set_nSnapshotsSimulation(const Foam::label& nss) { m_nSnapshotsSimulation = nss; }
    void set_initialTime(const double& it) { m_InitialTime = it; }
    void set_finalTime(const double& ft) { m_FinalTime = ft; }
    void set_finalTimeSimulation(const double& fts) { m_FinalTimeSimulation = fts; }
    void set_writeInterval(const double& wi) { m_writeInterval = wi; }
    void set_nSimu(const Foam::label& ns) { m_nSimu = ns; }

private:
    Foam::fileName m_casenameData;

    Foam::label m_startTime;
    Foam::label m_endTime;
    Foam::label m_nSnapshots;

    Foam::label m_endTimeSimulation;
    Foam::label m_nSnapshotsSimulation;

    double m_InitialTime;
    double m_FinalTime;
    double m_FinalTimeSimulation;
    double m_writeInterval;

    Foam::label m_nSimu; // Relationship between the simulation time step and the time step used for the ROM
};

#endif
