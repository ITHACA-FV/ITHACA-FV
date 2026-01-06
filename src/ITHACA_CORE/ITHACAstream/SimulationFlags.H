#pragma once
#include "word.H"
#ifndef SimulationFlags_H
#define SimulationFlags_H

/**
 * @brief Manages boolean flags for various simulation options
 */
class SimulationFlags
{
public:
    SimulationFlags():
        m_exportPython(false), m_exportMatlab(false), m_exportTxt(false),
        m_useDDES(false), m_useDNS(false), m_onLineReconstruct(false),
        m_forcing(false), m_symDiff(false), m_useDEIM(false) { }

    const bool& exportPython() const { return m_exportPython; }
    const bool& exportMatlab() const { return m_exportMatlab; }
    const bool& exportTxt() const { return m_exportTxt; }
    const bool& useDDES() const { return m_useDDES; }
    const bool& useDNS() const { return m_useDNS; }
    const bool& onLineReconstruct() const { return m_onLineReconstruct; }
    const bool& forcing() const { return m_forcing; }
    const bool& symDiff() const { return m_symDiff; }
    const bool& useDEIM() const { return m_useDEIM; }

    void set_exportPython(const bool& b) { m_exportPython = b; }
    void set_exportMatlab(const bool& b) { m_exportMatlab = b; }
    void set_exportTxt(const bool& b) { m_exportTxt = b; }
    void set_useDDES(const bool& b) { m_useDDES = b; }
    void set_useDNS(const bool& b) { m_useDNS = b; }
    void set_onLineReconstruct(const bool& b) { m_onLineReconstruct = b; }
    void set_forcingOrNot(const bool& b) { m_forcing = b; }
    void set_symDiff(const bool& b) { m_symDiff = b; }
    void set_useDEIM(const bool& b) { m_useDEIM = b; }

private:
    Foam::word m_HRMethod;
    bool m_exportPython;
    bool m_exportMatlab;
    bool m_exportTxt;
    bool m_useDDES;
    bool m_useDNS;
    bool m_onLineReconstruct;
    bool m_forcing;
    bool m_symDiff;
    bool m_useDEIM;
};

#endif
