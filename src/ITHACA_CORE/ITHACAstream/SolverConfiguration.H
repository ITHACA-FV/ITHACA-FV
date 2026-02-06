#ifndef SolverConfiguration_H
#define SolverConfiguration_H

#include <labelFwd.H>
#include <word.H>

/**
 * @brief Controls which eigensolver algorithm is used and output formatting.
 */
class SolverConfiguration
{
public:
    SolverConfiguration() : m_precision(0), m_nBlocks(0), m_centeredOrNot(false) {}

    const Foam::word& get_eigensolver() const { return m_eigensolver; }
    const int& get_precision() const { return m_precision; }
    const std::_Ios_Fmtflags& get_outytpe() const { return m_outytpe; }
    const Foam::label& get_nBlocks() const { return m_nBlocks; }
    const bool& get_centeredOrNot() const { return m_centeredOrNot; }
    
    void set_eigensolver(const Foam::word& es) { m_eigensolver = es; }
    void set_precision(const int& p) { m_precision = p; }
    void set_outytpe(const std::_Ios_Fmtflags& ot) { m_outytpe = ot; }
    void set_nBlocks(const Foam::label& nb) { m_nBlocks = nb; }
    void set_centeredOrNot(const bool& c) { m_centeredOrNot = c; }

private:
    Foam::word m_eigensolver;
    int m_precision;
    std::_Ios_Fmtflags m_outytpe;
    Foam::label m_nBlocks;
    bool m_centeredOrNot;
};

#endif  //SolverConfiguration_H
