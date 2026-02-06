#pragma once
#ifndef PODConfiguration_H
#define PODConfiguration_H
#include <Eigen/Eigen>
#include "HashTable.H"

/**
 * @brief Manages POD-specific parameters including modes, energy, and Hilbert spaces
 */
class PODConfiguration
{
public:
    PODConfiguration():
        m_weightH1(0.0), m_weightBC(0.0) { }

    const Foam::wordList& fieldlist() const { return m_fieldlist; }
    const Foam::HashTable<Foam::label, Foam::word>& nModes() const { return m_nModes; }
    const Foam::HashTable<Foam::word, Foam::word>& hilbertSpacePOD() const { return m_hilbertSpacePOD; }
    const Foam::List<Foam::word>& field_name() const { return m_field_name; }
    const Foam::List<Foam::word>& field_type() const { return m_field_type; }

    const Foam::HashTable<double, Foam::word>& resolvedVaryingEnergy() const { return m_resolvedVaryingEnergy; }
    const Foam::HashTable<double, Foam::word>& varyingEnergy() const { return m_varyingEnergy; }
    const Foam::HashTable<double, Foam::word>& meanEnergy() const { return m_meanEnergy; }

    const Eigen::VectorXd& eigenValues_U() const { return m_eigenValues_U; }
    const Eigen::VectorXd& eigenValues_p() const { return m_eigenValues_p; }
    const Eigen::VectorXd& lambda() const { return m_lambda; }

    const double& weightH1() const { return m_weightH1; }
    const double& weightBC() const { return m_weightBC; }
    const Foam::word& patchBC() const { return m_patchBC; }

    void set_fieldlist(const Foam::wordList& fl) { m_fieldlist = fl; }
    void insert_nModes(const Foam::word& field_name, const Foam::label& n) { m_nModes.insert(field_name, n); }
    void insert_hilbertSpacePOD(const Foam::word& field_name, const Foam::word& hilbertSp)
    {
        m_hilbertSpacePOD.insert(field_name, hilbertSp);
    }
    void set_field_name(const Foam::List<Foam::word>& fn) { m_field_name = fn; }
    void set_field_type(const Foam::List<Foam::word>& ft) { m_field_type = ft; }

    void set_resolvedVaryingEnergy(const Foam::word& field_name, const double& n)
    {
        m_resolvedVaryingEnergy.set(field_name, n);
    }

    void insert_resolvedVaryingEnergy(const Foam::word& field_name, const double& n)
    {
        m_resolvedVaryingEnergy.insert(field_name, n);
    }

    void set_varyingEnergy(const Foam::word& field_name, const double& n)
    {
        m_varyingEnergy.set(field_name, n);
    }

    void insert_varyingEnergy(const Foam::word& field_name, const double& n)
    {
        m_varyingEnergy.insert(field_name, n);
    }

    void insert_meanEnergy(const Foam::word& field_name, const double& n)
    {
        m_meanEnergy.insert(field_name, n);
    }

    void set_eigenValues_U(const Eigen::VectorXd& ev) { m_eigenValues_U = ev; }
    void set_eigenValues_p(const Eigen::VectorXd& ev) { m_eigenValues_p = ev; }
    void set_lambda(const Eigen::VectorXd& l) { m_lambda = l; }

    void set_weightH1(const double& w) { m_weightH1 = w; }
    void set_weightBC(const double& w) { m_weightBC = w; }
    void set_patchBC(const Foam::word& p) { m_patchBC = p; }

    void appendField(const Foam::word& fieldName, const Foam::word& fieldType)
    {
        m_field_name.append(fieldName);
        m_field_type.append(fieldType);
    }

private:
    Foam::wordList m_fieldlist;
    Foam::HashTable<Foam::label, Foam::word> m_nModes;
    Foam::HashTable<Foam::word, Foam::word> m_hilbertSpacePOD;
    Foam::List<Foam::word> m_field_name;
    Foam::List<Foam::word> m_field_type;

    Foam::HashTable<double, Foam::word> m_resolvedVaryingEnergy;
    Foam::HashTable<double, Foam::word> m_varyingEnergy;
    Foam::HashTable<double, Foam::word> m_meanEnergy;

    Eigen::VectorXd m_eigenValues_U;
    Eigen::VectorXd m_eigenValues_p;
    Eigen::VectorXd m_lambda;

    double m_weightH1;
    double m_weightBC;
    Foam::word m_patchBC;
};

#endif
