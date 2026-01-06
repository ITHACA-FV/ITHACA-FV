#ifndef FIELDPARAMETERS_H
#define FIELDPARAMETERS_H
#pragma once

#include "volFieldsFwd.H"

/**
 * @brief Stores template fields for velocity, pressure, and turbulence quantities
 */
class FieldTemplates
{
public:
    FieldTemplates() : m_U(nullptr), m_p(nullptr),
                       m_nut(nullptr), m_omega(nullptr),
                       m_k(nullptr), m_fullStressFunction(nullptr) {}

    const Foam::volVectorField& get_U() const { return *m_U; }
    const Foam::volScalarField& get_p() const { return *m_p; }
    Foam::volScalarField& get_nut() const { return *m_nut; }
    const Foam::volScalarField& get_omega() const { return *m_omega; }
    const Foam::volScalarField& get_k() const { return *m_k; }
    const Foam::volVectorField& get_fullStressFunction() const { return *m_fullStressFunction; }
    
    void set_U(Foam::volVectorField* u) { m_U = u; }
    void set_p(Foam::volScalarField* p) { m_p = p; }
    void set_nut(Foam::volScalarField* nut) { m_nut = nut; }
    void set_omega(Foam::volScalarField* omega) { m_omega = omega; }
    void set_k(Foam::volScalarField* k) { m_k = k; }
    void set_fullStressFunction(Foam::volVectorField& templateSmag); 


private:
    Foam::volVectorField* m_U;
    Foam::volScalarField* m_p;
    Foam::volScalarField* m_nut;
    Foam::volScalarField* m_omega;
    Foam::volScalarField* m_k;
    Foam::volVectorField* m_fullStressFunction;
};

#endif // MESHPARAMETERS_H
