#ifndef MESHPARAMETERS_H
#define MESHPARAMETERS_H
#pragma once

#include "fvMesh.H"
#include "volFieldsFwd.H"


/**
 * @brief Manages mesh-related parameters and geometric information
 */
class MeshConfiguration
{
public:
    MeshConfiguration() : m_mesh(nullptr), m_nCells(0), m_volume(nullptr),
                          m_totalVolume(0.0), m_delta(nullptr){}

    Foam::fvMesh& get_mesh() const { return *m_mesh; }
    const Foam::label& get_nCells() const { return m_nCells; }
    
    Foam::volScalarField& get_volume() const { return *m_volume; }
    const double& get_totalVolume() const { return m_totalVolume; }
    const Foam::volScalarField& get_delta() const { return *m_delta; }
    
    void set_mesh(Foam::fvMesh* m) { m_mesh = m; }
    void set_nCells(Foam::label n) { m_nCells = n; }
    void set_volume(Foam::volScalarField* v) { m_volume = v; }
    void set_totalVolume(double tv) { m_totalVolume = tv; }
    void set_delta(Foam::volScalarField* d) { m_delta = d; }

private:
    Foam::fvMesh* m_mesh;
    Foam::label m_nCells;
    Foam::volScalarField* m_volume;
    double m_totalVolume;
    Foam::volScalarField* m_delta;
};

#endif // MESHPARAMETERS_H
