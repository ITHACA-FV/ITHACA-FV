#pragma once
#ifndef DEIMConfiguration_H
#define DEIMConfiguration_H

#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>
#include "List.H"
#include "labelFwd.H"
#include "volFieldsFwd.H"


namespace ITHACAPOD
{

/**
 * @brief Discrete Empirical Interpolation Method settings
 *
 * Manages magic points, interpolation matrices, and DEIM-specific fields
 * for hyper-reduction of nonlinear terms.
 */
class DEIMConfiguration
{
public:
    DEIMConfiguration():
        m_nMagicPoints(0), m_Ck(0.0), m_Ce(0.0),
        m_meanVectorDEIM(nullptr), m_meanScalarDEIM(nullptr),
        m_meanVectorDEIMMagic(nullptr), m_meanScalarDEIMMagic(nullptr) { };

    // Magic points
    const Foam::label& nMagicPoints() const { return m_nMagicPoints; }
    const Foam::List<Foam::label>& magicPoints() const { return m_magicPoints; }
    const Foam::List<Foam::label>& localMagicPoints() const { return m_localMagicPoints; }

    void setNMagicPoints(Foam::label n) { m_nMagicPoints = n; }
    void setMagicPoints(const Foam::List<Foam::label>& points) { m_magicPoints = points; }
    void setLocalMagicPoints(const Foam::List<Foam::label>& points) { m_localMagicPoints = points; }

    // Turbulence model coefficients
    const float& Ck() const { return m_Ck; }
    const float& Ce() const { return m_Ce; }
    void setCk(float ck) { m_Ck = ck; }
    void setCe(float ce) { m_Ce = ce; }

    // DEIM folder
    const Foam::word& folderDEIM() const { return m_folderDEIM; }
    void setFolderDEIM(const Foam::word& folder) { m_folderDEIM = folder; }

    // Deformation tensor modes
    const Foam::PtrList<Foam::volTensorField>& deformationTensorOnMagicNeighborhoods() const
    {
        return m_deformationTensorOnMagicNeighborhoods;
    }
    const Foam::PtrList<Foam::volTensorField>& deformationTensorOnMagicPoints() const
    {
        return m_deformationTensorOnMagicPoints;
    }
    const Foam::PtrList<Foam::volVectorField>& tracerGradOnMagicNeighborhoods() const
    {
        return m_tracerGradOnMagicNeighborhoods;
    }
    const Foam::PtrList<Foam::volVectorField>& tracerGradOnMagicPoints() const
    {
        return m_tracerGradOnMagicPoints;
    }

    void setDeformationTensorOnMagicNeighborhoods(Foam::PtrList<Foam::volTensorField>& defT, Foam::label nModes);
    void setDeformationTensorOnMagicPoints(Foam::PtrList<Foam::volTensorField>& defT, Foam::label nModes);
    void setTracerGradOnMagicNeighborhoods(Foam::PtrList<Foam::volVectorField>& defT, Foam::label nModes);
    void setTracerGradOnMagicPoints(Foam::PtrList<Foam::volVectorField>& defT, Foam::label nModes);

    // DEIM matrices
    const Eigen::MatrixXd& K_DEIM() const { return m_K_DEIM; }
    const Eigen::MatrixXd& projectedK_DEIM() const { return m_projectedK_DEIM; }
    const Eigen::MatrixXd& projectedMK_DEIM() const { return m_projectedMK_DEIM; }

    void setK_DEIM(const Eigen::MatrixXd& K) { m_K_DEIM = K; }
    void setProjectedK_DEIM(const Eigen::MatrixXd& K) { m_projectedK_DEIM = K; }
    void setProjectedMK_DEIM(const Eigen::MatrixXd& MK) { m_projectedMK_DEIM = MK; }

    // DEIM tensors for turbulence
    const Eigen::Tensor<double, 3>& xiOnMagicPts() const { return m_xiOnMagicPts; }
    const Eigen::Tensor<double, 3>& M_xiOnMagicPts() const { return m_M_xiOnMagicPts; }

    void setXiOnMagicPts(const Eigen::Tensor<double, 3>& xi) { m_xiOnMagicPts = xi; }
    void setM_XiOnMagicPts(const Eigen::Tensor<double, 3>& M_xi) { m_M_xiOnMagicPts = M_xi; }

    // Mean fields for DEIM
    const Foam::volVectorField* meanVectorDEIM() const { return m_meanVectorDEIM; }
    const Foam::volScalarField* meanScalarDEIM() const { return m_meanScalarDEIM; }
    const Foam::volVectorField* meanVectorDEIMMagic() const { return m_meanVectorDEIMMagic; }
    const Foam::volScalarField* meanScalarDEIMMagic() const { return m_meanScalarDEIMMagic; }

    void setMeanVectorDEIM(const Foam::volVectorField& mean);
    void setMeanScalarDEIM(const Foam::volScalarField& mean);
    void setMeanVectorDEIMMagic(const Foam::volVectorField& mean);
    void setMeanScalarDEIMMagic(const Foam::volScalarField& mean);

private:
    Foam::label m_nMagicPoints;
    Foam::List<Foam::label> m_magicPoints;
    Foam::List<Foam::label> m_localMagicPoints;

    float m_Ck; // Smagorinsky coefficient
    float m_Ce; // Smagorinsky coefficient

    Foam::word m_folderDEIM;

    // Mode fields on magic points
    Foam::PtrList<Foam::volTensorField> m_deformationTensorOnMagicPoints;
    Foam::PtrList<Foam::volTensorField> m_deformationTensorOnMagicNeighborhoods;
    Foam::PtrList<Foam::volVectorField> m_tracerGradOnMagicPoints;
    Foam::PtrList<Foam::volVectorField> m_tracerGradOnMagicNeighborhoods;

    // DEIM matrices
    Eigen::MatrixXd m_K_DEIM;
    Eigen::MatrixXd m_projectedK_DEIM;
    Eigen::MatrixXd m_projectedMK_DEIM;

    // Turbulence tensors
    Eigen::Tensor<double, 3> m_xiOnMagicPts;
    Eigen::Tensor<double, 3> m_M_xiOnMagicPts;

    // Mean fields
    Foam::volVectorField* m_meanVectorDEIM;
    Foam::volScalarField* m_meanScalarDEIM;
    Foam::volVectorField* m_meanVectorDEIMMagic;
    Foam::volScalarField* m_meanScalarDEIMMagic;
};

} // namespace ITHACAPOD

#endif // DEIMConfiguration_H
