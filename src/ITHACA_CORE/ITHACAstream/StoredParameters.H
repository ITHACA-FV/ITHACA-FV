
#ifndef POD_PARAMETERS_H
#define POD_PARAMETERS_H

#include <FieldField.H>

#include "FieldTemplates.H"
#include "PODConfiguration.H"
#include "SnapshotConfiguration.H"
#include "MeshConfiguration.H"
#include "SolverConfiguration.H"
#include "HyperreductionConfiguration.H"
#include "ROMExecutionConfig.H"
#include "SimulationFlags.H"


#include "fvMesh.H"

#include "ITHACAparameters.H"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#include <Eigen/Eigen>
#pragma GCC diagnostic pop


class Parameters
{
public:
    virtual ~Parameters() { }
};
// only one derivate
/**
 *  @brief  Class that contains all parameters of the stochastic POD
 ***********************************************/
class StoredParameters : public Parameters
{
public:
    StoredParameters(int argc, char* argv[]);

    /// Getter to access all m_parameters members in the other classes
    ITHACAparameters* get_ITHACAparameters() const
    {
        return ithacaLibraryParameters;
    }

    const Foam::HashTable<Foam::label, Foam::word>& get_nModes() const
    {
        return PODConfiguration_->nModes();
    }

    const Foam::HashTable<Foam::word, Foam::word>& get_hilbertSpacePOD() const
    {
        return PODConfiguration_->hilbertSpacePOD(); // Contains both the field name and Hilbert space
    }

    Foam::word get_pathHilbertSpace_fromHS(Foam::word hilbertSp);
    const double& get_weightH1() const
    {
        return PODConfiguration_->weightH1();
    }
    const double& get_weightPOD() const
    {
        return PODConfiguration_->weightBC();
    }
    const Foam::word& get_patchBC() const
    {
        return PODConfiguration_->patchBC();
    }
    const Foam::fileName& get_casenameData() const
    {
        return snapshotConfiguration_->get_casenameData();
    }

    const Foam::label& get_startTime() const
    {
        return snapshotConfiguration_->get_startTime();
    }
    const Foam::label& get_endTime() const
    {
        return snapshotConfiguration_->get_endTime();
    }
    const Foam::label& get_nSnapshots() const
    {
        return snapshotConfiguration_->get_nSnapshots();
    }

    const Foam::label& get_endTimeSimulation() const
    {
        return snapshotConfiguration_->get_endTimeSimulation();
    }
    const Foam::label& get_nSnapshotsSimulation() const
    {
        return snapshotConfiguration_->get_nSnapshotsSimulation();
    }

    const Foam::label& get_nBlocks() const
    {
        return solverConfiguration_->get_nBlocks();
    }

    bool get_centeredOrNot() const
    {
        return solverConfiguration_->get_centeredOrNot();
    }
    int get_precision() const
    {
        return solverConfiguration_->get_precision();
    }
    const std::_Ios_Fmtflags& get_outytpe() const
    {
        return solverConfiguration_->get_outytpe();
    }

    const Foam::word& get_eigensolver() const
    {
        return solverConfiguration_->get_eigensolver();
    }

    Foam::fvMesh& get_mesh() const
    {
        return meshConfig_->get_mesh();
    }

    const Foam::label& get_nCells() const
    {
        return meshConfig_->get_nCells();
    }

    const Foam::volVectorField& get_template_field_U() const
    {
        return fieldTemplates_->get_U();
    }
    const Foam::volScalarField& get_template_field_p() const
    {
        return fieldTemplates_->get_p();
    }
    const Foam::volScalarField& get_template_field_nut() const
    {
        return fieldTemplates_->get_nut();
    }
    const Foam::volScalarField& get_template_field_omega() const
    {
        return fieldTemplates_->get_omega();
    }
    const Foam::volScalarField& get_template_field_k() const
    {
        return fieldTemplates_->get_k();
    }
    const Foam::volVectorField& get_template_field_fullStressFunction() const
    {
        return fieldTemplates_->get_fullStressFunction();
    }
    void set_template_field_fullStressFunction(volVectorField& templateSmag)
    {
        fieldTemplates_->set_fullStressFunction(templateSmag);
    }


    const double& get_saveTime() const
    {
        return snapshotConfiguration_->get_writeInterval();
    }

    const word& get_HRMethod() const
    {
        return HyperreductionConfiguration_->HRMethod();
    }
    const Foam::word& get_DEIMInterpolatedField() const
    {
        return HyperreductionConfiguration_->HRInterpolatedField();
    }
    const word& get_HRSnapshotsField() const
    {
        return HyperreductionConfiguration_->HRSnapshotsField();
    }
    const word& get_ECPAlgo() const
    {
        return HyperreductionConfiguration_->ECPAlgo();
    }
    Eigen::VectorXd get_deltaWeight() const
    {
        return HyperreductionConfiguration_->deltaWeight();
    }

    void set_centeredOrNot(const bool b)
    {
        solverConfiguration_->set_centeredOrNot(b);
    }

    /// Setters to the members of StoredParameters that cannot initialize at the beginning of the program
    void set_eigenValues_U(const Eigen::VectorXd& input_eigenValues_U)
    {
        PODConfiguration_->set_eigenValues_U(input_eigenValues_U);
    }
    void set_lambda(const Eigen::VectorXd& input_lambda)
    {
        PODConfiguration_->set_lambda(input_lambda);
    }
    /// Setters of POD and DEIM parameters
    void set_nModes(const Foam::word& field_name, const Foam::label& n)
    {
        PODConfiguration_->insert_nModes(field_name, n);
    }
    void set_resolvedVaryingEnergy(const Foam::word& field_name, const double& n)
    {
        PODConfiguration_->set_resolvedVaryingEnergy(field_name, n);
    }
    void set_varyingEnergy(const Foam::word& field_name, const double& n)
    {
        PODConfiguration_->set_varyingEnergy(field_name, n);
    }
    void set_meanEnergy(const Foam::word& field_name, const double& n)
    {
        PODConfiguration_->insert_meanEnergy(field_name, n);
    }
    void set_hilbertSpacePOD(const Foam::word& field_name,
        const Foam::word& hilbertSp)
    {
        PODConfiguration_->insert_hilbertSpacePOD(field_name, hilbertSp);
    }
    void set_weightH1(const double& c)
    {
        PODConfiguration_->set_weightH1(c);
    }

    const wordList& get_fieldlist() const { return PODConfiguration_->fieldlist(); }
    const HashTable<double, word>& get_resolvedVaryingEnergy() const { return PODConfiguration_->resolvedVaryingEnergy(); }
    const HashTable<double, word>& get_varyingEnergy() const { return PODConfiguration_->varyingEnergy(); }
    const HashTable<double, word>& get_meanEnergy() const { return PODConfiguration_->meanEnergy(); }
    Foam::word get_pathHilbertSpace(Foam::word fieldName);

    const List<word>& get_field_name() const { return PODConfiguration_->field_name(); }
    const List<word>& get_field_type() const { return PODConfiguration_->field_type(); }

    // getting current F_0 force vector and F_mask
    const Foam::Vector<double>& get_F_0() const;
    void compute_F_mask();

    const PressureResolutionKind& get_pressureResolutionKind() { return ROMExecutionConfig_->pressureResolutionKind(); }

    const word& get_ROMTemporalScheme() const { return ROMExecutionConfig_->ROMTemporalScheme(); }

    const double& get_initialTime() const { return snapshotConfiguration_->get_initialTime(); }
    const double& get_initialTimeSimulation() const
    {
        return snapshotConfiguration_->get_finalTime();
    }
    const double& get_finalTime() const
    {
        return snapshotConfiguration_->get_finalTime();
    }
    const double& get_finalTimeSimulation() const
    {
        return snapshotConfiguration_->get_finalTimeSimulation();
    }

    bool get_interpFieldCenteredOrNot() const { return HyperreductionConfiguration_->interpFieldCentered(); }

    bool get_exportPython() const { return SimulationFlags_->exportPython(); }
    bool get_exportMatlab() const { return SimulationFlags_->exportMatlab(); }
    bool get_exportTxt() const { return SimulationFlags_->exportTxt(); }

    const word& get_useSOTA() const { return ROMExecutionConfig_->useSOTA(); }

    const dimensionedScalar& get_nu() const { return ROMExecutionConfig_->nu(); }

    const Eigen::VectorXd& get_eigenValues_U() const { return PODConfiguration_->eigenValues_U(); }
    const Eigen::VectorXd& get_eigenValues_p() const { return PODConfiguration_->eigenValues_p(); }
    const Eigen::VectorXd& get_lambda() const { return PODConfiguration_->lambda(); }

    const label& get_nSimu() const { return snapshotConfiguration_->get_nSimu(); }
    bool get_forcingOrNot() const { return SimulationFlags_->forcing(); }
    bool get_symDiff() const { return SimulationFlags_->symDiff(); }

    /// on line velocity Reconstruction getter
    bool get_onLineReconstruct() const { return SimulationFlags_->onLineReconstruct(); }

    /// DDES getter
    bool get_useDDES() const { return SimulationFlags_->useDDES(); }

    /// DNS getter
    bool get_useDNS() const { return SimulationFlags_->useDNS(); }

    /// Getters of DEIM parameters
    bool get_useDEIM() const { return SimulationFlags_->useDEIM(); }
    const float& get_Ck() const { return HyperreductionConfiguration_->Ck(); }
    const float& get_Ce() const { return HyperreductionConfiguration_->Ce(); }
    const label& get_nMagicPoints() const { return HyperreductionConfiguration_->nMagicPoints(); }
    const word& get_folder_DEIM() const { return HyperreductionConfiguration_->folderDEIM(); }

    const PtrList<volVectorField>& get_tracerGradOfModesOnMagicNeighborhoods() const { return HyperreductionConfiguration_->tracerGradOnMagicNeighborhoods(); }
    const PtrList<volVectorField>& get_tracerGradOfModesOnMagicPoints() const { return HyperreductionConfiguration_->tracerGradOnMagicPoints(); }
    const PtrList<volTensorField>& get_deformationTensorOfModesOnMagicNeighborhoods() const { return HyperreductionConfiguration_->deformationTensorOnMagicNeighborhoods(); }
    const PtrList<volTensorField>& get_deformationTensorOfModesOnMagicPoints() const { return HyperreductionConfiguration_->deformationTensorOnMagicPoints(); }
    const List<label>& get_magicPoints() const { return HyperreductionConfiguration_->magicPoints(); }
    const List<label>& get_localMagicPoints() const { return HyperreductionConfiguration_->localMagicPoints(); }
    const Eigen::MatrixXd& get_projected_K_DEIM() const { return HyperreductionConfiguration_->projectedK_DEIM(); }
    const Eigen::MatrixXd& get_projected_MK_DEIM() const { return HyperreductionConfiguration_->projectedMK_DEIM(); }
    const volScalarField& get_delta() const { return meshConfig_->get_delta(); }
    const volScalarField& get_magicDelta() const { return HyperreductionConfiguration_->get_magicDelta(); }
    const fvMesh& get_submesh() const { return HyperreductionConfiguration_->get_submesh(); }
    const volScalarField& get_volume() const { return meshConfig_->get_volume(); }
    const double& get_totalVolume() const { return meshConfig_->get_totalVolume(); }
    const volVectorField& get_meanU() const { return *ROMExecutionConfig_->meanU(); }
    const volVectorField& get_meanVectorDEIM() const { return *HyperreductionConfiguration_->meanVectorDEIM(); }
    const volScalarField& get_meanScalarDEIM() const { return *HyperreductionConfiguration_->meanScalarDEIM(); }
    void get_meanDEIM(volVectorField& meanVectorDEIM_out) { meanVectorDEIM_out = get_meanVectorDEIM(); }
    void get_meanDEIM(volScalarField& meanVectorDEIM_out) { meanVectorDEIM_out = get_meanScalarDEIM(); }
    const volVectorField& get_meanVectorDEIMMagic() const { return *HyperreductionConfiguration_->meanVectorDEIMMagic(); }
    const volScalarField& get_meanScalarDEIMMagic() const { return *HyperreductionConfiguration_->meanScalarDEIMMagic(); }

    const Eigen::Tensor<double, 3>& get_xi_onMagicPts() const { return HyperreductionConfiguration_->xiOnMagicPts(); }
    const Eigen::Tensor<double, 3>& get_M_xi_onMagicPts() const { return HyperreductionConfiguration_->M_xiOnMagicPts(); }

    /// Getter and setter of DEIM matrix K
    const Eigen::MatrixXd& get_K_DEIM() const { return HyperreductionConfiguration_->K_DEIM(); }
    void set_K_DEIM(Eigen::MatrixXd& K) { HyperreductionConfiguration_->setK_DEIM(K); }

    /// Setters to the members of m_parameters that cannot initialize at the beginning of the program
    /// For instance, eigen values can be known only after POD was performed
    void set_nSnapshots(const label& input_nSnapshots) { snapshotConfiguration_->set_nSnapshots(input_nSnapshots); }
    void set_nSnapshotsSimulation(const label& input_nSnapshotsSimulation) { snapshotConfiguration_->set_nSnapshotsSimulation(input_nSnapshotsSimulation); }
    void set_eigenValues_p(const Eigen::VectorXd& input_eigenValues_p) { PODConfiguration_->set_eigenValues_p(input_eigenValues_p); }
    void set_saveTime(const double& input_saveTime) { snapshotConfiguration_->set_writeInterval(input_saveTime); }
    void set_xi_onMagicPts(const Eigen::Tensor<double, 3> xi) { HyperreductionConfiguration_->setXiOnMagicPts(xi); }
    void set_M_xi_onMagicPts(const Eigen::Tensor<double, 3> M_xi) { HyperreductionConfiguration_->setM_XiOnMagicPts(M_xi); }
    void set_useSOTA(word& input_useSOTA) { ROMExecutionConfig_->setUseSOTA(input_useSOTA); }

    /// set DDES
    void set_useDDES(const bool& a) { SimulationFlags_->set_useDDES(a); }

    /// set DNS
    void set_useDNS(const bool& c) { SimulationFlags_->set_useDNS(c); }

    /// Setters of POD and DEIM parameters
    void set_DEIMInterpolatedField(Foam::word name) { HyperreductionConfiguration_->setHRInterpolatedField(name); }
    void set_nMagicPoints(const label& n) { HyperreductionConfiguration_->setNMagicPoints(n); }
    void set_Ck(const float& c) { HyperreductionConfiguration_->setCk(c); }
    void set_Ce(const float& c) { HyperreductionConfiguration_->setCe(c); }
    void set_useDEIM(const bool& b) { SimulationFlags_->set_useDEIM(b); }

    void set_magicPoints(const List<label>& nn) { HyperreductionConfiguration_->setMagicPoints(nn); }
    void set_localMagicPoints(const List<label>& nn) { HyperreductionConfiguration_->setLocalMagicPoints(nn); }
    void set_projected_K_DEIM(const Eigen::MatrixXd& K_input) { HyperreductionConfiguration_->setProjectedK_DEIM(K_input); }
    void set_deformationTensorOfModesOnMagicNeighborhoods(PtrList<volVectorField>& defT)
    {
        HyperreductionConfiguration_->setTracerGradOnMagicNeighborhoods(defT, get_nModes()["U"] + 1);
    }

    void set_deformationTensorOfModesOnMagicPoints(PtrList<volVectorField>& defT)
    {
        HyperreductionConfiguration_->setTracerGradOnMagicPoints(defT, get_nModes()["U"] + 1);
    }

    void set_deformationTensorOfModesOnMagicNeighborhoods(PtrList<volTensorField>& defT)
    {
        HyperreductionConfiguration_->setDeformationTensorOnMagicNeighborhoods(defT, get_nModes()["U"] + 1);
    }

    void set_deformationTensorOfModesOnMagicPoints(PtrList<volTensorField>& defT)
    {
        HyperreductionConfiguration_->setDeformationTensorOnMagicPoints(defT, get_nModes()["U"] + 1);
    }

    void set_projected_MK_DEIM(const Eigen::MatrixXd& proj_MK_input) { HyperreductionConfiguration_->setProjectedMK_DEIM(proj_MK_input); }
    void set_magicDelta(const volScalarField& mD) { HyperreductionConfiguration_->set_magicDelta(mD); }
    void set_submesh(fvMesh& s) { HyperreductionConfiguration_->set_submesh(s); }


    void set_meanU(const volVectorField& mm) { ROMExecutionConfig_->setMeanU(mm); }
    void set_meanDEIM(const volVectorField& mm) { HyperreductionConfiguration_->setMeanVectorDEIM(mm); }
    void set_meanDEIM(const volScalarField& mm) { HyperreductionConfiguration_->setMeanScalarDEIM(mm); }
    void set_meanDEIMMagic(const volVectorField& mm)
    {
        HyperreductionConfiguration_->setMeanVectorDEIMMagic(mm);
    }
    void set_meanDEIMMagic(const volScalarField& mm)
    {
        HyperreductionConfiguration_->setMeanScalarDEIMMagic(mm);
    }
    void set_deltaWeight(Eigen::VectorXd dw) { HyperreductionConfiguration_->setDeltaWeight(dw); }

private:
    // configuration components
    std::unique_ptr<MeshConfiguration> meshConfig_;
    std::unique_ptr<FieldTemplates> fieldTemplates_;
    std::unique_ptr<SnapshotConfiguration> snapshotConfiguration_;
    std::unique_ptr<SolverConfiguration> solverConfiguration_;
    std::unique_ptr<HyperreductionConfiguration> HyperreductionConfiguration_;
    std::unique_ptr<ROMExecutionConfig> ROMExecutionConfig_;
    std::unique_ptr<SimulationFlags> SimulationFlags_;
    std::unique_ptr<PODConfiguration> PODConfiguration_;

    void initializeFieldConfiguration();

protected:
    // OPENFOAM runtime
    autoPtr<argList> _args;
    autoPtr<Time> runTime0;
    Foam::Time* runTimeData;


    // ITHACA integration
    ITHACAparameters* ithacaLibraryParameters;
    IOdictionary* ITHACAdict;
};


#endif
