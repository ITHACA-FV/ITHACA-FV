
#ifndef POD_PARAMETERS_H
#define POD_PARAMETERS_H

#include <FieldField.H>

#include "fvMesh.H"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#include <Eigen/Eigen>
#pragma GCC diagnostic pop

#include <unsupported/Eigen/CXX11/Tensor>

namespace ITHACAPOD
{

class Parameters
{
    public:
        virtual ~Parameters() { }
};
// only one derivate
/**
 *  @brief  Class that containing all parameters of the stochastic POD
 ***********************************************/
class PODParameters : public Parameters
{

    // ITHACA-FV 
    enum class PressureResolutionKind{
      FullOrder,
      ReducedOrder,
      Neglected,
      Undefined
    };

  private:
    autoPtr<argList> _args;
    // Dossier dans lequel on a le 0 cst system
    // utile pour ITHACA-FV ? un runtime est déjà utilisé pour lire ITHACAdict
    // pour l’instant on le garde ici
    // objet openfoam
    autoPtr<Time> runTime0;

    // ITHACA-FV, pour la DEIM
    // il est probablement déjà instancié quelque part dans ITHACA-FV
    fvMesh* submesh;//! ref
 
    /// Parameters from ITHACAparameters class from ITHACA-FV library

    ITHACAparameters* ithacaLibraryParameters;//! ref
    IOdictionary* ITHACAdict;
    wordList fieldlist;

    // Provient déjà d’ITHACA-FV
    bool exportPython;
    bool exportMatlab;
    bool exportTxt;

    /// boolean indicating if DEIM is used or not
    // ITHACA-FV
    bool useDEIM;

    /// boolean indicating use DDES or not
    // ça reste là
    bool useDDES;

    /// boolean indicating use of DNS or not
    // ITHACA-FV
    bool useDNS;


    /// boolean indicating on line reconstruction of velocity
    // ITHACA-FV
    bool onLineReconstruct;

    /// string indicating what temporal scheme to used
    // Les deux (ITHACA-FV & RedLUM)
    word ROMTemporalScheme;

    // ITHACA-FV
    bool interpFieldCenteredOrNot;

    // Forçage tour antigel, ça reste là
    bool forcingOrNot;

    // ITHACA-FV
    bool symDiff;

    // rapport entre le pas de temps de la simu et celui utilisé pour la résolution du ROM 
    label nSimu;


    // ITHACA-FV
    double FinalTime;
    double FinalTimeSimulation;

    // ITHACA-FV
    dimensionedScalar* nu;

    // ITHACA-FV
    PressureResolutionKind pressureResolutionKind;
    
    // ITHACA-FV
    /// Turbulent viscosity coefficient for Smagorinsky model
    /// Standard values used in OpenFOAM
    /// see : A. Yoshizawa, “Statistical theory for compressible turbulent shear flows, withthe application to subgrid modeling,” Phys. Fluids29, 2152 (1986).
    float Ck;
    float Ce;

    // ITHACA-FV
    /// Number of magic points for build_DEIM_U
    label nMagicPoints;

    // ITHACA-FV
    /// DEIM parameters
    PtrList<volTensorField> deformationTensorOfModesOnMagicPoints;
    PtrList<volTensorField> deformationTensorOfModesOnMagicNeighborhoods;
    PtrList<volVectorField> tracerGradOfModesOnMagicPoints;
    PtrList<volVectorField> tracerGradOfModesOnMagicNeighborhoods;
    List<label> magicPoints;
    List<label> localMagicPoints;

    // ITHACA-FV
    /// DEIM matrix (high dimensioned)
    Eigen::MatrixXd K_DEIM_matrix;
    Eigen::MatrixXd projected_K_DEIM_matrix;
    Eigen::MatrixXd projected_MK_DEIM_matrix;

    // ITHACA-FV
    /// Geometric Information
    volScalarField* volume;
    double totalVolume;
    volScalarField* delta;
    volScalarField* magicDelta;

    // ITHACA-FV
    volVectorField* meanU;
    volVectorField* meanVectorDEIM;
    volScalarField* meanScalarDEIM;
    volVectorField* meanVectorDEIMMagic;
    volScalarField* meanScalarDEIMMagic;

    // ITHACA-FV
    word folder_DEIM;

    /// for DEIM on nut
    Eigen::Tensor<double, 3> xi_onMagicPts;
    Eigen::Tensor<double, 3 > M_xi_onMagicPts;

    protected:
        // ITHACA-FV
        Foam::Time* runTimeData;
        Foam::fvMesh* mesh;
        Foam::label nCells;

        // ITHACA-FV
        Foam::word eigensolver;
        int precision;
        std::_Ios_Fmtflags outytpe;

        // ITHACA-FV
        Foam::label nBlocks;
        bool centeredOrNot;

        // Nécessaire pour la POD
        // ITHACA-FV
        Foam::volVectorField* template_field_U;
        Foam::volScalarField* template_field_p;
        Foam::volScalarField* template_field_nut;
        Foam::volScalarField* template_field_omega;
        Foam::volScalarField* template_field_k;

        // ITHACA-FV
        Foam::HashTable<Foam::label, Foam::word> nModes;
        Foam::HashTable<Foam::word, Foam::word> hilbertSpacePOD;
        Foam::HashTable<double, Foam::word> resolvedVaryingEnergy;
        Foam::HashTable<double, Foam::word> varyingEnergy;
        Foam::HashTable<double, Foam::word> meanEnergy;

        // ITHACA-FV
        // Sert pour les produits scalaires dans ITHACAFVTool
        double weightH1;
        double weightBC;
        Foam::word patchBC;
        Foam::List<Foam::word> field_name;
        List<word> field_type;

        // ITHACA-FV
        Foam::fileName casenameData;

        // ITHACA-FV
        Foam::label startTime;
        Foam::label endTime;
        Foam::label nSnapshots;

        // ITHACA-FV
        Foam::label endTimeSimulation;
        Foam::label nSnapshotsSimulation;

        // ITHACA-FV
        // À renommer en "writeInterval" 
        double saveTime;

        // ITHACA-FV
        Eigen::VectorXd eigenValues_U;
        Eigen::VectorXd eigenValues_p;
        Eigen::VectorXd lambda;

        // ITHACA-FV
        Foam::word DEIMInterpolatedField;
        Eigen::VectorXd* deltaWeight;

        // ITHACA-FV
        word useSOTA;

    public:
        PODParameters(int argc,char* argv[]);

        /// Getter to access all IthacaFVParameters members in the other classes
        ITHACAparameters* get_ITHACAparameters() const
        { 
            return ithacaLibraryParameters;
        }

        const Foam::HashTable<Foam::label, Foam::word>& get_nModes() const
        {
            return nModes;
        }

        const Foam::HashTable<Foam::word, Foam::word>& get_hilbertSpacePOD() const
        {
            return hilbertSpacePOD;    // Contains both the field name and Hilbert space
        }

        Foam::word get_pathHilbertSpace_fromHS(Foam::word hilbertSp);
        const double& get_weightH1() const
        {
            return weightH1;
        }
        const double& get_weightPOD() const
        {
            return weightBC;
        }
        const Foam::word& get_patchBC() const
        {
            return patchBC;
        }
        const Foam::fileName& get_casenameData() const
        {
            return casenameData;
        }

        const Foam::label& get_startTime() const
        {
            return startTime;
        }
        const Foam::label& get_endTime() const
        {
            return endTime;
        }
        const Foam::label& get_nSnapshots() const
        {
            return nSnapshots;
        }

        const Foam::label& get_endTimeSimulation() const
        {
            return endTimeSimulation;
        }
        const Foam::label& get_nSnapshotsSimulation() const
        {
            return nSnapshotsSimulation;
        }

        const Foam::label& get_nBlocks() const
        {
            return nBlocks;
        }

        const bool& get_centeredOrNot() const
        {
            return centeredOrNot;
        }
        const int& get_precision() const
        {
            return precision;
        }
        const std::_Ios_Fmtflags& get_outytpe() const
        {
            return outytpe;
        }

        const Foam::word& get_eigensolver() const
        {
            return eigensolver;
        }

        const Foam::fvMesh& get_mesh() const
        {
            return *mesh;
        }
        const Foam::label& get_nCells() const
        {
            return nCells;
        }

        const Foam::volVectorField& get_template_field_U() const
        {
            return *template_field_U;
        }
        const Foam::volScalarField& get_template_field_p() const
        {
            return *template_field_p;
        }
        const Foam::volScalarField& get_template_field_nut() const
        {
            return *template_field_nut;
        }
        const Foam::volScalarField& get_template_field_omega() const
        {
            return *template_field_omega;
        }
        const Foam::volScalarField& get_template_field_k() const
        {
            return *template_field_k;
        }

        const double& get_saveTime() const
        {
            return saveTime;
        }

        const Foam::word& get_DEIMInterpolatedField() const
        {
            return DEIMInterpolatedField;
        }
        Eigen::VectorXd get_deltaWeight() const
        {
            return *deltaWeight;
        }

        void set_centeredOrNot(const bool b)
        {
            centeredOrNot = b;
        }

        /// Setters to the members of PODParameters that cannot initialize at the beginning of the program
        void set_eigenValues_U(const Eigen::VectorXd& input_eigenValues_U)
        {
            eigenValues_U = input_eigenValues_U;
        }
        void set_lambda(const Eigen::VectorXd& input_lambda)
        {
            lambda = input_lambda;
        }
        /// Setters of POD and DEIM parameters
        void set_nModes(const Foam::word& field_name, const Foam::label& n)
        {
            nModes.insert(field_name, n);
        }
        void set_resolvedVaryingEnergy(const Foam::word& field_name, const double& n)
        {
            resolvedVaryingEnergy.set(field_name, n);
        }
        void set_varyingEnergy(const Foam::word& field_name, const double& n)
        {
            varyingEnergy.set(field_name, n);
        }
        void set_meanEnergy(const Foam::word& field_name, const double& n)
        {
            meanEnergy.insert(field_name, n);
        }
        void set_hilbertSpacePOD(const Foam::word& field_name,
                                 const Foam::word& hilbertSp)
        {
            hilbertSpacePOD.insert(field_name, hilbertSp);
        }
        void set_weightH1(const double& c)
        {
            weightH1 = c;
        }


           const wordList& get_fieldlist() const { return fieldlist; }
    const HashTable<double,word>& get_resolvedVaryingEnergy() const { return resolvedVaryingEnergy; }
    const HashTable<double,word>& get_varyingEnergy() const { return varyingEnergy; }
    const HashTable<double,word>& get_meanEnergy() const { return meanEnergy; }
    Foam::word get_pathHilbertSpace(Foam::word fieldName);

    const List<word>& get_field_name() const { return field_name; }
    const List<word>& get_field_type() const { return field_type; }

    // getting current F_0 force vector and F_mask
    const Foam::Vector<double>& get_F_0() const;
    void compute_F_mask();

    const PressureResolutionKind& get_pressureResolutionKind(){return pressureResolutionKind;}

    const word& get_ROMTemporalScheme() const { return ROMTemporalScheme; }

    const double& get_initialTimeSimulation() const {
      return FinalTime; }
    const double& get_finalTime() const {
      return FinalTime; }
    const double& get_finalTimeSimulation() const {
      return FinalTimeSimulation; }



    const bool& get_interpFieldCenteredOrNot() const { return interpFieldCenteredOrNot; }

    const bool& get_exportPython() const { return exportPython; }
    const bool& get_exportMatlab() const { return exportMatlab; }
    const bool& get_exportTxt() const { return exportTxt; }


    const word& get_useSOTA() const {return useSOTA;}


    const dimensionedScalar& get_nu() const { return *nu; }


    const Eigen::VectorXd& get_eigenValues_U() const { return eigenValues_U; }
    const Eigen::VectorXd& get_eigenValues_p() const { return eigenValues_p; }
    const Eigen::VectorXd& get_lambda() const { return lambda; }

    const label& get_nSimu() const { return nSimu; }
    const bool& get_forcingOrNot() const {return forcingOrNot;}
    const bool& get_symDiff() const {return symDiff;}

    /// on line velocity Reconstruction getter
    const bool& get_onLineReconstruct() const {return onLineReconstruct;}

    /// DDES getter
    const bool& get_useDDES() const {return useDDES;}

    /// DNS getter
    const bool& get_useDNS() const {return useDNS;}

    /// Getters of DEIM parameters
    const bool& get_useDEIM() const {return useDEIM;}
    const float& get_Ck() const {return Ck;}
    const float& get_Ce() const {return Ce;}
    const label& get_nMagicPoints() const {return nMagicPoints;}
    const word& get_folder_DEIM() const {return folder_DEIM;}

    const PtrList<volVectorField>& get_tracerGradOfModesOnMagicNeighborhoods() const {return tracerGradOfModesOnMagicNeighborhoods;}
    const PtrList<volVectorField>& get_tracerGradOfModesOnMagicPoints() const {return tracerGradOfModesOnMagicPoints;}
    const PtrList<volTensorField>& get_deformationTensorOfModesOnMagicNeighborhoods() const {return deformationTensorOfModesOnMagicNeighborhoods;}
    const PtrList<volTensorField>& get_deformationTensorOfModesOnMagicPoints() const {return deformationTensorOfModesOnMagicPoints;}
    const List<label>& get_magicPoints() const {return magicPoints;}
    const List<label>& get_localMagicPoints() const {return localMagicPoints;}
    const Eigen::MatrixXd& get_projected_K_DEIM() const {return projected_K_DEIM_matrix;}
    const Eigen::MatrixXd& get_projected_MK_DEIM() const {return projected_MK_DEIM_matrix;}
    const volScalarField& get_delta() const {return *delta;}
    const volScalarField& get_magicDelta() const {return *magicDelta;}
    const fvMesh& get_submesh() const {return *submesh;}
    const volScalarField& get_volume() const {return *volume;}
    const double& get_totalVolume() const {return totalVolume;}
    const volVectorField& get_meanU() const {return *meanU;}
    const volVectorField& get_meanVectorDEIM() const { return *meanVectorDEIM; }
    const volScalarField& get_meanScalarDEIM() const { return *meanScalarDEIM; }
    void get_meanDEIM(volVectorField& meanVectorDEIM_in) {meanVectorDEIM_in = *meanVectorDEIM;}
    void get_meanDEIM(volScalarField& meanScalarDEIM_in) {meanScalarDEIM_in = *meanScalarDEIM;}
    const volVectorField& get_meanVectorDEIMMagic() const { return *meanVectorDEIMMagic; }
    const volScalarField& get_meanScalarDEIMMagic() const { return *meanScalarDEIMMagic; }

    const Eigen::Tensor<double, 3>& get_xi_onMagicPts() const {return xi_onMagicPts;}
    const Eigen::Tensor<double, 3>& get_M_xi_onMagicPts() const {return M_xi_onMagicPts;}


    /// Getter and setter of DEIM matrix K
    const Eigen::MatrixXd& get_K_DEIM() const {return K_DEIM_matrix;}
    void set_K_DEIM(Eigen::MatrixXd& K_input) {K_DEIM_matrix=K_input;}


    /// Setters to the members of IthacaFVParameters that cannot initialize at the beginning of the program
    /// For instance, eigen values can be known only after POD was performed
    void set_nSnapshots(const label& input_nSnapshots) { nSnapshots = input_nSnapshots; }
    void set_nSnapshotsSimulation(const label& input_nSnapshotsSimulation) { nSnapshotsSimulation = input_nSnapshotsSimulation; }
    void set_eigenValues_p(const Eigen::VectorXd& input_eigenValues_p) { eigenValues_p = input_eigenValues_p; }
    void set_saveTime(const double& input_saveTime) { saveTime = input_saveTime; }
    void set_xi_onMagicPts(const Eigen::Tensor<double, 3> xi_input) {xi_onMagicPts=xi_input;}
    void set_M_xi_onMagicPts(const Eigen::Tensor<double, 3 > M_xi_input) {M_xi_onMagicPts=M_xi_input;}
    void set_useSOTA(word& input_useSOTA) {useSOTA=input_useSOTA;}

    /// set DDES
    void set_useDDES(const bool& a) {useDDES=a;}

    /// set DNS
    void set_useDNS(const bool & c) {useDNS= c;}

    /// Setters of POD and DEIM parameters
    void set_DEIMInterpolatedField(Foam::word name) {HRInterpolatedField=name;}
    void set_nMagicPoints(const label& n) {nMagicPoints=n;}
    void set_Ck(const float& c) {Ck=c;}
    void set_Ce(const float& c) {Ce=c;}
    void set_useDEIM(const bool& b) {useDEIM=b;}

    void set_deformationTensorOfModesOnMagicNeighborhoods(PtrList<volVectorField>& defT) {
      tracerGradOfModesOnMagicNeighborhoods.resize( nModes["U"]+1 );
      for (label k = 0 ; k < (nModes["U"]+1) ; k++)
      {
        tracerGradOfModesOnMagicNeighborhoods.set(k, new volVectorField(defT[k]));
      }
    }

    void set_deformationTensorOfModesOnMagicPoints(PtrList<volVectorField>& defT) {
      tracerGradOfModesOnMagicPoints.resize( nModes["U"]+1 );
      for (label k = 0 ; k < (nModes["U"]+1) ; k++)
      {
        tracerGradOfModesOnMagicPoints.set(k, new volVectorField(defT[k]));
      }
    }

    void set_deformationTensorOfModesOnMagicNeighborhoods(PtrList<volTensorField>& defT) {
      deformationTensorOfModesOnMagicNeighborhoods.resize( nModes["U"]+1 );
      for (label k = 0 ; k < (nModes["U"]+1) ; k++)
      {
        deformationTensorOfModesOnMagicNeighborhoods.set(k, new volTensorField(defT[k]));
      }
    }

    void set_deformationTensorOfModesOnMagicPoints(PtrList<volTensorField>& defT) {
      deformationTensorOfModesOnMagicPoints.resize( nModes["U"]+1 );
      for (label k = 0 ; k < (nModes["U"]+1) ; k++)
      {
        deformationTensorOfModesOnMagicPoints.set(k, new volTensorField(defT[k]));
      }
    }

    void set_magicPoints(const List<label>& nn) {magicPoints = nn;}
    void set_localMagicPoints(const List<label>& nn) {localMagicPoints = nn;}
    void set_projected_K_DEIM(const Eigen::MatrixXd& K_input) {projected_K_DEIM_matrix=K_input;}
    void set_projected_MK_DEIM(const Eigen::MatrixXd& proj_MK_input) {projected_MK_DEIM_matrix= proj_MK_input;}
    void set_magicDelta(const volScalarField& mD) {magicDelta = new volScalarField(mD);}
    void set_submesh(fvMesh& s) {submesh = &s;}
    void set_meanU(const volVectorField& mm) {meanU = new volVectorField(mm);}
    void set_meanDEIM(const volVectorField& mm) {meanVectorDEIM = new volVectorField(mm);}
    void set_meanDEIM(const volScalarField& mm) {meanScalarDEIM = new volScalarField(mm);}
    void set_meanDEIMMagic(const volVectorField& mm)
    {meanVectorDEIMMagic = new volVectorField(mm);}
    void set_meanDEIMMagic(const volScalarField& mm)
    {meanScalarDEIMMagic = new volScalarField(mm);}

    void set_deltaWeight(Eigen::VectorXd dw){deltaWeight = new Eigen::VectorXd(dw);}

};




}
#endif
