
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

protected:
    Foam::Time* runTimeData;
    Foam::fvMesh* mesh;
    Foam::label nCells;

    Foam::word eigensolver;
    int precision;
    std::_Ios_Fmtflags outytpe;


    /// boolean indicating if the large-scale advection is corrected (LU style)
    Foam::label nBlocks;
    bool centeredOrNot;

    Foam::volVectorField* template_field_U;
    Foam::volScalarField* template_field_p;
    Foam::volScalarField* template_field_nut;
    Foam::volScalarField* template_field_omega;
    Foam::volScalarField* template_field_k;

    Foam::HashTable<Foam::label, Foam::word> nModes;
    Foam::HashTable<Foam::word, Foam::word> hilbertSpacePOD;
    Foam::HashTable<double, Foam::word> resolvedVaryingEnergy;
    Foam::HashTable<double, Foam::word> varyingEnergy;
    Foam::HashTable<double, Foam::word> meanEnergy;
    double weightH1;
    double weightBC;
    Foam::word patchBC;
    Foam::List<Foam::word> field_name;

    Foam::fileName casenameData;

    Foam::label startTime;
    Foam::label endTime;
    Foam::label nSnapshots;

    Foam::label endTimeSimulation;
    Foam::label nSnapshotsSimulation;

    double saveTime;

    Eigen::VectorXd eigenValues_U;
    Eigen::VectorXd lambda;

    Foam::word DEIMInterpolatedField;
    Eigen::VectorXd* deltaWeight;

public:
    template <typename T>
    void read_snapshot(T& snapshot, const Foam::label& i_snap, Foam::word path = "default_path", Foam::word name = "default_name") const;

    const Foam::HashTable<Foam::label, Foam::word>& get_nModes() const { return nModes; }

    const Foam::HashTable<Foam::word, Foam::word>& get_hilbertSpacePOD() const { return hilbertSpacePOD; } // Contains both the field name and Hilbert space

    Foam::word get_pathHilbertSpace_fromHS(Foam::word hilbertSp);
    const double& get_weightH1() const { return weightH1; }
    const double& get_weightPOD() const { return weightBC; }
    const Foam::word& get_patchBC() const { return patchBC; }
    const Foam::fileName& get_casenameData() const { return casenameData; }

    const Foam::label& get_startTime() const { return startTime; }
    const Foam::label& get_endTime() const { return endTime; }
    const Foam::label& get_nSnapshots() const { return nSnapshots; }

    const Foam::label& get_endTimeSimulation() const { return endTimeSimulation; }
    const Foam::label& get_nSnapshotsSimulation() const { return nSnapshotsSimulation; }

    const Foam::label& get_nBlocks() const { return nBlocks; }

    const bool& get_centeredOrNot() const { return centeredOrNot; }
    const int& get_precision() const { return precision; }
    const std::_Ios_Fmtflags& get_outytpe() const { return outytpe; }

    const Foam::word& get_eigensolver() const { return eigensolver; }

    const Foam::fvMesh& get_mesh() const { return *mesh; }
    const Foam::label& get_nCells() const { return nCells; }

    const Foam::volVectorField& get_template_field_U() const { return *template_field_U; }
    const Foam::volScalarField& get_template_field_p() const { return *template_field_p; }
    const Foam::volScalarField& get_template_field_nut() const { return *template_field_nut; }
    const Foam::volScalarField& get_template_field_omega() const { return *template_field_omega; }
    const Foam::volScalarField& get_template_field_k() const { return *template_field_k; }

    const double& get_saveTime() const { return saveTime; }

    const Foam::word& get_DEIMInterpolatedField() const { return DEIMInterpolatedField; }
    Eigen::VectorXd get_deltaWeight() const { return *deltaWeight; }

    /// Setters to the members of PODParameters that cannot initialize at the beginning of the program
    void set_eigenValues_U(const Eigen::VectorXd& input_eigenValues_U) { eigenValues_U = input_eigenValues_U; }
    void set_lambda(const Eigen::VectorXd& input_lambda) { lambda = input_lambda; }
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
    void set_hilbertSpacePOD(const Foam::word& field_name, const Foam::word& hilbertSp)
    {
        hilbertSpacePOD.insert(field_name, hilbertSp);
    }
    void set_weightH1(const double& c) { weightH1 = c; }
};




}
#endif
