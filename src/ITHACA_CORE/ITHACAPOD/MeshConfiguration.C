#include "MeshConfiguration.H"
#include "GeometricFields.H"
#include "fvPatchFields.H"
#include "volMesh.H"

namespace ITHACAPOD
{
    void MeshConfiguration::set_magicDelta(const Foam::volScalarField& mD) 
    {
        m_magicDelta = new Foam::GeometricField<Foam::scalar, Foam::fvPatchField, Foam::volMesh>(mD); 
    }
}
