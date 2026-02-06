#include "FieldTemplates.H"

#include "GeometricFields.H"
#include "fvPatchFields.H"
#include "volMesh.H"

void FieldTemplates::set_fullStressFunction(Foam::volVectorField& templateSmag)
{ 
        m_fullStressFunction = new Foam::volVectorField(templateSmag);
}
