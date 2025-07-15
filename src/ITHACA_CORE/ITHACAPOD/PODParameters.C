
#include "ITHACAsystem.H"
#include "PODParameters.H"

namespace ITHACAPOD
{
Foam::word PODParameters::get_pathHilbertSpace_fromHS(Foam::word hilbertSp)
{
    Foam::word pathHilbertSpace = "";

    if (hilbertSp == "L2" || hilbertSp == "dL2")
    {
        pathHilbertSpace = "";
    }
    else if (hilbertSp == "L2wBC")
    {
        pathHilbertSpace = "_L2wBC";
    }
    else if (hilbertSp == "H1")
    {
        pathHilbertSpace = "_H1";
    }
    else if (hilbertSp == "wH1")
    {
        pathHilbertSpace = "_wH1";
    }
    else
    {
        Foam::Info << "Error: hilbertSpacePOD type " << hilbertSp
                   << " is not valid." << Foam::endl;
        Foam::Info << "dot_product_POD is available for L2, L2wBC, H1 and wH1 only." <<
                   Foam::endl;
        abort();
    }

    return pathHilbertSpace;
}


//specialisation:
template<typename T>
void PODParameters::read_snapshot(T& snapshot, const Foam::label& i_snap,
                                  Foam::word path, Foam::word name) const
{
    if (name == "default_name")
    {
        name = snapshot.name();
    }

    if (path == "default_path")
    {
        path = runTimeData->path() + "/" + runTimeData->times()[i_snap].name();
    }

    if (!ITHACAutilities::check_file(path))
    {
        Info << "Error: data not found at :" << endl;
        Info << path << endl;
        Info << endl;
        abort();
    }

    T snapshot_dummy(
        IOobject
        (
            name,
            path,
            *mesh,
            IOobject::MUST_READ
        ),
        *mesh
    );
    snapshot = snapshot_dummy;
}
template void PODParameters::read_snapshot(Foam::volScalarField& snapshot,
        const Foam::label& i_snap, Foam::word path, Foam::word name) const;
template void PODParameters::read_snapshot(Foam::volVectorField& snapshot,
        const Foam::label& i_snap, Foam::word path, Foam::word name) const;
template void PODParameters::read_snapshot(Foam::volTensorField& snapshot,
        const Foam::label& i_snap, Foam::word path, Foam::word name) const;
}
