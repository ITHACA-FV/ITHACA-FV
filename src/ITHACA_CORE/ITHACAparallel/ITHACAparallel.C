#include "ITHACAparallel.H"

List<label> ITHACAparallel::oldProcIDs_(0);
List<label> ITHACAparallel::newProcIDs_(0);

ITHACAparallel::ITHACAparallel(fvMesh& mesh)
{
    N_BF = 0;

    for (int i = 0; i < mesh.boundaryMesh().size(); i++)
    {
        if (mesh.boundaryMesh()[i].type() != "processor")
        {
            N_BF++;
        }
    }

    // Load cell addressing
    indices = new labelIOList(
        IOobject
        (
            "cellProcAddressing",
            mesh.facesInstance(),
            mesh.meshSubDir,
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    // Load face addressing
    indicesF = new labelIOList(
        IOobject
        (
            "faceProcAddressing",
            mesh.facesInstance(),
            mesh.meshSubDir,
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    // Calculate total number of cells
    N_IF_glob = (mesh.C().size());
    reduce(N_IF_glob, sumOp<int>());
    // BF construction
    Gsize_BF = new labelList (N_BF, 0);
    IndFaceLocal = new List<labelList> (N_BF, labelList(0, 0));

    for (int i = 0; i < N_BF; i++)
    {
        Gsize_BF()[i] = mesh.boundaryMesh()[i].size();
        reduce(Gsize_BF()[i], sumOp<label>());
        IndFaceLocal()[i].resize(mesh.boundaryMesh()[i].size());

        for (int k = 0; k < mesh.boundaryMesh()[i].size(); k++)
        {
            IndFaceLocal()[i][k] = indicesF()[mesh.boundaryMesh()[i].start() + k];
        }
    }

    Start = new labelList(N_BF, 0);

    for (int i = 0; i < N_BF; i++)
    {
        if (IndFaceLocal()[i].size() == 0)
        {
            Start()[i] = INT_MAX;
        }
        else if (IndFaceLocal()[i][0] < 0)
        {
            Start()[i] = INT_MAX;
        }

        if (IndFaceLocal()[i].size() != 0)
        {
            Start()[i] = IndFaceLocal()[i][0];
        }

        reduce(Start()[i], minOp<label>());
    }

    Info << Start() << endl;
}

void ITHACAparallel::suspendMPI()
{
    Pstream::parRun() = false;
    label comm        = Pstream::worldComm;
    oldProcIDs_       = Pstream::procID(comm);
    newProcIDs_       = List<label> (1, 0);
    Pstream::procID(comm) = newProcIDs_;
}

void ITHACAparallel::resumeMPI()
{
    label comm        = Pstream::worldComm;
    Pstream::procID(comm) = oldProcIDs_;
    Pstream::parRun() = true;
}