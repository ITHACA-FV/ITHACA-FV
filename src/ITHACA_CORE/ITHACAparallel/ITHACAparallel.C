#include "ITHACAparallel.H"

List<int> ITHACAparallel::oldProcIDs_(0);
List<int> ITHACAparallel::newProcIDs_(0);
ITHACAparallel* ITHACAparallel::instance = nullptr;

ITHACAparallel* ITHACAparallel::getInstance(fvMesh& mesh, Time& localTime)

{
    M_Assert(instance == nullptr,
             "ITHACAparallel is already initialized, call ITHACAparallel::getInstance() to return an instance of ITHACAparallel");
    instance = new ITHACAparallel(mesh, localTime);
    return instance;
}

ITHACAparallel* ITHACAparallel::getInstance()

{
    M_Assert(instance != nullptr,
             "ITHACAparallel needs to be initialized, call ITHACAparallel::getInstance(mesh, runTime) first");
    return instance;
}


ITHACAparallel::ITHACAparallel(fvMesh& mesh, Time& localTime)
    :
    runTime(localTime),
    mesh(mesh)
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
    indices = autoPtr<labelIOList>
              (
                  new labelIOList
                  (
                      IOobject
                      (
                          "cellProcAddressing",
                          mesh.facesInstance(),
                          mesh.meshSubDir,
                          mesh,
                          IOobject::MUST_READ,
                          IOobject::NO_WRITE
                      )
                  )
              );
    // Load face addressing
    indicesF = autoPtr<labelIOList>
               (
                   new labelIOList
                   (
                       IOobject
                       (
                           "faceProcAddressing",
                           mesh.facesInstance(),
                           mesh.meshSubDir,
                           mesh,
                           IOobject::MUST_READ,
                           IOobject::NO_WRITE
                       )
                   )
               );
    // Calculate total number of cells
    N_IF_glob = (mesh.C().size());
    reduce(N_IF_glob, sumOp<int>());
    // BF construction
    Gsize_BF = autoPtr<labelList>(new labelList (N_BF, label(0)));
    IndFaceLocal = autoPtr< List<labelList>> (new List<labelList> (N_BF,
                   labelList(label(0), label(0))));

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

    Start = autoPtr<labelList> (new labelList(N_BF, 0));

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
}

void ITHACAparallel::suspendMPI()
{
    Pstream::parRun() = false;
    label comm        = Pstream::worldComm;
    oldProcIDs_       = Pstream::procID(comm);
    newProcIDs_       = List<int> (1);
    Pstream::procID(comm) = newProcIDs_;
}

void ITHACAparallel::resumeMPI()
{
    label comm        = Pstream::worldComm;
    Pstream::procID(comm) = oldProcIDs_;
    Pstream::parRun() = true;
}

template<>
List<Field <scalar>> ITHACAparallel::combineFields(
                      GeometricField<scalar, fvPatchField, volMesh>& field)
{
    List<Field< scalar>> GlobField(field.boundaryFieldRef().size() + 1);
    GlobField[0].resize(N_IF_glob);
    GlobField[0] = GlobField[0] * 0;

    // Assemble internalField
    for (int i = 0; i < field.size(); i++)
    {
        GlobField[0][indices()[i]] = field[i];
    }

    reduce(GlobField[0], sumOp<Field<scalar>>());

    // Assemble BoundariField
    for (int i = 0; i < N_BF; i++)
    {
        GlobField[i + 1].resize(Gsize_BF()[i]);
        Field<scalar> zero(Gsize_BF()[i], 0.0);
        GlobField[i + 1] = zero;
    }

    for (int i = 0; i < N_BF; i++)
    {
        for (int k = 0; k < field.boundaryFieldRef()[i].size(); k++)
        {
            if (IndFaceLocal()[i].size() > 0
                    && field.boundaryFieldRef()[i].type() != "zeroGradient"
                    && field.boundaryFieldRef()[i].type() != "processor" )
            {
                GlobField[i + 1][abs(IndFaceLocal()[i][k]) - Start()[i]] =
                    field.boundaryFieldRef()[i][k];
            }
        }

        reduce(GlobField[i + 1], sumOp<Field<scalar>>());
    }

    return GlobField;
}

template<>
List<Field <vector>> ITHACAparallel::combineFields(
                      GeometricField<vector, fvPatchField, volMesh>& field)
{
    List<Field< vector>> GlobField(field.boundaryFieldRef().size() + 1);
    GlobField[0].resize(N_IF_glob);
    GlobField[0] = GlobField[0] * 0;

    // Assemble internalField
    for (int i = 0; i < field.size(); i++)
    {
        GlobField[0][indices()[i]] = field[i];
    }

    reduce(GlobField[0], sumOp<Field<vector>>());

    // Assemble BoundariField
    for (int i = 0; i < N_BF; i++)
    {
        GlobField[i + 1].resize(Gsize_BF()[i]);
        Field<vector> zero(Gsize_BF()[i], vector(0.0, 0.0, 0.0));
        GlobField[i + 1] = zero;
    }

    for (int i = 0; i < N_BF; i++)
    {
        for (int k = 0; k < field.boundaryFieldRef()[i].size(); k++)
        {
            if (IndFaceLocal()[i].size() > 0
                    && field.boundaryFieldRef()[i].type() != "zeroGradient"
                    && field.boundaryFieldRef()[i].type() != "processor" )
            {
                GlobField[i + 1][abs(IndFaceLocal()[i][k]) - Start()[i]] =
                    field.boundaryFieldRef()[i][k];
            }
        }

        reduce(GlobField[i + 1], sumOp<Field<vector>>());
    }

    return GlobField;
}
