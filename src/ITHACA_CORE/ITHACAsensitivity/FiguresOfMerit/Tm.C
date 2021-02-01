#include "Tm.H"

Tm::Tm() {}

Tm::Tm(int argc, char* argv[], label Nsampled)
{
    Npoints = Nsampled;
    modelOutput.resize(Npoints);
    modelOutput.setZero();
    name = {"Tm"};
#include "setRootCase.H"
#include "createTime.H"
#include "createMesh.H"
#include "createT.H"
}

void Tm::buildMO(std::string dir)
{
    Time& runTime = _runTime();
    fvMesh& mesh = _mesh();
    volScalarField& T = _T();
    std::string folder = dir;

    if (ITHACAutilities::check_folder(folder) == true)
    {
        for (label j = 0; j < Npoints; j++)
        {
            folder.append(std::to_string(j));
            folder.append("/");
            ITHACAstream::read_fields(ptrfield, T, folder);
            auto k = ptrfield.last();
            modelOutput(j) = k.weightedAverage(mesh.V()).value();
            folder = dir;
            ptrfield.clear();
        }

        MObuilt = true;
    }
    else
    {
        std::cout << "Outputs of the model are not computed yet, programm aborted" <<
                  std::endl;
        exit(0);
    }
}

