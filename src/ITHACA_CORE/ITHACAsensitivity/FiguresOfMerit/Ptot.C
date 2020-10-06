#include "Ptot.H"

Ptot::Ptot() {}

Ptot::Ptot(label argc, char* argv[], label Nsampled)
{
    Npoints = Nsampled;
    modelOutput.resize(Npoints);
    modelOutput.setZero();
    name = {"Ptot"};
#include "setRootCase.H"
#include "createTime.H"
#include "createMesh.H"
#include "createpowerDens.H"
}

void Ptot::buildMO(std::string dir)
{
    Time& runTime = _runTime();
    fvMesh& mesh = _mesh();
    volScalarField& powerDens = _powerDens();
    std::string folder = dir;

    if (ITHACAutilities::check_folder(folder) == true)
    {
        for (label j = 0; j < Npoints; j++)
        {
            folder.append(std::to_string(j));
            folder.append("/");
            ITHACAstream::read_fields(ptrfield, powerDens, folder);
            auto k = ptrfield.last();
            modelOutput(j) = fvc::domainIntegrate(k).value();
            folder = dir;
            ptrfield.clear();
        }

        MObuilt = true;
    }
    else
    {
        std::cout << "Outputs of the model are not computed yet, programm aborted" <<
                  std::endl;
    }
}
