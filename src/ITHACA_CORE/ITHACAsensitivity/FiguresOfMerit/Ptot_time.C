#include "Ptot_time.H"

Ptot_time::Ptot_time() {}

Ptot_time::Ptot_time(int argc, char* argv[], label Nsampled)
{
    Npoints = Nsampled;
    modelOutput.resize(Npoints);
    modelOutput.setZero();
    name = {"Ptot_time"};
#include "setRootCase.H"
#include "createTime.H"
#include "createMesh.H"
#include "createpowerDens.H"
}

void Ptot_time::buildMO(std::string dir, label t)
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

            for (label i = 0; i < ptrfield.size(); i++)
            {
                if (i == t)
                {
                    modelOutput(j) = fvc::domainIntegrate(ptrfield[i]).value();
                }
            }

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
