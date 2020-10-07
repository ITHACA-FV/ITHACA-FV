#include "FofM.H"

FofM::FofM() {}

FofM::FofM(int argc, char* argv[], label Nsampled)
{
    Npoints = Nsampled;
    modelOutput.resize(Npoints);
    modelOutput.setZero();
    name = {"Base class"};
}

void FofM::buildMO(std::string dir)
{
    Info << "The method FofM::createfOFm() in F.C is a virtual method"
         << endl;
    Info << "It must be overridden, exiting the code" << endl;
    exit(0);
}



