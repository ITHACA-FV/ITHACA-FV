#include "ITHACAstream.H"
#include "ITHACAsurfacetools.H"
#include "fvCFD.H"
#include "volPointInterpolation.H"
#include "fvOptions.H"
#include <complex>
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <tuple>

using namespace ITHACAutilities;
using namespace ITHACAsurfacetools;

int main(int argc, char **argv)
{
  #include "setRootCaseLists.H"
  #include "createTime.H"
  #include "createMesh.H"
  #include "createFields.H"

  IOdictionary ITHACAdict
  (
    IOobject
    (
      "ITHACAdict",
      runTime.system(),
      mesh,
      IOobject::MUST_READ,
      IOobject::NO_WRITE
    )
  );

  label patchInt = ITHACAdict.lookupOrDefault<label>("patchInt", 0);
  label patchExt = ITHACAdict.lookupOrDefault<label>("patchExt", 0);

  List<label> indexesInt = surfaceIndexInt(U, patchInt, patchExt);
  List<label> indexesExt = surfaceIndexExt(U, patchInt, patchExt);
  List<Foam::Vector<double>> patchValuesInt = *(new List<Foam::Vector<double>>);
  surfaceValuesInt(U, patchInt, patchExt, patchValuesInt);
  List<Foam::Vector<double>> patchValuesExt = *(new List<Foam::Vector<double>>);
  surfaceValuesExt(U, patchInt, patchExt, patchValuesExt);
  List<Foam::Vector<double>> average = *(new List<Foam::Vector<double>>);
  surfaceAverage(U, patchInt, patchExt, average);
  List<Foam::Vector<double>> jump = *(new List<Foam::Vector<double>>);
  surfaceJump(U, patchInt, patchExt, jump);

  Info << "surfaceIndexInt(U, 6) = " << endl;
  for (size_t i = 0; i < indexesInt.size(); i++) {
    Info << "   " << indexesInt[i] << endl;
  }

  Info << endl << "surfaceIndexExt(U, 6) = " << endl;
  for (size_t i = 0; i < indexesExt.size(); i++) {
    Info << "   " << indexesExt[i] << endl;
  }

  Info << endl << "surfaceValuesInt(U, 6) = " << endl;
  for (size_t i = 0; i < patchValuesInt.size(); i++) {
    Info << "   " << patchValuesInt[i] << endl;
  }

  Info << endl << "surfaceValuesExt(U, 6) = " << endl;
  for (size_t i = 0; i < patchValuesExt.size(); i++) {
    Info << "   " << patchValuesExt[i] << endl;
  }

  Info << endl << "surfaceAverage(U, 6) = " << endl;
  for (size_t i = 0; i < average.size(); i++) {
    Info << "   " << average[i] << endl;
  }

  Info << endl << "surfaceJump(U, 6) = " << endl;
  for (size_t i = 0; i < jump.size(); i++) {
    Info << i << " :   " << jump[i] << endl;
  }

  return 0;
}
