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

  Info << "surfaceIndexInt(U, 6) = " << indexesInt << endl;
  Info << "surfaceIndexExt(U, 6) = " << indexesExt << endl;
  Info << "surfaceValuesInt(U, 6) = " << patchValuesInt << endl;
  Info << "surfaceValuesExt(U, 6) = " << patchValuesExt << endl;
  Info << "surfaceAverage(U, 6) = " << average << endl;
  Info << "surfaceJump(U, 6) = " << jump << endl;

  return 0;
}
