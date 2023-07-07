/*---------------------------------------------------------------------------*\
     ██╗████████╗██╗  ██╗ █████╗  ██████╗ █████╗       ███████╗██╗   ██╗
     ██║╚══██╔══╝██║  ██║██╔══██╗██╔════╝██╔══██╗      ██╔════╝██║   ██║
     ██║   ██║   ███████║███████║██║     ███████║█████╗█████╗  ██║   ██║
     ██║   ██║   ██╔══██║██╔══██║██║     ██╔══██║╚════╝██╔══╝  ╚██╗ ██╔╝
     ██║   ██║   ██║  ██║██║  ██║╚██████╗██║  ██║      ██║      ╚████╔╝
     ╚═╝   ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝      ╚═╝       ╚═══╝

 * In real Time Highly Advanced Computational Applications for Finite Volumes
 * Copyright (C) 2017 by the ITHACA-FV authors
-------------------------------------------------------------------------------

License
    This file is part of ITHACA-FV

    ITHACA-FV is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    ITHACA-FV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Lesser General Public License for more details.


    You should have received a copy of the GNU Lesser General Public License
    along with ITHACA-FV. If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "ITHACAsurfacetools.H"


namespace ITHACAutilities
{
namespace ITHACAsurfacetools
{

  template<typename T>
  List<label> surfaceIndexInt(T& field, const label patchInt, const label patchExt)
  {
    if (field.mesh().boundaryMesh()[patchInt].type() != "cyclic") {
      Info << "PATCH TYPE MUST BE cyclic, not " << field.mesh().boundaryMesh()[patchInt].type() <<endl;
      return *(new List<label>);
    }

    List<label>* result = new List<label>;

    for (size_t i = 0; i < field.mesh().boundaryMesh()[patchInt].size(); i++) {
      result->append(field.mesh().boundaryMesh()[patchInt].faceCells()[i]);
    }

    return *result;
  }

  template List<label> surfaceIndexInt(volScalarField& field, const label patchInt, const label patchExt);
  template List<label> surfaceIndexInt(volVectorField& field, const label patchInt, const label patchExt);
  template List<label> surfaceIndexInt(volTensorField& field, const label patchInt, const label patchExt);

  template<typename T>
  List<label> surfaceIndexExt(T& field, const label patchInt, const label patchExt)
  {
    return surfaceIndexInt(field, patchExt, patchInt);
  }

  template List<label> surfaceIndexExt(volScalarField& field, const label patchInt, const label patchExt);
  template List<label> surfaceIndexExt(volVectorField& field, const label patchInt, const label patchExt);
  template List<label> surfaceIndexExt(volTensorField& field, const label patchInt, const label patchExt);

  template<typename T, typename V>
  void surfaceValuesInt(T& field, const label patchInt, const label patchExt, List<V>& result)
  {
    if (field.mesh().boundaryMesh()[patchInt].type() != "cyclic") {
      Info << "PATCH TYPE MUST BE cyclic, not " << field.mesh().boundaryMesh()[patchInt].type() <<endl;
    }
    else
    {
      List<label> indexes = surfaceIndexInt(field, patchInt, patchExt);

      for (size_t i = 0; i < indexes.size(); i++) {
        result.append(field[indexes[i]]);
      }
    }
  }

  template void surfaceValuesInt(volScalarField& field, const label patchInt, const label patchExt, List<scalar>& result);
  template void surfaceValuesInt(volVectorField& field, const label patchInt, const label patchExt, List<Foam::Vector<scalar>>& result);
  template void surfaceValuesInt(volTensorField& field, const label patchInt, const label patchExt, List<Foam::Tensor<scalar>>& result);

  template<typename T, typename V>
  void surfaceValuesExt(T& field, const label patchInt, const label patchExt, List<V>& result)
  {
    if (field.mesh().boundaryMesh()[patchInt].type() != "cyclic") {
      Info << "PATCH TYPE MUST BE cyclic, not " << field.mesh().boundaryMesh()[patchInt].type() <<endl;
    }
    else
    {
      surfaceValuesInt(field, patchExt, patchInt, result);
    }
  }

  template void surfaceValuesExt(volScalarField& field, const label patchInt, const label patchExt, List<scalar>& result);
  template void surfaceValuesExt(volVectorField& field, const label patchInt, const label patchExt, List<Foam::Vector<scalar>>& result);
  template void surfaceValuesExt(volTensorField& field, const label patchInt, const label patchExt, List<Foam::Tensor<scalar>>& result);

} // End namespace ITHACAsurfacetools
} // End namespace ITHACAutilities

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
