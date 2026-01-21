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
List<label> surfaceIndexInt(T& field, const label patchInt,
                            const label patchExt)
{
    List<label>* result = new List<label>;

    for (size_t i = 0; i < field.mesh().boundaryMesh()[patchInt].size(); i++)
    {
        result->append(field.mesh().boundaryMesh()[patchInt].faceCells()[i]);
    }

    return *result;
}

template List<label> surfaceIndexInt(volScalarField& field,
                                     const label patchInt, const label patchExt);
template List<label> surfaceIndexInt(volVectorField& field,
                                     const label patchInt, const label patchExt);
template List<label> surfaceIndexInt(volTensorField& field,
                                     const label patchInt, const label patchExt);

template<typename T>
List<label> surfaceIndexExt(T& field, const label patchInt,
                            const label patchExt)
{
    return surfaceIndexInt(field, patchExt, patchInt);
}

template List<label> surfaceIndexExt(volScalarField& field,
                                     const label patchInt, const label patchExt);
template List<label> surfaceIndexExt(volVectorField& field,
                                     const label patchInt, const label patchExt);
template List<label> surfaceIndexExt(volTensorField& field,
                                     const label patchInt, const label patchExt);

template<typename T, typename V>
void surfaceValuesInt(T& field, const label patchInt, const label patchExt,
                      List<V>& result)
{
    List<label> indexes = surfaceIndexInt(field, patchInt, patchExt);

    for (size_t i = 0; i < indexes.size(); i++)
    {
        result.append(field[indexes[i]]);
    }
}

template void surfaceValuesInt(volScalarField& field, const label patchInt,
                               const label patchExt, List<scalar>& result);
template void surfaceValuesInt(volVectorField& field, const label patchInt,
                               const label patchExt, List<Foam::Vector<scalar >>& result);
template void surfaceValuesInt(volTensorField& field, const label patchInt,
                               const label patchExt, List<Foam::Tensor<scalar >>& result);

template<typename T, typename V>
void surfaceValuesExt(T& field, const label patchInt, const label patchExt,
                      List<V>& result)
{
    surfaceValuesInt(field, patchExt, patchInt, result);
}

template void surfaceValuesExt(volScalarField& field, const label patchInt,
                               const label patchExt, List<scalar>& result);
template void surfaceValuesExt(volVectorField& field, const label patchInt,
                               const label patchExt, List<Foam::Vector<scalar >>& result);
template void surfaceValuesExt(volTensorField& field, const label patchInt,
                               const label patchExt, List<Foam::Tensor<scalar >>& result);

template<typename T>
Foam::Vector<scalar> surfaceFindMirrorPoint(T& field, const label patchInt,
        const label patchExt, const label cellID)
{
    Foam::Vector<scalar> result = * (new Foam::Vector<scalar>);
    List<label> indexesInt = surfaceIndexInt(field, patchInt, patchExt);
    result = 2.0 * field.mesh().boundaryMesh()[patchInt].faceCentres()[cellID] -
             field.mesh().C()[indexesInt[cellID]];
    return result;
}

template Foam::Vector<scalar> surfaceFindMirrorPoint(volScalarField& field,
        const label patchInt, const label patchExt, const label cellID);
template Foam::Vector<scalar> surfaceFindMirrorPoint(volVectorField& field,
        const label patchInt, const label patchExt, const label cellID);
template Foam::Vector<scalar> surfaceFindMirrorPoint(volTensorField& field,
        const label patchInt, const label patchExt, const label cellID);

template<typename T>
label surfaceFindClosest(T& field, const label patchInt, const label patchExt,
                         Foam::Vector<scalar> point)
{
    label result = 0;
    scalar dist = 0;
    List<label> indexesExt = surfaceIndexExt(field, patchInt, patchExt);
    dist = mag(point - field.mesh().C()[indexesExt[0]]);

    for (int i = 1; i < indexesExt.size(); i++)
    {
        scalar temp = mag(point - field.mesh().C()[indexesExt[i]]);

        if (temp <= dist)
        {
            dist = temp;
            result = i;
        }
    }

    result = indexesExt[result];
    return result;
}

template label surfaceFindClosest(volScalarField& field, const label patchInt,
                                  const label patchExt, Foam::Vector<scalar> point);
template label surfaceFindClosest(volVectorField& field, const label patchInt,
                                  const label patchExt, Foam::Vector<scalar> point);
template label surfaceFindClosest(volTensorField& field, const label patchInt,
                                  const label patchExt, Foam::Vector<scalar> point);

template<typename T, typename V>
void surfaceAverage(T& field, const label patchInt, const label patchExt,
                    List<V>& result)
{
    List<label> indexesInt = surfaceIndexInt(field, patchInt, patchExt);
    result.resize(0);

    for (int i = 0; i < indexesInt.size(); i++)
    {
        Foam::Vector<scalar> mirror = surfaceFindMirrorPoint(field, patchInt, patchExt,
                                      i);
        label closest = surfaceFindClosest(field, patchInt, patchExt, mirror);
        result.append(0.5 * field[closest] + 0.5 * field[indexesInt[i]]);
    }
}

template void surfaceAverage(volScalarField& field, const label patchInt,
                             const label patchExt, List<scalar>& result);
template void surfaceAverage(volVectorField& field, const label patchInt,
                             const label patchExt, List<Foam::Vector<scalar >>& result);
template void surfaceAverage(volTensorField& field, const label patchInt,
                             const label patchExt, List<Foam::Tensor<scalar >>& result);

template<typename T, typename V>
void surfaceJump(T& field, const label patchInt, const label patchExt,
                 List<V>& result)
{
    List<label> indexesInt = surfaceIndexInt(field, patchInt, patchExt);
    result.resize(0);

    for (int i = 0; i < indexesInt.size(); i++)
    {
        Foam::Vector<scalar> mirror = surfaceFindMirrorPoint(field, patchInt, patchExt,
                                      i);
        label closest = surfaceFindClosest(field, patchInt, patchExt, mirror);
        result.append(field[closest] - field[indexesInt[i]]);
    }
}

template void surfaceJump(volScalarField& field, const label patchInt,
                          const label patchExt, List<scalar>& result);
template void surfaceJump(volVectorField& field, const label patchInt,
                          const label patchExt, List<Foam::Vector<scalar >>& result);
template void surfaceJump(volTensorField& field, const label patchInt,
                          const label patchExt, List<Foam::Tensor<scalar >>& result);
} // End namespace ITHACAsurfacetools
} // End namespace ITHACAutilities

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
