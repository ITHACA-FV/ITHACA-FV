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
#include "ConvLayer.H"

namespace ITHACAtorch
{
template<class Type, template<class> class PatchField, class GeoMesh>
ConvLayer<Type, PatchField, GeoMesh>::ConvLayer(
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& snapshots):
    _snapshots(snapshots),
    mesh(snapshots[0].mesh()),
    convDict(autoPtr<IOdictionary>
             (
                 new IOdictionary
                 (
                     IOobject
                     (
                         "convDict",
                         mesh.time().system(),
                         mesh,
                         IOobject::MUST_READ,
                         IOobject::NO_WRITE
                     )
                 )
             ))
{
    flt = autoPtr<Filter>(Filter::New(word(convDict().lookup("Filter")),
                                      convDict()));
    domainDivision = Vector<label>(convDict().lookup("domainDivision"));
    filterSize = Vector<scalar>(convDict().lookup("filterSize"));
    domainSize = mesh.bounds().max() - mesh.bounds().min();
    setDomainDivision(domainDivision[0], domainDivision[1], domainDivision[2]);
    setFilterSize(filterSize[0], filterSize[1], filterSize[2]);
    weights = flt->apply(cellsInBoxes, convPoints, mesh);
}

template<class Type, template<class> class PatchField, class GeoMesh>
void ConvLayer<Type, PatchField, GeoMesh>::setDomainDivision(label Nx, label Ny,
        label Nz)
{
    M_Assert(((Nx != 1 && mesh.solutionD()[0] != -1 )  || (Nx == 1 &&
              mesh.solutionD()[0] == -1)),
             "The mesh has valid components only along the y and z directions, set Nx = 1");
    M_Assert(((Ny != 1 && mesh.solutionD()[1] != -1 )  || (Ny == 1 &&
              mesh.solutionD()[1] == -1)),
             "The mesh has valid components only along the x and z directions, set Ny = 1");
    M_Assert(((Nz != 1 && mesh.solutionD()[2] != -1 )  || (Nz == 1 &&
              mesh.solutionD()[2] == -1)),
             "The mesh has valid components only along the x and y directions, set Nz = 1");

    for (label i = 0; i < ds.size(); i++)
    {
        if (mesh.solutionD()[i] != -1 && domainDivision[i] != 1)
        {
            ds[i] = domainSize[i] / (domainDivision[i] - 1);
        }
        else
        {
            ds[i] = domainSize[i] / 2;
        }
    }

    convPoints = List<point>(Nx * Ny * Nz);
    label index = 0;

    for (label i = 0; i < Nx; i++)
    {
        for (label j = 0; j < Ny; j++)
        {
            for (label k = 0; k < Nz; k++)
            {
                if (i == 0 && Nx == 1)
                {
                    i = 1;
                }

                if (j == 0 && Ny == 1)
                {
                    j = 1;
                }

                if (k == 0 && Nz == 1)
                {
                    k = 1;
                }

                convPoints[index] = mesh.bounds().min() + cmptMultiply((ds * i), vector(1, 0,
                                    0)) + cmptMultiply((ds * j), vector(0, 1, 0)) + cmptMultiply((ds * k), vector(0,
                                            0, 1));
                index++;
            }
        }
    }

    isDomainDivisionSet = true;
}

template<class Type, template<class> class PatchField, class GeoMesh>
void ConvLayer<Type, PatchField, GeoMesh>::setFilterSize(double dx, double dy,
        double dz)
{
    M_Assert(isDomainDivisionSet, "You need to set the division before.");
    filterSize[0] = dx;
    filterSize[1] = dy;
    filterSize[2] = dz;
    cellsInBoxes.resize(convPoints.size());

    for (label i = 0; i < convPoints.size(); i++)
    {
        cellSet a(mesh, "set", 0);
        point mini = convPoints[i] - filterSize / 2;
        point maxi = convPoints[i] + filterSize / 2;
        treeBoundBox boxi(mini, maxi);
        List<treeBoundBox> l;
        l.append(boxi);
        boxToCell finding(mesh, l);
#if OPENFOAM >= 1812
        finding.verbose(false);
#endif
        finding.applyToSet(topoSetSource::ADD, a);
        cellsInBoxes[i] = a.toc();
    }

    isFilterSizeSet = true;
}

// template<class Type, template<class> class PatchField, class GeoMesh>
// torch::Tensor ConvLayer<Type, PatchField, GeoMesh>::filter()
// {
//     M_Assert(isDomainDivisionSet &&
//              isFilterSizeSet,
//              "You need to set domainDivision and filterSize before calling the filter funtion.");
//     torch::Tensor output;
//     return output;
// }

template<>
torch::Tensor ConvLayer<scalar, fvPatchField, volMesh>::filter()
{
    M_Assert(isDomainDivisionSet &&
             isFilterSizeSet,
             "You need to set domainDivision and filterSize before calling the filter funtion.");
    label Nx = domainDivision[0];
    label Ny = domainDivision[1];
    label Nz = domainDivision[2];
    torch::Tensor output = torch::zeros({_snapshots.size(), 1, Nx, Ny, Nz});
    auto foo_a = output.accessor<float, 5>();

    for (label i = 0; i < _snapshots.size(); i++)
    {
        label index = 0;

        for (label j = 0; j < Nx; j++)
        {
            for (label k = 0; k < Ny; k++)
            {
                for (label l = 0; l < Nz; l++)
                {
                    for (label p = 0; p < cellsInBoxes[index].size(); p++)
                    {
                        foo_a[i][0][j][k][l] += _snapshots[i][cellsInBoxes[index][p]] *
                                                weights[index][p];
                    }

                    index++;
                }
            }
        }
    }

    return output;
}


template<>
torch::Tensor ConvLayer<vector, fvPatchField, volMesh>::filter()
{
    M_Assert(isDomainDivisionSet &&
             isFilterSizeSet,
             "You need to set domainDivision and filterSize before calling the filter funtion.");
    label Nx = domainDivision[0];
    label Ny = domainDivision[1];
    label Nz = domainDivision[2];
    torch::Tensor output = torch::zeros({_snapshots.size(), 3, Nx, Ny, Nz});
    auto foo_a = output.accessor<float, 5>();

    for (label i = 0; i < _snapshots.size(); i++)
    {
        label index = 0;

        for (label j = 0; j < Nx; j++)
        {
            for (label k = 0; k < Ny; k++)
            {
                for (label l = 0; l < Nz; l++)
                {
                    for (label p = 0; p < cellsInBoxes[index].size(); p++)
                    {
                        foo_a[i][0][j][k][l] += _snapshots[i][cellsInBoxes[index][p]][0] *
                                                weights[index][p];
                        foo_a[i][1][j][k][l] += _snapshots[i][cellsInBoxes[index][p]][1] *
                                                weights[index][p];
                        foo_a[i][2][j][k][l] += _snapshots[i][cellsInBoxes[index][p]][2] *
                                                weights[index][p];
                    }

                    index++;
                }
            }
        }
    }

    return output;
}


template class ConvLayer<scalar, fvPatchField, volMesh>;
template class ConvLayer<vector, fvPatchField, volMesh>;
// template class ConvLayer<scalar, fvsPatchField, surfaceMesh>;

};
