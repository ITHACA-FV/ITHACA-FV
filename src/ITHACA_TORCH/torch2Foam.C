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
#include "torch2Foam.H"
using namespace ITHACAtorch::torch2Foam;

namespace ITHACAtorch
{
namespace torch2Foam
{

template<>
torch::Tensor field2Torch(Field<vector>& field)
{
    int rows = 1;
    int cols = field.size() * 3;
    double* dataPtr = &field[0][0];
    return torch::from_blob(dataPtr, {rows, cols}, {torch::kFloat64}).clone().to(
               torch::kFloat32);
}

template<>
torch::Tensor field2Torch(Field<scalar>& field)
{
    int rows = 1;
    int cols = field.size();
    double* dataPtr = &field[0];
    return torch::from_blob(dataPtr, {rows, cols}, {torch::kFloat64}).clone().to(
               torch::kFloat32);
}

template<>
Field<vector> torch2Field(torch::Tensor& torchTensor)
{
    std::string error_message("The provided tensor has " + std::to_string(
                                  torchTensor.dim()) +
                              " dimensions and with the current implementation only 1-D tensor can be casted in an OpenFOAM field.");
    M_Assert(torchTensor.dim() <= 2, error_message.c_str());
    M_Assert(torchTensor.dim() != 0, "The provided tensor has 0 dimension");
    Field<vector> a(torchTensor.numel() / 3);
    std::memcpy(&a[0][0], torchTensor.to(torch::kFloat64).data_ptr(),
                sizeof (double)*torchTensor.numel());
    return a;
}

template<>
Field<scalar> torch2Field(torch::Tensor& torchTensor)
{
    std::string error_message("The provided tensor has " + std::to_string(
                                  torchTensor.dim()) +
                              " dimensions and with the current implementation only 1-D tensor can be casted in an OpenFOAM field.");
    M_Assert(torchTensor.dim() <= 2, error_message.c_str());
    M_Assert(torchTensor.dim() != 0, "The provided tensor has 0 dimension");
    Field<scalar> a(torchTensor.numel());
    std::memcpy(&a[0], torchTensor.to(torch::kFloat64).data_ptr(),
                sizeof (double)*torchTensor.numel());
    return a;
}

template<>
torch::Tensor ptrList2Torch(PtrList<Field<vector>>& ptrList)
{
    int Nrows = ptrList.size();
    int Ncols = ptrList[0].size() * 3;
    torch::Tensor out = torch::randn({Nrows, Ncols});

    for (auto i = 0; i < ptrList.size(); i++)
    {
        out.slice(0, i, i + 1) = field2Torch(ptrList[i]);
    }

    return out;
}

template<>
torch::Tensor ptrList2Torch(PtrList<Field<scalar>>& ptrList)
{
    int Nrows = ptrList.size();
    int Ncols = ptrList[0].size();
    torch::Tensor out = torch::randn({Nrows, Ncols});

    for (auto i = 0; i < ptrList.size(); i++)
    {
        out.slice(0, i, i + 1) = field2Torch(ptrList[i]);
    }

    return out;
}

template<class type_f>
PtrList<Field<type_f>> torch2PtrList(torch::Tensor& tTensor)
{
    PtrList<Field<type_f>> out;

    for (auto i = 0; i < tTensor.size(0); i++)
    {
        torch::Tensor t = tTensor.slice(0, i, i + 1);
        out.append(tmp<Field<type_f>>(torch2Field<type_f>(t)));
    }

    return out;
}

template PtrList<Field<scalar>> torch2PtrList<scalar>(torch::Tensor& tTensor);
template PtrList<Field<vector>> torch2PtrList<vector>(torch::Tensor& tTensor);




}

}
