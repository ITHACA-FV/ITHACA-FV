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

#include "torchUTILITIES.H"

namespace ITHACAtorch
{
torch::Tensor removeConstValues(torch::Tensor input, std::vector<int>& indices,
                                std::vector<double>& constValues)
{
    torch::Tensor output;
    return output;
}

bool isConst(torch::Tensor& tTensor)
{
    const float* ptr = (float*)tTensor.data_ptr();
    bool all_equal = true;

    for (std::size_t i = 1, s = tTensor.numel(); i < s && all_equal; i++)
    {
        all_equal = *ptr == *(ptr + i);
    }

    return all_equal;
}

void save(const torch::Tensor& torchTensor, const std::string fname)
{
    unsigned int dim = torchTensor.dim();
    unsigned int shape[dim];
    float* data_p = torchTensor.data_ptr<float>();

    for (unsigned int i = 0; i < torchTensor.dim(); i++)
    {
        shape[i] = (unsigned int) torchTensor.size(i);
    }

    cnpy::npy_save(fname, data_p, shape, dim);
}

torch::Tensor load(const std::string fname)
{
    cnpy::NpyArray arr = cnpy::npy_load(fname);
    at::IntList shape[arr.shape.size()];
    std::vector<int64_t> dims(arr.shape.size());

    for (int i = 0; i < arr.shape.size(); i++)
    {
        dims[i] = (int64_t) arr.shape[i];
    }

    torch::Tensor tensor = torch::randn({2, 2, 3});
    return torch::from_blob(arr.data, dims).clone();
    //return tensor;
}


}
