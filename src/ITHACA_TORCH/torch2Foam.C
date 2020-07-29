#include "torch2Foam.H"
using namespace ITHACAtorch::torch2Foam;

namespace ITHACAtorch
{
namespace torch2Foam
{

template<>
torch::Tensor field2Torch(Field<vector>& field)
{
    int rows = field.size() * 3;
    int cols = 1;
    double* dataPtr = &field[0][0];
    return torch::from_blob(dataPtr, {rows * cols}, {torch::kFloat64}).clone().to(
               torch::kFloat32);
}

template<>
torch::Tensor field2Torch(Field<scalar>& field)
{
    int rows = field.size();
    int cols = 1;
    double* dataPtr = &field[0];
    return torch::from_blob(dataPtr, {rows * cols}, {torch::kFloat64}).clone().to(
               torch::kFloat32);
}

template<>
Field<vector> torch2Field(torch::Tensor& torchTensor)
{
    std::string error_message("The provided tensor has " + std::to_string(
                                  torchTensor.dim()) +
                              " dimensions and with the current implementation only 1-D tensor can be casted in an OpenFOAM field.");
    M_Assert(torchTensor.dim() == 1, error_message.c_str());
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
    M_Assert(torchTensor.dim() == 1, error_message.c_str());
    M_Assert(torchTensor.dim() != 0, "The provided tensor has 0 dimension");
    Field<scalar> a(torchTensor.numel());
    std::memcpy(&a[0], torchTensor.to(torch::kFloat64).data_ptr(),
                sizeof (double)*torchTensor.numel());
    return a;
}



}

}
