#include "torch2Eigen.H"
using namespace ITHACAtorch::torch2Eigen;

namespace ITHACAtorch
{
namespace torch2Eigen
{

torch::Tensor eigenMatrix2torchTensor(Eigen::MatrixXf eigenMatrix)
{
    int rows = eigenMatrix.rows();
    int cols = eigenMatrix.cols();

    if (!eigenMatrix.IsRowMajor)
    {
        eigenMatrix.transposeInPlace();
    }

    return torch::from_blob(eigenMatrix.data(), {rows, cols}).clone();
}

Eigen::MatrixXf torchTensor2eigenMatrix(torch::Tensor& torchTensor)
{
    std::string error_message("The provided tensor has " + std::to_string(
                                  torchTensor.dim()) + " dimensions and cannot be casted in a Matrix.");
    M_Assert(torchTensor.dim() <= 2, error_message.c_str());
    M_Assert(torchTensor.dim() != 0, "The provided tensor has 0 dimension");
    int rows;
    int cols;

    if (torchTensor.dim() == 1)
    {
        rows = torchTensor.size(0);
        cols = 1;
    }
    else
    {
        rows = torchTensor.size(0);
        cols = torchTensor.size(1);
    }

    typedef Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
    MatrixXf_rm;
    float* data = torchTensor.data_ptr<float>();
    Eigen::Map<MatrixXf_rm> eigenMatrix(data, rows,
                                        cols);
    return eigenMatrix;
}
}

}
