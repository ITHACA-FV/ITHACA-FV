#include "torch2Eigen.H"
using namespace ITHACAtorch::torch2Eigen;

namespace ITHACAtorch
{
namespace torch2Eigen
{

torch::Tensor eigenMatrix2torchTensor(Eigen::MatrixXf& eigenMatrix)
{
    torch::Tensor a;
    return a;
}

Eigen::MatrixXf torchTensor2eigenMatrix(torch::Tensor& torchTensor)
{
    typedef Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXf_rm;
    float* data = torchTensor.data_ptr<float>();
    Eigen::Map<MatrixXf_rm> eigenMatrix(data, torchTensor.size(0),
                        torchTensor.size(1));
    return eigenMatrix;
}
}

}
