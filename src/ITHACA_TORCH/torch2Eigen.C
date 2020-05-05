#include "torch2Eigen.H"
using namespace ITHACAtorch::torch2Eigen;

namespace ITHACAtorch
{
namespace torch2Eigen
{

template<class type>
torch::Tensor eigenMatrix2torchTensor(
    Eigen::Matrix<type, Eigen::Dynamic, Eigen::Dynamic> _eigenMatrix)
{
    Eigen::MatrixXf eigenMatrix = _eigenMatrix.template cast <float> ();
    int rows = eigenMatrix.rows();
    int cols = eigenMatrix.cols();

    if (!eigenMatrix.IsRowMajor)
    {
        eigenMatrix.transposeInPlace();
    }

    return torch::from_blob(eigenMatrix.data(), {rows, cols}).clone();
}

template<class type>
Eigen::Matrix<type, Eigen::Dynamic, Eigen::Dynamic> torchTensor2eigenMatrix(
    torch::Tensor& torchTensor)
{
    std::string error_message("The provided tensor has " + std::to_string(
                                  torchTensor.dim()) + " dimensions and cannot be casted in a Matrix.");
    M_Assert(torchTensor.dim() <= 2, error_message.c_str());
    M_Assert(torchTensor.dim() != 0, "The provided tensor has 0 dimension");
    int rows;
    int cols;
    int nElem = torchTensor.size(0);

    for (int i = 1; i < torchTensor.dim(); i++)
    {
        nElem *= torchTensor.size(i);
    }

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

    typedef Eigen::Matrix<type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
    MatrixXf_rm;
    float* data_p = torchTensor.data_ptr<float>();
    std::vector<type> raw(nElem);

    for (int i = 0; i < nElem; i++)
    {
        float d(*(data_p + i));
        type a = static_cast <type>(d);
        raw[i]=a;
    }

    Eigen::Map<MatrixXf_rm> eigenMatrix(&raw[0], rows,
                                        cols);
    return eigenMatrix;
}

template torch::Tensor eigenMatrix2torchTensor(
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> eigenMatrix);
template Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>
torchTensor2eigenMatrix<float>(torch::Tensor& torchTensor);
template torch::Tensor eigenMatrix2torchTensor(
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> eigenMatrix);
template Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
torchTensor2eigenMatrix<double>(torch::Tensor& torchTensor);
template torch::Tensor eigenMatrix2torchTensor(
    Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> eigenMatrix);
template Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>
torchTensor2eigenMatrix<int>(torch::Tensor& torchTensor);

}

}
