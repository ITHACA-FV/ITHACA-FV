#include <iostream>
#include <torch/script.h>
#include <torch/torch.h>
#include <Eigen/Eigen>

namespace ITHACAtorch
{
    namespace torch2Eigen
    {
        ///
        /// @brief      { function_description }
        ///
        /// @param      torchTensor  The torch tensor
        ///
        /// @return     { description_of_the_return_value }
        ///
        torch::Tensor eigenMatrix2torchTensor(Eigen::MatrixXf& eigenMatrix);
        

        ///
        /// @brief      { function_description }
        ///
        /// @param      torch  The torch
        ///
        /// @return     { description_of_the_return_value }
        ///
        Eigen::MatrixXf torchTensor2eigenMatrix(torch::Tensor& torchTensor);
    }
    
}
