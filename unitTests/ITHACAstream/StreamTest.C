#include "ITHACAstream.H"

bool ReadAndWriteTensor()
{
    bool esit = false;
    Eigen::Tensor<double, 3> output;
    Eigen::Tensor<double, 3> input;
    output.resize(3,4,5);
    output.setRandom();

    ITHACAstream::SaveDenseTensor(output, "./","output");
    ITHACAstream::ReadDenseTensor(input, "./","output");
    Eigen::Tensor<double,0> difference = ((output - input).abs().sum()); 
    Eigen::Tensor<double,0> comp;
    comp(0)=1e-18;
    if(difference(0) < comp(0))   
    {
         esit = true;
         std::cout << "> Read And Write Test for tensors succeeded!" << std::endl;
    }
    return esit;
}

int main(int argc, char **argv)
{
    ReadAndWriteTensor();
    return 0;
}
