#include"ensembleClass.H"

// Kabir: In the ensemble class, we have methods for calculating the mean, covariance, and cross-covariance of samples. 
namespace ITHACAmuq
{
ensemble::ensemble() {}

ensemble::ensemble(int _Nsamples, int _samplesSize)
{
    Nsamples = _Nsamples;
    samplesSize = _samplesSize;
    samples = Eigen::MatrixXd(samplesSize, Nsamples);
}

ensemble::ensemble(Eigen::MatrixXd _samples)
{
    samples = _samples;
    Nsamples = samples.cols();
    samplesSize = samples.rows();
}

int ensemble::getSize()
{
    return Nsamples;
}

Eigen::VectorXd ensemble::getSample(int sampleI)
{
    std::string message = "Sample index (" + std::to_string(sampleI) + 
        ")  is bigger than the number of samples ("  
        + std::to_string(samplesSize) + ")";

    M_Assert(sampleI < samples.cols(), message.c_str());
    return samples.col(sampleI);
}

Eigen::MatrixXd ensemble::getSamples()
{
    return samples;
}

void ensemble::assignSample(int sampleI, Eigen::VectorXd _sample)
{
    std::string message = "The input sample should have the right size. Here input size = " + std::to_string(_sample.size())
        + " while samplesSize = "
        + std::to_string(samplesSize);

    M_Assert(_sample.size() == samplesSize, message.c_str());
    M_Assert(sampleI < Nsamples, "The index of the input sample is bigger than the number of samples in the ensemble");
    samples.col(sampleI) = _sample;
}

void ensemble::assignSamples(Eigen::MatrixXd _samples)
{
    M_Assert(_samples.size() > 0, "Not valid input matrix");
    Nsamples = _samples.cols();
    if(samplesSize == 0)
    {    
        samplesSize = _samples.rows();
    }
    else
    {
        std::string message = "The input samples have a different size. Input size = " + std::to_string(_samples.rows())
            + " while samplesSize = "
            + std::to_string(samplesSize);
        M_Assert(_samples.rows() == samplesSize, message.c_str());
    }
    samples = _samples;
}

Eigen::VectorXd ensemble::mean()
{
    return samples.rowwise().mean();
}

Eigen::MatrixXd ensemble::cov()
{
    Eigen::MatrixXd centered = samples.colwise() - samples.rowwise().mean();
    return (centered * centered.adjoint()) / double(Nsamples - 1);
}

Eigen::MatrixXd ensemble::crossCov(Eigen::MatrixXd samples2)
{
    std::string message = "Input has wrong number of samples. Nsamples = " + 
        std::to_string(Nsamples) + "while samples2.cols() = " + 
        std::to_string(samples2.cols());
    M_Assert(samples2.cols() == Nsamples, message.c_str());
    Eigen::MatrixXd centered1 = samples.colwise() - samples.rowwise().mean();
    Eigen::MatrixXd centered2 = samples2.colwise() - samples2.rowwise().mean();
    return (centered1 * centered2.adjoint()) / double(Nsamples - 1);
}
}
