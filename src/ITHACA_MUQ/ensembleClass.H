#ifndef ensemble_H
#define ensemble_H

#include"muq2ithaca.H"

namespace ITHACAmuq
{
//--------------------------------------------------------------------------
/// @brief      Class for ensembles, each column is a sample
///
class ensemble
{
    int Nsamples = 0;
    int samplesSize = 0;
    Eigen::MatrixXd samples;

    public:
        // Constructors
        ensemble();
        ensemble(int _Nsamples, int _samplesSize);
        ensemble(Eigen::MatrixXd _samples);
        ~ensemble() {};

        // Functions

        //--------------------------------------------------------------------------
        /// Get number of samples 
        int getSize();

        //--------------------------------------------------------------------------
        /// Get matrix of samples
        Eigen::VectorXd getSample(int sampleI);

        //--------------------------------------------------------------------------
        /// Get matrix of samples
        Eigen::MatrixXd getSamples();

        //--------------------------------------------------------------------------
        /// Assign a sample
        void assignSample(int sampleI, Eigen::VectorXd sample);

        //--------------------------------------------------------------------------
        /// Assing the samples matrix 
        void assignSamples(Eigen::MatrixXd _samples);

        //--------------------------------------------------------------------------
        /// Samples mean
        Eigen::VectorXd mean();

        //--------------------------------------------------------------------------
        /// Samples autocovariance
        Eigen::MatrixXd cov();

        //--------------------------------------------------------------------------
        /// Samples cross covariance
        Eigen::MatrixXd crossCov(Eigen::MatrixXd samples2);
};
}
#endif
