#include "muq2ithaca.H"

namespace ITHACAmuq
{
namespace muq2ithaca
{
Eigen::MatrixXd EnsembleKalmanFilter(Eigen::MatrixXd prior,
                                     Eigen::VectorXd measurements,
                                     Eigen::MatrixXd measurementsCov,
                                     Eigen::MatrixXd observedState)
{
    M_Assert(measurements.rows() == observedState.rows(),
             "The observed state should have the same dimention of the measurements");
    M_Assert(observedState.cols() == prior.cols(),
             "The input matrices should all have the samples on the columns");
    M_Assert(measurementsCov.rows() == measurements.rows()
             && measurementsCov.cols() == measurements.rows(),
             "Wrong measurements covariance matrix");
    unsigned Nseeds = prior.cols();
    unsigned stateDim = prior.rows();
    unsigned measDim = measurements.size();
    Eigen::MatrixXd HA(measDim, Nseeds);
    Eigen::VectorXd observedStateMean = observedState.rowwise().mean();
    Eigen::MatrixXd A(stateDim, Nseeds);
    Eigen::VectorXd priorMean = prior.rowwise().mean();
    auto measNoise = std::make_shared<muq::Modeling::Gaussian>
                     (Eigen::VectorXd::Zero(measDim), measurementsCov);
    Eigen::MatrixXd D(measDim, Nseeds);

    for (unsigned i = 0; i < Nseeds; i++)
    {
        A.col(i) = prior.col(i) - priorMean;
        HA.col(i) = observedState.col(i) - observedStateMean;
        D.col(i) = measurements + measNoise->Sample();
    }

    Eigen::MatrixXd Y = D -
                        observedState; //diff measurement data and simulated data
    Eigen::MatrixXd P = HA * HA.transpose() / (Nseeds - 1.) + measurementsCov;
    Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(P);
    auto P_rank = lu_decomp.rank();

    if (P_rank < P.cols() || P_rank < P.rows())
    {
        std::cout << "Pseudo inverse of P should be implemented in the EnKF, exiting" <<
                  std::endl;
        std::cout << P << std::endl;
        exit(10);
    }
    else
    {
        P = P.inverse();
    }

    //KF update
    Eigen::MatrixXd M = P * Y;
    Eigen::MatrixXd Z = (1. / (Nseeds - 1.)) * A * HA.transpose();
    //Kalman Gain
    Eigen::MatrixXd K = (1. / (Nseeds - 1.)) * A * HA.transpose() * P;

    if ((K.array() < 0.0).any())
    {
        std::cout << "WARNING: Kalman Gain is negative.\nK = \n" << K << std::endl <<
                  std::endl;
    }

    return prior + Z * M;
}

Eigen::MatrixXd EnsembleKalmanFilter(PtrList<volScalarField>& prior,
                                     Eigen::VectorXd measurements,
                                     Eigen::MatrixXd measurementsCov,
                                     Eigen::MatrixXd observedState)
{
    M_Assert(measurements.rows() == observedState.rows(),
             "The observed state should have the same dimention of the measurements");
    M_Assert(observedState.cols() == prior.size(),
             "The input matrices should all have the samples on the columns");
    M_Assert(measurementsCov.rows() == measurements.rows()
             && measurementsCov.cols() == measurements.rows(),
             "Wrong measurements covariance matrix");
    unsigned Nseeds = prior.size();
    unsigned stateDim = prior[0].size();
    unsigned measDim = measurements.size();
    Eigen::MatrixXd priorMatrix = Foam2Eigen::PtrList2Eigen(prior);
    Eigen::MatrixXd HA(measDim, Nseeds);
    Eigen::VectorXd observedStateMean = observedState.rowwise().mean();
    Eigen::MatrixXd A(stateDim, Nseeds);
    Eigen::VectorXd priorMean = priorMatrix.rowwise().mean();
    auto measNoise = std::make_shared<muq::Modeling::Gaussian>
                     (Eigen::VectorXd::Zero(measDim), measurementsCov);
    Eigen::MatrixXd D(measDim, Nseeds);

    for (unsigned i = 0; i < Nseeds; i++)
    {
        A.col(i) = priorMatrix.col(i) - priorMean;
        HA.col(i) = observedState.col(i) - observedStateMean;
        D.col(i) = measurements + measNoise->Sample();
    }

    Eigen::MatrixXd Y = D -
                        observedState; //diff measurement data and simulated data
    Eigen::MatrixXd P = HA * HA.transpose() / (Nseeds - 1.) + measurementsCov;
    Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(P);
    auto P_rank = lu_decomp.rank();

    if (P_rank < P.cols() || P_rank < P.rows())
    {
        std::cout << "Pseudo inverse of P should be implemented in the EnKF, exiting" <<
                  std::endl;
        std::cout << P << std::endl;
        exit(10);
    }
    else
    {
        P = P.inverse();
    }

    //KF update
    Eigen::MatrixXd M = P * Y;
    Eigen::MatrixXd Z = (1. / (Nseeds - 1.)) * A * HA.transpose();
    //Kalman Gain
    Eigen::MatrixXd K = (1. / (Nseeds - 1.)) * A * HA.transpose() * P;

    if ((K.array() < 0.0).any())
    {
        std::cout << "WARNING: Kalman Gain is negative.\nK = \n" << K << std::endl <<
                  std::endl;
    }

    return priorMatrix + Z * M;
}

double quantile(Eigen::VectorXd samps, double p, int method)
{
    double m;
    int n = samps.size();
    int j;
    double gamma = 0.0;
    std::sort(samps.data(), samps.data() + samps.size());

    if (method == 1)
    {
        m = 0;
        j = std::floor(p * n + m);
        double g = p * n + m - j * 1.0;

        if (g > 0)
        {
            gamma = 1.0;
        }
    }
    else if (method == 2)
    {
        gamma = 0.5;
        m = 0;
        j = std::floor(p * n + m);
        double g = p * n + m - j * 1.0;

        if (g > 0)
        {
            gamma = 1.0;
        }
    }
    else if (method == 3)
    {
        gamma = 1.0;
        m = -0.5;
        j = std::floor(p * n + m);
        double g = p * n + m - j * 1.0;

        if (g == 0 && j % 2 == 0)
        {
            gamma = 0.0;
        }
    }
    else
    {
        std::cout << "Quantile method not implemented. Exiting." << std::endl;
        exit(10);
    }

    if (j >= samps.size() - 1)
    {
        j = samps.size() - 2;
    }

    return (1 - gamma) * samps(j) + gamma * samps(j + 1);
}

Eigen::VectorXd quantile(Eigen::MatrixXd samps, double p, int method)
{
    Eigen::VectorXd output(samps.rows());

    for (int i = 0; i < samps.rows(); i++)
    {
        Eigen::VectorXd sampsRow = samps.row(i);
        output(i) = quantile(sampsRow, p, method);
    }

    return output;
}

}
}
