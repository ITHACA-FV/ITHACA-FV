#include "RBFinterpolator.H"
#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include "dictionary.H"
#include "IFstream.H"

// Function to approximate: z = sin(x) * cos(y)
double targetFunction(double x, double y)
{
    return std::sin(x) * std::cos(y);
}

int main(int argc, char **argv)
{
    // 1. Generate Training Data
    int nSamples = 20;
    Eigen::MatrixXd X_train(2, nSamples * nSamples);
    Eigen::VectorXd y_train(nSamples * nSamples);

    double minVal = -3.14;
    double maxVal = 3.14;
    double step = (maxVal - minVal) / (nSamples - 1);

    int idx = 0;
    for (int i = 0; i < nSamples; ++i)
    {
        for (int j = 0; j < nSamples; ++j)
        {
            double x = minVal + i * step;
            double y = minVal + j * step;
            X_train(0, idx) = x;
            X_train(1, idx) = y;
            y_train(idx) = targetFunction(x, y);
            idx++;
        }
    }

    // 2. Initialize RBFinterpolator
    Foam::dictionary dict(Foam::IFstream("ITHACAdict")());
    Foam::dictionary viscDict = dict.subDict("viscRBFdict");
    RBFinterpolator rbfModel(viscDict);
    Foam::word package = viscDict.lookupOrDefault<Foam::word>("package", "mathtoolbox");
    std::cout << "Initializing RBFinterpolator with " << package << " package..." << std::endl;
    rbfModel.printInfo();

    // 3. Fit the model
    std::cout << "Fitting the model..." << std::endl;
    rbfModel.fit(X_train, y_train);

    // 4. Test the model
    std::cout << "Testing the model..." << std::endl;
    int nTest = 10;
    double testStep = (maxVal - minVal) / (nTest - 1);
    double errorSum = 0.0;
    int testCount = 0;

    for (int i = 0; i < nTest; ++i)
    {
        for (int j = 0; j < nTest; ++j)
        {
            // Shift slightly to avoid exact training points
            double x = minVal + i * testStep + step/1.9;
            double y = minVal + j * testStep + step/1.9;

            if (x > maxVal || y > maxVal) continue;

            Eigen::VectorXd query(2);
            query(0) = x;
            query(1) = y;

            double predicted = rbfModel.predict(query);
            double exact = targetFunction(x, y);
            double error = std::abs(predicted - exact);

            errorSum += error;
            testCount++;
        }
    }

    double mae = errorSum / testCount;
    std::cout << "Mean Absolute Error: " << mae << std::endl;

    bool test1Passed = (mae < 0.2);

    if (test1Passed)
    {
        std::cout << "Test 1 PASSED" << std::endl;
    }
    else
    {
        std::cout << "Test 1 FAILED" << std::endl;
    }
}
