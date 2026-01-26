/*---------------------------------------------------------------------------*\
     ██╗████████╗██╗  ██╗ █████╗  ██████╗ █████╗       ███████╗██╗   ██╗
     ██║╚══██╔══╝██║  ██║██╔══██╗██╔════╝██╔══██╗      ██╔════╝██║   ██║
     ██║   ██║   ███████║███████║██║     ███████║█████╗█████╗  ██║   ██║
     ██║   ██║   ██╔══██║██╔══██║██║     ██╔══██║╚════╝██╔══╝  ╚██╗ ██╔╝
     ██║   ██║   ██║  ██║██║  ██║╚██████╗██║  ██║      ██║      ╚████╔╝
     ╚═╝   ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝      ╚═╝       ╚═══╝

 * In real Time Highly Advanced Computational Applications for Finite Volumes
 * Copyright (C) 2017 by the ITHACA-FV authors
-------------------------------------------------------------------------------

License
    This file is part of ITHACA-FV

    ITHACA-FV is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    ITHACA-FV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with ITHACA-FV. If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "ithacaInterpolator.H"
#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <limits>
#include <iomanip>
#include "dictionary.H"
#include "IFstream.H"

// Function to approximate: z = sin(x) * cos(y)
double targetFunction(double x, double y)
{
    return std::sin(x) * std::cos(y);
}

// Function to create a dictionary for testing different configurations
Foam::dictionary createTestDict(const Foam::word& package, const Foam::word& algorithm, const Foam::word& kernel = "")
{
    Foam::dictionary dict;
    dict.add("package", package);
    dict.add("algorithm", algorithm);

    if (package == "mathtoolbox")
    {
        if (algorithm == "GPR")
        {
            dict.add("kernel", kernel);
            dict.add("normalize", true);
            dict.add("optimizeHyperparams", true);
            dict.add("kernelScale", 0.10);
            dict.add("lengthScale", 0.10);
            dict.add("noise", 1e-5);
        }
        else if (algorithm == "RBF")
        {
            dict.add("kernel", kernel);
            dict.add("regularization", false);
            dict.add("polynomial", false);
            dict.add("lambda", 1e-6);
            dict.add("epsilon", 0.5);
            dict.add("normalize", true);
        }
    }
    else if (package == "splinter")
    {
        dict.add("kernel", kernel);
        dict.add("epsilon", 0.5);
    }

    return dict;
}

// Function to evaluate model performance
double evaluateModel(ithacaInterpolator& model, int nTest, double minVal, double maxVal, double step)
{
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

            double predicted = model.predict(query);
            double exact = targetFunction(x, y);
            double error = std::abs(predicted - exact);

            errorSum += error;
            testCount++;
        }
    }

    return errorSum / testCount;
}

int main(int argc, char **argv)
{
    std::cout << "Comprehensive test of RBF and GPR kernels with different packages..." << std::endl;

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

    // 2. Define test configurations
    struct TestConfig {
        Foam::word package;
        Foam::word algorithm;
        Foam::word kernel;
        std::string description;
    };

    std::vector<TestConfig> configs = {
        // Mathtoolbox GPR kernels
        {"mathtoolbox", "GPR", "matern", "Mathtoolbox GPR with Matern kernel"},
        {"mathtoolbox", "GPR", "squared_exp", "Mathtoolbox GPR with Squared Exponential kernel"},

        // Mathtoolbox RBF kernels
        {"mathtoolbox", "RBF", "linear", "Mathtoolbox RBF with Linear kernel"},
        {"mathtoolbox", "RBF", "gaussian", "Mathtoolbox RBF with Gaussian kernel"},
        {"mathtoolbox", "RBF", "multiquadric", "Mathtoolbox RBF with Multiquadric kernel"},
        {"mathtoolbox", "RBF", "cubic", "Mathtoolbox RBF with Cubic kernel"},
        {"mathtoolbox", "RBF", "quintic", "Mathtoolbox RBF with Quintic kernel"},
        {"mathtoolbox", "RBF", "inverse_quadratic", "Mathtoolbox RBF with Inverse Quadratic kernel"},
        {"mathtoolbox", "RBF", "inverse_multiquadric", "Mathtoolbox RBF with Inverse Multiquadric kernel"},
        {"mathtoolbox", "RBF", "thin_plate_spline", "Mathtoolbox RBF with Thin Plate Spline kernel"},

        // Splinter RBF kernels
        {"splinter", "RBF", "linear", "Splinter RBF with Linear kernel"},
        {"splinter", "RBF", "gaussian", "Splinter RBF with Gaussian kernel"},
        {"splinter", "RBF", "multiquadric", "Splinter RBF with Multiquadric kernel"},
        {"splinter", "RBF", "cubic", "Splinter RBF with Cubic kernel"},
        {"splinter", "RBF", "quintic", "Splinter RBF with Quintic kernel"},
        {"splinter", "RBF", "inverse_quadratic", "Splinter RBF with Inverse Quadratic kernel"},
        {"splinter", "RBF", "inverse_multiquadric", "Splinter RBF with Inverse Multiquadric kernel"},
        {"splinter", "RBF", "thin_plate_spline", "Splinter RBF with Thin Plate Spline kernel"}
    };

    // 3. Test each configuration
    std::vector<double> mae_values;
    std::vector<std::string> descriptions;

    for (const auto& config : configs)
    {
        std::cout << "\n" << std::string(60, '=') << std::endl;
        std::cout << "Testing: " << config.description << std::endl;
        std::cout << std::string(60, '=') << std::endl;

        try
        {
            // Create dictionary for this configuration
            Foam::dictionary dict = createTestDict(config.package, config.algorithm, config.kernel);

            // Initialize interpolator
            ithacaInterpolator model(dict);
            model.printInfo();

            // Fit the model
            std::cout << "Fitting the model..." << std::endl;
            model.fit(X_train, y_train);

            // Evaluate performance
            int nTest = 10;
            double mae = evaluateModel(model, nTest, minVal, maxVal, step);
            mae_values.push_back(mae);
            descriptions.push_back(config.description);

            std::cout << "Mean Absolute Error: " << mae << std::endl;
            std::cout << "Test PASSED" << std::endl;
        }
        catch (const std::exception& e)
        {
            std::cout << "Test FAILED with exception: " << e.what() << std::endl;
            mae_values.push_back(std::numeric_limits<double>::max());
            descriptions.push_back(config.description + " (FAILED)");
        }
    }

    // 4. Summary and comparison
    std::cout << "\n" << std::string(80, '=') << std::endl;
    std::cout << "SUMMARY OF RESULTS" << std::endl;
    std::cout << std::string(80, '=') << std::endl;

    std::cout << std::left << std::setw(50) << "Configuration" << std::setw(15) << "MAE" << std::endl;
    std::cout << std::string(65, '-') << std::endl;

    for (size_t i = 0; i < descriptions.size(); ++i)
    {
        std::cout << std::left << std::setw(50) << descriptions[i]
                  << std::setw(15) << (mae_values[i] == std::numeric_limits<double>::max() ? "FAILED" : std::to_string(mae_values[i])) << std::endl;
    }

    // Find best performing configuration (excluding failed ones)
    double best_mae = std::numeric_limits<double>::max();
    size_t best_idx = 0;
    for (size_t i = 0; i < mae_values.size(); ++i)
    {
        if (mae_values[i] < best_mae && mae_values[i] != std::numeric_limits<double>::max())
        {
            best_mae = mae_values[i];
            best_idx = i;
        }
    }

    if (best_mae != std::numeric_limits<double>::max())
    {
        std::cout << "\nBest performing configuration:" << std::endl;
        std::cout << descriptions[best_idx] << " with MAE = " << best_mae << std::endl;
    }

    return 0;
}
