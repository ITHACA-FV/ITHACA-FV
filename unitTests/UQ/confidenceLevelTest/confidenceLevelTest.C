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
Description
    Test of the implementation of the confidence level function
SourceFiles
    confidenceLevelTest.C 
\*---------------------------------------------------------------------------*/


#include "MUQ/Modeling/Distributions/Gaussian.h"
#include "MUQ/Modeling/Distributions/DensityProduct.h"
#include "MUQ/Modeling/Distributions/InverseGamma.h"
#include "MUQ/Modeling/Distributions/Density.h"

#include "MUQ/Utilities/AnyHelpers.h"
#include "MUQ/Utilities/RandomGenerator.h"
#include "MUQ/Utilities/StringUtilities.h"

#include  <boost/property_tree/ptree.hpp>

#include <iostream>
#include "interpolation.H"
#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include "IOmanip.H"
#include "Time.H"
#include "laplacianProblem.H"
#include "inverseLaplacianProblem.H"
#include "reducedInverseLaplacian.H"
#include "ITHACAPOD.H"
#include "ITHACAutilities.H"
#include <Eigen/Dense>
#define _USE_MATH_DEFINES
#include <cmath>
#include "Foam2Eigen.H"
#include "muq2ithaca.H"
#include "mixedFvPatchFields.H"
#include "cellDistFuncs.H"


int main(int argc, char* argv[])
{
    word outputFolder = "./ITHACAoutput/";
    int Nseeds = 1000;

    unsigned dim = 1;
    double stdDev = 0.2;
    double mean = 0.0;
    auto dist = std::make_shared<muq::Modeling::Gaussian>(dim); 
    Eigen::VectorXd samps(Nseeds);


    for(int i = 0; i < Nseeds; i++)
    {
        samps(i) = dist->Sample()(0) * stdDev;
    }

    Eigen::VectorXd sampsNoMean = samps - Eigen::VectorXd::Ones(samps.size()) * samps.mean();
    double samps_stdDev = std::sqrt(sampsNoMean.dot(sampsNoMean) / (samps.size() - 1));


    auto [minConf, maxConf] = ITHACAmuq::muq2ithaca::confidenceLevel(samps, 0.99);
    std::cout << "Mean = " <<  mean << "   stdDev = " << stdDev << std::endl;
    std::cout << "Samples mean = " << samps.mean() << "\nSamples stdDev = " << samps_stdDev << std::endl;
    std::cout << "\nConfidence levels = " << minConf << " , " << maxConf << std::endl ;

    return 0;
}
