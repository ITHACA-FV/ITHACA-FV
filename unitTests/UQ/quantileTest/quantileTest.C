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
    quantileTest.C 
\*---------------------------------------------------------------------------*/


#include "MUQ/Modeling/Distributions/Gaussian.h"
#include "MUQ/Modeling/Distributions/Density.h"

#include <boost/math/special_functions/erf.hpp>

#include <iostream>
#include <stdio.h>
#include <Eigen/Dense>
#define _USE_MATH_DEFINES
#include <cmath>
#include "muq2ithaca.H"


int main(int argc, char* argv[])
{
    std::cout << "******************************************************" << std::endl;
    std::cout << "\nTEST of the function ITHACAmuq::muq2ithaca::quantile" << std::endl;
    std::cout << "We sample from a Gaussian distribution with mean mu = 0"<< std::endl;
    std::cout << "and variance sigma = 1.\n "<< std::endl;
    std::cout << "For this distribution, the p-quantile is given by\n\n" <<
         "mu + sigma * sqrt(2) * invErf(2 * p - 1)\n" << std::endl;
    std::cout << "Several methods to approximate the quantile have been implemented and are here tested." << std::endl;



    int nSeeds = 100000;
    auto density = std::make_shared<muq::Modeling::Gaussian>(nSeeds);
    Eigen::VectorXd samps = density->Sample();
    double p = 0.3;
    std::cout << "\nQuantile method 1 = " << ITHACAmuq::muq2ithaca::quantile(samps, p, 1) << std::endl;
    std::cout << "Quantile method 2 = " << ITHACAmuq::muq2ithaca::quantile(samps, p, 2) << std::endl;
    std::cout << "Quantile method 3 = " << ITHACAmuq::muq2ithaca::quantile(samps, p, 3) << std::endl;
    std::cout << "True quantile = " << std::sqrt(2) * boost::math::erf_inv(2 * p - 1) << std::endl;
    return 0;
}
