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
#ifndef ITHACAsampling_H
#define ITHACAsampling_H
#include <iostream>
#include <ctime>
#include <random>
#include <vector>
#include <string>
#include <Eigen/Eigen>
#include <sys/stat.h>
#include "fvCFD.H"
#include "ITHACAstream.H"

class ITHACAsampling
{
    public:

        //Methods
        ///-------------------------------------------------------------------------------------------
        /// method to sample with MC inverse transform according to the parameter distribution
        /// pdftype can be only one of the list present in the private members of the class.
        /// lowerE and upperE are the extremes of the desired range, values outside of it are rejected
        /// distpara1 and distpara2 are the parameters of the distibution, in case of distribution with
        /// only 1 parameter only the value in distpara1 is used
        static Eigen::VectorXd samplingMC(std::string pdftype, double& lowerE,
                                          double& upperE, double& distpara1, double& distpara2, label& Npoints);

    private:

        /// vector containing the name of the distribution managed in this class
        static std::vector<std::string> distributions;

};

#endif
