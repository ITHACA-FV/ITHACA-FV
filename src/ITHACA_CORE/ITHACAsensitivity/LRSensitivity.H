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
#ifndef LRSensitivity_H
#define LRSensitivity_H
#include <iostream>
#include <vector>
#include <string>
#include <Eigen/Eigen>
#include "FofM.H"
#include "ITHACAstream.H"
#include "ITHACAutilities.H"
#include "ITHACAsampling.H"


class LRSensitivity
{
    public:
        //Constructors
        LRSensitivity();
        explicit LRSensitivity(int Npara, int Np);
        ~LRSensitivity() {};

        //Members
        /// Matrix to store the range used in training/offline stage
        Eigen::MatrixXd trainingRange;
        /// Matrices storing the independent variables' values
        Eigen::MatrixXd MatX;
        /// Matrices storing the independent variables' normalized values
        Eigen::MatrixXd MatXn;
        /// Figure of merit object
        FofM M;
        /// Number of parameters used for the analysis
        label No_parameters;
        /// Number of points used for Monte Carlo sampling
        label Npoints;
        /// Vectors to store the output of the model
        Eigen::VectorXd y;
        ///Vector to store the normalized output of the model
        Eigen::VectorXd yn;
        ///Vector to store the result of linear approximation
        Eigen::VectorXd ylin;
        /// Mean value of the output
        double Ey;
        /// Variance of the output
        double Vy;
        /// Vector to store the meanvalues of the parameters
        Eigen::VectorXd EX;
        /// Vector to store the variances of the parameters
        Eigen::VectorXd VX;
        /// Vector to store the standardized regression coefficient (SRC)
        Eigen::VectorXd betas;
        /// double to quantify the goodness of linear regression from 0 (no linearity) to 1 (perfect linearity)
        /// also known as R^2
        double QI;
        /// boolean variable to check if EX and VX are computed
        bool Xdone = false;
        /// boolean variable to check if Ey and Vy are computed
        bool Ydone = false;
        /// boolean variable to check if  SRCs are computed
        bool bdone = false;

        //--------------------------------------------------------------------------------------------------------------------//
        //Methods

        /// Method to build MatX, it internally calls ITHACAsampling::samplingMC
        /// pdflist is the list with the names of the distributions of the parameters
        /// plist contains the parameters of the distributions, therefore its size must be (No_parameters,2)
        void buildSamplingSet(std::vector<std::string>& pdflist, Eigen::MatrixXd plist);
        ///----------------------------------------------------------------------------------------

        /// Method to load the output of FOM/ROM inside member LRsensitivity member y
        void load_output();
        ///----------------------------------------------------------------------------------------

        /// Method to compute EX and VX, it sets Xdone=true
        void getXstats();
        ///----------------------------------------------------------------------------------------

        /// Method to compute Ey and Vy, it sets Ydone=true
        void getYstat();
        ///----------------------------------------------------------------------------------------

        /// Method to compute SRCs, it sets bdone=true
        void getBetas();
        ///----------------------------------------------------------------------------------------

        /// Method to compute the quality of linear regression approximation
        void assessQuality();
        ///-----------------------------------------------------------------------------------------





    private:

        /// method used in the LRsensitivity object's constructor to set the size of all matrices and vectors of the class
        void setAll();

};

#endif
