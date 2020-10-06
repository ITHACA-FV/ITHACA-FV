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
#ifndef FofM_H
#define FofM_H
#include <string>
#include <iostream>
#include <math.h>
#include <Eigen/Eigen>
#include <sys/stat.h>
#include "fvCFD.H"
#include "ITHACAstream.H"
class FofM
{
    public:
        //constructors
        FofM();
        FofM(label argc, char* argv[], label Nsampled);
        ~FofM() {};

        //Members:
        /// Vector to store the desired output of the model, i.e. the figure of merit
        Eigen::VectorXd modelOutput;
        /// dummy variables to transform FofM into a class
        mutable autoPtr<fvMesh> _mesh;
        autoPtr<Time> _runTime;
        /// list of pointers to scalar and vector field
        PtrList<volScalarField> ptrfield;
        PtrList<volVectorField> Vptrfield;
        /// Number of points for which simulation is run
        label Npoints;
        /// Name of the figure of merit
        std::string name;
        /// boolean variable to check if the model output is computed (true) or not (false)
        bool MObuilt = false;

        //Methods

        /// virtual method that computes the figure of merit
        /// it must be overriden in derived classes of FofM
        /// dir is the folder where the output is sought to build the desired figure of merit
        virtual void buildMO(std::string dir);

};

#endif
