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
    Example of the hyperreduction of the Smagorinsky term in a ROM
SourceFiles
    27SmagorinskyHyperreduction.C
\*---------------------------------------------------------------------------*/

#include "ITHACAutilities.H"
#include <Eigen/Eigen>
#include "TurbDiffusionHyperreduction.H"
#include "hyperReduction.templates.H"


class tutorial27_offline
{
    public:

            StoredParameters* m_parameters;
            UnsteadyNSTurb* m_UnsteadyNSTurb;

            volVectorField* meanU;
            label nModesU;
            PtrList<volVectorField> spatialModesU;

            Eigen::MatrixXd covMatrixU;
            Eigen::MatrixXd temporalModesU;
            Eigen::MatrixXd temporalModesUSimulation;

            word interpolatedField;
            label nModesHR;
            label nMagicPoints;
            List<int> magicPoints;


    tutorial27_offline(int argc, char* argv[]);

    // Compute POD on velocity and on the hyperreduced term. Then select the magic points
    void decompose();

    // Compute the matrix for fast online prediction of the nonpolynomial term 
    void project();

};
