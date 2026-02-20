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
#include "ITHACAstream.H"
#include <Eigen/Eigen>
#include "TurbDiffusionHyperreduction.H"


class tutorial27_online
{
    public:

        StoredParameters* m_parameters;
        UnsteadyNSTurb* m_UnsteadyNSTurb;

        Foam::Time runTime2;

        volVectorField* meanU;
        label nModesU;
        PtrList<volVectorField> spatialModesU;
        Eigen::MatrixXd massMatrixInv;

        Eigen::MatrixXd temporalModesUSimulation;

        word interpolatedField;
        label nModesHR;
        label nMagicPoints;
        List<int> magicPoints;
        labelList localMagicPoints;
        Eigen::VectorXd weightAtMg;
        std::vector<volTensorField> defTensorOfModesAtMg;

        int d = 0;
        float Ck;
        float Ce;

        volScalarField* nut0=nullptr;
        volScalarField* delta=nullptr;
        volScalarField* aaa=nullptr;
        volScalarField* inv2aaa=nullptr;
        volTensorField* fullDefField=nullptr;
        volScalarField* bbb=nullptr;
        volScalarField* ccc=nullptr;
        volScalarField* sqrtk=nullptr;
        volScalarField* nut=nullptr;
        volVectorField* stressField=nullptr;
        volVectorField* meanSmagOnMagicNeighborhoods=nullptr;
        volScalarField* meanNutOnMagicPoints=nullptr;


    tutorial27_online(StoredParameters* parameters);

    // Evaluate nut at the magic points
    void evaluateApproxNut(const Eigen::VectorXd& reducedCoeffs);

    // Compute Smagorinsky term at the magic points
    Eigen::VectorXd computeApproxSmagMg(const Eigen::VectorXd& reducedCoeffs);

    // Reshape output for nut at the magic points
    Eigen::VectorXd computeApproxNutMg(const Eigen::VectorXd& reducedCoeffs);

    // Compute the hyperreduction prediction for Smagorinsky
    Eigen::VectorXd predictSmagROMCoeffs(const Eigen::VectorXd& reducedCoeffs);

    // Compute a vector field onto velocity modes from its reduced temporal modes. Homogenize to a stress if specified.
    volVectorField computeROMproj_fromCoeffs(const Eigen::VectorXd& reducedCoeffs, bool stressUnit=false);

    // Compute the coefficients of a vector field in the velocity basis (V^T f)
    Eigen::VectorXd computeROMcoeffs_fromFullDim(volVectorField& f_full);

    // Compute the reference Smagorinsky term (from reduced velocity)
    volVectorField computeSmagTerm_fromChronos(const Eigen::VectorXd& reducedCoeffs);

    // Perfom the prediction on the test time steps and compare the results to the reference
    void prediction();

};
