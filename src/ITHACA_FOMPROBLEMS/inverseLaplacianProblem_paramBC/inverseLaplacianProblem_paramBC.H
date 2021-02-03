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
Class
    inverseLaplacianProblem_paramBC
Description
    Implementation of the parameterization method for the solution of the
    inverse problem of estimating the boundary heat flux, given
    pointwise temperature measurements, in a Laplacian problem
SourceFiles
    inverseLaplacianProblem_paramBC.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the inverseLaplacianProblem_paramBC class.
/// \dir
/// Directory containing the header and source files for the inverseLaplacianProblem_paramBC class.

#ifndef inverseLaplacianProblem_paramBC_H
#define inverseLaplacianProblem_paramBC_H

#include "inverseLaplacianProblem.H"

#define _USE_MATH_DEFINES

struct metaData_offline
{
    int numberTC;
    int numberBasis;
    word basisType;
    double shapeParameter;


    metaData_offline() {}

    metaData_offline(int _numberTC, int _numberBasis, word _basisType,
                     double _shapeParameter)
        : numberTC(_numberTC), numberBasis(_numberBasis),
          basisType(_basisType), shapeParameter(_shapeParameter) {}
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                 Class inverseLaplacianProblem_paramBC Declaration
\*---------------------------------------------------------------------------*/

/// Implementation of a <b> parameterization method </b>to solve inverse Laplacian problems.
/** In this class, we implement the parameterization method for solving
the inverse problem of estimating the boundary heat flux, given pointwise
temperature measurements, in a Laplacian problem*/
class inverseLaplacianProblem_paramBC: public inverseLaplacianProblem
{

    public:
        // Constructors
        /// Null constructor
        inverseLaplacianProblem_paramBC();
        /// Construct with argc and argv. Reads the thermocouples dictionary
        inverseLaplacianProblem_paramBC(int argc, char* argv[]);
        /// Destructor
        virtual ~inverseLaplacianProblem_paramBC() {};

        /// Folder where the offline solutions and matrices are saved
        word folderOffline = "./ITHACAoutput/offlineParamBC/";

        /// Solutions of the direct problem with the basis of the parameterization as boundary heat flux
        PtrList<volScalarField> Tbasis;


        /// Solution of the additional problem
        PtrList<volScalarField> Tad_base;

        /// Residual of the linear system
        Eigen::VectorXd residual;

        /// Theta matrix of the lynear system
        Eigen::MatrixXd Theta;

        /// Modes of the POD
        Eigen::MatrixXd gPODmodes;

        /// Solution of the additional problem at the thermocouples positions
        Eigen::VectorXd addSol;

        /// Flag to know if the offline phase was computed
        bool offlineReady = 0;

        /// Basis of the heat flux parameterization
        List<List<scalar>> gBaseFunctions;

        /// Weights of the heat flux parameterization
        List<scalar> gWeights;

        /// RBF shape parameter
        double shapeParameter;

        /// Type of basis functions used for the parameterization of the heat flux
        word baseFuncType;


        // Functions

        //--------------------------------------------------------------------------
        /// Define the base functions used for the parametrization of the heat flux g
        /// The center of each base function is the projection of each thermocouple
        /// on the boundary hotSide.
        /// The type of basis function is defined by the value of baseFuncType
        /// (rbf = radial basis functions, pod = proper orthogonal decomposition).
        ///
        virtual void set_gBaseFunctions();

        //--------------------------------------------------------------------------
        /// Performs POD on the RBF basis functions
        ///
        /// @param[in]  Nmodes  Number of modes used for the POD
        ///
        void set_gBaseFunctionsPOD(label Nmodes);

        //--------------------------------------------------------------------------
        /// Set initial heat flux for the parameterized BC method
        ///
        /// @param[in]  baseFuncType        Type of basis function (rbf or pod)
        /// @param[in]  _shapeParameter     Value of the RBF shape parameter
        ///
        void set_gParametrized(word baseFuncType, scalar _shapeParameter = 1);

        //--------------------------------------------------------------------------
        /// Update the boundary heat flux
        ///
        /// @param[in]  weigths  New values of the weights of the basis functions
        ///
        void update_gParametrized(List<scalar> weigths);

        //--------------------------------------------------------------------------
        /// Performs offline computation for the parameterized BC method, if
        /// the offline directory "./ITHACAoutputs/offlineParamBC" exists,
        /// it reads the solution from there
        ///
        /// @param[in]  force   If 1, force the offline phase to be computed
        ///
        void parameterizedBCoffline(bool force = 0);

        //--------------------------------------------------------------------------
        /// Solve the online phase
        ///
        /// @param[in]  linSys_solver   Type of linear solver are fullPivLU, jacobiSvd, householderQr, ldlt, inverse, TSVD, Tikhonov
        /// @param[in]  regPar          Value of the regularization parameter
        ///
        Eigen::VectorXd parameterizedBC(word linSys_solver = "fullPivLU",
                                        double regPar = 0);

        //--------------------------------------------------------------------------
        /// Reconstruct the temperature field and compute the cost function J
        ///
        /// @param[in]  weigths Weights of the parameterization
        ///
        void parameterizedBCpostProcess(Eigen::VectorXd weigths);

        //--------------------------------------------------------------------------
        /// Set BC the additional problem and solves it
        ///
        virtual void solveAdditional();

        //--------------------------------------------------------------------------
        /// Reconstuct the field T using the offline computed fields
        ///
        void reconstructT();
};

#endif
