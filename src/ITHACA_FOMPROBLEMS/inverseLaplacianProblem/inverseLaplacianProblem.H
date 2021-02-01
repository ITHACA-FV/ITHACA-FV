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
    inverseLaplacianProblem
Description
    A general full order implementation of an inverse boundary condition estimation
    problem
SourceFiles
    inverseLaplacianProblem.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the inverseLaplacianProblem class.
/// \dir
/// Directory containing the header and source files for the inverseLaplacianProblem class.


#ifndef inverseLaplacianProblem_H
#define inverseLaplacianProblem_H
#include <iostream>
#include "fvCFD.H"
#include "interpolation.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include "IOmanip.H"
#include "Time.H"
#include "laplacianProblem.H"
#include "ITHACAPOD.H"
#include "ITHACAutilities.H"
#include "ITHACAregularization.H"
#include <Eigen/Dense>
#include <cmath>
#include "Foam2Eigen.H"
#include "cnpy.H"
#include "mixedFvPatchFields.H"
#include "cellDistFuncs.H"
#include "reductionProblem.H"
#include "Modes.H"
#define _USE_MATH_DEFINES


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                 Class inverseLaplacianProblem Declaration
\*---------------------------------------------------------------------------*/

/// Implementation of a <b> inverse Laplacian problem </b>.
/** In this class, a general inverse Laplacian problem is set up.
In this problem, the objective is to estimate the boundary heat flux at the
hotSide patch given pointwise temeprature measurements in the interior of the domain.
The different solution methodologies are implemented as child of this class */
class inverseLaplacianProblem: public laplacianProblem
{

    public:
        // Constructors
        /// Null constructor
        inverseLaplacianProblem();
        /// Construct with argc and argv. Reads the thermocouples dictionary.
        inverseLaplacianProblem(int argc, char* argv[]);
        /// Destructor
        virtual ~inverseLaplacianProblem() {};

        ITHACAparameters* para;

        /// Temperature field
        autoPtr<volScalarField> _T;

        /// Mesh
        mutable autoPtr<fvMesh> _mesh;

        /// simpleControl
        autoPtr<simpleControl> _simple;

        /// fvOptions
        autoPtr<fv::options> _fvOptions;

        /// Time
        autoPtr<Time> _runTime;

        /// Dummy thermal conductivity with unitary value
        dimensionedScalar DT;

        /// Thermal diffusivity [W/(m K)]
        double k;

        /// Heat transfer coefficient [W/(m2 K)]
        double H;

        /// Flag to know if thermocouples file was read
        bool thermocouplesRead = 0;

        /// Number of thermocouples
        int thermocouplesNum;

        /// Cost funtion  [K^2]
        double J;
        double L2norm;
        double LinfNorm;

        /// Homogenerous BC
        scalar homogeneousBC = 0.0;

        /// List of zeros of the size of coldSide patch
        List<scalar> homogeneousBCcoldSide;

        /// Temperature at coldSide [K]
        List<scalar> Tf;

        /// Reference gradient for the Robin BC
        List<scalar> refGrad;

        /// Value fraction for the Robin BC
        List<scalar> valueFraction;

        /// Index of the hotSide patch
        label hotSide_ind;

        /// Index of the coldSide patch
        label coldSide_ind;

        /// Heat flux at hotSide [W/m2]
        List<scalar> g;

        /// List of boundary heat fluxes
        List<List<scalar>> gList;

        /// True heat flux at hotSide [w/m2]
        List<scalar> gTrue;

        /// List of patch faces areas [m2]
        List<scalar> faceCellArea;

        /// List containing the positions of the thermocouples
        List<vector> thermocouplesPos;

        /// List of cells indices containing a thermocouple
        List<int> thermocouplesCellID;

        /// List of incedes of the processors contining each thermocouple
        List<int> thermocouplesCellProc;

        /// Vector of measured temperatures at the thermocouples locations [K]
        Eigen::VectorXd Tmeas;

        /// Vector of computed temperatures at the thermocouples locations [K]
        Eigen::VectorXd Tdirect;

        /// Difference between computed and measured temperatures at the thermocouples
        Eigen::VectorXd Tdiff;

        /// Number of processors
        label nProcs;


        // Functions

        //--------------------------------------------------------------------------
        /// Set the right g size and fills it with zeros
        ///
        void set_g();

        //--------------------------------------------------------------------------
        /// Create a field with the hotSide boundary heat flux at the hotSide
        /// bounday cells for visualization
        ///
        /// @param[in]  list        List of boundary heat flux
        /// @param[in]  innerField  Value of the inner field
        ///
        /// @return     Field with the heat flux value at the boundary cells
        ///
        volScalarField list2Field(List<scalar> list, scalar innerField = 0.0);

        //--------------------------------------------------------------------------
        /// Set valueFraction list values for Robin condition
        ///
        void set_valueFraction();

        //--------------------------------------------------------------------------
        /// Solve the direct problem with the true heat flux as boundary condition
        /// and fills the measured temeprature vector
        ///
        virtual void solveTrue() {};

        //--------------------------------------------------------------------------
        /// Set boundary condition of the direct problem
        virtual void assignDirectBC();

        //--------------------------------------------------------------------------
        /// Solve direct problem
        ///
        void solveDirect();

        //--------------------------------------------------------------------------
        /// Solve Laplacian problem without source term
        ///
        void solve(const char* problem);

        //--------------------------------------------------------------------------
        /// Identifies in the mesh the cells corresponding to the
        /// termocouples locations
        ///
        virtual void readThermocouples();

        //--------------------------------------------------------------------------
        /// Interpolates the field value at the thermocouples points
        ///
        /// @param[in]  field    Field to read the values
        ///
        /// @return  Vector of field values at thermocouples points
        ///
        Eigen::VectorXd fieldValueAtThermocouples(volScalarField& field);

        //--------------------------------------------------------------------------
        /// Computes the difference between direct problem solution and measures
        /// Saves the difference vector in Tdiff
        ///
        void differenceBetweenDirectAndMeasure();

        //--------------------------------------------------------------------------
        /// Compute maximum cell dimension in x, y and z
        ///
        /// @param[in]  ff
        /// @param[in]  pp
        /// @param[in]  cc
        /// @param[in]  pLabels
        /// @param[in]  pLocal
        ///
        /// @param[out] dim Maximum cell dimensions
        ///
        Foam::vector cellDim(const faceList& ff, const pointField& pp,
                             const cell& cc, labelList pLabels, pointField pLocal);

        //--------------------------------------------------------------------------
        /// Restart fields
        ///
        void restart();
};

#endif
