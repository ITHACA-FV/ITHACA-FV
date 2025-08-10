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
    sequentialIHTP
Description
    A general full order implementation of an inverse problem
SourceFiles
    sequentialIHTP.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the sequentialIHTP class.
/// \dir
/// Directory containing the header and source files for the sequentialIHTP class.
#ifndef sequentialIHTP_H
#define sequentialIHTP_H
#include <iostream>
#include "fvCFD.H"
#include "interpolation.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include "IOmanip.H"
#include "Time.H"
#include "laplacianProblem.H"
// #include "reducedLaplacian.H"
#include "ITHACAutilities.H"
#include <Eigen/Dense>
#include <cmath>
#include "Foam2Eigen.H"
#include "cnpy.H"
#include "mixedFvPatchFields.H"
#include "cellDistFuncs.H"
#include "ITHACAregularization.H"

#define _USE_MATH_DEFINES

using namespace SPLINTER;

/// Class to implement a full order inverse parameterized laplacian problem
class sequentialIHTP: public laplacianProblem
{

    public:
        // Constructors
        sequentialIHTP();
        /// Construct with argc and argv
        sequentialIHTP(int argc, char* argv[]);
        //sequentialIHTP(int argc, char* argv[], scalar _diffusivity);
        virtual ~sequentialIHTP() {};

        ITHACAparameters* para;

        // Dummy variables to transform inverseLaplacianFoam into a class
        /// Temperature field
        autoPtr<volScalarField> _T;
        PtrList<volScalarField> Ttime;
        List<PtrList<volScalarField>> Tbasis;

        /// Additional temperature field
        autoPtr<volScalarField> _Tad;
        PtrList<volScalarField> Tad_time;
        PtrList<volScalarField> T0_time;

        /// Mesh
        mutable autoPtr<fvMesh> _mesh;

        /// simpleControl
        autoPtr<simpleControl> _simple;

        /// fvOptions
        autoPtr<fv::options> _fvOptions;

        /// Time
        autoPtr<Time> _runTime;
        autoPtr<Time> _runTimeShort;

        /// Dummy thermal diffusivity with unitary value
        dimensionedScalar DT;

        /// Diffusivity value
	scalar diffusivity = 0.0; 

        /// List of snapshots for the T0 solutions
        PtrList<volScalarField> T0field;

        /// List of POD modes
        volScalarModes T0modes;

        /// T0 reduced explicit matrix
        Eigen::MatrixXd T0explicitMatrix_red;

        /// T0 reduced implicit matrix
        Eigen::MatrixXd T0implicitMatrix_red;

        /// Projection of the Tbasis on the reduced space
        Eigen::MatrixXd Tbasis_projectionMat;

        /// Projection of Tad on the reduced space
        Eigen::VectorXd Tad_projected;

        /// Magic points for the T0 projection error estimation
        List<label> magicPoints;

        /// Matrix that reconstruct some points into the full order space
        Eigen::MatrixXd pointsReconstructMatrix;

        /// Matrix that reconstructs the Tbasis at the magic points  
        Eigen::MatrixXd pointTbasis_reconstructionMat;

        /// Reconstruction of Tsd at the magic points 
        Eigen::VectorXd pointTad_reconstructed;

        /// L2 norm of the projection error of Tbasis on the T0 modes
        PtrList<volScalarField> projectionErrorTbasis;

        /// L2 norm of the projection error of Tad on the T0 modes
        //double projectionErrorTad;
        PtrList<volScalarField> projectionErrorTad;

        /// Number of POD modes
        int NmodesT0 = 0;

        /// Time discretization (filled in the constructor)
        scalar startTime;
        scalar deltaTime;
        scalar endTime;
        label Ntimes;
        List<scalar> timeSteps;

        /// List of times at which the measurements are acquired
        /// (this List is filled by readThermocouples())
        List<scalar> samplingTime;

        /// Time interval in between samples (read from thermocouplesDict)
        scalar timeSamplesDeltaT;

        /// First sampling time (read from thermocouplesDict)
        scalar timeSamplesT0;

        /// Number of time samples (filled in createThermocouples.H)
        label timeSamplesNum;

        /// Number of timesteps between two samples
	label NtimeStepsBetweenSamples;

        /// Time sample index
	label timeSampleI;

        /// Number of sampling steps to consider when computing offline phase
        label basisDeltaSample = 2; 

        /// List of timesteps at which measurements are available
        List<label> samplingSteps; 

        /// 1 if readThermocouples() was called, 0 elsewise
        bool thermocouplesRead = 0;

        /// Number of thermocouples
        int thermocouplesNum;

        /// Number of heat flux bases (both in time and space)
        int gBasisSize;

        /// Cost funtion  [K^2]
        double J;

        /// List of cost funtions  [K^2]
        Eigen::VectorXd Jlist;

        /// Thermal conductivity [W/(m K)]
        double thermalCond;

        /// Heat transfer coefficient [W/(m2 K)]
        double HTC;

        /// Density [kg /m3]
        double density;

        /// Specific heat capacity [J/kg/K]
        double specificHeat;

        /// Folder where the offline solutions are saved
        word folderOffline = "./ITHACAoutput/offlineParamBC/";

        /// Number of timestep to solve for during offline phase
        label offlineTimestepsSize = 0;

        /// End time for the ofline computation
        scalar offlineEndTime = 0.0;

        bool offlineFlag = 0;
        bool interpolationFlag = 0;

        /// Number of timesteps considered in each acquisition sequence
        label NtimestepsInSequence = 0;

        scalar homogeneousBC = 0.0;
        List<scalar> homogeneousBCcoldSide;
        List<scalar> Tf; //temperature at coldSide [K]
        List<scalar> refGrad;
        List<scalar> valueFraction;

        /// Index of the hotSide patch
        label hotSide_ind;

        /// Index of the coldSide patch
        label coldSide_ind;

        /// Heat flux at hotSide
        List<List<scalar>> g;

        /// Heat flux space basis
        List<List<scalar>> heatFluxSpaceBasis;

        /// Bases of the heat flux
        List<List<List<scalar>>> gBaseFunctions;

        /// Number of basis
        label Nbasis = 0;

        /// Number of basis in time
        label NbasisInTime = 0;

        /// Number of basis in space
        label NbasisInSpace = 0;

        /// Number of samples considered in the offline phase
        label NsamplesWindow = 0;

        /// Weights of the parameterization
        List<scalar> gWeights;

        /// Weights of the parameterization
        List<scalar> gWeightsOld;

        /// True heat flux at hotSide [W/m2]
        List<List<scalar>> gTrue;

        /// Parametrized BC
        Eigen::VectorXd residual;
        Eigen::VectorXd addSol;
        Eigen::VectorXd T0_vector;
        Eigen::VectorXd Tcomp;
        Eigen::MatrixXd Theta;
	word timeBasisType = "None";

        word linSys_solver;
        label TSVD_filter;
        scalar Tikhonov_filter;
        label CG_Nsteps;

        List<vector> thermocouplesPos;
        List<int> thermocouplesCellID;
        List<int> thermocouplesCellProc;
        List<Foam::vector> thermocouplesCellC;
        Eigen::VectorXd Tmeas;   //Temperature at the thermocouples locations [K]
	Eigen::VectorXd TmeasShort; //Temperature at the thermocouples locations in only two sampling times [K]
        Eigen::VectorXd Tdirect;
        Eigen::VectorXd Tdiff;

        // Functions

        //--------------------------------------------------------------------------

        /// Set diffusivity
        void setDiffusivity(scalar _diff);

        //--------------------------------------------------------------------------

        /// Define the base functions used for the parametrization of g
        void setSpaceBasis(word type, scalar shapeParameter, label Npod = 0);


        //--------------------------------------------------------------------------

        /// Set parameterized heat flux defining the basis
        void set_gParametrized(word spaceBaseFuncType,
            scalar shapeParameter_space);

        //--------------------------------------------------------------------------

        /// Convert list of boundary heat flux into field
        volScalarField list2Field(List<scalar> list, scalar innerField = 0.0);

        //--------------------------------------------------------------------------

        ///
        ///
        /// @param[in] 
        ///
        List<List<scalar>> interpolateWeights(List<scalar> Wold, List<scalar> Wnew);

        //--------------------------------------------------------------------------

        /// Update the boundary condition g when g is parameterized
        ///
        /// @param[in]  weigths  New values of the weights of the basis functions
        ///
        void update_gParametrized(List<scalar> weights);

        //--------------------------------------------------------------------------

        /// Performs offline computation for the parameterized BC method, if
        /// the offline directory ""./ITHACAoutputs/offlineParamBC" exists,
        /// it reads the solution from there
        ///
        /// @param[in]  force   If 1, force the offline phase to be computed
        ///
        void parameterizedBCoffline(bool force = 0);

        //--------------------------------------------------------------------------

        /// Reconstructs the temperature field using superposition of effects
        ///
        void reconstrucT(word outputFolder);

        //--------------------------------------------------------------------------

        /// Reconstructs the temperature field using superposition of effects at 
        /// some points at the last timestep of the offline phase
        ///
        /// @param[in]  cells  Cells IDs at which the field is reconstructed
        ///
        Eigen::VectorXd reconstrucT(Eigen::VectorXi cells);

        //--------------------------------------------------------------------------

        ///
        ///
        void parameterizedBC(word folder, volScalarField initialField);

        //--------------------------------------------------------------------------

        /// Set valueFraction list values for Robin condition
        void set_valueFraction();

        //--------------------------------------------------------------------------

        /// Set BC of the direct problem
        virtual void assignDirectBC(label timeI);

        //--------------------------------------------------------------------------

        /// Set IF of the T0 problem
        ///
        virtual void assignT0_IF(volScalarField& T0_init) = 0;

        //--------------------------------------------------------------------------

        /// Solve the T0 problem
        ///
        virtual void solveT0(volScalarField initialField);

        //--------------------------------------------------------------------------

        /// Compute T0 modes prome snapshots
        ///
        void getT0modes();

        //--------------------------------------------------------------------------

        /// Project T0 matrices onto the reduced spaced
        ///
        void projectT0();

        //--------------------------------------------------------------------------

        /// Assemble the matrices to go from the gWeights to the T0 reduced space
        ///
        void projectDirectOntoT0();

        //--------------------------------------------------------------------------

        /// Assemble the matrix pointsProjectionMatrix to project some points on the 
        /// reduced basis space
        ///
        void pointProjectionOffline(); 

        //--------------------------------------------------------------------------

        /// Compute the L2 norm of the projection error for each Tbasis and Tad
        ///
        void projectionErrorOffline();

        //--------------------------------------------------------------------------

        /// Assemble all the matrices required in the online phase
        ///
        void T0offline(int NmagicPoints);


        //--------------------------------------------------------------------------

        /// Set BC and IF of the additional problem for the parameterized heat flux
        ///
        void solveAdditional();

        //--------------------------------------------------------------------------

        /// Solve direct problem
        ///
        void solveDirect();

        //--------------------------------------------------------------------------

        /// Identifies in the mesh the cells corresponding to the
        /// termocouples locations.
        ///
        virtual void readThermocouples();

        //--------------------------------------------------------------------------

        /// Interpolates the field value at the thermocouples points
        /// NOTE: do NOT call whe field is an element of a PtrList
        ///
        /// @param[in]  field    Field to read the values
        ///
        /// @param[in]  Vector of field values at thermocouples points
        ///
        Eigen::VectorXd fieldValueAtThermocouples(volScalarField& field);

        //--------------------------------------------------------------------------

        /// Interpolates the field value at the thermocouples points
        /// Always use thisone for elementes of PtrList
        ///
        /// @param[in]  field    Field to read the values
        ///
        /// @param[in]  Vector of field values at thermocouples points
        ///
        Eigen::VectorXd fieldValueAtThermocouples(PtrList<volScalarField> fieldList,
                label fieldI);

        //--------------------------------------------------------------------------

        /// Interpolates the field value at the thermocouples points
        ///
        /// @param[in]  field    Field to read the values
        ///
        /// @param[in]  Vector of field values at thermocouples points
        ///
        Eigen::VectorXd fieldValueAtThermocouples(PtrList<volScalarField> fieldList);

        //--------------------------------------------------------------------------

        /// Restart temperature field
        ///
        void restart();

        //--------------------------------------------------------------------------

        /// Restart fields
        ///
        void restartOffline();

        //--------------------------------------------------------------------------

        /// Restart fields
        ///
        void restartT0();

        //--------------------------------------------------------------------------

        /// Fills the vector samplingSteps which contains the timesteps at which the measurements are taken
        ///
        void sampling2symulationTime();

        //--------------------------------------------------------------------------

        ///
        ///
        void parameterizedBC_postProcess(List<Eigen::MatrixXd> linSys,
                                         Eigen::VectorXd weigths, word outputFolder, label verbose = 0);

        //--------------------------------------------------------------------------

        /// Find the points at with the projection error is computed
        ///
        void findMagicPoints(int NmagicPoints);

};

#endif
