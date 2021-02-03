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
    reductionProblem
Description
    A general reduction problem class is implemented
SourceFiles
    reductionProblem.C
    reductionProblemTemplates.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the reductionProblem class.
/// \dir
/// Directory containing the header, source and template files for the reductionProblem class.


#ifndef reductionProblem_H
#define reductionProblem_H

#include <random>
#include "fvCFD.H"
#include "IOmanip.H"
#include "fixedFluxPressureFvPatchScalarField.H"
#include "freestreamFvPatchField.H"
#include <sys/stat.h>
#include "ITHACAutilities.H"
#include "ITHACAparallel.H"
#include "ITHACAstream.H"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wnon-virtual-dtor"
#pragma GCC diagnostic ignored "-Wreturn-type"
#include <datatable.h>
#include <bspline.h>
#include <bsplinebuilder.h>
#include <rbfspline.h>
#pragma GCC diagnostic pop

// #include <spline.h>


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


/*---------------------------------------------------------------------------*\
                        Class reductionProblem Declaration
\*---------------------------------------------------------------------------*/

/// A general class for the implementation of a full order parametrized problem
class reductionProblem
{
    private:

    public:
        // Constructors

        /// Construct Null
        reductionProblem();
        ~reductionProblem() {};

        // Members
        /// Number of parameters
        label Pnumber;
        /// Dimension of the training set (used only when gerating parameters without input)
        label Tnumber;
        /// Row matrix of parameters
        Eigen::MatrixXd mu;
        /// Range of the parameter spaces
        Eigen::MatrixXd mu_range;
        /// Matrix of parameters to be used for PODI, where each row corresponds to a sample point. In this matrix the time dimension is regarded as a parameter for unsteady problems
        Eigen::MatrixXd mu_samples;
        /// Current value of the parameter
        double mu_cur;
        /// Boolean variable, it is 1 if the POD has already been computed, else 0
        bool podex;
        /// Boolean variable, it is 1 if the Offline phase has already been computed, else 0
        bool offline;
        /// dictionary to store input output infos
        IOdictionary* ITHACAdict;
        /// argList
        autoPtr<argList> _args;
        /// parallel handling
        ITHACAparallel* paral;

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcomment"
        /// Matrix that contains informations about the inlet boundaries
        /// @details The dimension is: <br>
        /// Rows = Number of parametrized boundary conditions <br>
        /// Cols = 2 <br>
        /// Example: <br>
        /// example.inletIndex.resize(1, 2); <br>
        /// example.inletIndex(0, 0) = 0; <br>
        /// example.inletIndex(0, 1) = 0; <br>
        /// Means that there is one parametrized boundary conditions, the first col contains the
        /// index of the parametrized boundary, the second col contains the direction along which
        /// the BC is parametrized: <br>
        /// 0 for x <br>
        /// 1 for y <br>
        /// 2 for z <br>
        /// The Matrix should be implemented in the following way:
        /// \\[
        /// \\mbox{inletIndex}=\\begin{bmatrix}
        /// \\mbox{BC1_index} & \\mbox{Direction of BC1} \\\\ \
        /// \\mbox{BC2_index} & \\mbox{Direction of BC3} \\\\ \
        /// \\vdots & \\vdots \\\\ \
        /// \\mbox{BCN_index} & \\mbox{Direction of BCN} \\ \
        /// \\end{bmatrix}
        /// \\]
        Eigen::MatrixXi inletIndex;

        /// Matrix that contains informations about the inlet boundaries
        /// without specifing the direction
        /// Rows = Number of parametrized boundary conditions <br>
        /// Cols = 1 <br>
        /// Example: <br>
        /// example.inletIndex.resize(2, 1); <br>
        /// example.inletIndex(0, 0) = 0; <br>
        /// example.inletIndex(1, 0) = 1; <br>
        /// Means that there are two parametrized boundary conditions of which the
        /// first row is of patch 0 and the second row of patch 1.
        /// The Matrix should be implemented in the following way:
        /// \\[
        /// \\mbox{inletIndex}=\\begin{bmatrix}
        /// \\mbox{BC1_index}  \\\\ \
        /// \\mbox{BC2_index}  \\\\ \
        /// \\vdots & \\vdots \\\\ \
        /// \\mbox{BCN_index}  \\ \
        /// \\end{bmatrix}
        /// \\]
        Eigen::MatrixXi inletPatch;

#pragma GCC diagnostic pop

        // Functions

        //--------------------------------------------------------------------------
        /// Set Parameters Problems
        void setParameters();

        //--------------------------------------------------------------------------
        /// Generate Random Numbers
        void genRandPar();

        //--------------------------------------------------------------------------
        /// Generate Random Numbers given the dimension of the training set
        ///
        /// @param[in]  tsize  Dimension of the training set.
        ///
        void genRandPar(label tsize);

        //--------------------------------------------------------------------------
        /// Generate Equidistributed Numbers
        void genEquiPar();

        //--------------------------------------------------------------------------
        /// Perform a TruthSolve
        void truthSolve();

        Eigen::MatrixXi inletIndexT;

        //--------------------------------------------------------------------------
        /// Assign Boundary Condition to a volVectorField
        ///
        /// @param[in]  s       field where you want to assign the BC.
        /// @param[in]  BC_ind  The BC index.
        /// @param[in]  value   The value you want to assign (it must be an OpenFOAM vector).
        ///
        void assignBC(volVectorField& s, label BC_ind, Vector<double>& value);

        //--------------------------------------------------------------------------
        /// Assign Boundary Condition to a volScalarField
        ///
        /// @param[in]  s       field where you want to assign the BC.
        /// @param[in]  BC_ind  The BC index.
        /// @param[in]  value   The value you want to assign (it must be a double).
        ///
        void assignBC(volScalarField& s, label BC_ind, double& value);

        //--------------------------------------------------------------------------
        /// Exact reconstruction using a certain number of modes for vector list of
        /// fields and the projection coefficients (volVectorField)
        ///
        /// @param[out] rec_field2    The reconstructed field as PtrList of volVectorField.
        /// @param[in]  modes         The modes used for reconstruction as PtrList of volVectorField.
        /// @param[in]  Nmodes        The number of modes you want to use.
        /// @param[in]  coeff_matrix  The matrix of coefficients.
        ///
        void reconstructFromMatrix(PtrList<volVectorField>& rec_field2,
                                   PtrList<volVectorField>& modes, label Nmodes, Eigen::MatrixXd coeff_matrix);

        //--------------------------------------------------------------------------
        /// Exact reconstruction using a certain number of modes for vector list of
        /// fields and the projection coefficients (volScalarField)
        ///
        /// @param[out] rec_field2    The reconstructed field as PtrList of volScalarField.
        /// @param[in]  modes         The modes used for reconstruction as PtrList of volScalarField.
        /// @param[in]  Nmodes        The number of modes you want to use.
        /// @param[in]  coeff_matrix  The matrix of coefficients.
        ///
        /// Exact reconstruction using a certain number of modes for scalar list of fields
        void reconstructFromMatrix(PtrList<volScalarField>& rec_field2,
                                   PtrList<volScalarField>& modes, label Nmodes, Eigen::MatrixXd coeff_matrix);

        //--------------------------------------------------------------------------
        /// Assign internal field condition
        ///
        /// @param[in,out]  s      field where you want to assign the internal field condition
        /// @param          value  The value you want to assign
        ///
        /// @tparam         T      type of field (volVectorField or volScalarField)
        /// @tparam         G      type of value you want to assign (OpenFOAM vector or scalar)
        ///
        template<typename T, typename G>
        void assignIF(T& s, G& value);

        //--------------------------------------------------------------------------
        /// Homogenize the snapshot matrix, it works with PtrList of volVectorField and volScalarField
        ///
        /// @param[in]      Lfield     The list of snapshots to be homogenized.
        /// @param[in]      liftfield  The list containing.
        /// @param[out]     omfield    The homogenized snapshot matrix.
        ///
        /// @tparam     T          { description }
        ///
        template<typename T>
        void computeLift(T& Lfield, T& liftfield, T& omfield);

        //--------------------------------------------------------------------------
        /// Virtual function to compute the lifting function
        ///

        template<typename T>
        void computeLiftT(T& Lfield, T& liftfield, T& omfield);
        //--------------------------------------------------------------------------
        /// Virtual function to compute the lifting function for scalar field
        ///
        void liftSolve();
        void liftSolveT();
        //--------------------------------------------------------------------------
        /// General projection operation
        ///
        void project();

        //--------------------------------------------------------------------------
        /// Write out a list of scalar corresponding to the parameters used in the truthSolve
        ///
        /// @param[in]  mu_now  The list of scalars.
        ///
        void writeMu(List<scalar> mu_now);

        //--------------------------------------------------------------------------
        /// @brief      Constructs the parameters-coefficients manifold for vector fields, based on RBF-spline model
        /// @param[in]  snapshots   Snapshots vector fields, used to compute the coefficient matrix
        /// @param[in]  modes       POD modes vector fields, used to compute the coefficient matrix
        /// @param[in]  rbfBasis    The RBF basis type. Implemented bases are "GAUSSIAN", "THIN_PLATE",
        /// "MULTI_QUADRIC", "INVERSE_QUADRIC", and "INVERSE_MULTI_QUADRIC". Default basis is "Gaussian"
        ///
        /// @return     Vector of objects to the RBF splines corresponding to each mode
        ///
        std::vector<SPLINTER::RBFSpline> getCoeffManifoldRBF(PtrList<volVectorField>
                snapshots, PtrList<volVectorField>& modes, word rbfBasis = "GAUSSIAN");

        //--------------------------------------------------------------------------
        /// @brief      Constructs the parameters-coefficients manifold for scalar fields, based on RBF-spline model
        /// @param[in]  snapshots   Snapshots scalar fields, used to compute the coefficient matrix
        /// @param[in]  modes       POD modes scalar fields, used to compute the coefficient matrix
        /// @param[in]  rbfBasis    The RBF basis type. Implemented bases are "GAUSSIAN", "THIN_PLATE",
        /// "MULTI_QUADRIC", "INVERSE_QUADRIC", and "INVERSE_MULTI_QUADRIC". Default basis is "Gaussian"
        ///
        /// @return     Vector of objects to the RBF splines corresponding to each mode
        ///
        std::vector<SPLINTER::RBFSpline> getCoeffManifoldRBF(PtrList<volScalarField>
                snapshots, PtrList<volScalarField>& modes, word rbfBasis = "GAUSSIAN");

        //--------------------------------------------------------------------------
        /// @brief      Constructs the parameters-coefficients manifold for vector fields, based on the B-spline model
        /// @param[in]  snapshots   Snapshots vector fields, used to compute the coefficient matrix
        /// @param[in]  modes       POD modes vector fields, used to compute the coefficient matrix
        /// @param[in]  splDeg      The B-spline degree. Default value is 3
        ///
        /// @return     Vector of objects to the B-splines corresponding to each mode
        ///
        // std::vector<SPLINTER::BSpline*> getCoeffManifoldSPL(PtrList<volVectorField> snapshots, PtrList<volVectorField>& modes, label splDeg=3);
        std::vector<SPLINTER::BSpline> getCoeffManifoldSPL(PtrList<volVectorField>
                snapshots, PtrList<volVectorField>& modes, label splDeg = 3);

        //--------------------------------------------------------------------------
        /// @brief      Constructs the parameters-coefficients manifold for scalar fields, based on the B-spline model
        /// @param[in]  snapshots   Snapshots scalar fields, used to compute the coefficient matrix
        /// @param[in]  modes       POD modes scalar fields, used to compute the coefficient matrix
        /// @param[in]  splDeg      The B-spline degree. Default value is 3
        ///
        /// @return     Vector of objects to the B-splines corresponding to each mode
        ///
        std::vector<SPLINTER::BSpline> getCoeffManifoldSPL(PtrList<volScalarField>
                snapshots, PtrList<volScalarField>& modes, label splDeg = 3);

};

#ifdef NoRepository
#   include "reductionProblemTemplates.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// Additional Doxygen Documentation
/// \dir ../problems
/// Directory containing the implementation of the full order parametrized problems
/// and of the projection methods used to create the reduced matrices.

// Doxigen Documentation (Main Page)

/*! \mainpage In real Time Highly Advanced Computational Applications for Finite Volumes
 *
 * \section About About ITHACA-FV
 *
 * ITHACA-FV (In real Time Highly Advanced Computational Applications for
Finite Volumes) is C++ library based on the finite volume solver
[<b> OpenFOAM </b>](https://openfoam.org/). It consists of the implementation of several reduced order
modeling techniques for parametrized problems. ITHACA-FV is developed and maintained by
[<b> Dr. Giovanni Stabile </b>](http://people.sissa.it/~gstabile/) under the supervision of [<b> Prof. Gianluigi Rozza </b>](http://people.sissa.it/~grozza/) at [<b> SISSA mathLab </b>](https://mathlab.sissa.it/) an applied math laboratory of the
[<b> International School for advanced studies </b>](https://www.sissa.it/).

 *
 * \section structure ITHACA-FV directory structure
 * ITHACA-FV comprises of four main directories:

- \b src: the core ITHACA-FV library.
- \b applications: some utilities.
- \b tutorials: test-cases that demonstrate a range of ITHACA-FV functionalities.
- \b docs: documentation.
 *
 * \section step1 How to cite IHACA-FV.
 *
Most of the theoretical aspects behind ITHACA-FV are deeply explained in [<b> Stabile2017CAIM </b>](https://doi.org/10.1515/caim-2017-0011) and [<b> Stabile2017CAF </b>](https://doi.org/10.1016/j.compfluid.2018.01.035).
For this reason, if you use this software, please consider citing the mentioned works, reported in the following bibtex entries:


@Article{Stabile2017CAIM,
Title                    = {{POD-Galerkin reduced order methods for CFD using Finite Volume Discretisation: vortex shedding around a circular cylinder}},
Author                   = {Stabile, Giovanni and Hijazi, Saddam and Mola, Andrea and Lorenzi, Stefano and Rozza, Gianluigi},
Journal                  = {Communications in Applied and Industrial Mathematics},
Year                     = {(2017)},
Volume                   = {8},
Number                   = {1},
pages                    = {210-236},
Doi                      = {10.1515/caim-2017-0011}}

@Article{Stabile2017CAF,
Title                    = {{Finite volume POD-Galerkin stabilised reduced order methods for the parametrised incompressible Navier-Stokes equations}},
Author                   = {Stabile, Giovanni and Rozza, Gianluigi},
Journal                  = {Computers & Fluids},
Year                     = {2018},
Doi                      = {10.1016/j.compfluid.2018.01.035}}

And cite the ITHACA-FV [website](https://mathlab.sissa.it/ithaca-fv).

 */


