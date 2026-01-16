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
    ReducedUnsteadyNSTurb
Description
    A reduced problem for the unsteady turbulent NS equations
SourceFiles
    ReducedUnsteadyNSTurb.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the ReducedUnsteadyNSTurb class.
/// \dir
/// Directory containing the header and source files for the ReducedUnsteadyNSTurb class

#ifndef ReducedUnsteadyNSTurb_H
#define ReducedUnsteadyNSTurb_H

#include "fvCFD.H"
#include "IOmanip.H"
#include "ReducedSteadyNS.H"
#include "ReducedUnsteadyNS.H"
#include "UnsteadyNSTurb.H"
#include <Eigen/Dense>
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>

// Forward declarations for Splinter types we only keep pointers to
namespace SPLINTER { class RBFSpline; class DataTable; }

struct newtonUnsteadyNSTurbSUP: public newton_argument<double>
{
    public:
        newtonUnsteadyNSTurbSUP() {}

        newtonUnsteadyNSTurbSUP(int Nx, int Ny,
                                UnsteadyNSTurb& problem): newton_argument<double>(Nx, Ny),
            problem(& problem),
            Nphi_u(problem.NUmodes + problem.liftfield.size() + problem.NSUPmodes),
            nphiNut(problem.nNutModes),
            Nphi_p(problem.NPmodes),
            N_BC(problem.inletIndex.rows()),
            gNut(problem.nNutModes)
        {}

        int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const;
        int df(const Eigen::VectorXd& x,  Eigen::MatrixXd& fjac) const;

        UnsteadyNSTurb* problem;
        int Nphi_u;
        int nphiNut;
        int Nphi_p;
        int N_BC;
        scalar nu;
        scalar dt;
        Eigen::VectorXd y_old;
        Eigen::VectorXd yOldOld;
        Eigen::VectorXd bc;
        Eigen::MatrixXd tauU;
        Eigen::VectorXd gNut;
        std::vector<SPLINTER::RBFSpline*> SPLINES;
};



struct newtonUnsteadyNSTurbPPE: public newton_argument<double>
{
    public:
        newtonUnsteadyNSTurbPPE() {}

        newtonUnsteadyNSTurbPPE(int Nx, int Ny,
                                UnsteadyNSTurb& problem): newton_argument<double>(Nx, Ny),
            problem(& problem),
            Nphi_u(problem.NUmodes + problem.liftfield.size()),
            nphiNut(problem.nNutModes),
            Nphi_p(problem.NPmodes),
            N_BC(problem.inletIndex.rows()),
            gNut(problem.nNutModes)
        {}

        int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const;
        int df(const Eigen::VectorXd& x,  Eigen::MatrixXd& fjac) const;

        UnsteadyNSTurb* problem;
        int Nphi_u;
        int nphiNut;
        int Nphi_p;
        int N_BC;
        scalar nu;
        scalar dt;
        Eigen::VectorXd y_old;
        Eigen::VectorXd yOldOld;
        Eigen::VectorXd bc;
        Eigen::MatrixXd tauU;
        Eigen::VectorXd gNut;
        std::vector<SPLINTER::RBFSpline*> SPLINES;
};

struct newtonUnsteadyNSTurbSUPAve: public newton_argument<double>
{
    public:
        newtonUnsteadyNSTurbSUPAve() {}

        newtonUnsteadyNSTurbSUPAve(int Nx, int Ny,
                                   UnsteadyNSTurb& problem): newton_argument<double>(Nx, Ny),
            problem(& problem),
            Nphi_u(problem.NUmodes + problem.liftfield.size() + problem.NSUPmodes),
            nphiNut(problem.nNutModes),
            Nphi_p(problem.NPmodes),
            N_BC(problem.inletIndex.rows()),
            gNut(problem.nNutModes),
            gNutAve(problem.nutAve.size())
        {}

        int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const;
        int df(const Eigen::VectorXd& x,  Eigen::MatrixXd& fjac) const;

        UnsteadyNSTurb* problem;
        int Nphi_u;
        int nphiNut;
        int Nphi_p;
        int N_BC;
        scalar nu;
        scalar dt;
        Eigen::VectorXd y_old;
        Eigen::VectorXd yOldOld;
        Eigen::VectorXd bc;
        Eigen::MatrixXd tauU;
        Eigen::VectorXd gNut;
        Eigen::VectorXd gNutAve;
        // used in .C
        Eigen::VectorXd gNutFluct;
        std::vector<SPLINTER::RBFSpline*> SPLINES;
};

struct newtonUnsteadyNSTurbPPEAve: public newton_argument<double>
{
    public:
        newtonUnsteadyNSTurbPPEAve() {}

        newtonUnsteadyNSTurbPPEAve(int Nx, int Ny,
                                   UnsteadyNSTurb& problem): newton_argument<double>(Nx, Ny),
            problem(& problem),
            Nphi_u(problem.NUmodes + problem.liftfield.size()),
            nphiNut(problem.nNutModes),
            Nphi_p(problem.NPmodes),
            N_BC(problem.inletIndex.rows()),
            gNut(problem.nNutModes),
            gNutAve(problem.nutAve.size())

        {}

        int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const;
        int df(const Eigen::VectorXd& x,  Eigen::MatrixXd& fjac) const;

        UnsteadyNSTurb* problem;
        int Nphi_u;
        int nphiNut;
        int Nphi_p;
        int N_BC;
        scalar nu;
        scalar dt;
        Eigen::VectorXd y_old;
        Eigen::VectorXd yOldOld;
        Eigen::VectorXd bc;
        Eigen::MatrixXd tauU;
        Eigen::VectorXd gNut;
        Eigen::VectorXd gNutAve;
        // used in .C
        Eigen::VectorXd gNutFluct;
        std::vector<SPLINTER::RBFSpline*> SPLINES;
};



/*---------------------------------------------------------------------------*\
                        Class reducedProblem Declaration
\*---------------------------------------------------------------------------*/

class ReducedUnsteadyNSTurb: public reducedUnsteadyNS
{
    private:

    public:
        // Constructors
        ReducedUnsteadyNSTurb();
        explicit ReducedUnsteadyNSTurb(UnsteadyNSTurb& problem);
        ~ReducedUnsteadyNSTurb() {};

        /// Pointer to the FOM problem
        UnsteadyNSTurb* problem;

        /// Number of viscosity modes (combined)
        int nphiNut;

        /// Split counts for average/fluctuation ν_t
        int nphiNutAvg{0};
        int nphiNutFluct{0};

        /// Dimension of the interpolation independent variable
        int dimA;

        /// Legacy combined ν_t modes
        PtrList<volScalarField> nutModes;

        /// Split ν_t modes
        PtrList<volScalarField> nutAvgModes;
        PtrList<volScalarField> nutFluctModes;

        /// Coeff buffers for avg/fluct ν_t
        Eigen::VectorXd coeffNutAvg;
        Eigen::VectorXd coeffNutFluct;

        /// RBF sample tables (match FOM types)
        List<SPLINTER::DataTable*> samplesNutAvg;
        List<SPLINTER::DataTable*> samplesNutFluct;

        /// RBF splines (match FOM types)
        List<SPLINTER::RBFSpline*> rbfSplinesNutAvg;
        List<SPLINTER::RBFSpline*> rbfSplinesNutFluct;

        /// The matrix of the eddy viscosity RBF interpolated coefficients
        Eigen::MatrixXd rbfCoeffMat;

        /// The matrix of the initial velocity and pressure reduced coefficients
        Eigen::MatrixXd initCond;

        /// The initial eddy viscosity reduced coefficients
        Eigen::VectorXd nut0;

        /// Newton objects
        newtonUnsteadyNSTurbSUP     newtonObjectSUP;
        newtonUnsteadyNSTurbSUPAve  newtonObjectSUPAve;
        newtonUnsteadyNSTurbPPE     newtonObjectPPE;
        newtonUnsteadyNSTurbPPEAve  newtonObjectPPEAve;

        /// Online parameter value
        Eigen::VectorXd muStar;

        /// The reduced average vector for viscosity coefficients
        Eigen::VectorXd gNutAve;

        /// Interpolation independent variable choice
        int interChoice = 1;

        /// Interpolation boolean variable to skip lifting functions
        bool skipLift = false;

        /// Reconstructed eddy viscosity fields list
        PtrList<volScalarField> nutRecFields;

        /// Reconstructed U and p fields (used in reconstruct)
        PtrList<volVectorField> uRecFields;
        PtrList<volScalarField> pRecFields;

        /// Time series container (each entry is an (N+1)x1 column matrix [t; a; b])
        List<Eigen::MatrixXd> online_solution;

        /// Sizes used in constructor/reconstruct
        int Nphi_u{0};
        int Nphi_p{0};
        int N_BC  {0};

        // Functions
        void solveOnlinePPE(Eigen::MatrixXd velNow);
        void solveOnlineSUP(Eigen::MatrixXd velNow);
        void solveOnlineSUPAve(Eigen::MatrixXd velNow);
        void solveOnlinePPEAve(Eigen::MatrixXd velNow);

        void reconstruct(bool exportFields = false,
                         fileName folder = "./online_rec");

        Eigen::MatrixXd setOnlineVelocity(Eigen::MatrixXd vel);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
