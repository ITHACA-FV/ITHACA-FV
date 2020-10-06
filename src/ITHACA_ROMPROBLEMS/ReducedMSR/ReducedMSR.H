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
#ifndef ReducedMSR_H
#define ReducedMSR_H

#include "fvCFD.H"
#include "IOmanip.H"
#include "ReducedProblem.H"
#include "msrProblem.H"
#include "ITHACAutilities.H"
#include <Eigen/Dense>
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>

/// Structure to implement a newton object for a stationary MSR problem
struct newton_msr_fd: public newton_argument<double>
{
    public:
        newton_msr_fd() {}
        newton_msr_fd(label Nx, label Ny, msrProblem& problem) :
            newton_argument<double>(Nx, Ny),
            problem(&problem),
            Nphi_u(problem.NUmodes + problem.liftfield.size()),
            Nphi_p(problem.NPmodes),
            N_BC(problem.inletIndex.rows())
        {}
        label operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const;
        label df(const Eigen::VectorXd& x,  Eigen::MatrixXd& fjac) const;

        label Nphi_u;
        label Nphi_p;
        label N_BC;

        scalar nu;

        Eigen::VectorXd BC;
        msrProblem* problem;
};

struct newton_msr_n: public newton_argument<double>
{
    public:
        newton_msr_n() {}
        newton_msr_n(label Nx, label Ny, msrProblem& problem) :
            newton_argument<double>(Nx, Ny),
            problem(&problem),
            Nphi_flux(problem.NFluxmodes),
            Nphi_prec1(problem.NPrecmodes(0)),
            Nphi_prec2(problem.NPrecmodes(1)),
            Nphi_prec3(problem.NPrecmodes(2)),
            Nphi_prec4(problem.NPrecmodes(3)),
            Nphi_prec5(problem.NPrecmodes(4)),
            Nphi_prec6(problem.NPrecmodes(5)),
            Nphi_prec7(problem.NPrecmodes(6)),
            Nphi_prec8(problem.NPrecmodes(7)),
            Nphi_const(problem.NCmodes),
            d_c(problem.NCmodes),
            nsf_c(problem.NCmodes),
            a_c(problem.NCmodes)
        {}
        label operator()(const Eigen::VectorXd& n, Eigen::VectorXd& fvecn) const;
        label df(const Eigen::VectorXd& n,  Eigen::MatrixXd& fjacn) const;

        label Nphi_u;
        label Nphi_flux;
        label Nphi_prec1;
        label Nphi_prec2;
        label Nphi_prec3;
        label Nphi_prec4;
        label Nphi_prec5;
        label Nphi_prec6;
        label Nphi_prec7;
        label Nphi_prec8;
        label Nphi_const;
        //fluid-dynamics constants
        scalar nu;
        //neutronic constants
        scalar iv;
        scalar l1; //lambda-ith
        scalar l2;
        scalar l3;
        scalar l4;
        scalar l5;
        scalar l6;
        scalar l7;
        scalar l8;
        scalar b1; //beta-ith
        scalar b2;
        scalar b3;
        scalar b4;
        scalar b5;
        scalar b6;
        scalar b7;
        scalar b8;
        scalar btot;
        scalar Sc;

        Eigen::VectorXd d_c;
        Eigen::VectorXd nsf_c;
        Eigen::VectorXd a_c;
        Eigen::VectorXd a_tmp; //to store u time coeff
        std::vector<SPLINTER::RBFSpline*> SPLINES_d;
        std::vector<SPLINTER::RBFSpline*> SPLINES_nsf;
        std::vector<SPLINTER::RBFSpline*> SPLINES_a;
        msrProblem* problem;
};

struct newton_msr_t: public newton_argument<double>
{
    public:
        newton_msr_t() {}
        newton_msr_t(label Nx, label Ny, msrProblem& problem) :
            newton_argument<double>(Nx, Ny),
            problem(&problem),
            Nphi_T(problem.NTmodes),
            Nphi_dec1(problem.NDecmodes(0)),
            Nphi_dec2(problem.NDecmodes(1)),
            Nphi_dec3(problem.NDecmodes(2)),
            N_BCt(problem.inletIndexT.rows()),
            Nphi_const(problem.NCmodes),
            v_c(problem.NCmodes),
            sp_c(problem.NCmodes),
            txs_c(problem.NCmodes)
        {}
        label operator()(const Eigen::VectorXd& n, Eigen::VectorXd& fvecn) const;
        label df(const Eigen::VectorXd& n,  Eigen::MatrixXd& fjacn) const;

        label Nphi_u;
        label Nphi_flux;
        label Nphi_T;
        label Nphi_dec1;
        label Nphi_dec2;
        label Nphi_dec3;
        label N_BCt;
        label Nphi_const;
        //fluid-dynamics constants
        scalar nu;
        //thermal constants
        scalar cp;
        scalar dl1;
        scalar dl2;
        scalar dl3;
        scalar db1;
        scalar db2;
        scalar db3;
        scalar dbtot;
        scalar Sc;
        scalar Pr;

        Eigen::VectorXd a_tmp;
        Eigen::VectorXd c_tmp;
        Eigen::VectorXd BCt;
        Eigen::VectorXd v_c;
        Eigen::VectorXd sp_c;
        Eigen::VectorXd txs_c;
        std::vector<SPLINTER::RBFSpline*> SPLINES_v;
        std::vector<SPLINTER::RBFSpline*> SPLINES_sp;
        std::vector<SPLINTER::RBFSpline*> SPLINES_TXS;
        msrProblem* problem;
};


/*---------------------------------------------------------------------------*\
                        Class reducedProblem Declaration
\*---------------------------------------------------------------------------*/

class reducedMSR: public reducedProblem
{

    public:

        reducedMSR();
        explicit reducedMSR(msrProblem& problem);
        ~reducedMSR() {};

        /// Vector to store the previous solution during the Newton procedure
        Eigen::VectorXd y_old;
        Eigen::VectorXd w_old;
        Eigen::VectorXd z_old;

        /// Vector to store the solution during the Newton procedure
        Eigen::VectorXd y;
        Eigen::VectorXd w;
        Eigen::VectorXd z;


        /// characteristic constants of the problem
        //fluid-dynamics constants
        scalar nu;
        scalar rho;
        scalar Pr;
        scalar Sc;
        //neutronic constants
        scalar d;
        scalar m;
        scalar nsf;
        scalar a;
        scalar Keff;
        scalar iv;
        scalar l1; //lambda-ith
        scalar l2;
        scalar l3;
        scalar l4;
        scalar l5;
        scalar l6;
        scalar l7;
        scalar l8;
        scalar b1; //beta-ith
        scalar b2;
        scalar b3;
        scalar b4;
        scalar b5;
        scalar b6;
        scalar b7;
        scalar b8;
        scalar btot;
        //thermal constants
        scalar cp;
        scalar dl1;
        scalar dl2;
        scalar dl3;
        scalar db1;
        scalar db2;
        scalar db3;
        scalar dbtot;
        scalar sp;



        /// List of Eigen matrices to store the online solution
        List < Eigen::MatrixXd> online_solution_fd;
        List < Eigen::MatrixXd> online_solution_n;
        List < Eigen::MatrixXd> online_solution_t;
        List < Eigen::MatrixXd> online_solution_C;

        /// List of pointers to store the modes for each field
        PtrList<volVectorField> Umodes;
        PtrList<volScalarField> Pmodes;
        PtrList<volScalarField> Fluxmodes;
        PtrList<volScalarField> Prec1modes;
        PtrList<volScalarField> Prec2modes;
        PtrList<volScalarField> Prec3modes;
        PtrList<volScalarField> Prec4modes;
        PtrList<volScalarField> Prec5modes;
        PtrList<volScalarField> Prec6modes;
        PtrList<volScalarField> Prec7modes;
        PtrList<volScalarField> Prec8modes;
        PtrList<volScalarField> Tmodes;
        PtrList<volScalarField> Dec1modes;
        PtrList<volScalarField> Dec2modes;
        PtrList<volScalarField> Dec3modes;
        PtrList<volScalarField> vmodes;
        PtrList<volScalarField> Dmodes;
        PtrList<volScalarField> NSFmodes;
        PtrList<volScalarField> Amodes;
        PtrList<volScalarField> SPmodes;
        PtrList<volScalarField> TXSmodes;


        /// List of pointers to store the snapshots for each field
        PtrList<volVectorField> Usnapshots;
        PtrList<volScalarField> Psnapshots;
        PtrList<volScalarField> Fluxsnapshots;
        PtrList<volScalarField> Prec1snapshots;
        PtrList<volScalarField> Prec2snapshots;
        PtrList<volScalarField> Prec3snapshots;
        PtrList<volScalarField> Prec4snapshots;
        PtrList<volScalarField> Prec5snapshots;
        PtrList<volScalarField> Prec6snapshots;
        PtrList<volScalarField> Prec7snapshots;
        PtrList<volScalarField> Prec8snapshots;
        PtrList<volScalarField> Tsnapshots;
        PtrList<volScalarField> Dec1snapshots;
        PtrList<volScalarField> Dec2snapshots;
        PtrList<volScalarField> Dec3snapshots;
        PtrList<volScalarField> vsnapshots;
        PtrList<volScalarField> Dsnapshots;
        PtrList<volScalarField> NSFsnapshots;
        PtrList<volScalarField> Asnapshots;
        PtrList<volScalarField> SPsnapshots;
        PtrList<volScalarField> TXSsnapshots;

        /// Recontructed fields
        PtrList<volVectorField> UREC;
        PtrList<volScalarField> PREC;
        PtrList<volScalarField> FLUXREC;
        PtrList<volScalarField> PREC1REC;
        PtrList<volScalarField> PREC2REC;
        PtrList<volScalarField> PREC3REC;
        PtrList<volScalarField> PREC4REC;
        PtrList<volScalarField> PREC5REC;
        PtrList<volScalarField> PREC6REC;
        PtrList<volScalarField> PREC7REC;
        PtrList<volScalarField> PREC8REC;
        PtrList<volScalarField> TREC;
        PtrList<volScalarField> DEC1REC;
        PtrList<volScalarField> DEC2REC;
        PtrList<volScalarField> DEC3REC;
        PtrList<volScalarField> POWERDENSREC;
        PtrList<volScalarField> vREC;
        PtrList<volScalarField> DREC;
        PtrList<volScalarField> NSFREC;
        PtrList<volScalarField> AREC;
        PtrList<volScalarField> SPREC;
        PtrList<volScalarField> TXSREC;

        /// Newton object used to solve the non linear problem
        newton_msr_fd newton_object_fd;
        newton_msr_n  newton_object_n;
        newton_msr_t  newton_object_t;

        /// Pointer to the FOM problem
        msrProblem* problem;

        /// Number of modes for each field
        label Nphi_u;
        label Nphi_p;
        label Nphi_flux;
        label Nphi_prec1;
        label Nphi_prec2;
        label Nphi_prec3;
        label Nphi_prec4;
        label Nphi_prec5;
        label Nphi_prec6;
        label Nphi_prec7;
        label Nphi_prec8;
        label Nphi_T;
        label Nphi_dec1;
        label Nphi_dec2;
        label Nphi_dec3;
        label Nphi_const;

        ///boolean variable to check if the user wants to reconstruct
        /// all the three physics of the problem
        bool recall = false;

        /// Number of parametrized boundary conditions
        label N_BC;
        label N_BCt;

        /// Counter to count the online solutions
        label count_online_solve = 1;

        //Methods:



        /// Method to perform an online solve using a PPE stabilisation method
        void solveOnline(Eigen::MatrixXd vel_now, Eigen::MatrixXd temp_now,
                         Eigen::VectorXd mu_online);

        ///Method to clear all the fields of MSR (sets the size to zero)
        void clearFields();

        /// Methods to reconstruct a solution from an online solve with a PPE stabilisation technique.
        /// stabilisation method
        ///folder:      The folder where to output the solutions
        ///printevery=x:  save the solution every x timesteps
        ///reconstruct all the three physics together:
        void reconstructAP(fileName folder = "./ITHACAOutput/online_rec",
                           label printevery = 1);
        ///reconstruct fluid-dynamics
        void reconstruct_fd(fileName folder = "./ITHACAOutput/online_rec",
                            label printevery = 1);
        ///reconstruct neutronics
        void reconstruct_n(fileName folder = "./ITHACAOutput/online_rec",
                           label printevery = 1);
        ///reconstruct thermal fields
        void reconstruct_t(fileName folder = "./ITHACAOutput/online_rec",
                           label printevery = 1);
        ///reconstruct temperature dependent constants
        void reconstruct_C(fileName folder = "./ITHACAOutput/online_rec",
                           label printevery = 1);



    protected:

        //Methods:

        /// Method to load all the constants needed in the ROM from         ///the FOM
        void loadConstants(msrProblem* problem);



};


#endif
