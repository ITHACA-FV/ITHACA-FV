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
#ifndef ReducedusMSR_H
#define ReducedusMSR_H

#include "fvCFD.H"
#include "IOmanip.H"
#include "ReducedMSR.H"
#include "usmsrProblem.H"
#include <Eigen/Dense>
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>
#include <fstream>
#include <vector>
#include <stdlib.h>

struct newton_usmsr_fd: public newton_argument<double>
{
    public:
        newton_usmsr_fd() {}
        newton_usmsr_fd(int Nx, int Ny, usmsrProblem& problem):
            newton_argument<double>(Nx, Ny),
            problem(&problem),
            Nphi_u(problem.NUmodes + problem.liftfield.size()),
            Nphi_p(problem.NPmodes),
            N_BC(problem.inletIndex.rows())
        {}
        int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const;
        int df(const Eigen::VectorXd& x, Eigen::MatrixXd& fjac) const;

        int Nphi_u;
        int Nphi_p;
        int N_BC;

        scalar nu;
        scalar dt;
        Eigen::VectorXd y_old;
        Eigen::VectorXd BC;

        usmsrProblem* problem;
};

struct newton_usmsr_n: public newton_argument<double>
{
    public:
        newton_usmsr_n() {}
        newton_usmsr_n(int Nx, int Ny, usmsrProblem& problem):
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
        int operator()(const Eigen::VectorXd& n, Eigen::VectorXd& fvecn) const;
        int df(const Eigen::VectorXd& n, Eigen::MatrixXd& fjacn) const;
        int Nphi_flux;
        int Nphi_prec1;
        int Nphi_prec2;
        int Nphi_prec3;
        int Nphi_prec4;
        int Nphi_prec5;
        int Nphi_prec6;
        int Nphi_prec7;
        int Nphi_prec8;
        int Nphi_u;
        int Nphi_const;
        scalar nu;
        //neutronic constants
        scalar d;
        scalar m;
        scalar nsf;
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
        scalar Sc;

        scalar dt;
        Eigen::VectorXd d_c;
        Eigen::VectorXd nsf_c;
        Eigen::VectorXd a_c;
        Eigen::VectorXd a_tmp; //to store u time coeff
        Eigen::VectorXd w_old;
        std::vector<SPLINTER::RBFSpline*> SPLINES_d;
        std::vector<SPLINTER::RBFSpline*> SPLINES_nsf;
        std::vector<SPLINTER::RBFSpline*> SPLINES_a;
        usmsrProblem* problem;
};

struct newton_usmsr_t: public newton_argument<double>
{
    public:
        newton_usmsr_t() {}
        newton_usmsr_t(int Nx, int Ny, usmsrProblem& problem):
            newton_argument<double>(Nx, Ny),
            problem(&problem),
            Nphi_T(problem.NTmodes + problem.liftfieldT.size()),
            Nphi_dec1(problem.NDecmodes(0)),
            Nphi_dec2(problem.NDecmodes(1)),
            Nphi_dec3(problem.NDecmodes(2)),
            N_BCt(problem.inletIndexT.rows()),
            Nphi_const(problem.NCmodes),
            v_c(problem.NCmodes),
            sp_c(problem.NCmodes),
            txs_c(problem.NCmodes)
        {}
        int operator()(const Eigen::VectorXd& t, Eigen::VectorXd& fvect) const;
        int df(const Eigen::VectorXd& t, Eigen::MatrixXd& fjact) const;
        int Nphi_T;
        int Nphi_dec1;
        int Nphi_dec2;
        int Nphi_dec3;
        int Nphi_u;
        int Nphi_flux;
        int N_BCt;
        int Nphi_const;
        scalar nu;
        //neutronic constants
        scalar Keff;
        //thermal constants
        scalar sp;
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

        scalar dt;
        Eigen::VectorXd a_tmp; //to store u time coeff.
        Eigen::VectorXd c_tmp; //to store flux time coeff.
        Eigen::VectorXd z_old;
        Eigen::VectorXd BCt;
        Eigen::VectorXd v_c;
        Eigen::VectorXd sp_c;
        Eigen::VectorXd txs_c;
        std::vector<SPLINTER::RBFSpline*> SPLINES_v;
        std::vector<SPLINTER::RBFSpline*> SPLINES_sp;
        std::vector<SPLINTER::RBFSpline*> SPLINES_TXS;

        usmsrProblem* problem;
};

/*---------------------------------------------------------------------------*\
                        Class reducedProblem Declaration
\*---------------------------------------------------------------------------*/
class reducedusMSR: public reducedMSR
{
    public:
        reducedusMSR();
        explicit reducedusMSR(usmsrProblem& problem);
        ~reducedusMSR() {};

        newton_usmsr_fd newton_object_fd;
        newton_usmsr_n  newton_object_n;
        newton_usmsr_t  newton_object_t;
        scalar time;
        scalar dt;
        scalar finalTime;
        scalar tstart;
        usmsrProblem* problem;

        bool recall = false;
        void solveOnline(Eigen::MatrixXd vel_now, Eigen::MatrixXd temp_now,
                         Eigen::VectorXd mu_online, int startSnap = 0);
        void reconstructAP(fileName folder = "./ITHACAOutput/online_rec",
                           int printevery = 1);
        void reconstruct_fd(fileName folder = "./ITHACAOutput/online_rec",
                            int printevery = 1);
        void reconstruct_n(fileName folder = "./ITHACAOutput/online_rec",
                           int printevery = 1);
        void reconstruct_t(fileName folder = "./ITHACAOutput/online_rec",
                           int printevery = 1);
        void reconstruct_C(fileName folder = "./ITHACAOutput/online_rec",
                           int printevery = 1);


};




#endif

