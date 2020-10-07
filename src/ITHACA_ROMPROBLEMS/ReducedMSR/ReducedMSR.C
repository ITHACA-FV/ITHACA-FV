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
#include "ReducedMSR.H"

reducedMSR::reducedMSR()
{
}

reducedMSR::reducedMSR(msrProblem& FOMproblem)
{
    problem = &FOMproblem;
    N_BC = problem->inletIndex.rows();
    N_BCt = problem->inletIndexT.rows();
    Nphi_u = problem->B_matrix.rows();
    Nphi_p = problem->K_matrix.cols();
    Nphi_flux = problem->LF_matrix[0].cols();
    Nphi_prec1 = problem->PS1_matrix.cols();
    Nphi_prec2 = problem->PS2_matrix.cols();
    Nphi_prec3 = problem->PS3_matrix.cols();
    Nphi_prec4 = problem->PS4_matrix.cols();
    Nphi_prec5 = problem->PS5_matrix.cols();
    Nphi_prec6 = problem->PS6_matrix.cols();
    Nphi_prec7 = problem->PS7_matrix.cols();
    Nphi_prec8 = problem->PS8_matrix.cols();
    Nphi_T = problem->LT_matrix.rows();
    Nphi_dec1 = problem->LD1_matrix.rows();
    Nphi_dec2 = problem->LD2_matrix.rows();
    Nphi_dec3 = problem->LD3_matrix.rows();
    Nphi_const = problem->PF_matrix[0].rows();
    loadConstants(problem);

    //load the modes function
    for (int k = 0; k < problem->liftfield.size(); k++)
    {
        Umodes.append(problem->liftfield[k]);
    }

    for (int k = 0; k < problem->NUmodes; k++)
    {
        Umodes.append(problem->Umodes[k]);
    }

    for (int k = 0; k < problem->NPmodes; k++)
    {
        Pmodes.append(problem->Pmodes[k]);
    }

    for (int k = 0; k < problem->NFluxmodes; k++)
    {
        Fluxmodes.append(problem->Fluxmodes[k]);
    }

    for (int k = 0; k < Nphi_prec1; k++)
    {
        Prec1modes.append(problem->Prec1modes[k]);
    }

    for (int k = 0; k < Nphi_prec2; k++)
    {
        Prec2modes.append(problem->Prec2modes[k]);
    }

    for (int k = 0; k < Nphi_prec3; k++)
    {
        Prec3modes.append(problem->Prec3modes[k]);
    }

    for (int k = 0; k < Nphi_prec4; k++)
    {
        Prec4modes.append(problem->Prec4modes[k]);
    }

    for (int k = 0; k < Nphi_prec5; k++)
    {
        Prec5modes.append(problem->Prec5modes[k]);
    }

    for (int k = 0; k < Nphi_prec6; k++)
    {
        Prec6modes.append(problem->Prec6modes[k]);
    }

    for (int k = 0; k < Nphi_prec7; k++)
    {
        Prec7modes.append(problem->Prec7modes[k]);
    }

    for (int k = 0; k < Nphi_prec8; k++)
    {
        Prec8modes.append(problem->Prec8modes[k]);
    }

    for (int k = 0; k < problem->liftfieldT.size(); k++)
    {
        Tmodes.append(problem->liftfieldT[k]);
    }

    for (int k = 0; k < problem->NTmodes; k++)
    {
        Tmodes.append(problem->Tmodes[k]);
    }

    for (int k = 0; k < Nphi_dec1; k++)
    {
        Dec1modes.append(problem->Dec1modes[k]);
    }

    for (int k = 0; k < Nphi_dec2; k++)
    {
        Dec2modes.append(problem->Dec2modes[k]);
    }

    for (int k = 0; k < Nphi_dec3; k++)
    {
        Dec3modes.append(problem->Dec3modes[k]);
    }

    for (int k = 0; k < Nphi_const; k++)
    {
        vmodes.append(problem->vmodes[k]);
        Dmodes.append(problem->Dmodes[k]);
        NSFmodes.append(problem->NSFmodes[k]);
        Amodes.append(problem->Amodes[k]);
        SPmodes.append(problem->SPmodes[k]);
        TXSmodes.append(problem->TXSmodes[k]);
    }

    //load the snapshots
    for (int k = 0; k < problem->Ufield.size(); k++)
    {
        Usnapshots.append(problem->Ufield[k]);
        Psnapshots.append(problem->Pfield[k]);
        Fluxsnapshots.append(problem->Fluxfield[k]);
        Prec1snapshots.append(problem->Prec1field[k]);
        Prec2snapshots.append(problem->Prec2field[k]);
        Prec3snapshots.append(problem->Prec3field[k]);
        Prec4snapshots.append(problem->Prec4field[k]);
        Prec5snapshots.append(problem->Prec5field[k]);
        Prec6snapshots.append(problem->Prec6field[k]);
        Prec7snapshots.append(problem->Prec7field[k]);
        Prec8snapshots.append(problem->Prec8field[k]);
        Tsnapshots.append(problem->Tfield[k]);
        Dec1snapshots.append(problem->Dec1field[k]);
        Dec2snapshots.append(problem->Dec2field[k]);
        Dec3snapshots.append(problem->Dec3field[k]);
        vsnapshots.append(problem->vFields[k]);
        Dsnapshots.append(problem->DFields[k]);
        NSFsnapshots.append(problem->NSFFields[k]);
        Asnapshots.append(problem->AFields[k]);
        SPsnapshots.append(problem->SPFields[k]);
        TXSsnapshots.append(problem->TXSFields[k]);
    }

    newton_object_fd = newton_msr_fd(Nphi_u + Nphi_p, Nphi_u + Nphi_p, FOMproblem);
    newton_object_n = newton_msr_n(Nphi_flux + Nphi_prec1 + Nphi_prec2 + Nphi_prec3
                                   + Nphi_prec4 + Nphi_prec5 + Nphi_prec6 + Nphi_prec7 + Nphi_prec8,
                                   Nphi_flux + Nphi_prec1 + Nphi_prec2 + Nphi_prec3 + Nphi_prec4 + Nphi_prec5 +
                                   Nphi_prec6 + Nphi_prec7 + Nphi_prec8, FOMproblem);
    newton_object_t = newton_msr_t(Nphi_T + Nphi_dec1 + Nphi_dec2 + Nphi_dec3,
                                   Nphi_T + Nphi_dec1 + Nphi_dec2 + Nphi_dec3, FOMproblem);
}

int newton_msr_fd::operator()(const Eigen::VectorXd& x,
                              Eigen::VectorXd& fvec) const
{
    Eigen::VectorXd a_tmp(Nphi_u);
    Eigen::VectorXd b_tmp(Nphi_p);
    a_tmp = x.head(Nphi_u);
    b_tmp = x.tail(Nphi_p);
    /// Fluid-dynamics terms
    // Convective terms
    Eigen::MatrixXd cc(1, 1);
    Eigen::MatrixXd gg(1, 1);
    Eigen::MatrixXd bb(1, 1);
    // Mom Term
    Eigen::VectorXd M1 = problem->B_matrix * a_tmp * nu;
    // Gradient of pressure
    Eigen::VectorXd M2 = problem->K_matrix * b_tmp;
    // Pressure Term
    Eigen::VectorXd M3 = problem->D_matrix * b_tmp;
    // BC PPE
    Eigen::VectorXd M6 = problem->BC1_matrix * a_tmp * nu;
    // BC PPE
    Eigen::VectorXd M7 = problem->BC3_matrix * a_tmp * nu;

    for (int i = 0; i < Nphi_u; i++)
    {
        cc = a_tmp.transpose() * problem->C_matrix[i] * a_tmp;
        fvec(i) =  M1(i) - cc(0, 0) - M2(i);
    }

    int p_fvec = Nphi_u;

    for (int i = 0; i < Nphi_p; i++)
    {
        int k = i + p_fvec;
        gg = a_tmp.transpose() * problem->G_matrix[i] * a_tmp;
        //bb = a_tmp.transpose() * problem->BC2_matrix[i] * a_tmp;
        fvec(k) =  M3(i, 0) + gg(0, 0) - M7(i, 0);
    }

    for (int j = 0; j < N_BC; j++)
    {
        fvec(j) = x(j) - BC(j);
    }

    return 0;
}

int newton_msr_fd::df(const Eigen::VectorXd& x,  Eigen::MatrixXd& fjac) const
{
    Eigen::NumericalDiff<newton_msr_fd> numDiff(*this);
    numDiff.df(x, fjac);
    return 0;
}

int newton_msr_n::operator()(const Eigen::VectorXd& n,
                             Eigen::VectorXd& fvecn) const
{
    Eigen::VectorXd c_tmp(Nphi_flux);  //for flux
    Eigen::VectorXd d1_tmp(Nphi_prec1); //for prec1
    Eigen::VectorXd d2_tmp(Nphi_prec2); //for prec2
    Eigen::VectorXd d3_tmp(Nphi_prec3); //for prec3
    Eigen::VectorXd d4_tmp(Nphi_prec4); //for prec4
    Eigen::VectorXd d5_tmp(Nphi_prec5); //for prec5
    Eigen::VectorXd d6_tmp(Nphi_prec6); //for prec6
    Eigen::VectorXd d7_tmp(Nphi_prec7); //for prec7
    Eigen::VectorXd d8_tmp(Nphi_prec8); //for prec8
    c_tmp = n.head(Nphi_flux);
    int pos = Nphi_flux;
    d1_tmp = n.segment(pos, Nphi_prec1);
    pos = pos + Nphi_prec1;
    d2_tmp = n.segment(pos, Nphi_prec2);
    pos = pos + Nphi_prec2;
    d3_tmp = n.segment(pos, Nphi_prec3);
    pos = pos + Nphi_prec3;
    d4_tmp = n.segment(pos, Nphi_prec4);
    pos = pos + Nphi_prec4;
    d5_tmp = n.segment(pos, Nphi_prec5);
    pos = pos + Nphi_prec5;
    d6_tmp = n.segment(pos, Nphi_prec6);
    pos = pos + Nphi_prec6;
    d7_tmp = n.segment(pos, Nphi_prec7);
    pos = pos + Nphi_prec7;
    d8_tmp = n.segment(pos, Nphi_prec8);
    /// Neutronics terms
    // Laplacian flux term
    Eigen::MatrixXd lf(1, 1);
    // flux production term
    Eigen::MatrixXd pf(1, 1);
    // flux absorption term
    Eigen::MatrixXd af(1, 1);
    // precursor sources
    Eigen::VectorXd F3_1 = problem->PS1_matrix * d1_tmp * l1;
    Eigen::VectorXd F3_2 = problem->PS2_matrix * d2_tmp * l2;
    Eigen::VectorXd F3_3 = problem->PS3_matrix * d3_tmp * l3;
    Eigen::VectorXd F3_4 = problem->PS4_matrix * d4_tmp * l4;
    Eigen::VectorXd F3_5 = problem->PS5_matrix * d5_tmp * l5;
    Eigen::VectorXd F3_6 = problem->PS6_matrix * d6_tmp * l6;
    Eigen::VectorXd F3_7 = problem->PS7_matrix * d7_tmp * l7;
    Eigen::VectorXd F3_8 = problem->PS8_matrix * d8_tmp * l8;
    // Convective terms in prec-eq:
    Eigen::MatrixXd pp1(1, 1);
    Eigen::MatrixXd pp2(1, 1);
    Eigen::MatrixXd pp3(1, 1);
    Eigen::MatrixXd pp4(1, 1);
    Eigen::MatrixXd pp5(1, 1);
    Eigen::MatrixXd pp6(1, 1);
    Eigen::MatrixXd pp7(1, 1);
    Eigen::MatrixXd pp8(1, 1);
    // laplacian of precursor
    Eigen::VectorXd P1_1 = problem->LP1_matrix * d1_tmp * (nu / Sc);
    Eigen::VectorXd P1_2 = problem->LP2_matrix * d2_tmp * (nu / Sc);
    Eigen::VectorXd P1_3 = problem->LP3_matrix * d3_tmp * (nu / Sc);
    Eigen::VectorXd P1_4 = problem->LP4_matrix * d4_tmp * (nu / Sc);
    Eigen::VectorXd P1_5 = problem->LP5_matrix * d5_tmp * (nu / Sc);
    Eigen::VectorXd P1_6 = problem->LP6_matrix * d6_tmp * (nu / Sc);
    Eigen::VectorXd P1_7 = problem->LP7_matrix * d7_tmp * (nu / Sc);
    Eigen::VectorXd P1_8 = problem->LP8_matrix * d8_tmp * (nu / Sc);
    // algebric term
    Eigen::VectorXd P2_1 = problem->MP1_matrix * d1_tmp * l1;
    Eigen::VectorXd P2_2 = problem->MP2_matrix * d2_tmp * l2;
    Eigen::VectorXd P2_3 = problem->MP3_matrix * d3_tmp * l3;
    Eigen::VectorXd P2_4 = problem->MP4_matrix * d4_tmp * l4;
    Eigen::VectorXd P2_5 = problem->MP5_matrix * d5_tmp * l5;
    Eigen::VectorXd P2_6 = problem->MP6_matrix * d6_tmp * l6;
    Eigen::VectorXd P2_7 = problem->MP7_matrix * d7_tmp * l7;
    Eigen::VectorXd P2_8 = problem->MP8_matrix * d8_tmp * l8;
    // flux source term
    Eigen::MatrixXd fs1(1, 1);
    Eigen::MatrixXd fs2(1, 1);
    Eigen::MatrixXd fs3(1, 1);
    Eigen::MatrixXd fs4(1, 1);
    Eigen::MatrixXd fs5(1, 1);
    Eigen::MatrixXd fs6(1, 1);
    Eigen::MatrixXd fs7(1, 1);
    Eigen::MatrixXd fs8(1, 1);

    for (int i = 0; i < Nphi_flux; i++)
    {
        lf = d_c.transpose() * problem->LF_matrix[i] * c_tmp;
        pf = nsf_c.transpose() * problem->PF_matrix[i] * c_tmp * (1 - btot);
        af = a_c.transpose() * problem->AF_matrix[i] * c_tmp;
        fvecn(i) = lf(0, 0) + pf(0, 0) - af(0,
                                            0) + F3_1(i) + F3_2(i) + F3_3(i) + F3_4(i) + F3_5(i) + F3_6(i) + F3_7(i) + F3_8(
                       i);
    }

    int pfvecn = Nphi_flux;

    for (int i = 0; i < Nphi_prec1; i++)
    {
        int k = i + pfvecn;
        pp1 = a_tmp.transpose() * problem->ST1_matrix[i] * d1_tmp;
        fs1 = nsf_c.transpose() * problem->FS1_matrix[i] * c_tmp * b1;
        fvecn(k) = -pp1(0, 0) + P1_1(i) - P2_1(i) + fs1(0, 0);
    }

    pfvecn += Nphi_prec1;

    for (int i = 0; i < Nphi_prec2; i++)
    {
        int k = i + pfvecn;
        pp2 = a_tmp.transpose() * problem->ST2_matrix[i] * d2_tmp;
        fs2 = nsf_c.transpose() * problem->FS2_matrix[i] * c_tmp * b2;
        fvecn(k) = -pp2(0, 0) + P1_2(i) - P2_2(i) + fs2(0, 0);
    }

    pfvecn += Nphi_prec2;

    for (int i = 0; i < Nphi_prec3; i++)
    {
        int k = i + pfvecn;
        pp3 = a_tmp.transpose() * problem->ST3_matrix[i] * d3_tmp;
        fs3 = nsf_c.transpose() * problem->FS3_matrix[i] * c_tmp * b3;
        fvecn(k) = -pp3(0, 0) + P1_3(i) - P2_3(i) + fs3(0, 0);
    }

    pfvecn += Nphi_prec3;

    for (int i = 0; i < Nphi_prec4; i++)
    {
        int k = i + pfvecn;
        pp4 = a_tmp.transpose() * problem->ST4_matrix[i] * d4_tmp;
        fs4 = nsf_c.transpose() * problem->FS4_matrix[i] * c_tmp * b4;
        fvecn(k) = -pp4(0, 0) + P1_4(i) - P2_4(i) + fs4(0, 0);
    }

    pfvecn += Nphi_prec4;

    for (int i = 0; i < Nphi_prec5; i++)
    {
        int k = i + pfvecn;
        pp5 = a_tmp.transpose() * problem->ST5_matrix[i] * d5_tmp;
        fs5 = nsf_c.transpose() * problem->FS5_matrix[i] * c_tmp * b5;
        fvecn(k) = -pp5(0, 0) + P1_5(i) - P2_5(i) + fs5(0, 0);
    }

    pfvecn += Nphi_prec5;

    for (int i = 0; i < Nphi_prec6; i++)
    {
        int k = i + pfvecn;
        pp6 = a_tmp.transpose() * problem->ST6_matrix[i] * d6_tmp;
        fs6 = nsf_c.transpose() * problem->FS6_matrix[i] * c_tmp * b6;
        fvecn(k) = -pp6(0, 0) + P1_6(i) - P2_6(i) + fs6(0, 0);
    }

    pfvecn += Nphi_prec6;

    for (int i = 0; i < Nphi_prec7; i++)
    {
        int k = i + pfvecn;
        pp7 = a_tmp.transpose() * problem->ST7_matrix[i] * d7_tmp;
        fs7 = nsf_c.transpose() * problem->FS7_matrix[i] * c_tmp * b7;
        fvecn(k) = -pp7(0, 0) + P1_7(i) - P2_7(i) + fs7(0, 0);
    }

    pfvecn += Nphi_prec7;

    for (int i = 0; i < Nphi_prec8; i++)
    {
        int k = i + pfvecn;
        pp8 = a_tmp.transpose() * problem->ST8_matrix[i] * d8_tmp;
        fs8 = nsf_c.transpose() * problem->FS8_matrix[i] * c_tmp * b8;
        fvecn(k) = -pp8(0, 0) + P1_8(i) - P2_8(i) + fs8(0, 0);
    }

    return 0;
}

int newton_msr_n::df(const Eigen::VectorXd& n,
                     Eigen::MatrixXd& fjacn) const
{
    Eigen::NumericalDiff<newton_msr_n> numDiff(*this);
    numDiff.df(n, fjacn);
    return 0;
}

int newton_msr_t::operator()(const Eigen::VectorXd& t,
                             Eigen::VectorXd& fvect) const
{
    Eigen::VectorXd e_tmp(Nphi_T);  //for T
    Eigen::VectorXd f1_tmp(Nphi_dec1); //for dec1
    Eigen::VectorXd f2_tmp(Nphi_dec2); //for dec2
    Eigen::VectorXd f3_tmp(Nphi_dec3); //for dec3
    e_tmp = t.head(Nphi_T);
    int pos = Nphi_T;
    f1_tmp = t.segment(pos, Nphi_dec1);
    pos += Nphi_dec1;
    f2_tmp = t.segment(pos, Nphi_dec2);
    pos += Nphi_dec2;
    f3_tmp = t.segment(pos, Nphi_dec3);
    /// Thermal terms
    // convective term in T_eqn
    Eigen::MatrixXd tt(1, 1);
    // laplacian of T
    Eigen::VectorXd T1 = problem->LT_matrix * e_tmp * nu / Pr;
    // temp flux source
    Eigen::MatrixXd xsf(1, 1);
    // decay heat source term (*dli/cp)
    Eigen::MatrixXd dhs1(1, 1);
    Eigen::MatrixXd dhs2(1, 1);
    Eigen::MatrixXd dhs3(1, 1);
    // convective term in decay heat eq.
    Eigen::MatrixXd dh1(1, 1);
    Eigen::MatrixXd dh2(1, 1);
    Eigen::MatrixXd dh3(1, 1);
    //laplacian of dh
    Eigen::VectorXd DH1_1 = problem->LD1_matrix * f1_tmp * nu / Sc;
    Eigen::VectorXd DH1_2 = problem->LD2_matrix * f2_tmp * nu / Sc;
    Eigen::VectorXd DH1_3 = problem->LD3_matrix * f3_tmp * nu / Sc;
    //algebric term in dh eq
    Eigen::VectorXd DH2_1 = problem->MD1_matrix * f1_tmp * dl1;
    Eigen::VectorXd DH2_2 = problem->MD2_matrix * f2_tmp * dl2;
    Eigen::VectorXd DH2_3 = problem->MD3_matrix * f3_tmp * dl3;
    // flux source in term dh eq (*dbi)
    Eigen::MatrixXd dfs1(1, 1);
    Eigen::MatrixXd dfs2(1, 1);
    Eigen::MatrixXd dfs3(1, 1);

    for (int i = 0; i < Nphi_T; i++)
    {
        tt = a_tmp.transpose() * problem->TS_matrix[i] * e_tmp;
        xsf = txs_c.transpose() * problem->TXS_matrix[i] * c_tmp * ((1 - dbtot) / cp);
        dhs1 = v_c.transpose() * problem->THS1_matrix[i] * f1_tmp * (dl1 / cp);
        dhs2 = v_c.transpose() * problem->THS2_matrix[i] * f2_tmp * (dl2 / cp);
        dhs3 = v_c.transpose() * problem->THS3_matrix[i] * f3_tmp * (dl3 / cp);
        fvect(i) = -tt(0, 0) + T1(i) + xsf(0, 0) + dhs1(0, 0) + dhs2(0, 0) + dhs3(0, 0);
    }

    int pfvect = Nphi_T;

    for (int i = 0; i < Nphi_dec1; i++)
    {
        int k = i + pfvect;
        dh1 = a_tmp.transpose() * problem->SD1_matrix[i] * f1_tmp;
        dfs1 = sp_c.transpose() * problem->DFS1_matrix[i] * c_tmp * db1;
        fvect(k) = -dh1(0, 0) + DH1_1(i) - DH2_1(i) + dfs1(0, 0);
    }

    pfvect += Nphi_dec1;

    for (int i = 0; i < Nphi_dec2; i++)
    {
        int k = i + pfvect;
        dh2 = a_tmp.transpose() * problem->SD2_matrix[i] * f2_tmp;
        dfs2 = sp_c.transpose() * problem->DFS2_matrix[i] * c_tmp * db2;
        fvect(k) = -dh2(0, 0) + DH1_2(i) - DH2_2(i) + dfs2(0, 0);
    }

    pfvect += Nphi_dec2;

    for (int i = 0; i < Nphi_dec3; i++)
    {
        int k = i + pfvect;
        dh3 = a_tmp.transpose() * problem->SD3_matrix[i] * f3_tmp;
        dfs3 = sp_c.transpose() * problem->DFS3_matrix[i] * c_tmp * db3;
        fvect(k) = -dh3(0, 0) + DH1_3(i) - DH2_3(i) + dfs3(0, 0);
    }

    for (int i = 0; i < N_BCt; i++)
    {
        fvect(i) = t(i) - BCt(i);
    }

    return 0;
}

int newton_msr_t::df(const Eigen::VectorXd& t,
                     Eigen::MatrixXd& fjact) const
{
    Eigen::NumericalDiff<newton_msr_t> numDiff(*this);
    numDiff.df(t, fjact);
    return 0;
}


void reducedMSR::solveOnline(Eigen::MatrixXd vel_now, Eigen::MatrixXd temp_now,
                             Eigen::VectorXd mu_online)
{
    Info << "\n Starting online stage...\n" << endl;
    y.resize(Nphi_u + Nphi_p, 1);
    y.setZero();
    w.resize(Nphi_flux + Nphi_prec1 + Nphi_prec2 + Nphi_prec3 + Nphi_prec4 +
             Nphi_prec5 + Nphi_prec6 + Nphi_prec7 + Nphi_prec8, 1);
    w.setZero();
    z.resize(Nphi_T + Nphi_dec1 + Nphi_dec2 + Nphi_dec3, 1);
    z.setZero();

    for (int j = 0; j < N_BC; j++)
    {
        y(j) = vel_now(j, 0);
    }

    for (int j = 0; j < N_BCt; j++)
    {
        z(j) = temp_now(j, 0);
    }

    for (int i = 0; i < Nphi_const; i++)
    {
        newton_object_n.d_c(i) = problem->rbfsplines_D[i]->eval(mu_online);
        newton_object_n.nsf_c(i) = problem->rbfsplines_NSF[i]->eval(mu_online);
        newton_object_n.a_c(i) = problem->rbfsplines_A[i]->eval(mu_online);
        newton_object_t.v_c(i) = problem->rbfsplines_v[i]->eval(mu_online);
        newton_object_t.sp_c(i) = problem->rbfsplines_SP[i]->eval(mu_online);
        newton_object_t.txs_c(i) = problem->rbfsplines_TXS[i]->eval(mu_online);
    }

    online_solution_fd.resize(1);
    online_solution_n.resize(1);
    online_solution_t.resize(1);
    online_solution_C.resize(1);
    online_solution_fd[0].resize(y.rows() + 1, 1);
    online_solution_n[0].resize(w.rows() + 1, 1);
    online_solution_t[0].resize(z.rows() + 1, 1);
    online_solution_C[0].resize(6 * Nphi_const + 1, 1);
    Color::Modifier red(Color::FG_RED);
    Color::Modifier green(Color::FG_GREEN);
    Color::Modifier def(Color::FG_DEFAULT);
    Eigen::HybridNonLinearSolver<newton_msr_fd> hnls_fd(newton_object_fd);
    Eigen::HybridNonLinearSolver<newton_msr_n> hnls_n(newton_object_n);
    Eigen::HybridNonLinearSolver<newton_msr_t> hnls_t(newton_object_t);
    newton_object_fd.BC.resize(N_BC);

    for (int j = 0; j < N_BC; j++)
    {
        newton_object_fd.BC(j) = vel_now(j, 0);
    }

    newton_object_fd.nu = nu;
    newton_object_n.nu = nu;
    newton_object_n.iv = iv;
    newton_object_n.l1 = l1; //lambda-ith
    newton_object_n.l2 = l2;
    newton_object_n.l3 = l3;
    newton_object_n.l4 = l4;
    newton_object_n.l5 = l5;
    newton_object_n.l6 = l6;
    newton_object_n.l7 = l7;
    newton_object_n.l8 = l8;
    newton_object_n.b1 = b1; //beta-ith
    newton_object_n.b2 = b2;
    newton_object_n.b3 = b3;
    newton_object_n.b4 = b4;
    newton_object_n.b5 = b5;
    newton_object_n.b6 = b6;
    newton_object_n.b7 = b7;
    newton_object_n.b8 = b8;
    newton_object_n.btot = btot;
    newton_object_t.nu = nu;
    newton_object_t.cp = cp;
    newton_object_t.dl1 = dl1;
    newton_object_t.dl2 = dl2;
    newton_object_t.dl3 = dl3;
    newton_object_t.db1 = db1;
    newton_object_t.db2 = db2;
    newton_object_t.db3 = db3;
    newton_object_t.dbtot = dbtot;
    newton_object_t.BCt.resize(N_BCt);

    for (int j = 0; j < N_BCt; j++)
    {
        newton_object_t.BCt(j) = temp_now(j, 0);
    }

    hnls_fd.solve(y);
    Eigen::VectorXd res_fd(y);
    newton_object_n.a_tmp = y.head(Nphi_u);
    hnls_n.solve(w);
    Eigen::VectorXd res_n(w);
    newton_object_t.a_tmp = y.head(Nphi_u);
    newton_object_t.c_tmp = w.head(Nphi_flux);
    hnls_t.solve(z);
    Eigen::VectorXd res_t(z);
    newton_object_fd.operator()(y, res_fd);
    newton_object_n.operator()(w, res_n);
    newton_object_t.operator()(z, res_t);
    std::cout << "################## Online solve N° " <<  count_online_solve <<
              " ##################" << std::endl;

    if (res_fd.norm() / y.norm() < 1e-5)
    {
        std::cout << green << "|F_fd(x)| = " << res_fd.norm() / y.norm() <<
                  " - Minimun reached in " << hnls_fd.iter << " iterations " << def << std::endl
                  << std::endl;
    }
    else
    {
        std::cout << red << "|F_fd(x)| = " << res_fd.norm() / y.norm() <<
                  " - Minimun reached in " << hnls_fd.iter << " iterations " << def << std::endl
                  << std::endl;
    }

    if (res_n.norm() / w.norm() < 1e-5)
    {
        std::cout << green << "|F_n(x)| = " << res_n.norm() / w.norm() <<
                  " - Minimun reached in " << hnls_n.iter << " iterations " << def << std::endl <<
                  std::endl;
    }
    else
    {
        std::cout << red << "|F_n(x)| = " << res_n.norm() / w.norm() <<
                  " - Minimun reached in " << hnls_n.iter << " iterations " << def << std::endl <<
                  std::endl;
    }

    if (res_t.norm() / z.norm() < 1e-5)
    {
        std::cout << green << "|F_t(x)| = " << res_t.norm() / z.norm() <<
                  " - Minimun reached in " << hnls_t.iter << " iterations " << def << std::endl <<
                  std::endl;
    }
    else
    {
        std::cout << red << "|F_t(x)| = " << res_t.norm() / z.norm() <<
                  " - Minimun reached in " << hnls_t.iter << " iterations " << def << std::endl <<
                  std::endl;
    }

    online_solution_fd[0](0, 0) = count_online_solve;
    online_solution_fd[0].col(0).tail(y.rows()) = y;
    online_solution_n[0](0, 0) = count_online_solve;
    online_solution_n[0].col(0).tail(w.rows()) = w;
    online_solution_t[0](0, 0) = count_online_solve;
    online_solution_t[0].col(0).tail(z.rows()) = z;
    online_solution_C[0](0, 0) = count_online_solve;
    int pos_c = 1;
    online_solution_C[0].col(0).segment(pos_c, Nphi_const) = newton_object_t.v_c;
    pos_c += Nphi_const;
    online_solution_C[0].col(0).segment(pos_c, Nphi_const) = newton_object_n.d_c;
    pos_c += Nphi_const;
    online_solution_C[0].col(0).segment(pos_c, Nphi_const) = newton_object_n.nsf_c;
    pos_c += Nphi_const;
    online_solution_C[0].col(0).segment(pos_c, Nphi_const) = newton_object_n.a_c;
    pos_c += Nphi_const;
    online_solution_C[0].col(0).segment(pos_c, Nphi_const) = newton_object_t.sp_c;
    pos_c += Nphi_const;
    online_solution_C[0].col(0).segment(pos_c, Nphi_const) = newton_object_t.txs_c;
    ITHACAstream::exportMatrix(online_solution_fd, "red_coeff_fd", "matlab",
                               "./ITHACAoutput/red_coeff_fd");
    ITHACAstream::exportMatrix(online_solution_n, "red_coeff_n", "matlab",
                               "./ITHACAoutput/red_coeff_n");
    ITHACAstream::exportMatrix(online_solution_t, "red_coeff_t", "matlab",
                               "./ITHACAoutput/red_coeff_t");
    ITHACAstream::exportMatrix(online_solution_C, "red_coeff_C", "matlab",
                               "./ITHACAoutput/red_coeff_C");
    count_online_solve += 1;
}

void reducedMSR::reconstructAP(fileName folder, int printevery)
{
    recall = true;
    mkDir(folder);
    ITHACAutilities::createSymLink(folder);
    reconstruct_fd(folder, printevery);
    reconstruct_n(folder, printevery);
    reconstruct_C(folder, printevery);
    reconstruct_t(folder, printevery);
    return;
}

void reducedMSR::reconstruct_fd(fileName folder, int printevery)
{
    if (recall == false)
    {
        mkDir(folder);
        ITHACAutilities::createSymLink(folder);
    }

    Info << "Reconstructing online solution | fluid-dynamics" << endl;
    int counter = 0;
    int nextwrite = 0;
    int counter2 = 1;

    for (int i = 0; i < online_solution_fd.size(); i++)
    {
        if (counter == nextwrite)
        {
            volVectorField U_rec("U", Umodes[0] * 0);

            for (int j = 0; j < Nphi_u; j++)
            {
                U_rec += Umodes[j] * online_solution_fd[i](j + 1, 0);
            }

            ITHACAstream::exportSolution(U_rec,  name(counter2), folder);
            volScalarField P_rec("p", Pmodes[0] * 0);

            for (int j = 0; j < Nphi_p; j++)
            {
                P_rec += Pmodes[j] * online_solution_fd[i](j + Nphi_u + 1, 0);
            }

            ITHACAstream::exportSolution(P_rec, name(counter2), folder);
            nextwrite += printevery;
            counter2 ++;
            UREC.append(U_rec);
            PREC.append(P_rec);
        }

        counter++;
    }

    Info << "End" << endl;
}

void reducedMSR::reconstruct_n(fileName folder, int printevery)
{
    if (recall == false)
    {
        mkDir(folder);
        ITHACAutilities::createSymLink(folder);
    }

    Info << "Reconstructing online solution | neutronics" << endl;
    int counter = 0;
    int nextwrite = 0;
    int counter2 = 1;

    for (int i = 0; i < online_solution_fd.size(); i++)
    {
        if (counter == nextwrite)
        {
            volScalarField Flux_rec("flux", Fluxmodes[0] * 0);

            for (int j = 0; j < Nphi_flux; j++)
            {
                Flux_rec += Fluxmodes[j] * online_solution_n[i](j + 1, 0);
            }

            ITHACAstream::exportSolution(Flux_rec,  name(counter2), folder);
            int pos = Nphi_flux;
            volScalarField Prec1_rec("prec1", Prec1modes[0] * 0);

            for (int j = 0; j < Nphi_prec1; j++)
            {
                Prec1_rec += Prec1modes[j] * online_solution_n[i](j + pos + 1, 0);
            }

            ITHACAstream::exportSolution(Prec1_rec, name(counter2), folder);
            pos += Nphi_prec1;
            volScalarField Prec2_rec("prec2", Prec2modes[0] * 0);

            for (int j = 0; j < Nphi_prec2; j++)
            {
                Prec2_rec += Prec2modes[j] * online_solution_n[i](j + pos + 1, 0);
            }

            ITHACAstream::exportSolution(Prec2_rec, name(counter2), folder);
            pos += Nphi_prec2;
            volScalarField Prec3_rec("prec3", Prec3modes[0] * 0);

            for (int j = 0; j < Nphi_prec3; j++)
            {
                Prec3_rec += Prec3modes[j] * online_solution_n[i](j + pos + 1, 0);
            }

            ITHACAstream::exportSolution(Prec3_rec, name(counter2), folder);
            pos += Nphi_prec3;
            volScalarField Prec4_rec("prec4", Prec4modes[0] * 0);

            for (int j = 0; j < Nphi_prec4; j++)
            {
                Prec4_rec += Prec4modes[j] * online_solution_n[i](j + pos + 1, 0);
            }

            ITHACAstream::exportSolution(Prec4_rec, name(counter2), folder);
            pos += Nphi_prec4;
            volScalarField Prec5_rec("prec5", Prec5modes[0] * 0);

            for (int j = 0; j < Nphi_prec5; j++)
            {
                Prec5_rec += Prec5modes[j] * online_solution_n[i](j + pos + 1, 0);
            }

            ITHACAstream::exportSolution(Prec5_rec, name(counter2), folder);
            pos += Nphi_prec5;
            volScalarField Prec6_rec("prec6", Prec6modes[0] * 0);

            for (int j = 0; j < Nphi_prec6; j++)
            {
                Prec6_rec += Prec6modes[j] * online_solution_n[i](j + pos + 1, 0);
            }

            ITHACAstream::exportSolution(Prec6_rec, name(counter2), folder);
            pos += Nphi_prec6;
            volScalarField Prec7_rec("prec7", Prec7modes[0] * 0);

            for (int j = 0; j < Nphi_prec7; j++)
            {
                Prec7_rec += Prec7modes[j] * online_solution_n[i](j + pos + 1, 0);
            }

            ITHACAstream::exportSolution(Prec7_rec, name(counter2), folder);
            pos += Nphi_prec7;
            volScalarField Prec8_rec("prec8", Prec8modes[0] * 0);

            for (int j = 0; j < Nphi_prec8; j++)
            {
                Prec8_rec += Prec8modes[j] * online_solution_n[i](j + pos + 1, 0);
            }

            ITHACAstream::exportSolution(Prec8_rec, name(counter2), folder);
            nextwrite += printevery;
            counter2 ++;
            FLUXREC.append(Flux_rec);
            PREC1REC.append(Prec1_rec);
            PREC2REC.append(Prec2_rec);
            PREC3REC.append(Prec3_rec);
            PREC4REC.append(Prec4_rec);
            PREC5REC.append(Prec5_rec);
            PREC6REC.append(Prec6_rec);
            PREC7REC.append(Prec7_rec);
            PREC8REC.append(Prec8_rec);
        }

        counter++;
    }

    Info << "End" << endl;
}

void reducedMSR::reconstruct_t(fileName folder, int printevery)
{
    if (recall == false)
    {
        mkDir(folder);
        ITHACAutilities::createSymLink(folder);
    }

    Info << "Reconstructing online solution | thermal" << endl;
    int counter = 0;
    int nextwrite = 0;
    int counter2 = 1;
    dimensionedScalar decLam1("decLam1", dimensionSet(0, 0, -1, 0, 0, 0, 0), dl1);
    dimensionedScalar decLam2("decLam2", dimensionSet(0, 0, -1, 0, 0, 0, 0), dl2);
    dimensionedScalar decLam3("decLam3", dimensionSet(0, 0, -1, 0, 0, 0, 0), dl3);

    for (int i = 0; i < online_solution_t.size(); i++)
    {
        if (counter == nextwrite)
        {
            volScalarField T_rec("T", Tmodes[0] * 0);

            for (int j = 0; j < Nphi_T; j++)
            {
                T_rec += Tmodes[j] * online_solution_t[i](j + 1, 0);
            }

            ITHACAstream::exportSolution(T_rec,  name(counter2), folder);
            int pos = Nphi_T;
            volScalarField PowerDens_rec("powerDens", Dec1modes[0] * 0 * decLam1);
            volScalarField Dec1_rec("dec1", Dec1modes[0] * 0);

            for (int j = 0; j < Nphi_dec1; j++)
            {
                Dec1_rec += Dec1modes[j] * online_solution_t[i](j + pos + 1, 0);
                PowerDens_rec += Dec1modes[j] * online_solution_t[i](j + pos + 1, 0) * decLam1;
            }

            ITHACAstream::exportSolution(Dec1_rec, name(counter2), folder);
            pos += Nphi_dec1;
            volScalarField Dec2_rec("dec2", Dec2modes[0] * 0);

            for (int j = 0; j < Nphi_dec2; j++)
            {
                Dec2_rec += Dec2modes[j] * online_solution_t[i](j + pos + 1, 0);
                PowerDens_rec += Dec2modes[j] * online_solution_t[i](j + pos + 1, 0) * decLam2;
            }

            ITHACAstream::exportSolution(Dec2_rec, name(counter2), folder);
            pos += Nphi_dec2;
            volScalarField Dec3_rec("dec3", Dec3modes[0] * 0);

            for (int j = 0; j < Nphi_dec3; j++)
            {
                Dec3_rec += Dec3modes[j] * online_solution_t[i](j + pos + 1, 0);
                PowerDens_rec += Dec3modes[j] * online_solution_t[i](j + pos + 1, 0) * decLam3;
            }

            ITHACAstream::exportSolution(Dec3_rec, name(counter2), folder);
            PowerDens_rec += (1 - dbtot) * SPREC[counter2 - 1] * FLUXREC[counter2 - 1];
            ITHACAstream::exportSolution(PowerDens_rec, name(counter2), folder);
            nextwrite += printevery;
            counter2 ++;
            TREC.append(T_rec);
            DEC1REC.append(Dec1_rec);
            DEC2REC.append(Dec2_rec);
            DEC3REC.append(Dec3_rec);
            POWERDENSREC.append(PowerDens_rec);
        }

        counter++;
    }

    Info << "End" << endl;
}

void reducedMSR::reconstruct_C(fileName folder, int printevery)
{
    if (recall == false)
    {
        mkDir(folder);
        ITHACAutilities::createSymLink(folder);
    }

    Info << "Reconstructing temperature changing constants" << endl;
    int counter = 0;
    int nextwrite = 0;
    int counter2 = 1;

    for (int i = 0; i < online_solution_C.size(); i++)
    {
        if (counter == nextwrite)
        {
            volScalarField v_rec("v", vmodes[0] * 0);

            for (int j = 0; j < Nphi_const; j++)
            {
                v_rec += vmodes[j] * online_solution_C[i](j + 1, 0);
            }

            ITHACAstream::exportSolution(v_rec,  name(counter2), folder);
            int pos = Nphi_const;
            volScalarField D_rec("D", Dmodes[0] * 0);

            for (int j = 0; j < Nphi_const; j++)
            {
                D_rec += Dmodes[j] * online_solution_C[i](j + pos + 1, 0);
            }

            ITHACAstream::exportSolution(D_rec, name(counter2), folder);
            pos += Nphi_const;
            volScalarField NSF_rec("NSF", NSFmodes[0] * 0);

            for (int j = 0; j < Nphi_const; j++)
            {
                NSF_rec += NSFmodes[j] * online_solution_C[i](j + pos + 1, 0);
            }

            ITHACAstream::exportSolution(NSF_rec, name(counter2), folder);
            pos += Nphi_const;
            volScalarField A_rec("A", Amodes[0] * 0);

            for (int j = 0; j < Nphi_const; j++)
            {
                A_rec += Amodes[j] * online_solution_C[i](j + pos + 1, 0);
            }

            ITHACAstream::exportSolution(A_rec, name(counter2), folder);
            pos += Nphi_const;
            volScalarField SP_rec("SP", SPmodes[0] * 0);

            for (int j = 0; j < Nphi_const; j++)
            {
                SP_rec += SPmodes[j] * online_solution_C[i](j + pos + 1, 0);
            }

            ITHACAstream::exportSolution(SP_rec, name(counter2), folder);
            pos += Nphi_const;
            volScalarField TXS_rec("TXS", TXSmodes[0] * 0);

            for (int j = 0; j < Nphi_const; j++)
            {
                TXS_rec += TXSmodes[j] * online_solution_C[i](j + pos + 1, 0);
            }

            ITHACAstream::exportSolution(TXS_rec, name(counter2), folder);
            std::ofstream of(folder + "/" + name(counter2) + "/" + name(
                                 online_solution_C[i](0)));
            nextwrite += printevery;
            counter2 ++;
            vREC.append(v_rec);
            DREC.append(D_rec);
            NSFREC.append(NSF_rec);
            AREC.append(A_rec);
            SPREC.append(SP_rec);
            TXSREC.append(TXS_rec);
        }

        counter++;
    }

    Info << "End" << endl;
}


void reducedMSR::loadConstants(msrProblem* problem)
{
    nu = problem->_nu().value();
    iv = problem->_IV1().value();
    l1 = problem->_lam1().value();
    l2 = problem->_lam2().value();
    l3 = problem->_lam3().value();
    l4 = problem->_lam4().value();
    l5 = problem->_lam5().value();
    l6 = problem->_lam6().value();
    l7 = problem->_lam7().value();
    l8 = problem->_lam8().value();
    b1 = problem->_beta1().value();
    b2 = problem->_beta2().value();
    b3 = problem->_beta3().value();
    b4 = problem->_beta4().value();
    b5 = problem->_beta5().value();
    b6 = problem->_beta6().value();
    b7 = problem->_beta7().value();
    b8 = problem->_beta8().value();
    btot = problem->_betaTot().value();
    cp = problem->_CpRef().value();
    dl1 = problem->_decLam1().value();
    dl2 = problem->_decLam2().value();
    dl3 = problem->_decLam3().value();
    db1 = problem->_decBeta1().value();
    db2 = problem->_decBeta2().value();
    db3 = problem->_decBeta3().value();
    dbtot = problem->_decbetaTot().value();
    Pr = problem->_Pr().value();
    Sc = problem->_Sc().value();
}

void reducedMSR::clearFields()
{
    UREC.clear();
    PREC.clear();
    FLUXREC.clear();
    PREC1REC.clear();
    PREC2REC.clear();
    PREC3REC.clear();
    PREC4REC.clear();
    PREC5REC.clear();
    PREC6REC.clear();
    PREC7REC.clear();
    PREC8REC.clear();
    TREC.clear();
    DEC1REC.clear();
    DEC2REC.clear();
    DEC3REC.clear();
    POWERDENSREC.clear();
    vREC.clear();
    DREC.clear();
    NSFREC.clear();
    AREC.clear();
    SPREC.clear();
    TXSREC.clear();
}

