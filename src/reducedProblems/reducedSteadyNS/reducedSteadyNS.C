/*---------------------------------------------------------------------------*\
Copyright (C) 2017 by the ITHACA-FV authors

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

/// \file
/// Source file of the reducedSteadyNS class

#include "reducedSteadyNS.H"

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructor
reducedSteadyNS::reducedSteadyNS()
{

}

reducedSteadyNS::reducedSteadyNS(steadyNS& problem, word tipo)
{
	B_matrix = problem.B_matrix;
	C_matrix = problem.C_matrix;
	K_matrix = problem.K_matrix;

	N_BC = problem.inletIndex.rows();

	if (tipo == "PPE")
	{
		D_matrix = problem.D_matrix;
		G_matrix = problem.G_matrix;
	}
	if (tipo == "SUP")
	{
		P_matrix = problem.P_matrix;
	}

	Nphi_u = B_matrix.rows();
    Nphi_p = K_matrix.cols();

	for (label k = 0; k < problem.liftfield.size(); k++)
	{
		Umodes.append(problem.liftfield[k]);
	}
	for (label k = 0; k < problem.NUmodes; k++)
	{
		Umodes.append(problem.Umodes[k]);
	}
	for (label k = 0; k < problem.NSUPmodes; k++)
	{
		Umodes.append(problem.supmodes[k]);
	}

	for (label k = 0; k < problem.NPmodes; k++)
	{
		Pmodes.append(problem.Pmodes[k]);
	}

	newton_object = newton_steadyNS(Nphi_u + Nphi_p , Nphi_u + Nphi_p, problem, *this);
}

int newton_steadyNS::operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const
{
    Eigen::VectorXd a_tmp(Nphi_u);
    Eigen::VectorXd b_tmp(Nphi_p);
    a_tmp = x.head(Nphi_u);
    b_tmp = x.tail(Nphi_p);
    // Convective term
    Eigen::MatrixXd cc(1, 1);
    // Mom Term
    Eigen::VectorXd M1 = B_matrix * a_tmp;
    // Gradient of pressure
    Eigen::VectorXd M2 = K_matrix * b_tmp;
    // Pressure Term
    Eigen::VectorXd M3 = P_matrix * a_tmp;

    for (label i = 0; i < Nphi_u; i++)
    {
        cc = a_tmp.transpose() * C_matrix[i] * a_tmp;
        fvec(i) = M1(i) - cc(0, 0) - M2(i);
    }
    for (label j = 0; j < Nphi_p; j++)
    {
        label k = j + Nphi_u;
        fvec(k) = M3(j);
    }
    for (label j = 0; j < N_BC; j++)
    {
        fvec(j) = x(j) - BC(j);
    }
    return 0;
}


int newton_steadyNS::df(const Eigen::VectorXd &x,  Eigen::MatrixXd &fjac) const
{
    Eigen::NumericalDiff<newton_steadyNS> numDiff(*this);
    numDiff.df(x, fjac);

    return 0;
}


// * * * * * * * * * * * * * * * Solve Functions  * * * * * * * * * * * * * //

void reducedSteadyNS::solveOnline_PPE(Eigen::MatrixXd vel_now)
{
	y.resize(Nphi_u + Nphi_p, 1);
	y_old.resize(Nphi_u + Nphi_p, 1);
	y_old.setZero();

	B_matrix = B_matrix * nu;

	Eigen::VectorXd x(Nphi_u + Nphi_p);
	x.setZero();

	for (label j = 0; j < N_BC; j++)
	{
		x(j) = vel_now(j, 0);
	}
	for (label j = 0; j < N_BC; j++)
	{
		y_old(j) = vel_now(j, 0);
	}
	// y = NewtonRaphson_PPE(y_old);
}

void reducedSteadyNS::solveOnline_sup(Eigen::MatrixXd vel_now)
{

	y.resize(Nphi_u + Nphi_p, 1);
	B_matrix = B_matrix * nu;
	y.setZero();

	for (label j = 0; j < N_BC; j++)
	{
		y(j) = vel_now(j, 0);
	}
	//Eigen::LevenbergMarquardt<newton_steadyNS> lm(newton_object);
	Color::Modifier red(Color::FG_RED);
	Color::Modifier green(Color::FG_GREEN);
	Color::Modifier def(Color::FG_DEFAULT);
	Eigen::HybridNonLinearSolver<newton_steadyNS> hnls(newton_object);
	newton_object.BC.resize(N_BC);

	for (label j = 0; j < N_BC; j++)
    {
        newton_object.BC(j) = vel_now(j, 0);
    }

	hnls.solve(y);
	//lm.minimize(y);
	Eigen::VectorXd res(y);
	newton_object.operator()(y, res);

	std::cout << "################## Online solve N° " << count_online_solve << " ##################" << std::endl;
	std::cout << "Solving for the parameter: " << vel_now << std::endl;
	if (res.norm() < 1e-5)
	{
		std::cout << green << "|F(x)| = " << res.norm() << " - Minimun reached in " << hnls.iter << " iterations " << def << std::endl << std::endl;
	}
	else
	{
		std::cout << red << "|F(x)| = " << res.norm() << " - Minimun reached in " << hnls.iter << " iterations " << def << std::endl << std::endl;
	}
	count_online_solve += 1;

	//y = NewtonRaphson_sup(y_old, 1, 1e-7);
}


// * * * * * * * * * * * * * * * Jacobian Evaluation  * * * * * * * * * * * * * //

void reducedSteadyNS::reconstruct_PPE(steadyNS & problem, fileName folder, int printevery)
{
	mkDir(folder);
	system("ln -s ../../constant " + folder + "/constant");
	system("ln -s ../../0 " + folder + "/0");
	system("ln -s ../../system " + folder + "/system");

	int counter = 0;
	int nextwrite = 0;
	int counter2 = 1;

	for (label i = 0; i < online_solution.size(); i++)
	{
		if (counter == nextwrite)
		{
			volVectorField U_rec("U_rec", Umodes[0] * 0);
			for (label j = 0; j < Nphi_u; j++)
			{
				U_rec += Umodes[j] * online_solution[i](j + 1, 0);
			}
			problem.exportSolution(U_rec, name(online_solution[i](0, 0)), folder);
			//problem.exportSolution(U_rec,  name(counter2), folder);

			volScalarField P_rec("P_rec", Pmodes[0] * 0);
			for (label j = 0; j < Nphi_p; j++)
			{
				P_rec += Pmodes[j] * online_solution[i](j + Nphi_u + 1, 0);
			}
			problem.exportSolution(P_rec, name(online_solution[i](0, 0)), folder);
			//problem.exportSolution(P_rec, name(counter2), folder);
			nextwrite += printevery;
			counter2 ++;
		}
		counter++;
	}
}

void reducedSteadyNS::reconstruct_sup(steadyNS & problem, fileName folder, int printevery)
{
	mkDir(folder);
	system("ln -s ../../constant " + folder + "/constant");
	system("ln -s ../../0 " + folder + "/0");
	system("ln -s ../../system " + folder + "/system");

	int counter = 0;
	int nextwrite = 0;
	int counter2 = 1;

	for (label i = 0; i < online_solution.size(); i++)
	{
		if (counter == nextwrite)
		{
			volVectorField U_rec("U_rec", Umodes[0] * 0);
			for (label j = 0; j < Nphi_u; j++)
			{
				U_rec += Umodes[j] * online_solution[i](j + 1, 0);
			}
			problem.exportSolution(U_rec, name(online_solution[i](0, 0)), folder);
			//problem.exportSolution(U_rec,  name(counter2), folder);

			volScalarField P_rec("P_rec", Pmodes[0] * 0);
			for (label j = 0; j < Nphi_p; j++)
			{
				P_rec += Pmodes[j] * online_solution[i](j + Nphi_u + 1, 0);
			}
			problem.exportSolution(P_rec, name(online_solution[i](0, 0)), folder);
			//problem.exportSolution(P_rec, name(counter2), folder);
			nextwrite += printevery;
			counter2 ++;

			UREC.append(U_rec);
			PREC.append(P_rec);
		}
		counter++;
	}
}

double reducedSteadyNS::inf_sup_constant()
{
	double a;
	Eigen::VectorXd sup(Nphi_u);
	Eigen::VectorXd inf(Nphi_p);

	
	for (int i = 0; i < Nphi_p; i++)
	{
		for (int j = 0; j < Nphi_u; j++)
		{
			sup(j) = fvc::domainIntegrate(fvc::div(Umodes[j])*Pmodes[i]).value()/ITHACAutilities::H1seminorm(Umodes[j])/ITHACAutilities::L2norm(Pmodes[i]);
		}
		inf(i) = sup.maxCoeff();
	}
	Info << "output" << endl;
	a = inf.minCoeff();
	return a;
}


// ************************************************************************* //

