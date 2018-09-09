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

Description
    Example of NS-Stokes Reduction Problem

\*---------------------------------------------------------------------------*/
#include "steadyNS.H"
#include "steadyNSturb.H"
#include "Foam2Eigen.H"

#include "ITHACAutilities.H"
#include "reductionProblem.H"
#include "forces.H"
#include "forceCoeffs.H"
#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include "ITHACAstream.H"
#include "ITHACAPOD.H"
#include "EigenFunctions.H"
#include <chrono>
#include <Eigen/SVD>
#include <GenEigsSolver.h>
#include <Eigen/SparseLU>



#include "reducedSteadyNS.H"
#include "reducedSteadyNSturb.H"
#include <Eigen/Dense>
#include "forces.H"
#include "IOmanip.H"


class tutorial06 : public steadyNSturb
{
public:
	explicit tutorial06(int argc, char *argv[])
	:
	steadyNSturb(argc, argv),
	U(_U()),
	p(_p()),
	nut(_nut())
	{}

	// Relevant Fields
	volVectorField& U;
	volScalarField& p;
	volScalarField& nut;

	void offlineSolve()
	{
		Info << "here 11" << endl;
		Vector<double> inl(0, 0, 0);
		List<scalar> mu_now(2);
		if (offline)
		{
			ITHACAstream::read_fields(Ufield, U, "./ITHACAoutput/Offline/");
			ITHACAstream::read_fields(Pfield, p, "./ITHACAoutput/Offline/");
			ITHACAstream::read_fields(nutFields, nut, "./ITHACAoutput/Offline/");
			if(Ufield.size()<mu.rows())
			{

				Vector<double> Uinl(0, 0, 0);
				label BCind = 0;
				for (label i = Ufield.size(); i < mu.rows(); i++)
				{
					Uinl[0] = mu(i, 0);
					Uinl[1] = mu(i, 1);
					assignBC(U, BCind, Uinl);
					assignIF(U, Uinl);
					counter = Ufield.size() +1;
					truthSolve("./ITHACAoutput/Offline/");
				}
			}
		}
		else
		{
			Vector<double> Uinl(0, 0, 0);
			label BCind = 0;
			for (label i = 0; i < mu.rows(); i++)
			{
				Uinl[0] = mu(i, 0);
				Uinl[1] = mu(i, 1);
				assignBC(U, BCind, Uinl);
				assignIF(U, Uinl);
				truthSolve("./ITHACAoutput/Offline/");
			}

		}
	}


	void offlineSolve(Eigen::MatrixXd par, fileName folder)
	{
		Vector<double> inl(0, 0, 0);
		List<scalar> mu_now(2);
		
		
		Vector<double> Uinl(0, 0, 0);
		label BCind = 0;
		for (label i = 0; i < par.rows(); i++)
		{
			Uinl[0] = par(i, 0);
			Uinl[1] = par(i, 1);
			assignBC(U, BCind, Uinl);
			assignIF(U, Uinl);
			truthSolve(folder);
		}
	}

	void truthSolve(fileName folder)
	{
		Time& runTime = _runTime();
		fvMesh& mesh = _mesh();
		volScalarField& p = _p();
		volVectorField& U = _U();
		surfaceScalarField& phi = _phi();
		fv::options& fvOptions = _fvOptions();
		simpleControl& simple = _simple();
		IOMRFZoneList& MRF = _MRF();
		singlePhaseTransportModel& laminarTransport = _laminarTransport();
#include "NLsolve.H"
		exportSolution(U, name(counter), folder);
		exportSolution(p, name(counter), folder);
		volScalarField _nut(turbulence->nut());
	//volScalarField nuTilda = mesh.lookupObject<volScalarField>("nuTilda");
		exportSolution(_nut, name(counter), folder);
	//exportSolution(nuTilda, name(counter), "./ITHACAoutput/Offline/");
		Ufield.append(U);
		Pfield.append(p);
		nutFields.append(_nut);
		counter++;
		bool notconverged = 1;
	}
};

int main(int argc, char *argv[])
{
	//return 0;
	tutorial06 example(argc, argv);

	word par_offline("./par_offline"); // the offline samples.
	word par_new("./par_online"); // the samples which will be used in the online stage for checking the reduced order model.
	Eigen::MatrixXd par_online = ITHACAstream::readMatrix(par_new);
	
	example.mu = ITHACAstream::readMatrix(par_offline);
	example.inletIndex.resize(2, 2);
	example.inletIndex(0, 0) = 0;
	example.inletIndex(0, 1) = 0;
	example.inletIndex(1, 0) = 0;
	example.inletIndex(1, 1) = 1;

	ITHACAparameters para;
	int NmodesU = para.ITHACAdict->lookupOrDefault<int>("NmodesU", 5);
	int NmodesP = para.ITHACAdict->lookupOrDefault<int>("NmodesP", 5);
	int NmodesSUP = para.ITHACAdict->lookupOrDefault<int>("NmodesSUP", 5);
	int NmodesNUT = para.ITHACAdict->lookupOrDefault<int>("NmodesNUT", 5);
	int NmodesProject = para.ITHACAdict->lookupOrDefault<int>("NmodesProject", 5);
	int NmodesMatrixRec = para.ITHACAdict->lookupOrDefault<int>("NmodesMatrixRec", 5);

	example.offlineSolve();

	// Solve the supremizer problem
	example.solvesupremizer();
	//example.liftSolve();
	
	ITHACAstream::read_fields(example.liftfield, example.U, "./lift/");
	
	example.computeLift(example.Ufield, example.liftfield, example.Uomfield);


	ITHACAstream::exportFields(example.Uomfield,"./ITHACAoutput/Offline", "Uofield");

	ITHACAPOD::getModes(example.Uomfield, example.Umodes, example.podex,example.supex,0,NmodesProject);
	ITHACAPOD::getModes(example.Pfield, example.Pmodes, example.podex,example.supex,0,NmodesProject);
	ITHACAPOD::getModes(example.nutFields, example.nuTmodes, example.podex, example.supex,0,NmodesProject);
	ITHACAPOD::getModes(example.supfield, example.supmodes, example.podex, example.supex, 1,NmodesProject);
	

	PtrList<volVectorField> rec_field;
	PtrList<volVectorField> lift_and_modes;
	PtrList<volScalarField> rec_field_p;
	PtrList<volScalarField> rec_field_nu;

	Eigen::MatrixXd Mass_matrix_P = ITHACAutilities::get_mass_matrix(example.Pmodes);

	Eigen::MatrixXd Mass_matrix_nut = ITHACAutilities::get_mass_matrix(example.nuTmodes);


	ITHACAstream::exportMatrix(Mass_matrix_nut, "Mass_matrix_nut", "python", "./ITHACAoutput/Matrices_L2/");
	ITHACAstream::exportMatrix(Mass_matrix_P, "Mass_matrix_P", "python", "./ITHACAoutput/Matrices_L2/");


	example.projectSUP("./Matrices",NmodesU,NmodesP,NmodesSUP,NmodesNUT); // Read the matrices if they already exist and if not compute them, you have to delete the folder Matrices if you want to re-do the computations with different number of modes
	Info << "after projectSUP" << endl;
	reducedSteadyNSturb pod_rbf(example); //pod_rbf is an object of the class that used RBF in the online stage and thus takes into account the eddy viscosity in that stage
	Info << "after pod_rbf" << endl;

	pod_rbf.nu = 1e-3;
	for (label k = 0; k < par_online.rows(); k++)
	{	
		Eigen::MatrixXd vel_now(2, 1);
		vel_now(0, 0) = par_online(k, 0);
		vel_now(1, 0) = par_online(k, 1);
		pod_rbf.solveOnline_sup(vel_now);
		Eigen::MatrixXd tmp_sol(pod_rbf.y.rows() + 1, 1);
		tmp_sol(0) = k + 1;
		tmp_sol.col(0).tail(pod_rbf.y.rows()) = pod_rbf.y;
		pod_rbf.online_solution.append(tmp_sol);
		
	}
	ITHACAstream::exportMatrix(pod_rbf.online_solution, "red_coeff", "python", "./ITHACAoutput/red_coeff");
	ITHACAstream::exportMatrix(pod_rbf.online_solution, "red_coeff", "matlab", "./ITHACAoutput/red_coeff");
	ITHACAstream::exportMatrix(pod_rbf.online_solution, "red_coeff", "eigen", "./ITHACAoutput/red_coeff");

	pod_rbf.reconstruct_sup("./ITHACAoutput/Reconstruction/");
	ITHACAstream::exportFields(pod_rbf.nutREC,"./ITHACAoutput/Reconstruction","nuRec");





	reducedSteadyNS pod_normal(example); //pod_normal is an object of the classical POD model which does not use RBF in the online stage, here we want to make a comparison between the results of the two models
	pod_normal.nu = 1e-3;
	for (label k = 0; k < par_online.rows(); k++)
	{	
		Eigen::MatrixXd vel_now(2, 1);
		vel_now(0, 0) = par_online(k, 0);
		vel_now(1, 0) = par_online(k, 1);
		pod_normal.solveOnline_sup(vel_now);
		Eigen::MatrixXd tmp_sol(pod_normal.y.rows() + 1, 1);
		tmp_sol(0) = k + 1;
		tmp_sol.col(0).tail(pod_normal.y.rows()) = pod_normal.y;
		pod_normal.online_solution.append(tmp_sol);
		
	}
	ITHACAstream::exportMatrix(pod_normal.online_solution, "red_coeffnew", "python", "./ITHACAoutput/red_coeffnew");
	ITHACAstream::exportMatrix(pod_normal.online_solution, "red_coeffnew", "matlab", "./ITHACAoutput/red_coeffnew");
	ITHACAstream::exportMatrix(pod_normal.online_solution, "red_coeffnew", "eigen", "./ITHACAoutput/red_coeffnew");

	pod_normal.reconstruct_sup("./ITHACAoutput/Lam_Rec/");

	example.offlineSolve(par_online,"./ITHACAoutput/high_fidelity_online/");




	exit(0);
}








