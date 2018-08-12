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
#include "unsteadyNSturbPiso.H"
#include "unsteadyNSturb.H"

#include "ITHACAutilities.H"
#include "reductionProblem.H"

#include "reducedSteadyNS.H"
#include "reducedSteadyNSturb.H"
#include "reducedUnsteadyNSturb.H"

#include <Eigen/Dense>
#include "forces.H"
#include "IOmanip.H"
#include "forces.H"
#include "forceCoeffs.H"
#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include "Foam2Eigen.H"
#include <chrono>
#include "ITHACAstream.H"
#include "ITHACAPOD.H"
#include "EigenFunctions.H"
#include <chrono>
#include <Eigen/SVD>
#include <GenEigsSolver.h>
#include <Eigen/SparseLU>


class cylinder : public unsteadyNSturb
{
public:
	explicit cylinder(int argc, char *argv[])
	:
	unsteadyNSturb(argc, argv),
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
		Vector<double> inl(0, 0, 0);
		List<scalar> mu_now(1);
		mu_now[0] =  1.1386e-06;
		if (offline)
		{	
			Info << "Hey there" << endl;
			ITHACAstream::read_fields(Ufield, U, "./ITHACAoutput/Offline/");
			ITHACAstream::read_fields(Pfield, p, "./ITHACAoutput/Offline/");
			ITHACAstream::read_fields(nutFields, nut, "./ITHACAoutput/Offline/");

		}
		else
		{
			truthSolve(mu_now);
			//Vector<double> Uinl(0, 0, 0);
			//label BCind = 0;
			//for (label i = 0; i < mu.size(); i++)
			//{
				//change_viscosity( mu(i));
				//Info << "mu is " << mu(i) << endl;
				////assignIF(U, Uinl);
				//mu_now[0] =  mu(i);
				//if(i>0)
				//{
				//	startTime = 0.02;
				//	finalTime = 0.52;
				//}
				//truthSolve(mu_now);
			//}

		}
	}
	
};

int main(int argc, char *argv[])
{
	//return 0;
	cylinder example(argc, argv);
	word filename("./par");

	Eigen::VectorXd par;
	//par.setLinSpaced(100,1e-6,1e-4);
	// example.mu = Eigen::VectorXd::setLinSpaced(200,0.5e-4,1.5e-4);
	example.inletIndex.resize(1, 2);
	example.inletIndex << 1,0;



	ITHACAparameters para;
	int NmodesU = para.ITHACAdict->lookupOrDefault<int>("NmodesU", 5);
	int NmodesP = para.ITHACAdict->lookupOrDefault<int>("NmodesP", 5);
	int NmodesSUP = para.ITHACAdict->lookupOrDefault<int>("NmodesSUP", 5);
	int NmodesNUT = para.ITHACAdict->lookupOrDefault<int>("NmodesNUT", 5);
	int NmodesProject = para.ITHACAdict->lookupOrDefault<int>("NmodesProject", 5);
	int NmodesMatrixRec = para.ITHACAdict->lookupOrDefault<int>("NmodesMatrixRec", 5);


	//example.inletIndex(1, 0) = 0;/
	//example.inletIndex(1, 1) = 1;
 	// Time parameters
	example.startTime = 150;
	example.finalTime = 153.6;
	example.timeStep = 0.003;
	example.writeEvery = 0.006;
	

	scalar U_BC = 0.0896535433;

	//example.truthSolve(mu_now);
	example.offlineSolve();

	// volVectorField average("Ulift0",example.U);

	// average = example.Ufield[0];

	// 	for(unsigned int i = 1; i < example.Ufield.size(); i++ ){

	// 		average += example.Ufield[i]; 
	// 	}

	// 	average /= example.Ufield.size();
	// 	average = average/U_BC;
	// 	Vector<double> Uinl(1, 0, 0);
	// 	label BCind = 1;
	// 	ITHACAutilities::assignBC(average, BCind, Uinl);

	// 	example.liftfield.append(average);

	// 	ITHACAstream::exportFields(example.liftfield, "./ITHACAoutput/Offline", "lift");

	// Solve the supremizer problem
	//example.solvesupremizer();
	example.liftSolve();

	//ITHACAstream::read_fields(example.liftfield, example.U, "./lift/");

	example.computeLift(example.Ufield, example.liftfield, example.Uomfield);

	//ITHACAstream::exportFields(example.Uomfield,"./ITHACAoutput/Offline", "Uofield");

	ITHACAPOD::getModes(example.nutFields,example.nuTmodes, example.podex,0,0,NmodesProject); 
	ITHACAPOD::getModes(example.Uomfield, example.Umodes, example.podex,0,0,NmodesProject); 
	ITHACAPOD::getModes(example.Pfield, example.Pmodes, example.podex,0,0,NmodesProject); 
	//ITHACAPOD::getModes(example.supfield, example.supmodes, example.podex, example.supex, 1,NmodesProject);
	example.solvesupremizer(example.Pmodes);


	// PtrList<volVectorField> lift_and_modes;
	// PtrList<volVectorField> rec_field;
	// PtrList<volScalarField> rec_field_p;
	// PtrList<volScalarField> rec_field_nu; 

	// for (label k = 0; k<example.liftfield.size() ; k++)
	// {
	// 	lift_and_modes.append(example.liftfield[k]);
	// }	
	// for (label k = 0; k<example.Umodes.size() ; k++)
	// {
	// 	lift_and_modes.append(example.Umodes[k]);
	// }

 //    //std::cout << ITHACAutilities::get_mass_matrix(example.Umodes) << std::endl;

	// Eigen::MatrixXd coffNut = ITHACAutilities::get_coeffs_ortho(example.nutFields, example.nuTmodes);
	// //Eigen::MatrixXd coffT =  ITHACAutilities::get_coeffs(example.Ufield,lift_and_modes);
	// Eigen::MatrixXd coffU =  ITHACAutilities::get_coeffs(example.Uomfield,example.Umodes);

	// Eigen::MatrixXd coffT(coffU.rows()+1,coffU.cols());
	// coffT.row(0).fill(U_BC);
	// coffT.bottomRows(coffU.rows()) = coffU;
	// Eigen::MatrixXd coffP = ITHACAutilities::get_coeffs(example.Pfield, example.Pmodes);
	// Eigen::MatrixXd Mass_matrix_P = ITHACAutilities::get_mass_matrix(example.Pmodes);



	// ITHACAstream::exportMatrix(coffT, "coffT", "python", "./ITHACAoutput/Matrices_L2/");
	// ITHACAstream::exportMatrix(coffP, "coffP", "python", "./ITHACAoutput/Matrices_L2/");
	// ITHACAstream::exportMatrix(Mass_matrix_P, "Mass_matrix_P", "python", "./ITHACAoutput/Matrices_L2/");

	// example.reconstruct_from_matrix(rec_field_nu,example.nuTmodes,NmodesMatrixRec,coffNut);
	// example.reconstruct_from_matrix(rec_field,lift_and_modes,NmodesMatrixRec,coffT);
	// example.reconstruct_from_matrix(rec_field_p,example.Pmodes,NmodesMatrixRec,coffP);

	// ITHACAstream::exportFields(rec_field_nu, "./ITHACAoutput/Offline", "nuRec");
	// ITHACAstream::exportFields(rec_field, "./ITHACAoutput/Offline", "URec");
	// ITHACAstream::exportFields(rec_field_p, "./ITHACAoutput/Offline", "PRec");


	// Eigen::MatrixXd errorU;
	// Eigen::MatrixXd errorP;
	// Eigen::MatrixXd errorNut;

	// errorU = ITHACAutilities::error_listfields(example.Ufield, rec_field);
	// errorP = ITHACAutilities::error_listfields(example.Pfield, rec_field_p);
	// errorNut = ITHACAutilities::error_listfields(example.nutFields, rec_field_nu);

	// ITHACAstream::exportMatrix(errorU, "errorU", "python", "./ITHACAoutput/Offline/");
	// ITHACAstream::exportMatrix(errorP, "errorP", "python", "./ITHACAoutput/Offline/");
	// ITHACAstream::exportMatrix(errorNut, "errorNut", "python", "./ITHACAoutput/Offline/");

	example.mu.resize(example.Ufield.size(),1);
	Eigen::VectorXd temp;
	temp.setLinSpaced(example.Ufield.size(),example.startTime,example.finalTime);

	example.mu = temp;
	//example.mu(0) = example.timeStep;
	//example.mu.topLeftCorner(1,1) = 100;

	example.projectSUP("./Matrices",NmodesU,NmodesP,NmodesSUP,NmodesNUT);

		/*Eigen::JacobiSVD<Eigen::MatrixXd> svd(example.K_matrix);
		double cond = svd.singularValues()(0) 	/ svd.singularValues()(svd.singularValues().size()-1);*/





	//PtrList<volScalarField> rec_field_nu1;
	//Eigen::VectorXd t1(1);
	//t1(0) = temp(1);
	//for(label k=0; k<example.rbfsplines.size(); k++)
	//{	
	//std::cout << example.rbfsplines[k] ->eval(t1) << std::endl;
	//}
	//std::cout << coffNut.topLeftCorner(5,5) << std::endl;
	//exit(0);

	reducedUnsteadyNSturb ridotto(example, "SUP");

		// double Inf_sup = ridotto.inf_sup_constant();


		// std::ofstream cond_and_infsup;
		// cond_and_infsup.open("./ITHACAoutput/cond_and_infsup", std::ios_base::app);
		// cond_and_infsup << "Velocity modes " << NmodesU << std::endl;
		// cond_and_infsup << "Pressure modes " << NmodesP << std::endl;
		// cond_and_infsup << "Supremizer modes " << NmodesSUP << std::endl;
		// cond_and_infsup << "Condition Number " << cond << std::endl;
		// cond_and_infsup << "Inf_sup " << Inf_sup << std::endl;
		// cond_and_infsup.close();

	ridotto.nu = 1.1386e-06;
	ridotto.tstart = 150;
	ridotto.finalTime = 153.6;
	ridotto.dt = 0.003;

	Eigen::MatrixXd vel_now(1, 1);
	vel_now(0, 0) = U_BC;

	ridotto.solveOnline_sup(vel_now);
	ITHACAstream::exportMatrix(ridotto.online_solution, "red_coeff", "python", "./ITHACAoutput/red_coeff");
	ITHACAstream::exportMatrix(ridotto.online_solution, "red_coeff", "matlab", "./ITHACAoutput/red_coeff");
	ITHACAstream::exportMatrix(ridotto.online_solution, "red_coeff", "eigen", "./ITHACAoutput/red_coeff");
	ridotto.reconstruct_sup(example, "./ITHACAoutput/Reconstruction_SUP/",2);








	// reducedUnsteadyNSturb ridottoPPE(example, "PPE");

	// 	// double Inf_sup = ridotto.inf_sup_constant();


	// 	// std::ofstream cond_and_infsup;
	// 	// cond_and_infsup.open("./ITHACAoutput/cond_and_infsup", std::ios_base::app);
	// 	// cond_and_infsup << "Velocity modes " << NmodesU << std::endl;
	// 	// cond_and_infsup << "Pressure modes " << NmodesP << std::endl;
	// 	// cond_and_infsup << "Supremizer modes " << NmodesSUP << std::endl;
	// 	// cond_and_infsup << "Condition Number " << cond << std::endl;
	// 	// cond_and_infsup << "Inf_sup " << Inf_sup << std::endl;
	// 	// cond_and_infsup.close();

	// ridottoPPE.nu = 1.1386e-06;
	// ridottoPPE.tstart = 150;
	// ridottoPPE.finalTime = 153.6;
	// ridottoPPE.dt = 0.003;


	// ridottoPPE.solveOnline_PPE(vel_now);
	// ridottoPPE.reconstruct_PPE(example, "./ITHACAoutput/Reconstruction_PPE/",2);

	//exit(0);

	//for (label k = 0; k < example.mu.rows(); k++)
	//{
	//	//ridotto.nu = example.mu(k, 0);
	//	ridotto.solveOnline_sup(vel_now,1);
	//	Eigen::MatrixXd tmp_sol(ridotto.y.rows() + 1, 1);
	//	tmp_sol(0) = k + 1;
	//	tmp_sol.col(0).tail(ridotto.y.rows()) = ridotto.y;
	//	ridotto.online_solution.append(tmp_sol);
	//}
	
	//ITHACAstream::exportFields(example.Ufield, "./ITHACAoutput/Reconstruction", "U");




	reducedUnsteadyNS ridotto_normal(example, "SUP");
	ridotto_normal.nu = 1.1386e-06;
	ridotto_normal.tstart = 150;
	ridotto_normal.finalTime = 153.6;
	ridotto_normal.dt = 0.003;

	ridotto_normal.solveOnline_sup(vel_now);
	ridotto_normal.reconstruct_sup(example, "./ITHACAoutput/Reconstruction_normal_SUP/",2);


	exit(0);

}
