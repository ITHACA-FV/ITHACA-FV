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

#include "unsteadyNSturb.H"

/// \file
/// Source file of the unsteadyNS class.

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Construct Null
unsteadyNSturb::unsteadyNSturb() {};

// Construct from zero
unsteadyNSturb::unsteadyNSturb(int argc, char *argv[])
{
#include "setRootCase.H"
#include "createTime.H"
#include "createMesh.H"

	_pimple = autoPtr<pimpleControl>
	(
		new pimpleControl
		(
			mesh
			)
		);
#include "createFields.H"
#include "createFvOptions.H"

	supex = ITHACAutilities::check_sup();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void unsteadyNSturb::truthSolve(List<scalar> mu_now)
{
	#include "initContinuityErrs.H"
	Time& runTime = _runTime();
	surfaceScalarField& phi = _phi();
	fvMesh& mesh = _mesh();
	fv::options& fvOptions = _fvOptions();
	pimpleControl& pimple = _pimple();
	volScalarField p = _p();
	volVectorField U = _U();
	volScalarField nut = _nut();

	IOMRFZoneList& MRF = _MRF();
	singlePhaseTransportModel& laminarTransport = _laminarTransport();

	instantList Times = runTime.times();
	runTime.setEndTime(finalTime);
	// Perform a TruthSolve
	runTime.setTime(Times[1], 1);
	runTime.setDeltaT(timeStep);
	nextWrite = startTime;

	// Start the time loop
	while (runTime.run())
	{
#include "readTimeControls.H"
#include "CourantNo.H"
#include "setDeltaT.H"
		runTime.setEndTime(finalTime+timeStep);
		Info << "Time = " << runTime.timeName() << nl << endl;

		// --- Pressure-velocity PIMPLE corrector loop
		while (pimple.loop())
		{
#include "UEqn.H"
			// --- Pressure corrector loop
			while (pimple.correct())
			{
#include "pEqn.H"
			}
			if (pimple.turbCorr())
			{
				laminarTransport.correct();
				turbulence->correct();
			}
		}

		Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
		<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
		<< nl << endl;
		if (checkWrite(runTime))
		{
			volScalarField nut = turbulence->nut();
			exportSolution(U, name(counter), "./ITHACAoutput/Offline/");
			exportSolution(p, name(counter), "./ITHACAoutput/Offline/");
			exportSolution(nut, name(counter), "./ITHACAoutput/Offline/");

			std::ofstream of("./ITHACAoutput/Offline/"+name(counter)+"/"+runTime.timeName());
			Ufield.append(U);
			Pfield.append(p);
			nutFields.append(nut);

			counter++;
			nextWrite += writeEvery;
			writeMu(mu_now);
		}
		runTime++;
	}
}






List < Eigen::MatrixXd > unsteadyNSturb::turbulence_term1(label NUmodes, label NSUPmodes, label Nnutmodes)
{
	label Csize = NUmodes + NSUPmodes + liftfield.size();
	List < Eigen::MatrixXd > CT1_matrix;

	CT1_matrix.setSize(Csize);

	for (label j = 0; j < Csize; j++)
	{
		CT1_matrix[j].resize(Nnutmodes, Csize);
		CT1_matrix[j] = CT1_matrix[j] *0;
	}

	PtrList<volVectorField> Together(0);

	// Create PTRLIST with lift, velocities and supremizers
	
	if (liftfield.size() != 0)
	{
		for (label k = 0; k < liftfield.size(); k++)
		{
			Together.append(liftfield[k]);
		}
	}
	if (NUmodes != 0)
	{
		for (label k = 0; k < NUmodes; k++)
		{
			Together.append(Umodes[k]);
		}
	}
	if (NSUPmodes != 0)
	{
		for (label k = 0; k < NSUPmodes; k++)
		{
			Together.append(supmodes[k]);
		}
	}

	

	for (label i = 0; i < Csize; i++)
	{
		Info << "Filling layer number " << i+1 << " in the matrix CT1_matrix" << endl;

		for (label j = 0; j < Nnutmodes; j++)
		{
			for (label k = 0;k < Csize; k++)
			{
				CT1_matrix[i](j,k) = fvc::domainIntegrate(Together[i]&fvc::laplacian(nuTmodes[j],Together[k])).value();
			}
		}
	}
	// Export the matrix
	ITHACAstream::exportMatrix(CT1_matrix, "CT1_matrix", "python", "./ITHACAoutput/Matrices/");
	ITHACAstream::exportMatrix(CT1_matrix, "CT1_matrix", "matlab", "./ITHACAoutput/Matrices/");
	ITHACAstream::exportMatrix(CT1_matrix, "CT1_matrix", "eigen", "./ITHACAoutput/Matrices/CT1");
	return CT1_matrix;
}







List < Eigen::MatrixXd > unsteadyNSturb::turbulence_term2(label NUmodes, label NSUPmodes, label Nnutmodes)
{
	label Csize = NUmodes + NSUPmodes + liftfield.size();
	List < Eigen::MatrixXd > CT2_matrix;

	CT2_matrix.setSize(Csize);

	for (label j = 0; j < Csize; j++)
	{
		CT2_matrix[j].resize(Nnutmodes, Csize);
		CT2_matrix[j] = CT2_matrix[j] *0;

	}

	PtrList<volVectorField> Together(0);

	// Create PTRLIST with lift, velocities and supremizers
	
	if (liftfield.size() != 0)
	{
		for (label k = 0; k < liftfield.size(); k++)
		{
			Together.append(liftfield[k]);
		}
	}
	if (NUmodes != 0)
	{
		for (label k = 0; k < NUmodes; k++)
		{
			Together.append(Umodes[k]);
		}
	}
	if (NSUPmodes != 0)
	{
		for (label k = 0; k < NSUPmodes; k++)
		{
			Together.append(supmodes[k]);
		}
	}

	

	for (label i = 0; i < Csize; i++)
	{
		Info << "Filling layer number " << i+1 << " in the matrix CT2_matrix" << endl;

		for (label j = 0; j < Nnutmodes; j++)
		{
			for (label k = 0;k < Csize; k++)
			{
				CT2_matrix[i](j,k) = fvc::domainIntegrate(Together[i]&(fvc::div(nuTmodes[j]*dev((fvc::grad(Together[k]))().T())))).value();
			}
		}
	}
	// Export the matrix
	ITHACAstream::exportMatrix(CT2_matrix, "CT2_matrix", "python", "./ITHACAoutput/Matrices/");
	ITHACAstream::exportMatrix(CT2_matrix, "CT2_matrix", "matlab", "./ITHACAoutput/Matrices/");
	ITHACAstream::exportMatrix(CT2_matrix, "CT2_matrix", "eigen", "./ITHACAoutput/Matrices/CT2");
	return CT2_matrix;
}

Eigen::MatrixXd unsteadyNSturb::BT_turbulence(label NUmodes, label NSUPmodes)
{

	label BTsize = NUmodes + NSUPmodes + liftfield.size();

	Eigen::MatrixXd BT_matrix(BTsize, BTsize);

	BT_matrix = BT_matrix * 0;

    // Create PTRLIST with lift, velocities and supremizers
	PtrList<volVectorField> Together(0);


	if (liftfield.size() != 0)
	{
		for (label k = 0; k < liftfield.size(); k++)
		{
			Together.append(liftfield[k]);
		}
	}
	if (NUmodes != 0)
	{
		for (label k = 0; k < NUmodes; k++)
		{
			Together.append(Umodes[k]);
		}
	}
	if (NSUPmodes != 0)
	{
		for (label k = 0; k < NSUPmodes; k++)
		{
			Together.append(supmodes[k]);
		}
	}


	// Project everything
	for (label i = 0; i < BTsize; i++)
	{
		for (label j = 0; j < BTsize; j++)
		{
			BT_matrix(i , j) = fvc::domainIntegrate(Together[i]&(fvc::div(dev((T(fvc::grad(Together[j]))))))).value();


		}
	}

// Export the matrix
	ITHACAstream::exportMatrix(BT_matrix, "BT_matrix", "python", "./ITHACAoutput/Matrices/");
	ITHACAstream::exportMatrix(BT_matrix, "BT_matrix", "matlab", "./ITHACAoutput/Matrices/");
	ITHACAstream::exportMatrix(BT_matrix, "BT_matrix", "eigen", "./ITHACAoutput/Matrices/");
	return BT_matrix;
}

void unsteadyNSturb::projectSUP(fileName folder, label NU, label NP, label NSUP, label Nnut)
{

	if (ITHACAutilities::check_folder("./ITHACAoutput/Matrices/"))
	{
		B_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/B_mat.txt");
		C_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/C","C");
		K_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/K_mat.txt");
		P_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/P_mat.txt");
		M_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/M_mat.txt");
		BT_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/BT_matrix_mat.txt");
		CT1_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/CT1","CT1_matrix");
		CT2_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/CT2","CT2_matrix");

	}
	else
	{
		NUmodes = NU;
		NPmodes = NP;
		NSUPmodes = NSUP;
		Nnutmodes = Nnut;

		B_matrix = diffusive_term(NUmodes, NPmodes, NSUPmodes);
		C_matrix = convective_term(NUmodes, NPmodes, NSUPmodes);
		K_matrix = pressure_gradient_term(NUmodes, NPmodes, NSUPmodes);
		P_matrix = divergence_term(NUmodes, NPmodes, NSUPmodes);
		M_matrix = mass_term(NUmodes, NPmodes, NSUPmodes);
		BT_matrix = BT_turbulence(NUmodes,NSUPmodes);
		CT1_matrix = turbulence_term1(NUmodes,NSUPmodes,Nnutmodes);
		CT2_matrix = turbulence_term2(NUmodes,NSUPmodes,Nnutmodes);
	}	
	
	B_total_matrix = B_matrix + BT_matrix;
	
	C_total_matrix.setSize(C_matrix.size());
	
	for (label i = 0; i < C_matrix.size(); i++)
	{
		C_total_matrix[i] =  CT2_matrix[i] + CT1_matrix[i];
	}
    // Get the coeffs for interpolation (the orthonormal one is used because basis are orthogonal)
	Eigen::MatrixXd Ncoeff = ITHACAutilities::get_coeffs_ortho(nutFields, nuTmodes);
	ITHACAstream::exportMatrix(Ncoeff, "Ncoeff", "python", "./ITHACAoutput/Matrices/");

	NUmodes = NU;
	NPmodes = NP;
	NSUPmodes = NSUP;
	Nnutmodes = Nnut;
	SAMPLES.resize(Nnutmodes);
	rbfsplines.resize(Nnutmodes);

	for (int i = 0; i < Nnutmodes; i++)
	{

		SAMPLES[i] = new SPLINTER::DataTable(1,1);
		for(int j = 0; j < Ncoeff.cols(); j++)
		{
			SAMPLES[i]->addSample(mu.row(j),Ncoeff(i,j));

		}
		rbfsplines[i] = new SPLINTER::RBFSpline(*SAMPLES[i], SPLINTER::RadialBasisFunctionType::GAUSSIAN);
		std::cout << "Constructing RadialBasisFunction for mode " << i+1 << std::endl;
	}
}

void unsteadyNSturb::projectPPE(fileName folder, label NU, label NP, label NSUP, label Nnut)
{
	if (ITHACAutilities::check_folder("./ITHACAoutput/Matrices/"))
	{
		B_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/B_mat.txt");
		C_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/C","C");
		K_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/K_mat.txt");
		M_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/M_mat.txt");
		D_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/D_mat.txt");
		G_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/G","G");
		BC1_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/BC1_mat.txt");
		BC2_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/BC2/","BC2");
		BC3_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/BC3_mat.txt");
		BT_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/BT_matrix_mat.txt");
		CT1_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/CT1","CT1_matrix");
		CT2_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/CT2","CT2_matrix");

	}
	else
	{
		NUmodes = NU;
		NPmodes = NP;
		NSUPmodes = 0;
		Nnutmodes = Nnut;

		B_matrix = diffusive_term(NUmodes, NPmodes, NSUPmodes);
		C_matrix = convective_term(NUmodes, NPmodes, NSUPmodes);
		M_matrix = mass_term(NUmodes, NPmodes, NSUPmodes);
		K_matrix = pressure_gradient_term(NUmodes, NPmodes, NSUPmodes);
		D_matrix = laplacian_pressure(NPmodes);
		G_matrix = div_momentum(NUmodes, NPmodes);
		BC1_matrix = pressure_BC1(NUmodes, NPmodes);
		BC2_matrix = pressure_BC2(NUmodes, NPmodes);
		BC3_matrix = pressure_BC3(NUmodes, NPmodes);
		BT_matrix = BT_turbulence(NUmodes,NSUPmodes);
		CT1_matrix = turbulence_term1(NUmodes,NSUPmodes,Nnutmodes);
		CT2_matrix = turbulence_term2(NUmodes,NSUPmodes,Nnutmodes);
	}	

	B_total_matrix = B_matrix + BT_matrix;
	
	C_total_matrix.setSize(C_matrix.size());
	
	for (label i = 0; i < C_matrix.size(); i++)
	{
		C_total_matrix[i] =  CT2_matrix[i] + CT1_matrix[i];
	}
    // Get the coeffs for interpolation (the orthonormal one is used because basis are orthogonal)
	Eigen::MatrixXd Ncoeff = ITHACAutilities::get_coeffs_ortho(nutFields, nuTmodes);
	ITHACAstream::exportMatrix(Ncoeff, "Ncoeff", "python", "./ITHACAoutput/Matrices/");

	NUmodes = NU;
	NPmodes = NP;
	NSUPmodes = 0;
	Nnutmodes = Nnut;
	SAMPLES.resize(Nnutmodes);
	rbfsplines.resize(Nnutmodes);

	for (int i = 0; i < Nnutmodes; i++)
	{

		SAMPLES[i] = new SPLINTER::DataTable(1,1);
		for(int j = 0; j < Ncoeff.cols(); j++)
		{
			SAMPLES[i]->addSample(mu.row(j),Ncoeff(i,j));

		}
		rbfsplines[i] = new SPLINTER::RBFSpline(*SAMPLES[i], SPLINTER::RadialBasisFunctionType::GAUSSIAN);
		std::cout << "Constructing RadialBasisFunction for mode " << i+1 << std::endl;
	}


	for (int i = 0; i < Nnutmodes; i++)
	{

		for(int j = 0; j < Ncoeff.cols(); j++)
		{
			Info << rbfsplines[i]->eval(mu.row(j)) -Ncoeff(i,j) << endl;
		}
		
	}

}