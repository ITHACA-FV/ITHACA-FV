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


/// \file
/// Source file of the steadyNS class.

#include "steadyNS.H"
#include "viscosityModel.H"

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructor
steadyNS::steadyNS() {}
steadyNS::steadyNS(int argc, char *argv[])
{
#include "setRootCase.H"
#include "createTime.H"
#include "createMesh.H"
	_simple = autoPtr<simpleControl>
	          (
	              new simpleControl
	              (
	                  mesh
	              )
	          );
	simpleControl& simple = _simple();
#include "createFields.H"
#include "createFvOptions.H"
	supex = ITHACAutilities::check_sup();
	turbulence->validate();
	ITHACAdict = new IOdictionary
	(
	    IOobject
	    (
	        "ITHACAdict",
	        runTime.system(),
	        mesh,
	        IOobject::MUST_READ,
	        IOobject::NO_WRITE
	    )
	);
	tolerance = ITHACAdict->lookupOrDefault<scalar>("tolerance", 1e-5);
}


// * * * * * * * * * * * * * * Full Order Methods * * * * * * * * * * * * * * //

// Method to performa a truthSolve
void steadyNS::truthSolve()
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
	exportSolution(U, name(counter), "./ITHACAoutput/Offline/");
	exportSolution(p, name(counter), "./ITHACAoutput/Offline/");
	Ufield.append(U);
	Pfield.append(p);
	counter++;
}

// Method to solve the supremizer problem
void steadyNS::solvesupremizer()
{
	if (supex == 1)
	{
		volVectorField U = _U();

		volVectorField Usup
		(
		    IOobject
		    (
		        "Usup",
		        U.time().timeName(),
		        U.mesh(),
		        IOobject::NO_READ,
		        IOobject::AUTO_WRITE
		    ),
		    U.mesh(),
		    dimensionedVector("zero", U.dimensions(), vector::zero)
		);
		ITHACAstream::read_fields(supfield, Usup, "./ITHACAoutput/supfield/");
	}
	else
	{
		volVectorField U = _U();

		volVectorField Usup
		(
		    IOobject
		    (
		        "Usup",
		        U.time().timeName(),
		        U.mesh(),
		        IOobject::NO_READ,
		        IOobject::AUTO_WRITE
		    ),
		    U.mesh(),
		    dimensionedVector("zero", U.dimensions(), vector::zero)
		);

		dimensionedScalar nu_fake
		(
		    "nu_fake",
		    dimensionSet(0, 2, -1, 0, 0, 0, 0),
		    scalar(1)
		);

		Vector<double> v(0, 0, 0);
		for (label i = 0; i < Usup.boundaryField().size(); i++)
		{
			changeBCtype(Usup, "fixedValue", i);
			assignBC(Usup, i, v);
			assignIF(Usup, v);
		}

		for (label i = 0; i < Pfield.size(); i++)
		{

			fvVectorMatrix u_sup_eqn
			(
			    - fvm::laplacian(nu_fake, Usup)
			);
			solve
			(
			    u_sup_eqn == fvc::grad(Pfield[i])
			);
			supfield.append(Usup);
			exportSolution(Usup, name(i+1), "./ITHACAoutput/supfield/");
		}
		int systemRet = system("ln -s ../../constant ./ITHACAoutput/supfield/constant");
		systemRet += system("ln -s ../../0 ./ITHACAoutput/supfield/0");
		systemRet += system("ln -s ../../system ./ITHACAoutput/supfield/system");
		if (systemRet < 0)
		{
			Info << "System Command Failed in steadyNS.C" << endl;
			exit(0);
		}
	}
}

// Method to compute the lifting function
void steadyNS::liftSolve()
{
	for (label k = 0; k < inletIndex.rows(); k++)
	{
		Time& runTime = _runTime();
		surfaceScalarField& phi = _phi();
		fvMesh& mesh = _mesh();
		volScalarField p = _p();
		volVectorField U = _U();
		IOMRFZoneList& MRF = _MRF();
		label BCind = inletIndex(k, 0);
		volVectorField Ulift("Ulift" + name(k), U);
		instantList Times = runTime.times();
		runTime.setTime(Times[1], 1);

		pisoControl potentialFlow(mesh, "potentialFlow");

		Info << "Solving a lifting Problem" << endl;

		Vector<double> v1(0, 0, 0);
		v1[inletIndex(k, 1)] = 1;
		Vector<double> v0(0, 0, 0);

		for (label j = 0; j < U.boundaryField().size(); j++)
		{
			if (j == BCind)
			{
				assignBC(Ulift, j, v1);
			}
			else if (U.boundaryField()[BCind].type() == "fixedValue")
			{
				assignBC(Ulift, j, v0);
			}
			else
			{

			}
			assignIF(Ulift, v0);
			phi = linearInterpolate(Ulift) & mesh.Sf();
		}
		Info << "Constructing velocity potential field Phi\n" << endl;
		volScalarField Phi
		(
		    IOobject
		    (
		        "Phi",
		        runTime.timeName(),
		        mesh,
		        IOobject::READ_IF_PRESENT,
		        IOobject::NO_WRITE
		    ),
		    mesh,
		    dimensionedScalar("Phi", dimLength * dimVelocity, 0),
		    p.boundaryField().types()
		);

		label PhiRefCell = 0;
		scalar PhiRefValue = 0;
		setRefCell
		(
		    Phi,
		    potentialFlow.dict(),
		    PhiRefCell,
		    PhiRefValue
		);

		mesh.setFluxRequired(Phi.name());

		runTime.functionObjects().start();

		MRF.makeRelative(phi);
		adjustPhi(phi, Ulift, p);


		while (potentialFlow.correctNonOrthogonal())
		{
			fvScalarMatrix PhiEqn
			(
			    fvm::laplacian(dimensionedScalar("1", dimless, 1), Phi)
			    ==
			    fvc::div(phi)
			);

			PhiEqn.setReference(PhiRefCell, PhiRefValue);
			PhiEqn.solve();

			if (potentialFlow.finalNonOrthogonalIter())
			{
				phi -= PhiEqn.flux();
			}
		}

		MRF.makeAbsolute(phi);

		Info << "Continuity error = "
		     << mag(fvc::div(phi))().weightedAverage(mesh.V()).value()
		     << endl;

		Ulift = fvc::reconstruct(phi);
		Ulift.correctBoundaryConditions();

		Info << "Interpolated velocity error = "
		     << (sqrt(sum(sqr((fvc::interpolate(U) & mesh.Sf()) - phi)))
		         / sum(mesh.magSf())).value()
		     << endl;
		Ulift.write();
		liftfield.append(Ulift);
	}
}

// * * * * * * * * * * * * * * Projection Methods * * * * * * * * * * * * * * //

void steadyNS::projectPPE(fileName folder, label NU, label NP, label NSUP)
{
	NUmodes = NU;
	NPmodes = NP;
	NSUPmodes = 0;

	B_matrix = diffusive_term(NUmodes, NPmodes, NSUPmodes);
	C_matrix = convective_term(NUmodes, NPmodes, NSUPmodes);
	M_matrix = mass_term(NUmodes, NPmodes, NSUPmodes);
	K_matrix = pressure_gradient_term(NUmodes, NPmodes, NSUPmodes);

	D_matrix = laplacian_pressure(NPmodes);
	G_matrix = div_momentum(NUmodes, NPmodes);

	BC1_matrix = pressure_BC1(NUmodes, NPmodes);
	BC2_matrix = pressure_BC2(NUmodes, NPmodes);
	BC3_matrix = pressure_BC3(NUmodes, NPmodes);
}

void steadyNS::projectSUP(fileName folder, label NU, label NP, label NSUP)
{
	NUmodes = NU;
	NPmodes = NP;
	NSUPmodes = NSUP;

	B_matrix = diffusive_term(NUmodes, NPmodes, NSUPmodes);
	C_matrix = convective_term(NUmodes, NPmodes, NSUPmodes);
	K_matrix = pressure_gradient_term(NUmodes, NPmodes, NSUPmodes);
	P_matrix = divergence_term(NUmodes, NPmodes, NSUPmodes);
	M_matrix = mass_term(NUmodes, NPmodes, NSUPmodes);
}

// * * * * * * * * * * * * * * Momentum Eq. Methods * * * * * * * * * * * * * //

Eigen::MatrixXd steadyNS::diffusive_term(label NUmodes, label NPmodes, label NSUPmodes)
{
	label Bsize = NUmodes + NSUPmodes + liftfield.size();
	Eigen::MatrixXd B_matrix;
	B_matrix.resize(Bsize, Bsize);

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
	for (label i = 0; i < Bsize; i++)
	{
		for (label j = 0; j < Bsize; j++)
		{
			B_matrix(i, j) = fvc::domainIntegrate(Together[i] & fvc::laplacian(dimensionedScalar("1", dimless, 1), Together[j])).value();
		}
	}
	// Export the matrix
	ITHACAstream::exportMatrix(B_matrix, "B", "python", "./ITHACAoutput/Matrices/");
	ITHACAstream::exportMatrix(B_matrix, "B", "matlab", "./ITHACAoutput/Matrices/");
	ITHACAstream::exportMatrix(B_matrix, "B", "eigen", "./ITHACAoutput/Matrices/");
	return B_matrix;
}

Eigen::MatrixXd steadyNS::pressure_gradient_term(label NUmodes, label NPmodes, label NSUPmodes)
{

	label K1size = NUmodes + NSUPmodes + liftfield.size();
	label K2size = NPmodes;

	Eigen::MatrixXd K_matrix(K1size, K2size);

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
	for (label i = 0; i < K1size; i++)
	{
		for (label j = 0; j < K2size; j++)
		{
			K_matrix(i, j) = fvc::domainIntegrate(Together[i] & fvc::grad(Pmodes[j])).value();
		}
	}

// Export the matrix
	ITHACAstream::exportMatrix(K_matrix, "K", "python", "./ITHACAoutput/Matrices/");
	ITHACAstream::exportMatrix(K_matrix, "K", "matlab", "./ITHACAoutput/Matrices/");
	ITHACAstream::exportMatrix(K_matrix, "K", "eigen", "./ITHACAoutput/Matrices/");
	return K_matrix;
}

List < Eigen::MatrixXd > steadyNS::convective_term(label NUmodes, label NPmodes, label NSUPmodes)
{
	label Csize = NUmodes + NSUPmodes + liftfield.size();
	List < Eigen::MatrixXd > C_matrix;

	C_matrix.setSize(Csize);

	for (label j = 0; j < Csize; j++)
	{
		C_matrix[j].resize(Csize, Csize);
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
		for (label j = 0; j < Csize; j++)
		{
			for (label k = 0; k < Csize; k++)
			{
				C_matrix[i](j, k) = fvc::domainIntegrate(Together[i] & fvc::div(linearInterpolate(Together[j]) & Together[j].mesh().Sf(), Together[k])).value();
			}
		}
	}
	// Export the matrix
	ITHACAstream::exportMatrix(C_matrix, "C", "python", "./ITHACAoutput/Matrices/");
	ITHACAstream::exportMatrix(C_matrix, "C", "matlab", "./ITHACAoutput/Matrices/");
	ITHACAstream::exportMatrix(C_matrix, "C", "eigen", "./ITHACAoutput/Matrices/C");
	return C_matrix;
}

Eigen::MatrixXd steadyNS::mass_term(label NUmodes, label NPmodes, label NSUPmodes)
{
	label Msize = NUmodes + NSUPmodes + liftfield.size();

	Eigen::MatrixXd M_matrix(Msize, Msize);

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
	for (label i = 0; i < Msize; i++)
	{
		for (label j = 0; j < Msize; j++)
		{
			M_matrix(i, j) = fvc::domainIntegrate(Together[i] & Together[j]).value();
		}
	}
	// Export the matrix
	ITHACAstream::exportMatrix(M_matrix, "M", "python", "./ITHACAoutput/Matrices/");
	ITHACAstream::exportMatrix(M_matrix, "M", "matlab", "./ITHACAoutput/Matrices/");
	ITHACAstream::exportMatrix(M_matrix, "M", "eigen", "./ITHACAoutput/Matrices/");
	return M_matrix;
}

// * * * * * * * * * * * * * * Continuity Eq. Methods * * * * * * * * * * * * * //

Eigen::MatrixXd steadyNS::divergence_term(label NUmodes, label NPmodes, label NSUPmodes)
{
	label P1size = NPmodes;
	label P2size = NUmodes + NSUPmodes + liftfield.size();

	Eigen::MatrixXd P_matrix(P1size, P2size);

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
	for (label i = 0; i < P1size; i++)
	{
		for (label j = 0; j < P2size; j++)
		{
			P_matrix(i, j) = fvc::domainIntegrate(Pmodes[i] * fvc::div (Together[j])).value();
		}
	}
	//Export the matrix
	ITHACAstream::exportMatrix(P_matrix, "P", "python", "./ITHACAoutput/Matrices/");
	ITHACAstream::exportMatrix(P_matrix, "P", "matlab", "./ITHACAoutput/Matrices/");
	ITHACAstream::exportMatrix(P_matrix, "P", "eigen", "./ITHACAoutput/Matrices/");
	return P_matrix;
}


List < Eigen::MatrixXd > steadyNS::div_momentum(label NUmodes, label NPmodes)
{
	label G1size = NPmodes;
	label G2size = NUmodes + NSUPmodes + liftfield.size();
	List < Eigen::MatrixXd > G_matrix;

	G_matrix.setSize(G1size);

	for (label j = 0; j < G1size; j++)
	{
		G_matrix[j].resize(G2size, G2size);
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

	for (label i = 0; i < G1size; i++)
	{
		for (label j = 0; j < G2size; j++)
		{
			for (label k = 0; k < G2size; k++)
			{
				G_matrix[i](j, k) = fvc::domainIntegrate(fvc::grad(Pmodes[i]) & (fvc::div(fvc::interpolate(Together[j]) & Together[j].mesh().Sf(), Together[k]))).value();
			}
		}
	}
	// Export the matrix
	ITHACAstream::exportMatrix(G_matrix, "G", "python", "./ITHACAoutput/Matrices/");
	ITHACAstream::exportMatrix(G_matrix, "G", "matlab", "./ITHACAoutput/Matrices/");
	ITHACAstream::exportMatrix(G_matrix, "G", "eigen", "./ITHACAoutput/Matrices/G");
	return G_matrix;
}

Eigen::MatrixXd steadyNS::laplacian_pressure(label NPmodes)
{

	label Dsize = NPmodes;

	Eigen::MatrixXd D_matrix(Dsize, Dsize);

	// Project everything
	for (label i = 0; i < Dsize; i++)
	{
		for (label j = 0; j < Dsize; j++)
		{
			D_matrix(i, j) = fvc::domainIntegrate(fvc::grad(Pmodes[i])&fvc::grad(Pmodes[j])).value();
		}
	}

	//Export the matrix
	ITHACAstream::exportMatrix(D_matrix, "D", "python", "./ITHACAoutput/Matrices/");
	ITHACAstream::exportMatrix(D_matrix, "D", "matlab", "./ITHACAoutput/Matrices/");
	ITHACAstream::exportMatrix(D_matrix, "D", "eigen", "./ITHACAoutput/Matrices/");
	return D_matrix;
}

Eigen::MatrixXd steadyNS::pressure_BC1(label NUmodes, label NPmodes)
{

	label P_BC1size = NPmodes;
	label P_BC2size = NUmodes + liftfield.size();

	Eigen::MatrixXd BC1_matrix(P_BC1size, P_BC2size);
	fvMesh& mesh = _mesh();


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


	for (label i = 0; i < P_BC1size; i++)
	{
		for (label j = 0; j < P_BC2size; j++)
		{
			surfaceScalarField lpl((fvc::interpolate(fvc::laplacian(Together[j]))&mesh.Sf())*fvc::interpolate(Pmodes[i]));
			double s = 0;
			for (label k = 0; k < lpl.boundaryField().size(); k++)
			{
				s += gSum(lpl.boundaryField()[k]);
			}
			BC1_matrix(i, j) = s;
		}
	}
	return BC1_matrix;
}


List < Eigen::MatrixXd > steadyNS::pressure_BC2(label NUmodes, label NPmodes)
{
	label P2_BC1size = NPmodes;
	label P2_BC2size = NUmodes + NSUPmodes + liftfield.size();
	List < Eigen::MatrixXd > BC2_matrix;

	fvMesh& mesh = _mesh();

	BC2_matrix.setSize(P2_BC1size);

	for (label j = 0; j < P2_BC1size; j++)
	{
		BC2_matrix[j].resize(P2_BC2size, P2_BC2size);
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

	for (label i = 0; i < P2_BC1size; i++)
	{
		for (label j = 0; j < P2_BC2size; j++)
		{
			for (label k = 0; k < P2_BC2size; k++)
			{
				surfaceScalarField div_m(fvc::interpolate(fvc::div(fvc::interpolate(Together[j]) & mesh.Sf(), Together[k]))&mesh.Sf()*fvc::interpolate(Pmodes[i]));
				double s = 0;
				for (label k = 0; k < div_m.boundaryField().size(); k++)
				{
					s += gSum(div_m.boundaryField()[k]);
				}

				BC2_matrix[i](j, k) = s;
			}
		}
	}
	// Export the matrix
	return BC2_matrix;
}

Eigen::MatrixXd steadyNS::pressure_BC3(label NUmodes, label NPmodes)
{
	label P3_BC1size = NPmodes;
	label P3_BC2size = NUmodes + liftfield.size();

	Eigen::MatrixXd BC3_matrix(P3_BC1size, P3_BC2size);
	fvMesh& mesh = _mesh();

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

	surfaceVectorField n(mesh.Sf() / mesh.magSf());
	for (label i = 0; i < P3_BC1size; i++)
	{
		for (label j = 0; j < P3_BC2size; j++)
		{
			surfaceVectorField BC3 = fvc::interpolate(fvc::curl(Together[j]));
			surfaceVectorField BC4 = n ^ fvc::interpolate(fvc::grad(Pmodes[i]));
			surfaceScalarField BC5 = (BC3 & BC4) * mesh.magSf();
			double s = 0;
			for (label k = 0; k < BC5.boundaryField().size(); k++)
			{
				s += gSum(BC5.boundaryField()[k]);
			}
			BC3_matrix(i, j) = s;
		}
	}
	return BC3_matrix;
}


void steadyNS::change_viscosity(double mu)
{
	const volScalarField& nu = 	_laminarTransport().nu();
	volScalarField& ciao = const_cast<volScalarField&>(nu);

	this->assignIF(ciao, mu);
	for (int i = 0; i < ciao.boundaryFieldRef().size(); i++)
	{
		this->assignBC(ciao, i, mu);
	}
}


void steadyNS::Forces_matrices(label NUmodes, label NPmodes, label NSUPmodes)
{
	PtrList<volVectorField> Together(0);
	for (label k = 0; k < liftfield.size(); k++)
	{
		Together.append(liftfield[k]);
	}
	for (label k = 0; k < NUmodes; k++)
	{
		Together.append(Umodes[k]);
	}
	for (label k = 0; k < NSUPmodes; k++)
	{
		Together.append(supmodes[k]);
	}

	tau_matrix.resize(Together.size(), 3);
	n_matrix.resize(NPmodes, 3);

	tau_matrix = tau_matrix * 0;
	n_matrix = n_matrix * 0;

	Time& runTime = _runTime();

	instantList Times = runTime.times();

	fvMesh& mesh = _mesh();
	volScalarField& p = _p();
	volVectorField& U = _U();

	//Read FORCESdict
	IOdictionary FORCESdict
	(
	    IOobject
	    (
	        "FORCESdict",
	        runTime.system(),
	        mesh,
	        IOobject::MUST_READ,
	        IOobject::NO_WRITE
	    )
	);

	IOdictionary transportProperties
	(
	    IOobject
	    (
	        "transportProperties",
	        runTime.constant(),
	        mesh,
	        IOobject::MUST_READ,
	        IOobject::NO_WRITE
	    )
	);

	word pName = FORCESdict.lookup("pName");
	word UName = FORCESdict.lookup("UName");

	functionObjects::ITHACAforces f("Forces", mesh, FORCESdict);
	for (label i = 0; i < Together.size(); i++)
	{
		U = Together[i];
		p = Pmodes[0];
		mesh.readUpdate();
		f.write();
		f.calcForcesMoment();

		for (label j = 0; j < 3; j++)
		{

			tau_matrix(i, j) = f.force_tau()[j];
		}
	}




	for (label i = 0; i < NPmodes; i++)
	{
		U = Together[0];
		p = Pmodes[i];
		mesh.readUpdate();
		f.write();
		f.calcForcesMoment();

		for (label j = 0; j < 3; j++)
		{
			n_matrix(i, j) = f.force_pressure()[j];
		}

	}

	ITHACAstream::exportMatrix(tau_matrix, "tau", "python", "./ITHACAoutput/Matrices/");
	ITHACAstream::exportMatrix(tau_matrix, "tau", "matlab", "./ITHACAoutput/Matrices/");
	ITHACAstream::exportMatrix(tau_matrix, "tau", "eigen", "./ITHACAoutput/Matrices/");

	ITHACAstream::exportMatrix(n_matrix, "n", "python", "./ITHACAoutput/Matrices/");
	ITHACAstream::exportMatrix(n_matrix, "n", "matlab", "./ITHACAoutput/Matrices/");
	ITHACAstream::exportMatrix(n_matrix, "n", "eigen", "./ITHACAoutput/Matrices/");



}



