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
#include "ITHACAstream.H"
#include "ITHACAPOD.H"
#include "reducedSteadyNS.H"
#include "forces.H"
#include "IOmanip.H"


class tutorial03 : public steadyNS
{
public:
	/// Constructor
	explicit tutorial03(int argc, char *argv[])
		:
		steadyNS(argc, argv),
		U(_U()),
		p(_p())
	{}

	/// Velocity field
	volVectorField& U;
	/// Pressure field
	volScalarField& p;

	/// Perform an Offline solve
	void offlineSolve()
	{
		Vector<double> inl(0, 0, 0);
		List<scalar> mu_now(1);
		// if the offline solution is already performed read the fields
		if (offline)
		{
			ITHACAstream::read_fields(Ufield, U, "./ITHACAoutput/Offline/");
			ITHACAstream::read_fields(Pfield, p, "./ITHACAoutput/Offline/");
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
				truthSolve();
			}

		}
	}

};

int main(int argc, char *argv[])
{
	// Construct the tutorial object
	tutorial03 example(argc, argv);

	// Read the par file where the parameters are stored
	word filename("./par");
	example.mu = ITHACAstream::readMatrix(filename);

	// Set the inlet boundaries patch 0 directions x and y
	example.inletIndex.resize(2, 2);
	example.inletIndex(0, 0) = 0;
	example.inletIndex(0, 1) = 0;
	example.inletIndex(1, 0) = 0;
	example.inletIndex(1, 1) = 1;

	// Perform the offline solve
	example.offlineSolve();

	// Solve the supremizer problem
	example.solvesupremizer();
	
	//example.liftSolve();
	// Read the lifting functions from file
	ITHACAstream::read_fields(example.liftfield, example.U, "./lift/");
	
	// Homogenize the snapshots
	example.computeLift(example.Ufield, example.liftfield, example.Uomfield);

	// Perform POD on velocity pressure and supremizers and store the first 50 modes
	ITHACAPOD::getModes(example.Uomfield, example.Umodes, example.podex, 0, 20);
	ITHACAPOD::getModes(example.Pfield, example.Pmodes, example.podex, 0, 20);
	ITHACAPOD::getModes(example.supfield, example.supmodes, example.podex, example.supex, 1, 20);

	// Perform the Galerkin Projection
	example.projectSUP("./Matrices", 5, 5 , 5);

	// Create the reduced object
	reducedSteadyNS ridotto(example, "SUP");

	// Set the reduced viscosity
	ridotto.nu = 1;

	// Perform an online solve for the new values of inlet velocities
	for (label k = 0; k < example.mu.rows(); k++)
	{
		Eigen::MatrixXd vel_now(2, 1);
		vel_now(0, 0) = example.mu(k, 0);
		vel_now(1, 0) = example.mu(k, 1);
		ridotto.solveOnline_sup(vel_now);
		Eigen::MatrixXd tmp_sol(ridotto.y.rows() + 1, 1);
		tmp_sol(0) = k + 1;
		tmp_sol.col(0).tail(ridotto.y.rows()) = ridotto.y;
		ridotto.online_solution.append(tmp_sol);
	}

	// Save the online solution
	ITHACAstream::exportMatrix(ridotto.online_solution, "red_coeff", "python", "./ITHACAoutput/red_coeff");
	ITHACAstream::exportMatrix(ridotto.online_solution, "red_coeff", "matlab", "./ITHACAoutput/red_coeff");
	ITHACAstream::exportMatrix(ridotto.online_solution, "red_coeff", "eigen", "./ITHACAoutput/red_coeff");

	// Reconstruct and export the solution
	ridotto.reconstruct_sup(example, "./ITHACAoutput/Reconstruction/");
	exit(0);
}








