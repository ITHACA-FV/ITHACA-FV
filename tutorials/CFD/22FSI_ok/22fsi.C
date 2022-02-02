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
Description
    Example of an unsteady NS Reduction Problem
SourceFiles
    04unsteadyNS.C
\*---------------------------------------------------------------------------*/

#include "fsiBasic.H"
#include "ITHACAPOD.H"
#include "ReducedUnsteadyNS.H"
#include "ReducedSimpleSteadyNS.H"
#include "ITHACAstream.H"
#include "dynamicFvMesh.H"
#include <chrono>
#include<math.h>
#include<iomanip>

class tutorial22: public fsiBasic
{
    public:
        explicit tutorial22(int argc, char* argv[])
            :
            fsiBasic(argc, argv),
            U(_U()),
            p(_p()),
            phi(_phi())
        {}

        // Fields To Perform
       /// Velocity field
        volVectorField& U;
        /// Pressure field
        volScalarField& p;
        ///
        surfaceScalarField& phi;
        /// Dynamic mesh field
        //Foam::autoPtr<Foam::dynamicFvMesh> meshPtr;

        void offlineSolve()
        {
            Vector<double> inl(1, 0, 0);
            List<scalar> mu_now(1);

            // if (offline)
            // {
            //     ITHACAstream::read_fields(Ufield, U, "./ITHACAoutput/Offline/");
            //     ITHACAstream::read_fields(Pfield, p, "./ITHACAoutput/Offline/");
            //      // mu_samples =
            //      //    ITHACAstream::readMatrix("./ITHACAoutput/Offline/mu_samples_mat.txt");
            // }
            // else
            //{
                //for (label i = 0; i < mu.cols(); i++)
                //{
                    //inl[0] = mu(0, i);
                    mu_now[0] = mu(0, 0); //mu.cols()=50
                    //mu_now[0] = mu(0, i);
                    //std::cout << "////////////////////////////mu is :///////////////////"<< mu.cols() <<  std::endl;                 
                    //std::cout << "////////////////////////////mu_now[0] is :///////////////////"<< mu_now[0] <<  std::endl;                

                    //assignBC(U, BCind, inl);
                    //assignIF(U, inl);
                    //change_viscosity(mu(0, i));
                    //std::cout << "////////////////////////////i is :///////////////////"<< i <<  std::endl;

                    truthSolve3(mu_now);
                //}
            //}
        }
};

// reduced class problem

class reducedFSI : public reducedSimpleSteadyNS
{
    public:
        /// Constructor
        reducedFSI() {}
        explicit reducedFSI(tutorial22& FOMproblem)
            :
            problem(&FOMproblem)
        {


            // Create a new Umodes set where the first ones are the lift functions
            for (int i = 0; i < problem->inletIndex.rows(); i++)
            {
                ULmodes.append((problem->liftfield[i]).clone());
            }

            for (int i = 0; i < problem->Umodes.size(); i++)
            {
                ULmodes.append((problem->Umodes.toPtrList()[i]).clone());
            }
            std::cout << "////////////////////////////////////reduced problem Constructor calling////////////////////////////////////////////"<< std::endl;
        }
        // Variables

        /// Lifted velocity modes.
        volVectorModes ULmodes;
        /// Full problem.
        tutorial22* problem;
       
        //Eigen::MatrixXd projGradModP;
        /// Projected gradient of the pressure modes.
        Eigen::MatrixXd projGradModP;

        /// Imposed boundary conditions.
        Eigen::MatrixXd vel_now;

        /// Maximum iterations number for the online step
        int maxIterOn = 1000;

        /// Counter.
        int counter = 0;

        int UprojN;
        int PprojN;

        // Function to perform the online phase
        void solveOnline_Pimple(scalar mu_now,
                                int NmodesUproj, int NmodesPproj, int NmodesSup = 0,
                                fileName  folder = "./ITHACAoutput/Reconstruct/")
        {

            ULmodes.resize(0);

            for (int i = 0; i < problem->inletIndex.rows(); i++)
            {
                ULmodes.append((problem->liftfield[i]).clone());
            }

            for (int i = 0; i < NmodesUproj; i++)
            {
                ULmodes.append((problem->Umodes.toPtrList()[i]).clone());
            }

            // for (int i = 0; i < NmodesSup; i++)
            // {
            //     ULmodes.append((problem->supmodes.toPtrList()[i]).clone());
            // }

            counter++;

            if (NmodesUproj == 0)
            {
                UprojN = ULmodes.size();
            }
            else
            {
                UprojN = NmodesUproj + NmodesSup;
            }

            if (NmodesPproj == 0)
            {
                PprojN = problem->Pmodes.size();
            }
            else
            {
            	PprojN = NmodesPproj;
            }

           /* if (NmodesNut == 0)
            {
                NmodesNut = problem->nutModes.size();
            }*/

            // Initializations
            Eigen::VectorXd uresidualOld = Eigen::VectorXd::Zero(UprojN);
            Eigen::VectorXd presidualOld = Eigen::VectorXd::Zero(PprojN);
            Eigen::VectorXd uresidual = Eigen::VectorXd::Zero(UprojN);
            Eigen::VectorXd presidual = Eigen::VectorXd::Zero(PprojN);
            scalar U_norm_res(1);
            scalar P_norm_res(1);
            // coefficients declaration
            Eigen::MatrixXd a = Eigen::VectorXd::Zero(UprojN);
            //Eigen::MatrixXd a0 = Eigen::VectorXd::Zero(NmodesUproj);
            Eigen::MatrixXd b = Eigen::VectorXd::Zero(PprojN);
            a(0) = vel_now(0, 0);
            float residualJumpLim =
                problem->para->ITHACAdict->lookupOrDefault<float>("residualJumpLim", 1e-5);
            float normalizedResidualLim =
                problem->para->ITHACAdict->lookupOrDefault<float>("normalizedResidualLim",
                        1e-5);
            scalar residual_jump(1 + residualJumpLim);
            volScalarField& P = problem->_p();
            volVectorField& U = problem->_U();
            dynamicFvMesh& mesh = problem->meshPtr();
            Time& runTime = problem->_runTime();
            P.rename("p");
            surfaceScalarField& phi(problem->_phi());

            problem->Umodes.reconstruct(U, a, "U");
            problem->Pmodes.reconstruct(P, b, "p");
            phi = fvc::flux(U);
            int iter = 0;
            label pRefCell = 0;
            scalar pRefValue = 0.0;
            pimpleControl& pimple = problem->_pimple();
           
            IOMRFZoneList& MRF = problem->_MRF();
            fv::options& fvOptions = problem->_fvOptions();
            // std::ofstream res_os_U;
            // std::ofstream res_os_P;
            // res_os_U.open(folder + name(counter) + "/residualsU", std::ios_base::app);
            // res_os_P.open(folder + name(counter) + "/residualsP", std::ios_base::app);

            std::cout << "///////////////////////////solveOnline_Pimple method line 229//////////////////////////////////"<< std::endl;

            PtrList<volVectorField> gradModP;

		    for (int i = 0; i < NmodesPproj; i++)
		    {
		        gradModP.append(fvc::grad(problem->Pmodes[i]));
		    }

		    projGradModP = ULmodes.project(gradModP, NmodesUproj);
            // Pimple algorithm starts here
            while ((residual_jump > residualJumpLim
                    || std::max(U_norm_res, P_norm_res) > normalizedResidualLim) &&
                    iter < maxIterOn)
            {
                iter++;
                std::cout << "///////////////////////////solveOnline_Pimple method line 249: iter value is:  //////////////////////////////////"<< iter <<  std::endl;

                //std::cout << iter << std::endl;

#if defined(OFVER) && (OFVER == 6)
                pimple.loop(runTime);
                 std::cout << "///////////////////////////solveOnline_Pimple method line 251 after pimple.loop(runTime) //////////////////////////////////"<< std::endl;
#else
                pimple.loop();
                 std::cout << "///////////////////////////solveOnline_Pimple method line 254 after pimple.loop() //////////////////////////////////"<< std::endl;
#endif

                //volScalarField nueff = problem->turbulence->nuEff();
                //vector v(1, 0, 0);
                //ITHACAutilities::assignBC(U, 0, v);
                /*// Momentum equation
                fvVectorMatrix UEqn
                (
                    fvm::div(phi, U)
                    - fvm::laplacian(nueff, U)
                    - fvc::div(nueff * dev2(T(fvc::grad(U))))
                );
                UEqn.relax();*/

                // Solve the Momentum equation

                MRF.correctBoundaryVelocity(U);

                tmp<fvVectorMatrix> tUEqn
                (
                fvm::ddt(U) + fvm::div(phi, U)
                + MRF.DDt(U)
                + problem->turbulence->divDevReff(U)
                ==
                fvOptions(U)
                );
                fvVectorMatrix& UEqn = tUEqn.ref();

                UEqn.relax();
                std::cout << "///////////////////////////solveOnline_Pimple method line 285 //////////////////////////////////"<< std::endl;

                //////////////////////////////reduced the linear system for the velocity////////////////////////////////////////////
                List<Eigen::MatrixXd> RedLinSysU = ULmodes.project(UEqn, UprojN);
                std::cout << "///////////////////////////solveOnline_Pimple line 288: RedLinSysU  //////////////////////////////////"<< RedLinSysU.size() << std::endl;
                std::cout << "///////////////////////////solveOnline_Pimple line 289: b //////////////////////////////////"<< b << std::endl;

                RedLinSysU[1] = RedLinSysU[1] - projGradModP * b;
                std::cout << "///////////////////////////solveOnline_Pimple method line 292 //////////////////////////////////"<< std::endl;

                a = reducedProblem::solveLinearSys(RedLinSysU, a, uresidual, vel_now);
                std::cout << "///////////////////////////solveOnline_Pimple method line 295//////////////////////////////////"<< std::endl;

                ULmodes.reconstruct(U, a, "U");
        
                // fvOptions.constrain(UEqn);

                // if (pimple.momentumPredictor())
                // {
                //    solve(UEqn == -fvc::grad(P));
        
                //    fvOptions.correct(U);//??
                // }
                //solve(UEqn == - fvc::grad(P));
                //ITHACAutilities::assignBC(U, 0, v);
                std::cout << "///////////////////////////solveOnline_Pimple method line 309 //////////////////////////////////"<< std::endl;

                
              ///////////////////////////////////////// begining of solving the pEqn.H//////////////////////////////////////   

                volScalarField rAU(1.0/UEqn.A());
                volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, P));
                surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));

                if (pimple.ddtCorr())
                {
                    #include "createUfIfPresent.H"
                    phiHbyA += MRF.zeroFilter(fvc::interpolate(rAU)*fvc::ddtCorr(U, phi, Uf));
                }
                else
                {
                    phiHbyA += MRF.zeroFilter(fvc::interpolate(rAU));
                }

                MRF.makeRelative(phiHbyA);

                if (P.needReference())
                {
                    fvc::makeRelative(phiHbyA, U);
                    adjustPhi(phiHbyA, U, P);
                    fvc::makeAbsolute(phiHbyA, U);
                }

                tmp<volScalarField> rAtU(rAU);

                if (pimple.consistent())
                {
                    rAtU = 1.0/max(1.0/rAU - UEqn.H1(), 0.1/rAU);
                    phiHbyA +=
                        fvc::interpolate(rAtU() - rAU)*fvc::snGrad(P)*mesh.magSf();
                    HbyA -= (rAU - rAtU())*fvc::grad(P);
                }
                List<Eigen::MatrixXd> RedLinSysP;
                //bOld = b;

                if (pimple.nCorrPISO() <= 1)
                {
                    tUEqn.clear();
                }

                // Update the pressure BCs to ensure flux consistency
                constrainPressure(P, U, phiHbyA, rAtU(), MRF);

                // Non-orthogonal pressure corrector loop
                while (pimple.correctNonOrthogonal())
                {
                    // Continuity equation
                    fvScalarMatrix pEqn
                    (
                        fvm::laplacian(rAtU(), P) == fvc::div(phiHbyA)
                    );
                    ////////////////////////////////// Added for reduced problem to project the pressure ///////////////////////////////////////////////
                    RedLinSysP = problem->Pmodes.project(pEqn, NmodesPproj);
                    //pEqn.solve();
                    b = reducedProblem::solveLinearSys(RedLinSysP, b, presidual);
                    problem->Pmodes.reconstruct(P, b, "p");
                    pEqn.setReference(pRefCell, pRefValue);

                    pEqn.solve(mesh.solver(P.select(pimple.finalInnerIter())));

                    if (pimple.finalNonOrthogonalIter())
                    {
                        phi = phiHbyA - pEqn.flux();
                    }
                }
                scalar cumulativeContErr = 0.0;

                #include "continuityErrs.H"

                // Explicitly relax pressure for momentum corrector
                //p.relax();

                /*U = HbyA - rAtU*fvc::grad(p);
                U.correctBoundaryConditions();*/
               

                // Correct Uf if the mesh is moving
                //#include "correctPhi.H"
                #include "createUfIfPresent.H"
                fvc::correctUf(Uf, U, phi);

                // Make the fluxes relative to the mesh motion
                fvc::makeRelative(phi, U);

                //b = bOld + mesh.fieldRelaxationFactor("p") * (b - bOld);
                problem->Pmodes.reconstruct(P, b, "p");
                //nutCoeffOld = nutCoeff;
                // P.relax();
                U = HbyA - rAtU() * fvc::grad(P);
                U.correctBoundaryConditions();
                fvOptions.correct(U);

                uresidualOld = uresidualOld - uresidual;
                presidualOld = presidualOld - presidual;
                uresidualOld = uresidualOld.cwiseAbs();
                presidualOld = presidualOld.cwiseAbs();
                residual_jump = std::max(uresidualOld.sum(), presidualOld.sum());
                uresidualOld = uresidual;
                presidualOld = presidual;
                uresidual = uresidual.cwiseAbs();
                presidual = presidual.cwiseAbs();
                U_norm_res = uresidual.sum() / (RedLinSysU[1].cwiseAbs()).sum();
                P_norm_res = presidual.sum() / (RedLinSysP[1].cwiseAbs()).sum();

                if (problem->para->debug)
                {
                    std::cout << "Residual jump = " << residual_jump << std::endl;
                    std::cout << "Normalized residual = " << std::max(U_norm_res,
                              P_norm_res) << std::endl;
                    std::cout << "Final normalized residual for velocity: " << U_norm_res <<
                              std::endl;
                    std::cout << "Final normalized residual for pressure: " << P_norm_res <<
                              std::endl;
                }

                // res_os_U << U_norm_res << std::endl;
                // res_os_P << P_norm_res << std::endl;
            }

            // res_os_U.close();
            // res_os_P.close();
            std::cout << "Solution " << counter << " converged in " << iter <<
                      " iterations." << std::endl;
            std::cout << "Final normalized residual for velocity: " << U_norm_res <<
                      std::endl;
            std::cout << "Final normalized residual for pressure: " << P_norm_res <<
                      std::endl;
            problem->Umodes.reconstruct(U, a, "Uaux");
            problem->Pmodes.reconstruct(P, b, "Paux");

            /*if (ITHACAutilities::isTurbulent())
            {
                volScalarField& nut = const_cast<volScalarField&>
                                      (problem->_mesh().lookupObject<volScalarField>("nut"));
                nut.rename("nutAux");
        ITHACAstream::exportSolution(nut, name(counter), Folder);
            }*/

            ITHACAstream::exportSolution(U, name(counter), folder);
            ITHACAstream::exportSolution(P, name(counter), folder);
            runTime.setTime(runTime.startTime(), 0);
            std::cout << "////////////////////////////solveOnline_Pimple method end///////////////////////////"<< std::endl;                 

        }


	void OnlineVelocity(Eigen::MatrixXd vel)
	{
	        M_Assert(problem->inletIndex.rows() == vel.size(),
	                 "Imposed boundary conditions dimensions do not match given values matrix dimensions");
	        Eigen::MatrixXd vel_scal;
	        vel_scal.resize(vel.rows(), vel.cols());

	        for (int k = 0; k < problem->inletIndex.rows(); k++)
	        {
	            int p = problem->inletIndex(k, 0);
	            //int p = problem->inletIndex(0, 0);

	            //std::cout << "/////////////////////////// Value of p is:  //////////////////////////////////"<< p << std::endl;
	            int l = problem->inletIndex(k, 1);
	            //int l = problem->inletIndex(0, 1);
	            //std::cout << "/////////////////////////// Value of l is ://////////////////////////////////" << l << std::endl;
	            scalar area = gSum(problem->liftfield[0].mesh().magSf().boundaryField()[p]);
	            //std::cout << "/////////////////////////// the value of the area is : //////////////////////////////////"<< area << std::endl;
	            scalar u_lf = gSum(problem->liftfield[k].mesh().magSf().boundaryField()[p] *
	                              problem->liftfield[k].boundaryField()[p]).component(l) / area;
	            //std::cout << "///////////////////////////the value of u_lf is : //////////////////////////////////"<< u_lf << std::endl;
	            vel_scal(k, 0) = vel(k, 0) / u_lf;
	            //vel_scal(0, 0) = vel(0, 0) / u_lf;
	            //std::cout << "///////////////////////////vel_scal //////////////////////////////////"<< vel_scal(k,0) << std::endl;

	        }
            //std::cout << "///////////////////////////vel_scal dimension //////////////////////////////////"<< vel_scal << std::endl;
	        vel_now = vel_scal;
	        std::cout << "/////////////////////////// the OnlineVelocity method //////////////////////////////////"<< std::endl;
	}
};

/*---------------------------------------------------------------------------*\
                               Starting the MAIN
\*---------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
    // Construct the tutorial object
    tutorial22 example(argc, argv);
    // Read some parameters from file
    ITHACAparameters* para = ITHACAparameters::getInstance(example.meshPtr(),
                             example._runTime());

    int NmodesUout = para->ITHACAdict->lookupOrDefault<int>("NmodesUout", 20);
    int NmodesPout = para->ITHACAdict->lookupOrDefault<int>("NmodesPout", 20);
    int NmodesSUPout = para->ITHACAdict->lookupOrDefault<int>("NmodesSUPout", 15);
    int NmodesUproj = para->ITHACAdict->lookupOrDefault<int>("NmodesUproj", 10);
    int NmodesPproj = para->ITHACAdict->lookupOrDefault<int>("NmodesPproj", 15);
    int NmodesSUPproj = para->ITHACAdict->lookupOrDefault<int>("NmodesSUPproj", 10);
   
    // Read the par file where the parameters are stored
    word filename("./par");
    example.mu = ITHACAstream::readMatrix(filename);
    // example.Pnumber = 1; // Number of parameters.
    // example.Tnumber = 1; //Dimension of the training set (used only when gerating parameters without input)
    // example.setParameters();
    // // Set the parameter ranges: Range of the parameter spaces. 
    // example.mu_range(0, 0) = 0.005;
    // example.mu_range(0, 1) = 0.005;
    //   // Generate equispaced samples inside the parameter range
    // example.genEquiPar();
    //Set the inlet boundaries patch 0 directions x and y
    example.inletIndex.resize(1, 2);
    example.inletIndex(0, 0) = 0;
    example.inletIndex(0, 1) = 0;
    // Time parameters: W e can use Ioodictionnary to access time parameters
    example.startTime = 0;
    example.finalTime = 2;
    example.timeStep = 0.005;
    example.writeEvery = 0.1;

    //Perform the offline solve
    example.offlineSolve();
    //Search the lift function
    example.liftSolve3();
    //Info << example.liftfield << endl;
    ITHACAutilities::normalizeFields(example.liftfield);
    // Homogenize the snapshots
    example.computeLift(example.Ufield, example.liftfield, example.Uomfield);
    //ITHACAstream::read_fields(example.liftfield, example.U, "./lift/");
    ITHACAPOD::getModes(example.Uomfield, example.Umodes, example._U().name(),
                            example.podex, 0, 0, NmodesUout);


    // ITHACAPOD::getModes(example.Ufield, example.Umodes, example._U().name(),
    //                         example.podex, 0, 0, NmodesUout);

    // //Perform POD on velocity pressure and supremizers and store the first 10 modes
    ITHACAPOD::getModes(example.Pfield, example.Pmodes, example._p().name(),
                        example.podex, 0, 0,
                        NmodesPout);



    reducedFSI reduced(example);

 //    scalar mu_now = example.mu(0, 0);
 //    std::cout << "///////////////////////////value of mu_now //////////////////////////////////"<< mu_now<< std::endl;
 //    example.change_viscosity(mu_now);
 //    std::cout << "///////////////////////////value of mu_now after //////////////////////////////////"<< mu_now<< std::endl;
    // Reads inlet volocities boundary conditions.
    word vel_file(para->ITHACAdict->lookup("online_velocities"));
    std::cout << "///////////////////////////size  of vel_file //////////////////////////////////"<<vel_file.size() << std::endl;
    Eigen::MatrixXd vel = ITHACAstream::readMatrix(vel_file);
   
 //    //reduced.OnlineVelocity(vel);
	// //   reducedFSI.nu = 0.005;
 //    //   reducedFSI.tstart = 50;
 //   //    reducedFSI.finalTime = 70;
 //   //    reduced.dt = 0.005;
 //   //    reducedFSI.storeEvery = 0.005;
 //   //    reducedFSI.exportEvery = 0.1;
 //   // Set the online velocity
 //   Eigen::MatrixXd vel(1, 1);
 //   vel(0,0) = 1;
 //   //reduced.solveOnline_sup(vel_now, 1);
 //    //std::cout << "///////////////////////////value of mu_now after OnlineVelocity call //////////////////////////////////"<< mu_now<< std::endl;

 //    reduced.solveOnline_Pimple(vel(0,0), NmodesUproj, NmodesPproj);
    
    //Perform the online solutions
    //for (label k = 0; k < (example.mu).size(); k++)
    //{
        //scalar mu_now = example.mu(0, k);
        scalar mu_now = example.mu(0, 0);
        example.change_viscosity(mu_now);
        reduced.OnlineVelocity(vel);
        reduced.solveOnline_Pimple(mu_now, NmodesUproj, NmodesPproj);
    //}
   
    exit(0);
}



