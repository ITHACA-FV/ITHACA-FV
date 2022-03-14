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
    22fsi.C
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
        : fsiBasic(argc, argv), U(_U()), p(_p()), phi(_phi())
    {}

    // Fields To Perform
    /// Velocity field
    volVectorField& U;
    /// Pressure field
    volScalarField& p;
    /// flux field
    surfaceScalarField& phi;
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

        //assignBC(U, BCind, inl);
        //assignIF(U, inl);
        //change_viscosity(mu(0, i));

        truthSolve3(mu_now);
        //}
        //}
    }
};


class reducedBasicFsi: public reducedSimpleSteadyNS
{
public:
    explicit reducedBasicFsi(tutorial22& FoamPb)
        : problem(&FoamPb)
    {}

    tutorial22* problem;

    //volVectorModes ULmodes;
    ///////// time control variables
    scalar startTime;
    scalar finalTime;
    scalar timeStep;
    scalar writeEvery = timeStep;
    scalar nextWrite;

    label counter = 1;


    void solveOnline_Pimple(scalar mu_now, int NmodesUproj, int NmodesPproj,
                            fileName folder = "./ITHACAoutput/Reconstruct/")
    {

        // Initializations
        Eigen::VectorXd uresidual = Eigen::VectorXd::Zero(NmodesUproj);
        Eigen::VectorXd presidual = Eigen::VectorXd::Zero(NmodesPproj);
        Eigen::MatrixXd a = Eigen::VectorXd::Zero(NmodesUproj);
        Eigen::MatrixXd b = Eigen::VectorXd::Zero(NmodesPproj);

        Time& runTime = problem->_runTime();
        surfaceScalarField& phi = problem->_phi();
        dynamicFvMesh& mesh = problem->meshPtr();
#include "initContinuityErrs.H"
        fv::options& fvOptions = problem->_fvOptions();
        pimpleControl& pimple = problem->_pimple();
     
        volScalarField& p = problem->_p();
        volVectorField& U = problem->_U();
        IOMRFZoneList& MRF = problem->_MRF();
        singlePhaseTransportModel& laminarTransport = problem->_laminarTransport();
        autoPtr<incompressible::turbulenceModel> turbulence;
        turbulence = autoPtr<incompressible::turbulenceModel>(incompressible::turbulenceModel::New(U, phi, laminarTransport));

        instantList Times = runTime.times();
        runTime.setEndTime(finalTime);
        runTime.setTime(Times[1], 1);
        runTime.setDeltaT(timeStep);
        nextWrite = startTime;
        label pRefCell = 0;
        scalar pRefValue = 0.0;

          /// first reconstruction  of velocity and pressure
        a = ITHACAutilities::getCoeffs(U, problem->Umodes, NmodesUproj, true);
        b = ITHACAutilities::getCoeffs(p, problem->Pmodes, NmodesPproj, true);

#include "createDyMControls.H"
#include "createUfIfPresent.H"
#include "CourantNo.H"
#include "setInitialDeltaT.H"

        turbulence->validate();
        problem->Umodes.reconstruct(U, a, "U");
        problem->Pmodes.reconstruct(p, b, "p");

        ITHACAstream::exportSolution(U, name(counter), folder);
        ITHACAstream::exportSolution(p, name(counter), folder);
        ITHACAstream::writePoints(mesh.points(), folder, name(counter) + "/polyMesh/");
        std::ofstream of(folder + name(counter) + "/" + runTime.timeName());
        counter++;
        nextWrite += writeEvery;

        // PIMPLE algorithm starts here
        Info<< "\nStarting time loop\n" << endl;
        while (runTime.run())
        {
            
#include "readDyMControls.H"
#include "CourantNo.H"
#include "setDeltaT.H"
            runTime.setEndTime(finalTime);
            runTime++;

            Info << "Time = " << runTime.timeName() << nl << endl;

            while (pimple.loop())
            {
                if (pimple.firstIter() || moveMeshOuterCorrectors)
                {
                    // Do any mesh changes
                    mesh.controlledUpdate();

                    if (mesh.changing())
                    {
                        MRF.update();

                        if (correctPhi)
                        {
                            // Calculate absolute flux
                            // from the mapped surface velocity
                            phi = mesh.Sf() & Uf();

#include "correctPhi.H"

                            // Make the flux relative to the mesh motion
                            fvc::makeRelative(phi, U);
                        }

                        if (checkMeshCourantNo)
                        {
#include "meshCourantNo.H"
                        }
                    }
                }

                //#include "UEqn.H"

                // Solve the Momentum equation
                MRF.correctBoundaryVelocity(U);

                tmp<fvVectorMatrix> tUEqn
                (
                    fvm::ddt(U) + fvm::div(phi, U)
                    + MRF.DDt(U)
                    + turbulence->divDevReff(U)
                    ==
                    fvOptions(U)
                );
                fvVectorMatrix& UEqn = tUEqn.ref();

                UEqn.relax();

                fvOptions.constrain(UEqn);
                // #############Galerkin projection for the velocity ###########################
                List<Eigen::MatrixXd> RedLinSysU;
                if (pimple.momentumPredictor())
                {
                    //solve(UEqn == -fvc::grad(p));
                    RedLinSysU = problem->Umodes.project(UEqn, NmodesUproj);
                    volVectorField gradpfull = -fvc::grad(p);
                    Eigen::MatrixXd projGrad = problem->Umodes.project(gradpfull, NmodesUproj);
                    RedLinSysU[1] = RedLinSysU[1] + projGrad;
                    a = reducedProblem::solveLinearSys(RedLinSysU, a, uresidual);
                    problem->Umodes.reconstruct(U, a, "U");
                    fvOptions.correct(U);
                }

                /////////////////////////////////////////////////////////

                // --- Pressure corrector loop
                while (pimple.correct())
                {
                    //#include "pEqn.H"

                    volScalarField rAU(1.0 / UEqn.A());
                    volVectorField HbyA(constrainHbyA(rAU * UEqn.H(), U, p)); //p
                    surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
                    if (pimple.ddtCorr())
                    {
                        phiHbyA += MRF.zeroFilter(fvc::interpolate(rAU) * fvc::ddtCorr(U, phi, Uf));
                    }
                    else
                    {
                        phiHbyA += MRF.zeroFilter(fvc::interpolate(rAU));
                    }

                    MRF.makeRelative(phiHbyA);

                    if (p.needReference())
                    {
                        fvc::makeRelative(phiHbyA, U);
                        adjustPhi(phiHbyA, U, p);
                        fvc::makeAbsolute(phiHbyA, U);
                    }

                    tmp<volScalarField> rAtU(rAU);

                    if (pimple.consistent())
                    {
                        rAtU = 1.0 / max(1.0 / rAU - UEqn.H1(), 0.1 / rAU);
                        phiHbyA +=
                            fvc::interpolate(rAtU() - rAU) * fvc::snGrad(p) * mesh.magSf(); // p
                        HbyA -= (rAU - rAtU()) * fvc::grad(p); //p
                    }

                    if (pimple.nCorrPISO() <= 1)
                    {
                        tUEqn.clear();
                    }

                    // Update the pressure BCs to ensure flux consistency
                    constrainPressure(p, U, phiHbyA, rAtU(), MRF); //p

                    // ### Reduced linear system for Pressure
                    List<Eigen::MatrixXd> RedLinSysP;

                    // Non-orthogonal pressure corrector loop
                    while (pimple.correctNonOrthogonal())
                    {
                        fvScalarMatrix pEqn
                        (
                            fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA) //p
                        );

                        pEqn.setReference(pRefCell, pRefValue);

                        //pEqn.solve(mesh.solver(P.select(pimple.finalInnerIter()))); //p
                        //// Added for reduced problem to project the pressure /////////
                        RedLinSysP = problem->Pmodes.project(pEqn, NmodesPproj);
                        b = reducedProblem::solveLinearSys(RedLinSysP, b, presidual);
                        problem->Pmodes.reconstruct(p, b, "p");

                        if (pimple.finalNonOrthogonalIter())
                        {
                            phi = phiHbyA - pEqn.flux();
                        }
                    }

#include "continuityErrs.H"

                    // Explicitly relax pressure for momentum corrector
                    p.relax();
                    U = HbyA - rAtU * fvc::grad(p); //p
                    U.correctBoundaryConditions();
                    fvOptions.correct(U);
                    // Correct Uf if the mesh is moving
                    fvc::correctUf(Uf, U, phi);

                    // Make the fluxes relative to the mesh motion
                    fvc::makeRelative(phi, U);

                }// end of the pimple.correct()

            }// end of the pimple.loop()

            ITHACAstream::exportSolution(U, name(counter), folder);
            ITHACAstream::exportSolution(p, name(counter), folder);
            ITHACAstream::writePoints(mesh.points(), folder, name(counter) + "/polyMesh/");
            std::ofstream of(folder + name(counter) + "/" + runTime.timeName());
            counter++;
        } // end of the runTime.run()

    } // end of the method

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

    int NmodesUout = para->ITHACAdict->lookupOrDefault<int>("NmodesUout", 50);
    int NmodesPout = para->ITHACAdict->lookupOrDefault<int>("NmodesPout", 50);
    //int NmodesSUPout = para->ITHACAdict->lookupOrDefault<int>("NmodesSUPout", 15);
    int NmodesUproj = para->ITHACAdict->lookupOrDefault<int>("NmodesUproj", 10);
    int NmodesPproj = para->ITHACAdict->lookupOrDefault<int>("NmodesPproj", 10);
    //int NmodesSUPproj = para->ITHACAdict->lookupOrDefault<int>("NmodesSUPproj", 10);

    // Read the par file where the parameters are stored
    word filename("./par");
    example.mu = ITHACAstream::readMatrix(filename);
    // // example.Pnumber = 1; // Number of parameters.
    // // example.Tnumber = 1; //Dimension of the training set (used only when gerating parameters without input)
    // // example.setParameters();
    // // // Set the parameter ranges: Range of the parameter spaces.
    // // example.mu_range(0, 0) = 0.005;
    // // example.mu_range(0, 1) = 0.005;
    // //   // Generate equispaced samples inside the parameter range
    // // example.genEquiPar();

    //Set the inlet boundaries patch 0 directions x and y
    example.inletIndex.resize(1, 2);
    example.inletIndex(0, 0) = 0;
    example.inletIndex(0, 1) = 0;
    // Time parameters: We can use Ioodictionnary to access time parameters
    example.startTime = 0;
    example.finalTime = 5;
    example.timeStep = 0.005; //0.01;
    example.writeEvery = 0.005;

    // //Perform the offline solve
    example.offlineSolve();
    // //Search the lift function
    // //example.liftSolve3();
    // // Normalize the lifting function
    // //ITHACAutilities::normalizeFields(example.liftfield);
    // //Create homogeneous basis functions for velocity
    // //example.computeLift(example.Ufield, example.liftfield, example.Uomfield);
    // // Perform a POD decomposition for velocity and pressure
    // // ITHACAPOD::getModes(example.Uomfield, example.Umodes, example._U().name(),
    // //                         example.podex, 0, 0, NmodesUout);

    //Perform POD on velocity pressure store the first 20 modes

    ITHACAPOD::getModes(example.Ufield, example.Umodes, example._U().name(),
                        example.podex, 0, 0, NmodesUout);


    ITHACAPOD::getModes(example.Pfield, example.Pmodes, example._p().name(),
                        example.podex, 0, 0,
                        NmodesPout);
    // ################ Restart method is necessary #############################
    example.restart();

    // ############### contruct the reduced the class object ###################
    reducedBasicFsi reduced(example);
    // Reads inlet volocities boundary conditions.
    word vel_file(para->ITHACAdict->lookup("online_velocities"));
    Eigen::MatrixXd vel = ITHACAstream::readMatrix(vel_file);

    reduced.startTime = 0;
    reduced.finalTime = 5;
    reduced.timeStep = 0.005;
    reduced.writeEvery = 0.005;
    //reduced.nextStore = 0.1;
    //reduced.exportEvery = 0.005;

    //Perform the online solutions
    //for (label k = 0; k < (example.mu).size(); k++)
    //{
    //scalar mu_now = example.mu(0, k);
    scalar mu_now = example.mu(0, 0);
    //example.change_viscosity(mu_now);
    //reduced.OnlineVelocity(vel);
    reduced.solveOnline_Pimple(mu_now, NmodesUproj, NmodesPproj);
    //}

    // Eigen::MatrixXd errFrobU = ITHACAutilities::errorFrobRel(example.Ufield, reduced.U);
    // Eigen::MatrixXd errFrobP =  ITHACAutilities::errorFrobRel(example.Pfield,reduced.p);
    
    // ITHACAstream::exportMatrix(errFrobU, "errFrobU", "matlab","./ITHACAoutput/ErrorsFrob/");
    // ITHACAstream::exportMatrix(errFrobP, "errFrobP", "matlab","./ITHACAoutput/ErrorsFrob/");
   
    // Eigen::MatrixXd errL2U = ITHACAutilities::errorL2Rel(example.Ufield,reduced.U);
    // Eigen::MatrixXd errL2P =  ITHACAutilities::errorL2Rel(example.Pfield,reduced.p);
    
    // ITHACAstream::exportMatrix(errL2U, "errL2U", "matlab","./ITHACAoutput/ErrorsL2/");
    // ITHACAstream::exportMatrix(errL2P, "errL2P", "matlab","./ITHACAoutput/ErrorsL2/");

    exit(0);
}


//////////////////////////////////////////////////////////////////////
