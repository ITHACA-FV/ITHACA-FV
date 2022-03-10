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
        : fsiBasic(argc, argv), U(_U()), p(_p()), phi(_phi())
    {}

    // Fields To Perform
    /// Velocity field
    volVectorField& U;
    /// Pressure field
    volScalarField& p;
    /// flux field
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
        restart();
        //}
        //}
    }
};


class reducedBasicFsi: public reducedSimpleSteadyNS
{
public:
    explicit reducedBasicFsi(tutorial22& FoamPb)
        : problem(&FoamPb)
    {
        //std::cout << " #############" << "problem->Umodes.size() is: " << problem->Umodes.size() << std::endl;
        //startTime = problem->startTime;
        finalTime = problem->finalTime;
        timeStep = problem->timeStep;
        writeEvery = timeStep;

        //#include "setRootCase.H"

        for (int i = 0; i < problem->Umodes.size(); i++)
        {

            ULmodes.append(problem->Umodes.toPtrList()[i].clone());
            //std::cout << "############# the i ULmodes is " << ULmodes(i) << std::endl;
        }

        //std::cout << timeStep << std::endl;
        //std::cout << "The size of P is: " <<  p.size()<<  "##################" << std::endl;

        std::cout << "################ ctor of reducedBasicFsi ##################" << std::endl;
    }

    tutorial22* problem;

    //turbulence = problem->turbulence();

    volVectorModes ULmodes;
    ///////// time control variables
    scalar startTime;
    scalar finalTime;
    scalar timeStep;
    scalar writeEvery;
    scalar nextWrite;
    // label pRefCell = 0;
    // scalar pRefValue = 0.0;
    // scalar cumulativeContErr = 0.0;
    label counter = 1;


    void solveOnline_Pimple(scalar mu_now, int NmodesUproj, int NmodesPproj,
                            fileName folder = "./ITHACAoutput/Reconstruct/")
    {

        //ULmodes.resize(0);

        for (int i = 0; i < NmodesUproj; i++)
        {
            ULmodes.append((problem->Umodes.toPtrList()[i]).clone());
        }


        // Initializations
        Eigen::VectorXd uresidualOld = Eigen::VectorXd::Zero(NmodesUproj);
        Eigen::VectorXd presidualOld = Eigen::VectorXd::Zero(NmodesPproj);
        Eigen::VectorXd uresidual = Eigen::VectorXd::Zero(NmodesUproj);
        Eigen::VectorXd presidual = Eigen::VectorXd::Zero(NmodesPproj);
        scalar U_norm_res(1);
        scalar P_norm_res(1);
        Eigen::MatrixXd a = Eigen::VectorXd::Zero(NmodesUproj);
        Eigen::MatrixXd a0 = Eigen::VectorXd::Zero(NmodesUproj);
        Eigen::MatrixXd b = Eigen::VectorXd::Zero(NmodesPproj);
        Eigen::MatrixXd b0 = Eigen::VectorXd::Zero(NmodesPproj);
        Eigen::MatrixXd bOld = b;

        pimpleControl& pimple = problem->_pimple();
        IOMRFZoneList& MRF = problem->_MRF();
        fv::options& fvOptions = problem->_fvOptions();
        surfaceScalarField& phi = problem->_phi();
        volScalarField& p = problem->_p();
        volVectorField& U = problem->_U();
        dynamicFvMesh& mesh = problem->meshPtr();
        Time& runTime = problem->_runTime();
        singlePhaseTransportModel& laminarTransport = problem->_laminarTransport();
        //autoPtr<incompressible::turbulenceModel> turbulence;
        //volScalarField& p = problem->_p();
        //volVectorField& U = problem->_U();
        //surfaceScalarField& phi(problem->_phi());
        //volVectorField u2 = U;
        a = ITHACAutilities::getCoeffs(U, problem->Umodes, NmodesUproj, true);
        b = ITHACAutilities::getCoeffs(p, problem->Pmodes, NmodesPproj, true);
        //a(0) = a0(0); // Do not remove: it is not working without this condition
        //b = b0;
        float residualJumpLim = problem->para->ITHACAdict->lookupOrDefault<float>("residualJumpLim", 1e-5);
        float normalizedResidualLim = problem->para->ITHACAdict->lookupOrDefault<float>("normalizedResidualLim", 1e-5);
        scalar residual_jump(1 + residualJumpLim);

        //Time& runTime = problem->_runTime();
        instantList Times = runTime.times();
        runTime.setEndTime(finalTime);
        // runTime.setTime(Times[1], 1);
        const label startTime = 1;
        const label endTime = Times.size();
        Info << "endTime" << endTime << endl;
        Info << "runtime.run(): " << runTime.run() << endl;
        runTime.setDeltaT(timeStep);
        // runTime.setTime(Times[startTime], startTime);

        autoPtr<incompressible::turbulenceModel> turbulence;
        //singlePhaseTransportModel& laminarTransport = problem->_laminarTransport();
        turbulence = autoPtr<incompressible::turbulenceModel>
                     (
                         incompressible::turbulenceModel::New(U, phi, laminarTransport)
                     );

        label pRefCell = 0;
        scalar pRefValue = 0.0;
        scalar cumulativeContErr = 0.0;
        /// first reconstruction  of velocity and pressure
        problem->Umodes.reconstruct(U, a, "U");
        problem->Pmodes.reconstruct(p, b, "p");

        ITHACAstream::exportSolution(U, name(counter), folder);
        ITHACAstream::exportSolution(p, name(counter), folder);
        ITHACAstream::writePoints(mesh.points(), folder, name(counter) + "/polyMesh/");
        std::ofstream of(folder + name(counter) + "/" + runTime.timeName());
        counter++;

        //Project method return A list of Eigen Matrices of dimension 2.
        //The first element of the list is the reduced matrix of the linear system,
        //the second element is the reduced source term of the linear system.
        //#include "addCheckCaseOptions.H"
#include "createDyMControls.H"
#include "createUfIfPresent.H"
#include "CourantNo.H"
#include "setInitialDeltaT.H"
        // PIMPLE algorithm starts here
        while (runTime.loop())
        {
            std::cerr << "File: 22fsi.C, Line: 230"<< std::endl;
            //runTime.setTime(Times[i], i);

            //Info<< "Time = " << runTime.timeName() << endl;

            //while (runTime.run())
            //{
            //#include "readDyMControls.H"
            bool adjustTimeStep = runTime.controlDict().getOrDefault("adjustTimeStep", false);
            scalar maxCo = runTime.controlDict().getOrDefault<scalar>("maxCo", 1);
            scalar maxDeltaT = runTime.controlDict().getOrDefault<scalar>("maxDeltaT", GREAT);

            bool correctPhi(pimple.dict().getOrDefault("correctPhi", mesh.dynamic()));
            bool checkMeshCourantNo(pimple.dict().getOrDefault("checkMeshCourantNo", false));
            bool moveMeshOuterCorrectors(pimple.dict().getOrDefault("moveMeshOuterCorrectors", false));
#include "CourantNo.H"
#include "setDeltaT.H"
            //runTime.setEndTime(finalTime);
            //++runTime;
            //runTime++;
            // runTime.setTime(Times[i], i);

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
                            //CorrectPhi(U, phi, P, dimensionedScalar("rAUf", dimTime, 1), geometricZeroField(),pimple);

                            // Make the flux relative to the mesh motion
                            fvc::makeRelative(phi, U);
                        }

                        if (checkMeshCourantNo)
                        {
#include "meshCourantNo.H"
                        }
                    }
                }

                //vector v(1, 0, 0);
                // ITHACAutilities::assignBC(U, 0, v);
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
                    //std::cout << "#############"<< U_norm_res << "###################"<< std::endl;
                    //solve(UEqn == -fvc::grad(p));
                    RedLinSysU = problem->Umodes.project(UEqn, NmodesUproj);
                    volVectorField gradpfull = -fvc::grad(p);
                    Eigen::MatrixXd projGrad = problem->Umodes.project(gradpfull, NmodesUproj);
                    RedLinSysU[1] = RedLinSysU[1] + projGrad;
                    //Project method return A list of Eigen Matrices of dimension 2.
                    //The first element of the list is the reduced matrix of the linear system,
                    //the second element is the reduced source term of the linear system.
                    //RedLinSysU[1] = RedLinSysU[1] - projGradModP * b; //?

                    a = reducedProblem::solveLinearSys(RedLinSysU, a, uresidual);

                    problem->Umodes.reconstruct(U, a, "U");
                    fvOptions.correct(U);
                }

                /////////////////////////////////////////////////////////

                // --- Pressure corrector loop
                while (pimple.correct())
                {
                    //#include "pEqn.H"
                    //std::cout << "#############"<< P_norm_res << "###################"<< std::endl;
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

            // uresidualOld = uresidualOld - uresidual;
            // presidualOld = presidualOld - presidual;
            // uresidualOld = uresidualOld.cwiseAbs();
            // presidualOld = presidualOld.cwiseAbs();
            // residual_jump = std::max(uresidualOld.sum(), presidualOld.sum());
            // uresidualOld = uresidual;
            // presidualOld = presidual;
            // uresidual = uresidual.cwiseAbs();
            // presidual = presidual.cwiseAbs();
            //runTime.write();
            problem->Umodes.reconstruct(U, a, "U");
            problem->Pmodes.reconstruct(p, b, "p");

            ITHACAstream::exportSolution(U, name(counter), folder);
            ITHACAstream::exportSolution(p, name(counter), folder);
            ITHACAstream::writePoints(mesh.points(), folder, name(counter) + "/polyMesh/");
            std::ofstream of(folder + name(counter) + "/" + runTime.timeName());
            //runTime.setTime(runTime.startTime(), 0);
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
    example.finalTime = 0.5;
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

    example.restart();
    // ############### contruct the reduced the class object ##########
    reducedBasicFsi reduced(example);
    //std::cout << "########################" << reduced.startTime << "#####################" << std::endl;
    // Reads inlet volocities boundary conditions.
    word vel_file(para->ITHACAdict->lookup("online_velocities"));
    Eigen::MatrixXd vel = ITHACAstream::readMatrix(vel_file);

    reduced.startTime = example.startTime;
    reduced.finalTime = example.finalTime;
    reduced.timeStep = example.timeStep;
    reduced.writeEvery = example.writeEvery;
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
    //reduced.solveOnline_Pimple(mu_now, NmodesUproj, NmodesPproj);
    //}

    exit(0);
}


//////////////////////////////////////////////////////////////////////
