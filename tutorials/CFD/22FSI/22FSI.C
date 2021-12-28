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
    Example of steady NS Reduction Problem solved by the use of the SIMPLE algorithm
SourceFiles
    22FSI.C
\*---------------------------------------------------------------------------*/

//#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H"

#include "unsteadyNS.H"
#include "ITHACAstream.H"
#include "ITHACAutilities.H"
#include "ITHACAPOD.H"
#include "forces.H"
#include "IOmanip.H"
#include "IOstreams.H"
//#include "Switch.H"
//#include "objectRegistry.H"
#include "localEulerDdtScheme.H"
#include "ReducedUnsteadyNS.H"
#include "SteadyNSSimple.H"
#include "ReducedSimpleSteadyNS.H"
#include "reductionProblem.H"

class tutorial22 : public unsteadyNS
{
public:
    /// Constructors
    tutorial22();
    explicit tutorial22(int argc, char* argv[])
        :
        unsteadyNS(argc, argv),
        U(_U()),
        p(_p()),
        phi(_phi())
    {
        _args = autoPtr<argList>
                (
                    new argList(argc, argv)
                );


        if (!_args->checkRootCase())
        {
            Foam::FatalError.exit();
        }

        argList& args = _args();

        #include "createTime.H"
        //#include "createMesh.H"
        #include "createDynamicFvMesh.H"
        _pimple = autoPtr<pimpleControl>
                  (
                      new pimpleControl
                      (
                          mesh
                      )
                  );
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

        // #include "createFields.H"
        // #include "createFvOptions.H"
        para = ITHACAparameters::getInstance(mesh, runTime);
        
    }// end of the Constructor.


    // members data
    //para = ITHACAparameters::getInstance(mesh, runTime);
    /// Velocity field
    volVectorField& U;
    /// Pressure field
    volScalarField& p;
    ///
    surfaceScalarField& phi;
    //surfaceVectorField& Uf;
    //pointDisplacement &pd;

    //Foam::autoPtr<Foam::fvMesh> _mesh;
    Foam::autoPtr<Foam::dynamicFvMesh> meshPtr;

    int folderN;
    int saver;
    int middleStep;
    int middleExport;

    /// Perform an Offline solve
    /// Reimplement this starting from pimplefoam
    void offlineSolve()
    {
        //std::cerr << "File: 22FSI.C, Line: 135"<< std::endl;
        Vector<double> inl(0, 0, 0); //inl(0, 0, 0)
        List<scalar> mu_now(1);

        // if the offline solution is already performed read the fields
        // if (offline)
        //  {
        //     ITHACAstream::read_fields(Ufield, U, "./ITHACAoutput/Offline/");
        //     ITHACAstream::read_fields(Pfield, p, "./ITHACAoutput/Offline/");
        //     // mu_samples = ITHACAstream::readMatrix("./ITHACAoutput/Offline/mu_samples_mat.txt");
        // }
        // else
        // {
            Vector<double> Uinl(1, 0, 0);
            label BCind = 0;

            for (label i = 0; i < mu.cols(); i++)
            {
                mu_now[0] = mu(0, i);
                // change_viscosity(mu(0, i));
                // assignIF(U, Uinl);
                truthSolve3(mu_now);
                //restart();
            }
        //}
    }

    void truthSolve3(List<scalar> mu_now, word Folder = "/ITHACAoutput/Offline")
    {
        Time& runTime = _runTime();
        volScalarField& p = _p();
        volVectorField& U = _U();
        dynamicFvMesh& mesh = meshPtr();
        surfaceScalarField& phi = _phi();
        #include "initContinuityErrs.H"
        fv::options& fvOptions = _fvOptions();
        pimpleControl& pimple = _pimple();
        IOMRFZoneList& MRF = _MRF();
        singlePhaseTransportModel& laminarTransport = _laminarTransport();
        /*autoPtr<incompressible::turbulenceModel> _turbulence
        (
             incompressible::turbulenceModel::New(U, phi, _laminarTransport());
        );*/
        instantList Times = runTime.times();
        runTime.setEndTime(finalTime);
        // Perform a TruthSolve
        runTime.setTime(Times[1], 1);
        runTime.setDeltaT(timeStep);
        nextWrite = startTime;

        scalar residual = 1;
        scalar uresidual = 1;
        Vector<double> uresidual_v(0, 0, 0);
        scalar presidual = 1;
        scalar csolve = 0;
        turbulence->read();
        std::ofstream res_os;
        std::ofstream snaps_os;
        std::ofstream iters;
        std::ofstream res_U;
        std::ofstream res_P;
        res_os.open(Folder + "/residuals", std::ios_base::app);
        snaps_os.open(Folder + "/snaps", std::ios_base::app);
        iters.open(Folder + "/iters", std::ios_base::app);
        res_U.open(Folder + name(counter) + "/res_U", std::ios_base::app);
        res_P.open(Folder + name(counter) + "/res_P", std::ios_base::app);

        folderN = 0;
        saver = 0;

        // middleStep = para->ITHACAdict->lookupOrDefault<label>("middleStep", 20);
        // middleExport = para->ITHACAdict->lookupOrDefault<bool>("middleExport", true);

        //***************pimpleFoam algorithm***************
       
        //#include "postProcess.H"

        #include "addCheckCaseOptions.H"
        //#include "setRootCaseLists.H"
        //#include "createTime.H" //already done
        
        #include "initContinuityErrs.H"
        #include "createDyMControls.H"
        //#include "createFields.H"
        #include "createUfIfPresent.H"
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"

        turbulence->validate();
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        Info << "\nStarting time loop\n" << endl;

        while (runTime.run())
        {
            #include "readDyMControls.H"
            #include "CourantNo.H"
            #include "setDeltaT.H"

            ++runTime;

            Info << "Time = " << runTime.timeName() << nl << endl;

            // --- Pressure-velocity PIMPLE corrector loop
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

                            // ************ end of include correctPhi.H **********************

                            // Make the flux relative to the mesh motion
                            fvc::makeRelative(phi, U);
                        }

                        if (checkMeshCourantNo)
                        {
                            #include "meshCourantNo.H"
                        }
                    }
                }

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
            //*****************************************end of pimpleFoam*****************************************
            snaps_os << folderN + 1 << std::endl;
            iters << csolve << std::endl;
            res_os << residual << std::endl;
            res_os.close();
            res_U.close();
            res_P.close();
            snaps_os.close();
            iters.close();
            runTime.setTime(runTime.startTime(), 0);

            // if (middleExport)
            // {
            //     ITHACAstream::exportSolution(U, name(folderN + 1), Folder + name(counter));
            //     ITHACAstream::exportSolution(p, name(folderN + 1), Folder + name(counter));
            // }
            // else
            // {
            //     ITHACAstream::exportSolution(U, name(counter), Folder);
            //     ITHACAstream::exportSolution(p, name(counter), Folder);
            // }

            /*if (ITHACAutilities::isTurbulent())
            {
                auto nut = mesh.lookupObject<volScalarField>("nut");
                ITHACAstream::exportSolution(nut, name(folderN + 1), Folder + name(counter));
                nutFields.append(nut.clone());
            }*/

            Ufield.append(U.clone());
            Pfield.append(p.clone());
            counter++;
            writeMu(mu_now);
            // --- Fill in the mu_samples with parameters (mu) to be used for the POD sample points
            // mu_samples.conservativeResize(mu_samples.rows() + 1, mu_now.size());

            // for (label i = 0; i < mu_now.size(); i++)
            // {
            //     mu_samples(mu_samples.rows() - 1, i) = mu_now[i];
            // }

            // Resize to Unitary if not initialized by user (i.e. non-parametric problem)
            if (mu.cols() == 0)
            {
                mu.resize(1, 1);
            }

            // if (mu_samples.rows() == mu.cols())
            // {
            //     ITHACAstream::exportMatrix(mu_samples, "mu_samples", "eigen",
            //                                Folder);
            // }

        }
    }

};

// reduced class problem

// class reducedFSI : public reducedSimpleSteadyNS
// {
//     public:
//         /// Constructor
//         //reducedFSI();
//         explicit reducedFSI(tutorial22& FOMproblem)
//             :
//             problem(&FOMproblem)
//         {}

//         /// Full problem.
//         tutorial22* problem;

//         // Function to perform the online phase
//         void solveOnline_Pimple(scalar mu_now,
//                                 int NmodesUproj, int NmodesPproj, int NmodesNut = 0,
//                                 word Folder = "./ITHACAoutput/Reconstruct/")
//         {
//             counter++;

//             // For all the variables, in case one wants to use all the available modes it is just necessary to set
//             // the requested number into the ITHACAdict to zero
//             if (NmodesUproj == 0)
//             {
//                 NmodesUproj = problem->Umodes.size();
//             }

//             if (NmodesPproj == 0)
//             {
//                 NmodesPproj = problem->Pmodes.size();
//             }

//            /* if (NmodesNut == 0)
//             {
//                 NmodesNut = problem->nutModes.size();
//             }*/

//             // Initializations
//             Eigen::VectorXd uresidualOld = Eigen::VectorXd::Zero(NmodesUproj);
//             Eigen::VectorXd presidualOld = Eigen::VectorXd::Zero(NmodesPproj);
//             Eigen::VectorXd uresidual = Eigen::VectorXd::Zero(NmodesUproj);
//             Eigen::VectorXd presidual = Eigen::VectorXd::Zero(NmodesPproj);
//             scalar U_norm_res(1);
//             scalar P_norm_res(1);
//             Eigen::MatrixXd a = Eigen::VectorXd::Zero(NmodesUproj);
//             Eigen::MatrixXd a0 = Eigen::VectorXd::Zero(NmodesUproj);
//             Eigen::MatrixXd b = Eigen::VectorXd::Zero(NmodesPproj);
//             Eigen::MatrixXd b0 = Eigen::VectorXd::Zero(NmodesPproj);
//             Eigen::MatrixXd bOld = b;
//             //Eigen::MatrixXd nutCoeff = Eigen::VectorXd::Zero(NmodesNut);
//             //Eigen::MatrixXd nutCoeffOld = Eigen::VectorXd::Zero(NmodesNut);
//             float residualJumpLim =
//                 problem->para->ITHACAdict->lookupOrDefault<float>("residualJumpLim", 1e-5);
//             float normalizedResidualLim =
//                 problem->para->ITHACAdict->lookupOrDefault<float>("normalizedResidualLim",
//                         1e-5);
//             scalar residual_jump(1 + residualJumpLim);
//             volScalarField& P = problem->_p();
//             volVectorField& U = problem->_U();
//             //volScalarField& nut = const_cast<volScalarField&>
//             //                     (problem->_mesh().lookupObject<volScalarField>("nut"));
//             volVectorField u2 = U;

//             // Getting coefficients
//             a0 = ITHACAutilities::getCoeffs(U, problem->Umodes, NmodesUproj, true);
//             b = ITHACAutilities::getCoeffs(P, problem->Pmodes, NmodesPproj, true);
//             //nutCoeff = ITHACAutilities::getCoeffs(nut, problem->nutModes, NmodesNut, true);


//             a(0) = a0(0); // Do not remove: it is not working without this condition
//             b = b0;
//             fvMesh& mesh = problem->_mesh();
//             Time& runTime = problem->_runTime();
//             P.rename("p");
//             surfaceScalarField& phi(problem->_phi());

//             problem->Umodes.reconstruct(U, a, "U");
//             problem->Pmodes.reconstruct(P, b, "p");
//             vector v(1, 0, 0);
//             ITHACAutilities::assignBC(U, 0, v);
//             //problem->nutModes.reconstruct(nut, nutCoeff, "nut");
//             phi = fvc::flux(U);
//             int iter = 0;
//             label pRefCell = 0;
//             scalar pRefValue = 0.0;
//             pimpleControl& pimple = problem->_pimple();
//             /*_turbulence autoPtr<incompressible::turbulenceModel> turbulence
//             (
//                    incompressible::turbulenceModel::New(U, phi, laminarTransport)
//             );*/
//             IOMRFZoneList& MRF = problem->_MRF();
//             fv::options& fvOptions = problem->_fvOptions();
//             std::ofstream res_os_U;
//             std::ofstream res_os_P;
//             res_os_U.open(Folder + name(counter) + "/residualsU", std::ios_base::app);
//             res_os_P.open(Folder + name(counter) + "/residualsP", std::ios_base::app);

//             // Pimple algorithm starts here
//             while ((residual_jump > residualJumpLim
//                     || std::max(U_norm_res, P_norm_res) > normalizedResidualLim) &&
//                     iter < maxIterOn)
//             {
//                 iter++;
//                 std::cout << iter << std::endl;

// #if defined(OFVER) && (OFVER == 6)
//                 pimple.loop(runTime);
// #else
//                 problem->_pimple().loop();
// #endif

//                 // If the case is turbulent, then the network is evaluated
//                 /*if (ITHACAutilities::isTurbulent())
//                 {
//                     nutCoeff = problem->evalNet(a, mu_now);
//                     //nutCoeff = nutCoeffOld + 0.7*(nutCoeff - nutCoeffOld);
//                     volScalarField& nut = const_cast<volScalarField&>
//                                           (problem->_mesh().lookupObject<volScalarField>("nut"));
//                     problem->nutModes.reconstruct(nut, nutCoeff, "nut");
//                     //ITHACAstream::exportSolution(nut, name(counter), Folder);
//                 }*/

//                 //volScalarField nueff = problem->turbulence->nuEff();
//                 vector v(1, 0, 0);
//                 ITHACAutilities::assignBC(U, 0, v);
//                 /*// Momentum equation
//                 fvVectorMatrix UEqn
//                 (
//                     fvm::div(phi, U)
//                     - fvm::laplacian(nueff, U)
//                     - fvc::div(nueff * dev2(T(fvc::grad(U))))
//                 );
//                 UEqn.relax();*/

//                 // Solve the Momentum equation

//                 MRF.correctBoundaryVelocity(U);

//                 tmp<fvVectorMatrix> tUEqn
//                 (
//                 fvm::ddt(U) + fvm::div(phi, U)
//                 + MRF.DDt(U)
//                 + problem->turbulence->divDevReff(U)
//                 ==
//                 fvOptions(U)
//                 );
//                 fvVectorMatrix& UEqn = tUEqn.ref();

//                 UEqn.relax();
//                 List<Eigen::MatrixXd> RedLinSysU = ULmodes.project(UEqn, UprojN);
//                 RedLinSysU[1] = RedLinSysU[1] - projGradModP * b;
//                 a = reducedProblem::solveLinearSys(RedLinSysU, a, uresidual, vel_now);
//                 ULmodes.reconstruct(U, a, "U");
        
//                 fvOptions.constrain(UEqn);

//                 if (pimple.momentumPredictor())
//                 {
//                    solve(UEqn == -fvc::grad(P));
        
//                    fvOptions.correct(U);//??
//                 }
//                 //solve(UEqn == - fvc::grad(P));
//                 ITHACAutilities::assignBC(U, 0, v);
                
//               // begining of solving the pEqn.H   

//                 volScalarField rAU(1.0/UEqn.A());
//                 volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, P));
//                 surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));

//                 if (pimple.ddtCorr())
//                 {
//                     #include "createUfIfPresent.H"
//                     phiHbyA += MRF.zeroFilter(fvc::interpolate(rAU)*fvc::ddtCorr(U, phi, Uf));
//                 }
//                 else
//                 {
//                     phiHbyA += MRF.zeroFilter(fvc::interpolate(rAU));
//                 }

//                 MRF.makeRelative(phiHbyA);

//                 if (P.needReference())
//                 {
//                     fvc::makeRelative(phiHbyA, U);
//                     adjustPhi(phiHbyA, U, P);
//                     fvc::makeAbsolute(phiHbyA, U);
//                 }

//                 tmp<volScalarField> rAtU(rAU);

//                 if (pimple.consistent())
//                 {
//                     rAtU = 1.0/max(1.0/rAU - UEqn.H1(), 0.1/rAU);
//                     phiHbyA +=
//                         fvc::interpolate(rAtU() - rAU)*fvc::snGrad(P)*mesh.magSf();
//                     HbyA -= (rAU - rAtU())*fvc::grad(P);
//                 }
//                 List<Eigen::MatrixXd> RedLinSysP;
//                 bOld = b;

//                 if (pimple.nCorrPISO() <= 1)
//                 {
//                     tUEqn.clear();
//                 }

//                 // Update the pressure BCs to ensure flux consistency
//                 constrainPressure(P, U, phiHbyA, rAtU(), MRF);

//                 // Non-orthogonal pressure corrector loop
//                 while (pimple.correctNonOrthogonal())
//                 {
//                     // Continuity equation
//                     fvScalarMatrix pEqn
//                     (
//                         fvm::laplacian(rAtU(), P) == fvc::div(phiHbyA)
//                     );
//                     // Added for reduced problem
//                     RedLinSysP = problem->Pmodes.project(pEqn, NmodesPproj);
//                     //pEqn.solve();
//                     b = reducedProblem::solveLinearSys(RedLinSysP, b, presidual);
//                     problem->Pmodes.reconstruct(P, b, "p");
//                     pEqn.setReference(pRefCell, pRefValue);

//                     pEqn.solve(mesh.solver(P.select(pimple.finalInnerIter())));

//                     if (pimple.finalNonOrthogonalIter())
//                     {
//                         phi = phiHbyA - pEqn.flux();
//                     }
//                 }
//                 scalar cumulativeContErr = 0.0;

//                 #include "continuityErrs.H"

//                 // Explicitly relax pressure for momentum corrector
//                 //p.relax();

//                 /*U = HbyA - rAtU*fvc::grad(p);
//                 U.correctBoundaryConditions();*/
               

//                 // Correct Uf if the mesh is moving
//                 //#include "correctPhi.H"
//                 #include "createUfIfPresent.H"
//                 fvc::correctUf(Uf, U, phi);

//                 // Make the fluxes relative to the mesh motion
//                 fvc::makeRelative(phi, U);

//                 b = bOld + mesh.fieldRelaxationFactor("p") * (b - bOld);
//                 problem->Pmodes.reconstruct(P, b, "p");
//                 //nutCoeffOld = nutCoeff;
//                 // P.relax();
//                 U = HbyA - rAtU() * fvc::grad(P);
//                 U.correctBoundaryConditions();
//                 fvOptions.correct(U);

//                 uresidualOld = uresidualOld - uresidual;
//                 presidualOld = presidualOld - presidual;
//                 uresidualOld = uresidualOld.cwiseAbs();
//                 presidualOld = presidualOld.cwiseAbs();
//                 residual_jump = std::max(uresidualOld.sum(), presidualOld.sum());
//                 uresidualOld = uresidual;
//                 presidualOld = presidual;
//                 uresidual = uresidual.cwiseAbs();
//                 presidual = presidual.cwiseAbs();
//                 U_norm_res = uresidual.sum() / (RedLinSysU[1].cwiseAbs()).sum();
//                 P_norm_res = presidual.sum() / (RedLinSysP[1].cwiseAbs()).sum();

//                 if (problem->para->debug)
//                 {
//                     std::cout << "Residual jump = " << residual_jump << std::endl;
//                     std::cout << "Normalized residual = " << std::max(U_norm_res,
//                               P_norm_res) << std::endl;
//                     std::cout << "Final normalized residual for velocity: " << U_norm_res <<
//                               std::endl;
//                     std::cout << "Final normalized residual for pressure: " << P_norm_res <<
//                               std::endl;
//                 }

//                 res_os_U << U_norm_res << std::endl;
//                 res_os_P << P_norm_res << std::endl;
//             }

//             res_os_U.close();
//             res_os_P.close();
//             std::cout << "Solution " << counter << " converged in " << iter <<
//                       " iterations." << std::endl;
//             std::cout << "Final normalized residual for velocity: " << U_norm_res <<
//                       std::endl;
//             std::cout << "Final normalized residual for pressure: " << P_norm_res <<
//                       std::endl;
//             problem->Umodes.reconstruct(U, a, "Uaux");
//             problem->Pmodes.reconstruct(P, b, "Paux");

//             /*if (ITHACAutilities::isTurbulent())
//             {
//                 volScalarField& nut = const_cast<volScalarField&>
//                                       (problem->_mesh().lookupObject<volScalarField>("nut"));
//                 nut.rename("nutAux");
//         ITHACAstream::exportSolution(nut, name(counter), Folder);
//             }*/

//             ITHACAstream::exportSolution(U, name(counter), Folder);
//             ITHACAstream::exportSolution(P, name(counter), Folder);
//             runTime.setTime(runTime.startTime(), 0);
//         }
// };


int main(int argc, char* argv[])
{
    // Construct the tutorial object
    tutorial22 example(argc, argv);
    // Read some parameters from file
    ITHACAparameters* para = ITHACAparameters::getInstance(example._mesh(),
                            example._runTime());

    int NmodesUout = para->ITHACAdict->lookupOrDefault<int>("NmodesUout", 15);
    int NmodesPout = para->ITHACAdict->lookupOrDefault<int>("NmodesPout", 15);
    int NmodesSUPout = para->ITHACAdict->lookupOrDefault<int>("NmodesSUPout", 15);
    int NmodesUproj = para->ITHACAdict->lookupOrDefault<int>("NmodesUproj", 10);
    int NmodesPproj = para->ITHACAdict->lookupOrDefault<int>("NmodesPproj", 10);
    int NmodesSUPproj = para->ITHACAdict->lookupOrDefault<int>("NmodesSUPproj", 10);


    // Read the par file where the parameters are stored
    //word filename("./par");
    //example.mu = ITHACAstream::readMatrix(filename);
    // Set the inlet boundaries patch 0 directions x and y
    example.inletIndex.resize(1, 2);
    example.inletIndex(0, 0) = 0;
    example.inletIndex(0, 1) = 0;
    // Perform the offline solve
    example.offlineSolve();

    if (example.bcMethod == "lift")
    {
        ITHACAstream::read_fields(example.liftfield, example.U, "./lift/");
        ITHACAutilities::normalizeFields(example.liftfield);
        // Homogenize the snapshots
        example.computeLift(example.Ufield, example.liftfield, example.Uomfield);
        // Perform POD on the velocity snapshots
        ITHACAPOD::getModes(example.Uomfield, example.Umodes, example._U().name(),
                            example.podex, 0, 0, NmodesUout);
    }
    else
    {
        ITHACAPOD::getModes(example.Ufield, example.Umodes, example._U().name(),
                            example.podex, 0, 0, NmodesUout);
    }

    // Perform POD on velocity pressure and supremizers and store the first 10 modes
    ITHACAPOD::getModes(example.Pfield, example.Pmodes, example._p().name(),
                        example.podex, 0, 0,
                        NmodesPout);
    // ITHACAstream::read_fields(example.liftfield, example.U, "./lift/");
    // // // Homogenize the snapshots
    // //example.computeLift(example.Ufield, example.liftfield, example.Uomfield);

    // // Perform POD on velocity and pressure and store the first 10 modes
    // ITHACAPOD::getModes(example.Uomfield, example.Umodes, example._U().name(),
    //                     example.podex, 0, 0,
    //                      NmodesUout);
    //  ITHACAPOD::getModes(example.Pfield, example.Pmodes, example._p().name(),
    //                     example.podex, 0, 0,
    //                      NmodesPout);

    //  // Create the reduced object
    // reducedFSI reduced(example);
    
    // // Create the reduced object
    //reducedUnsteadyNS reduced(example);
    //reducedSimpleSteadyNS reduced(example);
    PtrList<volVectorField> U_rec_list;
    PtrList<volScalarField> P_rec_list;
    // // Reads inlet volocities boundary conditions.
    // word vel_file(para->ITHACAdict->lookup("online_velocities"));
    // Eigen::MatrixXd vel = ITHACAstream::readMatrix(vel_file);

    // //Perform the online solutions
    // for (label k = 0; k < (example.mu).size(); k++)
    // {
    //     scalar mu_now = example.mu(0, k);
    //     example.change_viscosity(mu_now);
    //     reduced.setOnlineVelocity(vel);
    //     reduced.solveOnline_Simple(mu_now, NmodesUproj, NmodesPproj);
    // }

    exit(0);
}
