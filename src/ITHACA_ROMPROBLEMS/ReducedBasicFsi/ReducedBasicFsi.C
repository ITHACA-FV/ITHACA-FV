#include "ReducedBasicFsi.H"


/// Constructor
ReducedBasicFsi() {}
explicit ReducedBasicFsi(fsiBasic& FOMproblem)
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
    
}


// Function to perform the online phase
void ReducedBasicFsi::solveOnline_Pimple(scalar mu_now, int NmodesUproj, int NmodesPproj, int NmodesSup = 0, fileName  folder = "./ITHACAoutput/Reconstruct/")
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

    // Initializations
    Eigen::VectorXd uresidualOld = Eigen::VectorXd::Zero(UprojN);
    Eigen::VectorXd presidualOld = Eigen::VectorXd::Zero(PprojN);
    Eigen::VectorXd uresidual = Eigen::VectorXd::Zero(UprojN);
    Eigen::VectorXd presidual = Eigen::VectorXd::Zero(PprojN);
    scalar U_norm_res(1);
    scalar P_norm_res(1);
    // coefficients declaration
    Eigen::MatrixXd a = Eigen::VectorXd::Zero(UprojN);
    Eigen::MatrixXd b = Eigen::VectorXd::Zero(PprojN);

    a(0) = vel_now(0, 0);
    float residualJumpLim = problem->para->ITHACAdict->lookupOrDefault<float>("residualJumpLim", 1e-5);
    float normalizedResidualLim = problem->para->ITHACAdict->lookupOrDefault<float>("normalizedResidualLim", 1e-5);
    scalar residual_jump(1 + residualJumpLim);
    volScalarField& P = problem->_p();
    volVectorField& U = problem->_U();
    dynamicFvMesh& mesh = problem->meshPtr();
    Time& runTime = problem->_runTime();
    P.rename("p");
    surfaceScalarField& phi(problem->_phi());
    problem->Umodes.reconstruct(U, a, "U");
    problem->Pmodes.reconstruct(P, b, "p");

    //phi = fvc::flux(U);
    // ### Create Uf field ###########
    autoPtr<surfaceVectorField> Uf;

    if (mesh.dynamic())
    {
        //Info<< "Constructing face velocity Uf\n" << endl;

        Uf.reset
        (
            new surfaceVectorField
            (
                IOobject
                (
                    "Uf",
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                fvc::interpolate(U)
            )
        );
    }

    //surfaceScalarField Uf = fvc::interpolate(U);
    //surfaceScalarField Uf = fvc::flux(U);
    int iter = 0;
    label pRefCell = 0;
    scalar pRefValue = 0.0;
    scalar cumulativeContErr = 0.0;
    scalar meshCoNum = 0.0;
    scalar meanMeshCoNum = 0.0;
    pimpleControl& pimple = problem->_pimple();
    IOMRFZoneList& MRF = problem->_MRF();
    fv::options& fvOptions = problem->_fvOptions();
    PtrList<volVectorField> gradModP;

    for (int i = 0; i < NmodesPproj; i++)
    {
        gradModP.append(fvc::grad(problem->Pmodes[i]));
    }

    projGradModP = ULmodes.project(gradModP, NmodesUproj);
    // ####### Starting of the reduced algorithm ##############
    time = tstart;

    while (time < finalTime ) //runTime.run()
    {
            time = time + dt;
            //runTime++;
            while ((residual_jump > residualJumpLim || std::max(U_norm_res, P_norm_res) > normalizedResidualLim) && iter < maxIterOn)
            {
                iter++;
                //std::cout << "///////////////////////////solveOnline_Pimple method line 249: iter value is:  //////////////////////////////////"<< iter <<  std::endl;

                //std::cout << iter << std::endl;
                // OFVER = Openfoam Version
                #if defined(OFVER) && (OFVER == 6)
                   pimple.loop(runTime);
                 //std::cout << "///////////////////////////solveOnline_Pimple method line 251 after pimple.loop(runTime) //////////////////////////////////"<< std::endl;
                #else
                   pimple.loop();
                 //std::cout << "///////////////////////////solveOnline_Pimple method line 254 after pimple.loop() //////////////////////////////////"<< std::endl;
                #endif

                bool correctPhi(pimple.dict().getOrDefault("correctPhi", mesh.dynamic()));

                bool checkMeshCourantNo(pimple.dict().getOrDefault("checkMeshCourantNo", false));

                bool moveMeshOuterCorrectors(pimple.dict().getOrDefault("moveMeshOuterCorrectors", false));

                // #################### Check the deformation of the mesh properties. #################

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

                            //#include "correctPhi.H"
                            CorrectPhi(U,phi,P,dimensionedScalar("rAUf", dimTime, 1), geometricZeroField(),pimple);

                            //#include "continuityErrs.H"
                            {
                                volScalarField contErr(fvc::div(phi));

                                scalar sumLocalContErr = runTime.deltaTValue()*
                                    mag(contErr)().weightedAverage(mesh.V()).value();

                                scalar globalContErr = runTime.deltaTValue()*
                                    contErr.weightedAverage(mesh.V()).value();
                                cumulativeContErr += globalContErr;

                                // Info<< "time step continuity errors : sum local = " << sumLocalContErr
                                //     << ", global = " << globalContErr
                                //     << ", cumulative = " << cumulativeContErr
                                //     << endl;
                            }

                            // Make the flux relative to the mesh motion
                            fvc::makeRelative(phi, U);
                        }

                         if (checkMeshCourantNo)
                        {
                             //#include "meshCourantNo.H"
                            {
                                scalarField sumPhi
                                (
                                    fvc::surfaceSum(mag(mesh.phi()))().primitiveField()
                                );

                                meshCoNum = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();

                                meanMeshCoNum =
                                    0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();
                            }

                            // Info<< "Mesh Courant Number mean: " << meanMeshCoNum
                            //     << " max: " << meshCoNum << endl;
                        }
                   }
                }

                // Assemble the Momentum equation

                MRF.correctBoundaryVelocity(U);
                tmp<fvVectorMatrix> tUEqn
                (
                  fvm::div(phi, U)
                + MRF.DDt(U)
                + problem->turbulence->divDevReff(U)
                ==
                fvOptions(U)
                );
                fvVectorMatrix& UEqn = tUEqn.ref();

                UEqn.relax();
                //std::cout << "///////////////////////////solveOnline_Pimple method line 285 //////////////////////////////////"<< std::endl;

                // #############Galerkin projection for the velocity ###########################
                List<Eigen::MatrixXd> RedLinSysU = ULmodes.project(UEqn, UprojN);
                //std::cout << "///////////////////////////solveOnline_Pimple line 288: RedLinSysU  //////////////////////////////////"<< RedLinSysU.size() << std::endl;
                //std::cout << "///////////////////////////solveOnline_Pimple line 289: b //////////////////////////////////"<< b << std::endl;
                RedLinSysU[1] = RedLinSysU[1] - projGradModP * b;
                //std::cout << "///////////////////////////solveOnline_Pimple method line 292 //////////////////////////////////"<< std::endl;

                a = reducedProblem::solveLinearSys(RedLinSysU, a, uresidual, vel_now);
                //std::cout << "///////////////////////////solveOnline_Pimple method line 295//"<< std::endl;
                ULmodes.reconstruct(U, a, "U");
        
                // fvOptions.constrain(UEqn);

                // if (pimple.momentumPredictor())
                // {
                //    solve(UEqn == -fvc::grad(P));
        
                //    fvOptions.correct(U);//??
                // }
                //solve(UEqn == - fvc::grad(P));
                //ITHACAutilities::assignBC(U, 0, v);
                //std::cout << "///////////////////////////solveOnline_Pimple method line 309"<< std::endl;

                
              // ########### Assemble matrix for the  Poisson pressure equation ############### 

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
                problem->Pmodes.reconstruct(P, b, "p");
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
            }

            std::cout << "Solution " << counter << " converged in " << iter <<
                      " iterations." << std::endl;
            std::cout << "Final normalized residual for velocity: " << U_norm_res <<
                      std::endl;
            std::cout << "Final normalized residual for pressure: " << P_norm_res <<
                      std::endl;
            // ###### Reconstruction of the solution ###############         
            ULmodes.reconstruct(U, a, "Uaux");
            P.rename("Paux");
            problem->Pmodes.reconstruct(P, b, "Paux");

            ITHACAstream::exportSolution(U, name(counter), folder);
            ITHACAstream::exportSolution(P, name(counter), folder);
            runTime.setTime(runTime.startTime(), 0);
            std::cout << "################### solveOnline_Pimple end ##################"<< std::endl;       
    }                  

}


void ReducedBasicFsi::OnlineVelocity(Eigen::MatrixXd vel)
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
