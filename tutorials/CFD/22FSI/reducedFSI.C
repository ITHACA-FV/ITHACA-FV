class reducedFSI : public reducedSimpleSteadyNS
{
    public:
        /// Constructor
        reducedFSI();
        explicit reducedFSI(FSI& FOMproblem)
            :
            problem(&FOMproblem)
        {}

        /// Full problem.
        FSI* problem;

        // Function to perform the online phase
        void solveOnline_Pimple(scalar mu_now,
                                int NmodesUproj, int NmodesPproj, int NmodesNut = 0,
                                word Folder = "./ITHACAoutput/Reconstruct/")
        {
            counter++;

            // For all the variables, in case one wants to use all the available modes it is just necessary to set
            // the requested number into the ITHACAdict to zero
            if (NmodesUproj == 0)
            {
                NmodesUproj = problem->Umodes.size();
            }

            if (NmodesPproj == 0)
            {
                NmodesPproj = problem->Pmodes.size();
            }

            if (NmodesNut == 0)
            {
                NmodesNut = problem->nutModes.size();
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
            //Eigen::MatrixXd nutCoeff = Eigen::VectorXd::Zero(NmodesNut);
            //Eigen::MatrixXd nutCoeffOld = Eigen::VectorXd::Zero(NmodesNut);
            float residualJumpLim =
                problem->para->ITHACAdict->lookupOrDefault<float>("residualJumpLim", 1e-5);
            float normalizedResidualLim =
                problem->para->ITHACAdict->lookupOrDefault<float>("normalizedResidualLim",
                        1e-5);
            scalar residual_jump(1 + residualJumpLim);
            volScalarField& P = problem->_p();
            volVectorField& U = problem->_U();
            //volScalarField& nut = const_cast<volScalarField&>
            //                     (problem->_mesh().lookupObject<volScalarField>("nut"));
            volVectorField u2 = U;

            // Getting coefficients
            a0 = ITHACAutilities::getCoeffs(U, problem->Umodes, NmodesUproj, true);
            b = ITHACAutilities::getCoeffs(P, problem->Pmodes, NmodesPproj, true);
            //nutCoeff = ITHACAutilities::getCoeffs(nut, problem->nutModes, NmodesNut, true);


            a(0) = a0(0); // Do not remove: it is not working without this condition
            b = b0;
            fvMesh& mesh = problem->_mesh();
            Time& runTime = problem->_runTime();
            P.rename("p");
            surfaceScalarField& phi(problem->_phi());
            problem->Umodes.reconstruct(U, a, "U");
            problem->Pmodes.reconstruct(P, b, "p");
            vector v(1, 0, 0);
            ITHACAutilities::assignBC(U, 0, v);
            //problem->nutModes.reconstruct(nut, nutCoeff, "nut");
            phi = fvc::flux(U);
            int iter = 0;
            pimpleControl& pimple = _pimple();
            //IOMRFZoneList& MRF = _MRF();
            //singlePhaseTransportModel& laminarTransport = _laminarTransport(); 
            //simpleControl& simple = problem->_simple();
            std::ofstream res_os_U;
            std::ofstream res_os_P;
            res_os_U.open(Folder + name(counter) + "/residualsU", std::ios_base::app);
            res_os_P.open(Folder + name(counter) + "/residualsP", std::ios_base::app);

            // Pimple algorithm starts here
            while ((residual_jump > residualJumpLim
                    || std::max(U_norm_res, P_norm_res) > normalizedResidualLim) &&
                    iter < maxIterOn)
            {
                iter++;
                std::cout << iter << std::endl;

#if defined(OFVER) && (OFVER == 6)
                pimple.loop(runTime);
#else
                problem->_pimple().loop();
#endif

                // If the case is turbulent, then the network is evaluated
                /*if (ITHACAutilities::isTurbulent())
                {
                    nutCoeff = problem->evalNet(a, mu_now);
                    //nutCoeff = nutCoeffOld + 0.7*(nutCoeff - nutCoeffOld);
                    volScalarField& nut = const_cast<volScalarField&>
                                          (problem->_mesh().lookupObject<volScalarField>("nut"));
                    problem->nutModes.reconstruct(nut, nutCoeff, "nut");
                    //ITHACAstream::exportSolution(nut, name(counter), Folder);
                }*/

                //volScalarField nueff = problem->turbulence->nuEff();
                vector v(1, 0, 0);
                ITHACAutilities::assignBC(U, 0, v);
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
                + turbulence->divDevReff(U)
                ==
                fvOptions(U)
                );
                fvVectorMatrix& UEqn = tUEqn.ref();

                UEqn.relax();
                /*List<Eigen::MatrixXd> RedLinSysU = ULmodes.project(UEqn, UprojN);
                RedLinSysU[1] = RedLinSysU[1] - projGradModP * b;
                a = reducedProblem::solveLinearSys(RedLinSysU, a, uresidual, vel_now);
                ULmodes.reconstruct(U, a, "U");
        */
                fvOptions.constrain(UEqn);

                if (pimple.momentumPredictor())
                {
                   //solve(UEqn == -fvc::grad(p));
                    List<Eigen::MatrixXd> RedLinSysU = ULmodes.project(UEqn, UprojN);
                    RedLinSysU[1] = RedLinSysU[1] - projGradModP * b;
                    a = reducedProblem::solveLinearSys(RedLinSysU, a, uresidual, vel_now);
                    ULmodes.reconstruct(U, a, "U");
        
                   //fvOptions.correct(U);??
                }
                //solve(UEqn == - fvc::grad(P));
                ITHACAutilities::assignBC(U, 0, v);
                
              // begining of solving the pEqn.H   

                volScalarField rAU(1.0/UEqn.A());
                volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, P));
                surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));

                if (pimple.ddtCorr())
                {
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
                        fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh.magSf();
                    HbyA -= (rAU - rAtU())*fvc::grad(P);
                }
                List<Eigen::MatrixXd> RedLinSysP;
                bOld = b;

                if (pimple.nCorrPISO() <= 1)
                {
                    tUEqn.clear();
                }

                // Update the pressure BCs to ensure flux consistency
                constrainPressure(p, U, phiHbyA, rAtU(), MRF);

                // Non-orthogonal pressure corrector loop
                while (pimple.correctNonOrthogonal())
                {
                    // Continuity equation
                    fvScalarMatrix pEqn
                    (
                        fvm::laplacian(rAtU(), P) == fvc::div(phiHbyA)
                    );
                    // Added for reduced problem
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

                #include "continuityErrs.H"

                // Explicitly relax pressure for momentum corrector
                P.relax();

                U = HbyA - rAtU*fvc::grad(p);
                U.correctBoundaryConditions();
                fvOptions.correct(U);

                // Correct Uf if the mesh is moving
                fvc::correctUf(Uf, U, phi);

                // Make the fluxes relative to the mesh motion
                fvc::makeRelative(phi, U);

                b = bOld + mesh.fieldRelaxationFactor("p") * (b - bOld);
                problem->Pmodes.reconstruct(P, b, "p");
                //nutCoeffOld = nutCoeff;
                // P.relax();
                U = HbyA - rAtU() * fvc::grad(P);
                U.correctBoundaryConditions();
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

                res_os_U << U_norm_res << std::endl;
                res_os_P << P_norm_res << std::endl;
            }

            res_os_U.close();
            res_os_P.close();
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

            ITHACAstream::exportSolution(U, name(counter), Folder);
            ITHACAstream::exportSolution(P, name(counter), Folder);
            runTime.setTime(runTime.startTime(), 0);
        }
};