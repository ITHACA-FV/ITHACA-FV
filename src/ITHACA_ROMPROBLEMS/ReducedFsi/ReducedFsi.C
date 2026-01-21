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
/// Source file of the ReducedFsi class
#include "ReducedFsi.H"

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Null Constructor
ReducedFsi::ReducedFsi() {}


ReducedFsi::ReducedFsi(fsiBasic& FoamPb): problem(&FoamPb)

{
    for (int i = 0; i < problem->Umodes.size(); i++)
    {
        Umodes.append((problem->Umodes.toPtrList()[i]).clone());
    }

    for (int i = 0; i < problem->Pmodes.size(); i++)
    {
        Pmodes.append((problem->Pmodes.toPtrList()[i]).clone());
    }

    for (int i = 0; i < problem->Dmodes.size(); i++)
    {
        Dmodes.append((problem->Dmodes.toPtrList()[i]).clone());
    }

    //problem->restart();
    //std::cout << "################ ctor of POD-I Fsi ##################" << std::endl;
}
void ReducedFsi::PODI(Eigen::MatrixXd coeffL2, Eigen::MatrixXd muu,
                      label NPdModes)
{
    if (NPdModes == 0)
    {
        NPdModes = Dmodes.size();
    }

    problem->samples.resize(NPdModes);
    problem->rbfSplines.resize(NPdModes);
    Eigen::MatrixXd weights;

    for (label i = 0; i < NPdModes; i++) // i is the nnumber of th mode
    {
        word weightName = "wRBF_M" + name(i + 1);

        if (ITHACAutilities::check_file("./ITHACAoutput/weights/" + weightName))
        {
            problem->samples[i] = new SPLINTER::DataTable(true, true);

            for (label j = 0; j < coeffL2.cols();
                    j++) // j is the number of the nut snapshot
            {
                problem->samples[i]->addSample(muu.row(j), coeffL2(i, j));
            }

            ITHACAstream::ReadDenseMatrix(weights, "./ITHACAoutput/weights/", weightName);
            problem->rbfSplines[i] = new SPLINTER::RBFSpline(*(problem->samples)[i],
                SPLINTER::RadialBasisFunctionType::THIN_PLATE_SPLINE, weights);
            std::cout << "Constructing RadialBasisFunction for mode " << i + 1 << std::endl;
        }
        else
        {
            problem->samples[i] = new SPLINTER::DataTable(true, true);

            for (label j = 0; j < coeffL2.cols();
                    j++) // j is the number of the nut snapshot
            {
                problem->samples[i]->addSample(muu.row(j), coeffL2(i, j));
            }

            problem->rbfSplines[i] = new SPLINTER::RBFSpline(*(problem->samples)[i],
                SPLINTER::RadialBasisFunctionType::THIN_PLATE_SPLINE);
            ITHACAstream::SaveDenseMatrix(problem->rbfSplines[i]->weights,
                                          "./ITHACAoutput/weights/", weightName);
            std::cout << "Constructing RadialBasisFunction for mode " << i + 1 << std::endl;
        }
    }
}



void ReducedFsi::solveOnline_Pimple(int NmodesUproj,
                                    int NmodesPproj,
                                    int NmodesDproj,
                                    fileName folder)
{
    Eigen::MatrixXd a = Eigen::VectorXd::Zero(NmodesUproj);
    Eigen::MatrixXd b = Eigen::VectorXd::Zero(NmodesPproj);
    Time& runTime = problem->_runTime();
    dynamicFvMesh& mesh = problem->meshPtr();
    //const pointMesh& pMesh = pointMesh::New(mesh);
    fv::options& fvOptions = problem->_fvOptions();
    pimpleControl& pimple = problem->_pimple();
    volScalarField& p = problem->_p();
    volVectorField& U = problem->_U();
    surfaceScalarField& phi = problem->_phi();
    pointVectorField& pointDisplacement = problem->_pointDisplacement();//????
    IOMRFZoneList& MRF = problem->_MRF();
    singlePhaseTransportModel& laminarTransport = problem->_laminarTransport();
    //autoPtr<incompressible::turbulenceModel> turbulence = problem->turbulence;
    autoPtr<incompressible::turbulenceModel> turbulence(
        incompressible::turbulenceModel::New(U,
                phi, laminarTransport));
    instantList Times = runTime.times();
    runTime.setEndTime(finalTime);
    PODI(problem->coeffL2,  problem->CylDispl,  NmodesDproj);
    runTime.setTime(Times[1], 1);
    runTime.setDeltaT(timeStep);
    nextWrite = startTime;
    label pRefCell = 0;
    scalar pRefValue = 0.0;
    bool    correctPhi = problem->correctPhi;
    bool    checkMeshCourantNo = problem->checkMeshCourantNo;
    bool    moveMeshOuterCorrectors = problem->moveMeshOuterCorrectors;
    scalar  cumulativeContErr = problem->cumulativeContErr;
#include "createUfIfPresent.H"
    turbulence->validate();
    dictionary dictCoeffs(
        problem->dyndict->findDict("sixDoFRigidBodyMotionCoeffs"));
    Foam::functionObjects::forces romforces("romforces", mesh, dictCoeffs);
    sixDoFRigidBodyMotion sDRBM(dictCoeffs, dictCoeffs, runTime );
    Foam::dimensionedVector g("g", dimAcceleration, Zero);
    dictCoeffs.readIfPresent("g", g);
    // Eigen::VectorXd pdCoeff;
    // pdCoeff.resize(NmodesDproj);
    bool firstIter = false;
    Eigen::MatrixXd pdCoeff(NmodesDproj, 1);
    //- Current time index (used for updating)
    label curTimeIndex_ = -1;
    // pointField points0 = mesh.points();
    Eigen::MatrixXd muEval(1, 1);
    //Eigen::VectorXd b_ref = Pmodes.project(p, NmodesPproj);
    //Eigen::VectorXd a_ref = Umodes.project(U, NmodesUproj);
    //std::cout << " ==== b_ref === " << b_ref << std::endl;
    //std::cout << " ==== a_ref === " << a_ref << std::endl;
    //exit(0);
    // double lambda0 = 1.0; // Initial weight
    // double lambda_t = 0.0;
    // double alpha = 0.125;

    // PIMPLE algorithm starts here
    //Info<< "\nStarting time loop\n" << endl;
    //std::ofstream res_p, res_u;
    //res_u.open("./res_u", std::ios_base::app);
    //res_p.open("./res_p", std::ios_base::app);
    //Errors << "Time, res_u, res_p" << endl;
    while (runTime.run())
    {
        runTime.setEndTime(finalTime);
        runTime++;
        //p.storePrevIter();
        Info << "Time = " << runTime.timeName() << nl << endl;
        //res_u << runTime.timeName() << std::endl;
        //res_p << runTime.timeName() << std::endl;

        while (pimple.loop())
        {
#include "CourantNo.H"

            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
#include"CylinderMotion.H"

                if (mesh.changing())
                {
                    //MRF.update();
                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();
                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);
                    }

                    if (checkMeshCourantNo)
                    {
#include "meshCourantNo.H"
                    }
                }
            }

            // Solve the Momentum equation
            //MRF.correctBoundaryVelocity(U);
            fvVectorMatrix UEqn
            (
                fvm::ddt(U)
                + fvm::div(phi, U)
                //+ MRF.DDt(U)
                + turbulence->divDevReff(U)
                == fvOptions(U)
            );
            //fvVectorMatrix& UEqn = tUEqn.ref();
            UEqn.relax();
            fvOptions.constrain(UEqn);
            List<Eigen::MatrixXd> RedLinSysU;

            if (pimple.momentumPredictor())
            {
                //solve(UEqn == -fvc::grad(p));
                RedLinSysU = Umodes.project(UEqn, NmodesUproj, "G");
                volVectorField gradpfull = -fvc::grad(p);
                Eigen::MatrixXd projGrad = Umodes.project(gradpfull, NmodesUproj);
                RedLinSysU[1] = RedLinSysU[1] + projGrad;
                //a = RedLinSysU[0].householderQr().solve(RedLinSysU[1]);
                //Eigen::MatrixXd I_u = 1e-6*Eigen::MatrixXd::Identity(NmodesUproj,NmodesUproj );
                //a = (RedLinSysU[0] - I_u).ldlt().solve(RedLinSysU[1]);
                //a=RedLinSysU[0].llt().solve(RedLinSysU[1]);
                a = RedLinSysU[0].colPivHouseholderQr().solve(RedLinSysU[1]);
                //a = RedLinSysU[0].completeOrthogonalDecomposition().solve(RedLinSysU[1]);
                //a=RedLinSysU[0].fullPivLu().solve(RedLinSysU[1]);
                //Eigen::ConjugateGradient<Eigen::MatrixXd> cg;
                //cg.setTolerance(1e-6);       // Tighter than default (1e-2)
                //cg.setMaxIterations(300);    // Prevent excessive iterations
                //cg.compute(RedLinSysU[0]);
                //Eigen::JacobiSVD<Eigen::MatrixXd> svd(RedLinSysU[0]);
                /*Eigen::JacobiSVD<Eigen::MatrixXd> svd(
                RedLinSysU[0],
                Eigen::ComputeThinU | Eigen::ComputeThinV);
                svd.setThreshold(1e-6);  // Truncate small singular values
                if (svd.rank() < NmodesUproj)
                {
                    Warning << "Rank-deficient matrix! Rank = " << svd.rank() << endl;
                }
                P    a = svd.solve(RedLinSysU[1]);*/
                //Eigen::MatrixXd A_reg = RedLinSysU[0] + 1e-6 *Eigen::MatrixXd::Identity(NmodesUproj, NmodesUproj);
                //a = A_reg.ldlt().solve(RedLinSysU[1]);
                //Eigen::MatrixXd aNew = cg.solve(RedLinSysU[1]);
                //std::cout << "res_u = " << (RedLinSysU[0] * a - RedLinSysU[1]).norm() << std::endl;
                //res_u << (RedLinSysU[0] * a - RedLinSysU[1]).norm() << std::endl;
                /*
                if (cg.info() != Eigen::Success)
                {
                Warning << "ROM velocity solve failed. Using previous coefficients." << endl;
                }*/
                //lambda_t = lambda0 * std::exp(-alpha * problem->_runTime().value());
                //Eigen::VectorXd a_cal = (a + lambda_t * a_ref) / (1.0 + lambda_t);
                //a_ref = a;
                Umodes.reconstruct(U, a, "U");
                //U.correctBoundaryConditions();
                fvOptions.correct(U);
                //volVectorField errorU = U -problem->Ufield[0];
                //Info << "errorU= " << errorU.boundaryFieldRef() << endl;
                //exit(0);
            }

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                volScalarField rAU(1.0 / UEqn.A());
                volVectorField HbyA(constrainHbyA(rAU * UEqn.H(), U, p)); //p
                surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
                tmp<volScalarField> rAtU(rAU);
                // Update the pressure BCs to ensure flux consistency
                constrainPressure(p, U, phiHbyA, rAtU(), MRF); //p
                List<Eigen::MatrixXd> RedLinSysP;

                // Non-orthogonal pressure corrector loop
                while (pimple.correctNonOrthogonal())
                {
                    fvScalarMatrix pEqn
                    (
                        fvm::laplacian(rAtU(), p)
                        ==
                        fvc::div(phiHbyA)
                    );
                    pEqn.setReference(pRefCell, pRefValue);
                    //pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter()))); //p
                    RedLinSysP = Pmodes.project(pEqn, NmodesPproj, "G");
                    /// Solve for the reduced coefficient for pressure
                    //b = RedLinSysP[0].householderQr().solve(RedLinSysP[1]);
                    b = RedLinSysP[0].colPivHouseholderQr().solve(RedLinSysP[1]);
                    //b = RedLinSysP[0].ldlt().solve(RedLinSysP[1]);
                    //Eigen::VectorXd b_cal = (b + lambda_t * b_ref) / (1.0 + lambda_t);
                    //b_ref = b;
                    //b =RedLinSysP[0].fullPivLu().solve(RedLinSysP[1]);
                    //Eigen::MatrixXd I_p = 1e-6*Eigen::MatrixXd::Identity(NmodesPproj,NmodesPproj );
                    //b = (RedLinSysP[0] - I_p).ldlt().solve(RedLinSysP[1]);
                    //b = RedLinSysP[0].ldlt().solve(RedLinSysP[1]);
                    // Solve pressure ROM system (with iterative solver)
                    /*
                    Eigen::ConjugateGradient<Eigen::MatrixXd> cgP;
                    cgP.setTolerance(1e-6);
                    cgP.compute(RedLinSysP[0]);
                    Eigen::MatrixXd bNew = cgP.solve(RedLinSysP[1]);*/
                    //std::cout << "res_p = " << (RedLinSysP[0] * b - RedLinSysP[1]).norm() << std::endl;
                    //res_p << (RedLinSysP[0] * b - RedLinSysP[1]).norm() << std::endl;
                    Pmodes.reconstruct(p, b, "p");

                    //p.correctBoundaryConditions();
                    //volScalarField errorp = p -problem->Pfield[0] ;
                    //Info << "errorp= " << errorp.internalField() << endl;
                    //exit(0);
                    if (pimple.finalNonOrthogonalIter())
                    {
                        phi = phiHbyA - pEqn.flux();
                    }
                }

                // Explicitly relax pressure for momentum corrector
                p.relax();
                /// Correct the velocity
                U = HbyA - rAtU * fvc::grad(p); //p
                U.correctBoundaryConditions();
                fvOptions.correct(U);
                //volVectorField error = U -problem->Ufield[0] ;
                //Info << "Initial error: " << error.internalField() << endl;
                //exit(0);
                // Correct Uf if the mesh is moving
                fvc::correctUf(Uf, U, phi);
                // Make the fluxes relative to the mesh motion
                fvc::makeRelative(phi, U);
            }// end of the pimple.correct()
        }// end of the pimple.loop()

        if (checkWrite(runTime))
        {
            centerofmassy.append(sDRBM.centreOfMass().y());
            pdcoeffrbf.append(pdCoeff(0, 0));
            CoeffP.append(b);
            CoeffU.append(a);
            romforcey.append(romforces.forceEff().y());
            romforcex.append(romforces.forceEff().x());

            if (runTime.time().value() + runTime.deltaT().value() >= finalTime)
            {
                Info << "===== Storing final mesh =====" << endl;
                OnlineMeshes.append(problem->meshPtr.ptr());  // Transfer ownership
                // NOTE: meshPtr is now empty! Handle accordingly.
            }

            ListOfpoints.append(mesh.points().clone());
            //std::ofstream of(folder + name(counter) + "/" + runTime.timeName());
            UredFields.append(U.clone());
            PredFields.append(p.clone());
            Dfield.append(pointDisplacement.clone());
            counter++;
            nextWrite += writeEvery;
        }
    } // end of the runTime.run() loop

    //res_u.close();
    //res_p.close();
} // end of the method Solve


bool ReducedFsi::checkWrite(Time& timeObject)
{
    scalar diffnow = mag(nextWrite - atof(timeObject.timeName().c_str()));
    scalar diffnext = mag(nextWrite - atof(timeObject.timeName().c_str()) -
                          timeObject.deltaTValue());

    if ( diffnow < diffnext)
    {
        return true;
    }
    else
    {
        return false;
    }
}

void  ReducedFsi::prepareRomData(const word& outputPath)
{
    word fullPath = "./" + outputPath;

    if (!ITHACAutilities::check_folder(fullPath))
    {
        mkDir(fullPath);
        problem->exportFoamFieldToNpy(fullPath, "romforcex", romforcex);
        problem->exportFoamFieldToNpy(fullPath, "romforcey", romforcey);
        problem->exportFoamFieldToNpy(fullPath, "CentreOfMassY", centerofmassy);
    }
}
