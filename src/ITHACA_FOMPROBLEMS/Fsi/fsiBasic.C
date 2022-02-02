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
    fsiBasic.C
\*---------------------------------------------------------------------------*/

#include "fsiBasic.H"

// Construct Null
fsiBasic::fsiBasic() {}
/// Construct from zero
fsiBasic::fsiBasic(int argc, char* argv[])
:unsteadyNS(argc, argv)
{
        // to create argument list
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
        //#include "createDynamicFvMesh.H"
        Info << "Create a dynamic mesh for time = "
             << runTime.timeName() << nl << endl;

        meshPtr = autoPtr<dynamicFvMesh> (dynamicFvMesh::New(args, runTime));

        dynamicFvMesh& mesh = meshPtr();
        _pimple = autoPtr<pimpleControl>
                   (
                       new pimpleControl
                       (
                           mesh
                       )
                   );
        //turbulence->validate();

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
         #include "createFields.H" 
        para = ITHACAparameters::getInstance(mesh, runTime);
        //pointVectorField & PointDisplacement = const_cast<pointVectorField&>(mesh.objectRegistry::lookupObject<pointVectorField>("pointDisplacement"));               
}

void fsiBasic::truthSolve3(List<scalar> mu_now)
{

    Time& runTime = _runTime();
    surfaceScalarField& phi = _phi();
    dynamicFvMesh& mesh = meshPtr();
    #include "initContinuityErrs.H"
    fv::options& fvOptions = _fvOptions();
    pimpleControl& pimple = _pimple();
    volScalarField& p = _p();
    volVectorField& U = _U();
    IOMRFZoneList& MRF = _MRF();
    singlePhaseTransportModel& laminarTransport = _laminarTransport();

    instantList Times = runTime.times();
    runTime.setEndTime(finalTime);

    // Perform a TruthSolve
    runTime.setTime(Times[1], 1);
    runTime.setDeltaT(timeStep);
    nextWrite = startTime; // timeStep initialization

    //****************************pimpleFoam algorithm******************************************
    #include "addCheckCaseOptions.H"
    #include "createDyMControls.H"
    #include "createUfIfPresent.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    turbulence->validate();
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    // Export and store the initial conditions for velocity and pressure
    ITHACAstream::exportSolution(U, name(counter), folder);
    ITHACAstream::exportSolution(p, name(counter), folder);
    std::ofstream of(folder + name(counter) + "/" +
    	     runTime.timeName());
    Ufield.append(U.clone());
    Pfield.append(p.clone());
    counter++;
    nextWrite += writeEvery;

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"
        runTime.setEndTime(finalTime);
        runTime++;

        //++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

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
        runTime.write();

        runTime.printExecutionTime(Info);
        
    //*****************************************end of pimpleFoam*****************************************
        if (checkWrite(runTime))
        {
            ITHACAstream::exportSolution(U, name(counter), folder);
            ITHACAstream::exportSolution(p, name(counter), folder);
            std::ofstream of(folder + name(counter) + "/" +
                             runTime.timeName());
            Ufield.append(U.clone());
            Pfield.append(p.clone());
            counter++;
            nextWrite += writeEvery;

            writeMu(mu_now);
            // --- Fill in the mu_samples with parameters (time, mu) to be used for the PODI sample points
            mu_samples.conservativeResize(mu_samples.rows() + 1, mu_now.size() + 1);
            mu_samples(mu_samples.rows() - 1, 0) = atof(runTime.timeName().c_str());

            for (label i = 0; i < mu_now.size(); i++)
            {
                mu_samples(mu_samples.rows() - 1, i + 1) = mu_now[i];
            }
        }

    }
    // Resize to Unitary if not initialized by user (i.e. non-parametric problem)
    if (mu.cols() == 0)
    {
        mu.resize(1, 1);
    }

    if (mu_samples.rows() == mu.cols())
    {
        ITHACAstream::exportMatrix(mu_samples, "mu_samples", "eigen",
                                   folder);
    }

} 


void fsiBasic::liftSolve3()
{
    //std::cout << "342 my lift solve "<< std::endl;
    for (label k = 0; k < inletIndex.rows(); k++)
    {
        //std::cout << "345 my lift solve "<< std::endl;
        Time& runTime = _runTime();
        std::cout << "???????????? 234 my lift solve ???????????????????? "<< std::endl;
        Foam::dynamicFvMesh& mesh = meshPtr();
        std::cout << "***********236 my lift solve********* "<< std::endl;
        volScalarField& p = _p();
        volVectorField& U = _U();
        std::cout << "239 my lift solve "<< std::endl;
        //Foam::dynamicFvMesh& mesh = meshPtr();
        std::cout << "351 my lift solve "<< std::endl;
        surfaceScalarField& phi = _phi();
        //#include "initContinuityErrs.H"
        fv::options& fvOptions = _fvOptions();
        pimpleControl& pimple = _pimple();
        IOMRFZoneList& MRF = _MRF();
        std::cout << "354 my lift solve "<< std::endl;

        label BCind = inletIndex(k, 0);
        volVectorField Ulift("Ulift" + name(k), U);
        instantList Times = runTime.times();
        runTime.setTime(Times[1], 1);
        pisoControl piso(mesh);
        Info << "Solving a lifting Problem" << endl;
        Vector<double> v1(0, 0, 0);
        v1[inletIndex(k, 1)] = 1;
        Vector<double> v0(0, 0, 0);
        std::cout << "365 my lift solve "<< std::endl;

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
         std::cout << "277 my lift solve "<< std::endl;

        Info << "Constructing velocity pimple field Phi\n" << endl;
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
      
        label pRefCell = 0;
        scalar pRefValue = 0.0;
        setRefCell
        (
            Phi, 
            mesh.solutionDict().subDict("PIMPLE"), 
            pRefCell, 
            pRefValue
        );
        mesh.setFluxRequired(Phi.name());
        runTime.functionObjects().start();
        MRF.makeRelative(phi);
        adjustPhi(phi, Ulift, p);
        #include "UEqn.H"
        #include "createUfIfPresent.H"
        volScalarField rAU(1.0/UEqn.A());
        volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
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

        if (p.needReference())
        {
            fvc::makeRelative(phiHbyA, U);
            adjustPhi(phiHbyA, U, p);
            fvc::makeAbsolute(phiHbyA, U);
        }

        tmp<volScalarField> rAtU(rAU);

        if (pimple.consistent())
        {
            rAtU = 1.0/max(1.0/rAU - UEqn.H1(), 0.1/rAU);
            phiHbyA +=
                fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh.magSf();
            HbyA -= (rAU - rAtU())*fvc::grad(p);
        }

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
                fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA)
            );
            // Added for reduced problem
            //RedLinSysP = problem->Pmodes.project(pEqn, NmodesPproj);
            //pEqn.solve();
            //b = reducedProblem::solveLinearSys(RedLinSysP, b, presidual);
            //problem->Pmodes.reconstruct(P, b, "p");
            pEqn.setReference(pRefCell, pRefValue);

            pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

            if (pimple.finalNonOrthogonalIter())
            {
                phi = phiHbyA - pEqn.flux();
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
        liftfield.append(Ulift.clone());
    }
}
