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
    Example of Fsi simulation for PIMPLE algorithm
SourceFiles
    fsiBasic.C
\*---------------------------------------------------------------------------*/

#include "fsiBasic.H"

// Construct Null
fsiBasic::fsiBasic() {}
/// Construct from zero
fsiBasic::fsiBasic(int argc, char* argv[])
:UnsteadyProblem()
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
        //pimpleControl& pimple = _pimple();
#include "createFields.H" 
#include "initContinuityErrs.H"
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
        ITHACAparameters* para = ITHACAparameters::getInstance(mesh,_runTime());
        //para = ITHACAparameters::getInstance(mesh, runTime); 
        offline = ITHACAutilities::check_off();
        podex = ITHACAutilities::check_pod();
        setTimes(runTime);
        point0 = mesh.points();  
        Info << offline << endl;
    /// Number of velocity modes to be calculated
        NUmodesOut = para->ITHACAdict->lookupOrDefault<label>("NmodesUout", 15);
        /// Number of pressure modes to be calculated
        NPmodesOut = para->ITHACAdict->lookupOrDefault<label>("NmodesPout", 15);
        /// Number of nut modes to be calculated
        NNutModesOut = para->ITHACAdict->lookupOrDefault<label>("NmodesNutOut", 15);
        /// Number of velocity modes used for the projection
        NUmodes = para->ITHACAdict->lookupOrDefault<label>("NmodesUproj", 10);
        /// Number of pressure modes used for the projection
        NPmodes = para->ITHACAdict->lookupOrDefault<label>("NmodesPproj", 10);
        /// Number of nut modes used for the projection
        NNutModes = para->ITHACAdict->lookupOrDefault<label>("NmodesNutProj", 0);
              
}

void fsiBasic::truthSolve(List<scalar> mu_now, fileName folder)
{

    Time& runTime = _runTime();
    surfaceScalarField& phi = _phi();
    dynamicFvMesh& mesh = meshPtr();
    fv::options& fvOptions = _fvOptions();
    pimpleControl& pimple = _pimple();
    volScalarField& p = _p();
    volVectorField& U = _U();
    IOMRFZoneList& MRF = _MRF();
    //surfaceVectorField& Uf = _Uf();
    singlePhaseTransportModel& laminarTransport = _laminarTransport();
    instantList Times = runTime.times();
    runTime.setEndTime(finalTime);

    // Perform a TruthSolve
    runTime.setTime(Times[1], 1);
    runTime.setDeltaT(timeStep);
    nextWrite = startTime; // timeStep initialization

    dictionary dictCoeffs(dyndict->findDict("sixDoFRigidBodyMotionCoeffs"));
    Foam::functionObjects::forces fomforces("fomforces", mesh, dictCoeffs);
   

    turbulence->validate();
#include "createUfIfPresent.H"

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {

#include "CourantNo.H"
        
        runTime.setEndTime(finalTime);
        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                fomforces.calcForcesMoment();

                // Do any mesh changes
                //mesh.controlledUpdate();
                // The following line remplace the above controlledUpdate() method
                sDRBMS().solve();
                mesh.movePoints(sDRBMS().curPoints());
                // std::cerr << "################"<< "Before six dof motion solver" << "#############"<< std::endl;
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
        //std::cout << "/////////////" << runTime.deltaTValue() << "////////////" << std::endl;
        scalar alffa = sDRBMS().motion().omega().z() / runTime.deltaTValue();

        if (checkWrite(runTime))
        {

            fomforcex.append(fomforces.forceEff().x());
            fomforcey.append(fomforces.forceEff().y());
            centerofmassx.append(sDRBMS().motion().centreOfMass().x());
            centerofmassy.append(sDRBMS().motion().centreOfMass().y());
            centerofmassz.append(quaternion(sDRBMS().motion().orientation()).eulerAngles(quaternion::XYZ).z());
      
            ITHACAstream::exportSolution(U, name(counter), folder);
            ITHACAstream::exportSolution(p, name(counter), folder);
            ITHACAstream::exportSolution(sDRBMS().pointDisplacement(), name(counter), folder);
            ITHACAstream::writePoints(meshPtr().points(), folder, name(counter) + "/polyMesh/");

            std::ofstream of(folder + name(counter) + "/" + runTime.timeName());
            Ufield.append(U.clone());
            Pfield.append(p.clone());
            Dfield.append(sDRBMS().pointDisplacement().clone());
            
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
        ITHACAstream::exportMatrix(mu_samples, "mu_samples", "eigen",folder);
    }

} 




void fsiBasic::restart()
{

    turbulence.clear();
    _fvOptions.clear();
    _laminarTransport.clear();
    _p.clear();
    _U.clear();
    _phi.clear();
    _Uf.clear();
    _pointDisplacement.clear();
    sDRBMS.clear();
    argList& args = _args();
    Time& runTime = _runTime();
    runTime.setTime(0, 1);
    // meshPtr().resetMotion();
    meshPtr().movePoints(point0);    
    pointField& pointOld = const_cast<pointField&> (meshPtr().oldPoints());
    pointOld = point0;
    _pimple.clear();
    Foam::dynamicFvMesh& mesh = meshPtr();
    _pimple = autoPtr<pimpleControl>
                   (
                       new pimpleControl
                       (
                           mesh
                       )
               );
    
#include "createFields.H" 
}

