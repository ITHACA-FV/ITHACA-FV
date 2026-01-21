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
    : UnsteadyProblem()
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
    // oMesh.reset
    // (
    //     new IOobject
    //     (
    //         "OriginalMesh",
    //         runTime.timeName(),
    //         runTime,
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     )
    // );
    /// Create a deep copy
    ITHACAparameters* para = ITHACAparameters::getInstance(mesh, _runTime());
    //para = ITHACAparameters::getInstance(mesh, runTime);
    offline = ITHACAutilities::check_off();
    podex = ITHACAutilities::check_pod();
    //setTimes(runTime);
    // point0 = mesh.points().clone();
    // faces0 =  mesh.faces().clone();
    // celllist0 = mesh.cells().clone();
    // ///
    // //const polyBoundaryMesh& boundary = mesh.boundaryMesh();
    // /// Construct the initial mesh
    // Mesh0.reset
    // (new fvMesh
    //     (
    //       oMesh,
    //       std::move(point0),
    //       std::move(faces0),
    //       std::move(celllist0),
    //       true
    //      )
    // );
    // /*
    // Mesh0
    // (
    //     oMesh,
    //     std::move(point0),
    //     std::move(faces0),
    //     std::move(celllist0),
    //     true  // syncPar
    // );*/
    // PtrList<polyPatch> patches(meshPtr->boundaryMesh().size());
    // forAll(meshPtr->boundaryMesh(), patchI)
    // {
    //     patches.set
    //     (
    //         patchI,
    //         meshPtr->boundaryMesh()[patchI].clone
    //         (
    //             Mesh0->boundaryMesh(),
    //             patchI,
    //             meshPtr->boundaryMesh()[patchI].size(),
    //             meshPtr->boundaryMesh()[patchI].start()
    //         )
    //     );
    // }
    // Mesh0->addFvPatches(patches);
    // // Now meshCopyPtr is fully independent
    // Mesh0->write();
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

void fsiBasic::truthSolve(label folderN, fileName folder)
{
    Time& runTime = _runTime();
    dynamicFvMesh& mesh = meshPtr();
    //fvMesh& mesh0 =  Mesh0();
    // Create a new independent copy (if dynamicFvMesh supports copying)
    //autoPtr<dynamicFvMesh> mesh2Ptr(meshPtr->clone());  // Requires clone() method
    //dynamicFvMesh& mesh = mesh2Ptr();  // Now 'mesh' is independent of 'meshPtr'
    fv::options& fvOptions = _fvOptions();
    pimpleControl& pimple = _pimple();
    volScalarField& p = _p();
    volVectorField& U = _U();
    surfaceScalarField& phi = _phi();
    IOMRFZoneList& MRF = _MRF();
    //surfaceVectorField& Uf = _Uf();
    singlePhaseTransportModel& laminarTransport = _laminarTransport();
    instantList Times = runTime.times();
    runTime.setEndTime(finalTime);
    // Perform a TruthSolve
    runTime.setTime(Times[1], 1);
    runTime.setDeltaT(timeStep);
    nextWrite = startTime; // timeStep initialization
    //const fvMesh& toMeshInit = meshPtr();
    //meshToMesh0::order mapOrder = meshToMesh0::INTERPOLATE;
    //meshToMesh0::order mapOrder = meshToMesh0::CELL_VOLUME_WEIGHT;
    //meshToMesh0::order mapOrder = meshToMesh0::CELL_POINT_INTERPOLATE;
    meshToMesh0::order mapOrder = meshToMesh0::MAP;
    dictionary dictCoeffs(dyndict().findDict("sixDoFRigidBodyMotionCoeffs"));
    Foam::functionObjects::forces fomforces("fomforces", mesh, dictCoeffs);
    surfaceVectorField N = mesh.Sf() / mesh.magSf();
    turbulence->validate();
#include "createUfIfPresent.H"
    Info << "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
#include "CourantNo.H"
        runTime.setEndTime(finalTime);
        runTime++;
        Info << "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                //#if (DOFVER == 2412 || DOFVER == 2506) //DOPENFOAM
                fomforces.execute();
                // #else
                //fomforces.calcForcesMoment();
                // #endif
                // #if defined(OFVERSION) && (OFVERSION > 2106)Do any mesh changes
                //mesh.controlledUpdate();
                // The following line remplace the above controlledUpdate() method
                sDRBMS().solve();
                mesh.movePoints(sDRBMS().curPoints());

                //meshToMesh0 mapper(mesh, meshPtr());
                //Info << mapper.toMesh().points()[1000] << endl;
                //Info << mapper.fromMesh().points()[1000] << endl;
                if (mesh.changing())
                {
                    //MRF.update();
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

        scalar alffa = sDRBMS().motion().omega().z() / runTime.deltaTValue();
        // Face area normal vectors
        //surfaceVectorField Sf = mesh.Sf();
        // Face areas
        //surfaceScalarField magSf = mesh.magSf();
        // Face normals
        //N = mesh.Sf()/mesh.magSf();
        //N.rename("Uf");

        // Access the result
        //volVectorField& old_U = t_old_U.ref();

        //Info << "old_U =" << old_U <<

        if (checkWrite(runTime))
        {
            //folderN++;
            fomforcex.append(fomforces.forceEff().x());
            fomforcey.append(fomforces.forceEff().y());
            centerofmassx.append(sDRBMS().motion().centreOfMass().x());
            centerofmassy.append(sDRBMS().motion().centreOfMass().y());
            centerofmassz.append(quaternion(sDRBMS().motion().orientation()).eulerAngles(
                                     quaternion::XYZ).z());
            /*
             Foam::meshToMesh0 mapper(mesh,mesh0);
            /// Interpolate from U on current mesh to U_mapped on original mesh
             //tmp<volVectorField> U_mapped(U);
              //tmp<volScalarField> p_mapped(p);
              volVectorField U_mapped
             (
                 IOobject("U", runTime.timeName(), mesh0),
                 mesh0,
                 dimensionedVector("U", U.dimensions(), Zero)
             );

             volScalarField p_mapped
             (
                 IOobject("p", runTime.timeName(), mesh0),
                 mesh0,
                 dimensionedScalar("p", p.dimensions(), Zero)
             );
            if (&mesh == &mesh0)
             {
                 Info << "[mapFieldToMesh0] Meshes are identical. Copying field..." << endl;
                 U_mapped = U;
                 p_mapped = p;
             }
             else
             {
                U_mapped = mapper.interpolate<Foam::vector,                    Foam::plusEqOp<Foam::vector>>(U, mapOrder, Foam::plusEqOp<Foam::vector>() );

            p_mapped = mapper.interpolate<scalar,Foam::plusEqOp<scalar>>(p, mapOrder, Foam::plusEqOp<scalar>() );

             }
             U_mapped.ref().rename("U");
             p_mapped.ref().rename("p");
             U_mapped.correctBoundaryConditions();
             p_mapped.correctBoundaryConditions();*/
            //ITHACAstream::exportSolution(N, name(counter), folder);
            word localFolder = folder +  name(folderN + 1);
            //old_U.rename("old_U");
            ITHACAstream::exportSolution(U, name(counter), localFolder );
            //ITHACAstream::exportSolution(N, name(counter), localFolder );
            ITHACAstream::exportSolution(p, name(counter),  localFolder );
            ITHACAstream::exportSolution(sDRBMS().pointDisplacement(),
                                         name(counter), localFolder );
            ITHACAstream::writePoints(mesh.points(),
                                      localFolder,  name(counter) + "/polyMesh/");
            //word saveDir = name(counter) + "/polyMesh/";
            //OFstream os(localFolder/saveDir/"points");
            // Write points
            //OFstream os(localFolder/saveDir);
            //List<vector> list(mesh.points());
            //list.write(os);
            //os << mesh.points();
            // Copy files to save directory
            //mesh.write();
            //word saveDir = name(counter) + "/polyMesh/";
            //cp(runTime.path()/runTime.timeName(), runTime.path()/saveDir);
            //std::ofstream of(folder + name(counter) + "/" + runTime.timeName());
            Ufield.append(U.clone() );
            Pfield.append(p.clone() );
            Dfield.append(sDRBMS().pointDisplacement().clone());
            NormalFields.append(N.clone());
            // Check if this is the last time step

            if (runTime.time().value() + runTime.deltaT().value() >= finalTime)
            {
                Info << "===== Storing final mesh =====" << endl;
                meshes.append(meshPtr.ptr());  // Transfer ownership
                // NOTE: meshPtr is now empty! Handle accordingly.
            }

            counter++;
            nextWrite += writeEvery;
        }
    }

    //exit(0);
    //const pointField& points2 = mesh.points();
    //Foam::meshToMesh0 mapper(U.mesh(), meshPtr());
    //Foam::MapConsistentMesh(U.mesh(), meshPtr(),
    //meshToMesh0::order::MAP);
    //Info<< nl
    //    << "Consistently creating and mapping fields for time "
    //    << mesh.time().timeIndex() << nl << endl;
    /*
    bool meshesDiffer = false;
    const pointField& points1 = mesh0.points();  // e.g., original mesh
    const pointField& points2 = mesh.points();  // e.g., moved mesh
    forAll(points1, i)
    {
        if (mesh0.V()[i] != mesh.V()[i])  // or a tolerance like 1e-10
        {
            meshesDiffer = true;
            break;
        }
    }
    if (meshesDiffer)
        Info << "Meshes differ in point locations." << endl;
    else
    Info << "Meshes are geometrically identical." << endl;
    exit(0);*/
    /// Store the mesh (transfer ownership)
    //meshes.append(meshPtr.ptr());
    //meshes.append(meshPtr.release()); // Releases ownership from autoPtr to PtrList
}


void fsiBasic::restart()
{
    turbulence.clear();
    _laminarTransport.clear();
    _p.clear();
    _U.clear();
    _phi.clear();
    _Uf.clear();
    _pointDisplacement.clear();
    sDRBMS.clear();
    _fvOptions.clear();
    _pimple.clear();
    argList& args = _args();
    Time& runTime = _runTime();
    instantList Times = runTime.times();
    runTime.setTime(0, 1);
    //runTime.stopAt(runTime.stopAtControls::saEndTime);
    //meshPtr().resetMotion();
    //meshPtr().movePoints(point0);
    //pointField& pointOld = const_cast<pointField&> (meshPtr().oldPoints());
    //pointOld = point0;
    /// Recreating the mesh
    Info << "ReCreating dynamic mesh for time = "
         << runTime.timeName() << nl << endl;
    meshPtr = autoPtr<dynamicFvMesh> (dynamicFvMesh::New(args, runTime));
    //meshPtr.reset(newMesh);
    // Take ownership from another autoPtr
    //meshPtr.reset(newMesh.ptr()); // Transfers ownership
    Foam::dynamicFvMesh& mesh = meshPtr();
    //exit(0);
    _pimple = autoPtr<pimpleControl>
              (
                  new pimpleControl
                  (
                      mesh
                  )
              );
#include "createFields.H"
    /// Reset the counter to zero
    counter = 1;
    /// clear list data members
}

void fsiBasic::change_viscosity(double mu)
{
    const volScalarField& nu =  _laminarTransport().nu();
    volScalarField& mu_new = const_cast<volScalarField&>(nu);
    this->assignIF(mu_new, mu);

    for (label i = 0; i < mu_new.boundaryFieldRef().size(); i++)
    {
        this->assignBC(mu_new, i, mu);
    }
}

void fsiBasic::change_stiffness(scalar& mu)
{
    dictionary& dictCoeffs = dyndict->subDict("sixDoFRigidBodyMotionCoeffs");
    dictionary& restraints = dictCoeffs.subDict("restraints");
    dictionary& spring = restraints.subDict("verticalSpring1");
    scalar stiffness = spring.get<scalar>("stiffness");
    Info << "==== stiffness ==== " << stiffness << endl;
    // Set new stiffness value (e.g., 0.1)
    spring.set("stiffness", mu);
    stiffness = spring.get<scalar>("stiffness");
    // Verify the change
    //scalar newStiffness = spring.lookup<scalar>("stiffness");
    //Info << "==== New stiffness ==== " << newStiffness << endl;
    //scalar stiffness = spring.lookupOrDefault<scalar>("stiffness", 0.0);
    Info << "==== New stiffness ==== " << stiffness << endl;
    //dyndict->write();
    dyndict->regIOobject::write();
}

void fsiBasic::exportFoamFieldToNpy(const word& outputDir,
                                    const word& fileName,
                                    const List<scalar>& foamField)
{
    Eigen::VectorXd  eigenData = Foam2Eigen::field2Eigen(foamField);
    cnpy::save(eigenData, outputDir + "/" + fileName + ".npy");
}




void fsiBasic::prepareFoamData(const word& outputPath)
{
    word fullPath = "./" + outputPath;

    if (!ITHACAutilities::check_folder(fullPath))
    {
        mkDir(fullPath);
        exportFoamFieldToNpy(fullPath, "fomforcex",     this->fomforcex);
        exportFoamFieldToNpy(fullPath, "fomforcey",     this->fomforcey);
        exportFoamFieldToNpy(fullPath, "CentreOfMassY", this->centerofmassy);
    }
}


void fsiBasic::loadCentreOfMassY(const fileName& baseDir)
{
    fileNameList dirs = readDir(baseDir, fileName::DIRECTORY);
    Eigen::VectorXd vecOfCentreOfMasses;
    // First, filter and sort the directories
    wordList sortedDirs;
    forAll(dirs, i)
    {
        if (dirs[i].find("DataFromFoam_") == 0)
        {
            sortedDirs.append(dirs[i]);
        }
    }
    // Custom sorting function to sort numerically
    std::sort(sortedDirs.begin(), sortedDirs.end(),
              [](const word & a, const word & b)
    {
        int numA = std::stoi(a.substr(std::string("DataFromFoam_").length()));
        int numB = std::stoi(b.substr(std::string("DataFromFoam_").length()));
        return numA < numB;
    });
    // Now process in sorted order
    forAll(sortedDirs, i)
    {
        fileName targetFile = baseDir / sortedDirs[i] / "CentreOfMassY.npy";

        if (exists(targetFile))
        {
            Info << "Found CentreOfMassY.npy in " << targetFile << endl;
            Eigen::MatrixXd CentreOfMassY;
            cnpy::load(CentreOfMassY, targetFile);
            Eigen::VectorXd currentVec = CentreOfMassY.col(0);
            // std::cout << "CentreOfMassY dimensions: "
            //          << CentreOfMassY.rows() << " x " << CentreOfMassY.cols() << std::endl;
            vecOfCentreOfMasses.conservativeResize(vecOfCentreOfMasses.size() +
                                                   currentVec.size());
            vecOfCentreOfMasses.tail(currentVec.size()) = currentVec;
            Info << "Loaded: " << targetFile << endl;
        }
        else
        {
            Info << "CentreOfMassY.npy not found in " << targetFile << endl;
        }
    }
    this->CylDispl = vecOfCentreOfMasses;
    // std::cout << "=== vecOfCentereOfMasses size ===" << vecOfCentreOfMasses.size() << std::endl;
}


void fsiBasic::updateStiffnessAndRebuildSolver(scalar& newMu, word param_name)
{
    dictionary& dictCoeffs = dyndict().subDict("sixDoFRigidBodyMotionCoeffs");
    dictionary& restraints = dictCoeffs.subDict("restraints");
    dictionary& spring = restraints.subDict("verticalSpring1");
    scalar oldMu = spring.get<scalar>(param_name);
    //Info << ">>> Replacing damping: " << oldMu << " -> " << newMu << endl;
    spring.set(param_name, newMu);
    /// Optional: Write back to disk
    dyndict().regIOobject::write(true);
    /// Clear and recreate solver and dictionary
    sDRBMS.clear();
    dyndict.clear();
    /// Recreate dictionary
    dyndict = autoPtr<IOdictionary>
              (
                  new IOdictionary
                  (
                      IOobject
                      (
                          "dynamicMeshDict",
                          meshPtr().time().constant(),
                          meshPtr(),
                          IOobject::MUST_READ,
                          IOobject::NO_WRITE,
                          false
                      )
                  )
              );
    /// Recreate solver
    sDRBMS = autoPtr<sixDoFRigidBodyMotionSolver>
             (
                 new sixDoFRigidBodyMotionSolver(meshPtr(), dyndict())
             );
    Info << ">>> Motion solver rebuilt with new stiffness.\n";
}
/*
template<class Type, class PatchField, class GeoMesh>
void mapFieldToMesh0(const Foam::GeometricField<Type, PatchField, GeoMesh>& sourceField,
const Foam::fvMesh& mesh0,
Foam::GeometricField<Type, PatchField, GeoMesh>& mappedField)
{
    const Foam::fvMesh& fromMesh = sourceField.mesh();

    if (&fromMesh == &mesh0)
    {
        Info << "[mapFieldToMesh0] Meshes are identical. Copying field..." << endl;
        mappedField = sourceField;
        return;
    }

    Info << "[mapFieldToMesh0] Interpolating field " << sourceField.name()
         << " to mesh0..." << endl;

    Foam::meshToMesh0 mapper(fromMesh, mesh0);

    Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>> interpField =
        mapper.interpolate<Type, Foam::plusEqOp<Type>>(
            Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>(sourceField),
            Foam::meshToMesh0::order::LINEAR,
            Foam::plusEqOp<Type>());

    mappedField = interpField;
}
mapFieldToMesh0( const volScalarField& sourceField,
                 const Foam::fvMesh& mesh0,
                 volScalarField& mappedField);

mapFieldToMesh0( const volVectorField& sourceField,
                 const Foam::fvMesh& mesh0,
                 volVectorField& mappedField); */



