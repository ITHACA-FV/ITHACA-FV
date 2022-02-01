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

#include "fvCFD.H"
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
#include "Switch.H"
#include "objectRegistry.H"
#include "localEulerDdtScheme.H"
#include "ReducedSimpleSteadyNS.H"
#include "reductionProblem.H"

class tutorial22 : public unsteadyNS
{
public:
    /// Constructors
    explicit tutorial22(int argc, char* argv[])
        :
        unsteadyNS(argc, argv),
        U(_U()),
        p(_p()),
        phi(_phi())

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
         #include "createFields.H" //Ok
        //#include "createFvOptions.H" // could be included in createFields.H
        //turbulence->read();
        para = ITHACAparameters::getInstance(mesh, runTime);
        // //pointVectorField & PointDisplacement = const_cast<pointVectorField&>(mesh.objectRegistry::lookupObject<pointVectorField>("pointDisplacement"));
        std::cout << "*******************copy constructor of the tutorial22***************"<< std::endl;                 


    }// end of the Constructor.
    
    // members data
    /// Velocity field
    volVectorField& U;
    /// Pressure field
    volScalarField& p;
    ///
    surfaceScalarField& phi;
    //surfaceVectorField& Uf;
    /// Dynamic mesh field
    Foam::autoPtr<Foam::dynamicFvMesh> meshPtr;
    // Dummy variable to transform pimplefoam into a class
    /// pimpleControl
    //autoPtr<pimpleControl> _pimple;

    /// Perform an Offline solve
    /// Reimplement this starting from pimplefoam
    void offlineSolve()
    {
        //std::cerr << "File: 22FSI.C, Line: 135"<< std::endl;
        Vector<double> inl(1, 0, 0); //inl(0, 0, 0)
        List<scalar> mu_now(1);

        // if the offline solution is already performed read the fields
          //if (offline)
          //{
             //ITHACAstream::read_fields(Ufield, U, "./ITHACAoutput/Offline/");
             //ITHACAstream::read_fields(Pfield, p, "./ITHACAoutput/Offline/");
             //mu_samples = ITHACAstream::readMatrix("./ITHACAoutput/Offline/mu_samples_mat.txt");
          //}
          //else
          //{
            //Vector<double> Uinl(1, 0, 0);
            //label BCind = 0;

            //for (label i = 0; i < mu.cols(); i++)
            //{
                mu_now[0] = mu(0, 0);
                // change_viscosity(mu(0, i));
                //assignIF(U, Uinl);
                truthSolve3(mu_now);
            //}
          //}
                        
        std::cout << "**************************offlineSolve method end****************************\n"<< std::endl; 

    }

    void truthSolve3(List<scalar> mu_now, word folder = "./ITHACAoutput/Offline")
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


        std::cout <<"************************The finalTime is ****************"<< finalTime << std::endl;           
        std::cout <<"************************The startTime is ****************"<<startTime << std::endl;                 
      
//****************************pimpleFoam algorithm******************************************
       
        //#include "postProcess.H"

        #include "addCheckCaseOptions.H"
        //#include "setRootCaseLists.H"
        //#include "createTime.H" //already done

        //#include "initContinuityErrs.H"
        #include "createDyMControls.H"
        //#include "createFields.H"
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
        std::cout<< "///////////////////////////////////////////beginning of the while for runTime////////////////////////"<< std::endl;

            #include "readDyMControls.H"
            #include "CourantNo.H"
            #include "setDeltaT.H"
            runTime.setEndTime(finalTime);
            std::cout<< "/////////////////////////////////Time = " <<finalTime << nl << std::endl;
            runTime++;

            //++runTime;
            //std::cout << runTime

            Info<< "Time = " << runTime.timeName() << nl << endl;
            std::cout<< "*************Time = " <<runTime.timeName() << nl << std::endl;

            // --- Pressure-velocity PIMPLE corrector loop
            while (pimple.loop())
            {
                std::cout<< "///////////////////////////////////////////beginning of the while for loop////////////////////////"<< std::endl;

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
                            std::cout<< "///////////////////////////////// inside the checkMeshCourantNo if cond/////////////////////////////"<< std::endl;
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
        std::cout << "*******************truthSolve3 method in the tutorial22***************\n"<< std::endl; 
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

    int NmodesUout = para->ITHACAdict->lookupOrDefault<int>("NmodesUout", 10);
    int NmodesPout = para->ITHACAdict->lookupOrDefault<int>("NmodesPout", 10);
    int NmodesSUPout = para->ITHACAdict->lookupOrDefault<int>("NmodesSUPout", 15);
    int NmodesUproj = para->ITHACAdict->lookupOrDefault<int>("NmodesUproj", 10);
    int NmodesPproj = para->ITHACAdict->lookupOrDefault<int>("NmodesPproj", 15);
    int NmodesSUPproj = para->ITHACAdict->lookupOrDefault<int>("NmodesSUPproj", 10);
   
    // Read the par file where the parameters are stored
    word filename("./par");
    example.mu = ITHACAstream::readMatrix(filename);
    //Set the inlet boundaries patch 0 directions x and y
    example.inletIndex.resize(1, 2);
    example.inletIndex(0, 0) = 0;
    example.inletIndex(0, 1) = 0;
    // Time parameters: W e can use Ioodictionnary to access time parameters
    example.startTime = 0;
    example.finalTime = 0.5;
    example.timeStep = 0.025;
    example.writeEvery = 0.025;

    //Perform the offline solve
    example.offlineSolve();
    // ITHACAPOD::getModes(example.Uomfield, example.Umodes, example._U().name(),
    //                         example.podex, 0, 0, NmodesUout);

    // ITHACAPOD::getModes(example.Pfield, example.Pmodes, example._p().name(),
    //                     example.podex, 0, 0,
    //                     NmodesPout);
    //Info<< example.Umodes;
    //Info << "example field: " <<  example.Uomfield << endl;

    //std::cout<< example.Ufield;
    // ITHACAstream::read_fields(example.liftfield, example.U, "./lift/");
    // // // Homogenize the snapshots
    // //example.computeLift(example.Ufield, example.liftfield, example.Uomfield);

    //Perform POD on velocity and pressure and store the first 10 modes
    //ITHACAPOD::getModes(example.Ufield, example.Umodes, example._U().name(),
    //                    example.podex, 0, 0,
    //                   NmodesUout);
    //ITHACAPOD::getModes(example.Pfield, example.Pmodes, example._p().name(),
    //                    example.podex, 0, 0,
    //                    NmodesPout);

    // Create the reduced object
    // reducedFSI reduced(example);

    // // Create the reduced object
    //reducedUnsteadyNS reduced(example);
    //reducedSimpleSteadyNS reduced(example);
    //PtrList<volVectorField> U_rec_list;
    //PtrList<volScalarField> P_rec_list;
    // // Reads inlet volocities boundary conditions.
    // word vel_file(para->ITHACAdict->lookup("online_velocities"));
    // Eigen::MatrixXd vel = ITHACAstream::readMatrix(vel_file);
    // //std::cout << "the vel size is: " << vel.size() << std::endl;
    // std::cout << " example size is: " << (example.mu).size() << std::endl;

    // //Perform the online solutions
    // for (label k = 0; k < (example.mu).size(); k++)
    // {
    //     scalar mu_now = example.mu(0, k);
    //     example.change_viscosity(mu_now);
    //     reduced.OnlineVelocity(vel);
    //     reduced.solveOnline_Pimple(mu_now, NmodesUproj, NmodesPproj);
    // }
    std::cout << "///////////////////////////////////////////////////////////////////////////////////////";
    // word temp(couplingDict.lookup("cylinder1"));
    // word interface = temp;
    // label fluidPatchID = example.meshPtr().boundaryMesh().findPatchID(interface);

    // std::cout<< "the label is: " << fluidPatchID;

    Info << "Hello there, the most recent time folder found is " << example._runTime().timeName() << nl
         << "The mesh has " << example.meshPtr().C().size() << " cells and " << example.meshPtr().Cf().size()
         << " internal faces in it. Wubalubadubdub!" << nl << endl;

    // // It's possible to iterate over every cell in a standard C++ for loop
    // for (label cellI = 0; cellI < example.meshPtr().C().size(); cellI++)
    //     if (cellI%20 == 0) // only show every twentieth cell not to spam the screen too much
    //         Info << "Cell " << cellI << " with centre at " << example.meshPtr().C()[cellI] << endl;
    // Info << endl; // spacer

    // Each cell is constructed of faces - these may either be internal or constitute a
    // boundary, or a patch in OpenFOAM terms; internal faces have an owner cell
    // // and a neighbour.
    // for (label faceI = 0; faceI < example.meshPtr().owner().size(); faceI++)
    //     if (faceI%40 == 0)
    //         Info << "Internal face " << faceI << " with centre at " << example.meshPtr().Cf()[faceI]
    //              << " with owner cell " << example.meshPtr().owner()[faceI]
    //              << " and neighbour " << example.meshPtr().neighbour()[faceI] << endl;
    // Info << endl;  

    const faceList& fcs = example.meshPtr().faces();
    const List<point>& pts = example.meshPtr().points();
    const List<point>& cents = example.meshPtr().faceCentres();

    // forAll(fcs,faceI)
    //     if (faceI%80==0)
    //     {
    //         if (faceI<example.meshPtr().Cf().size())
    //             Info << "Internal face ";
    //         else
    //         {
    //             forAll(example.meshPtr().boundary(),patchI)
    //                 if ((example.meshPtr().boundary()[patchI].start()<= faceI) &&
    //                     (faceI < example.meshPtr().boundary()[patchI].start()+example.meshPtr().boundary()[patchI].Cf().size()))
    //                 {
    //                     Info << "Face on patch " << patchI << ", faceI ";
    //                     break; // exit the forAll loop prematurely
    //                 }
    //         }

    //         Info << faceI << " with centre at " << cents[faceI]
    //              << " has " << fcs[faceI].size() << " vertices:";
    //         forAll(fcs[faceI],vertexI)
    //             // Note how fcs[faceI] holds the indices of points whose coordinates
    //             // are stored in the pts list.
    //             Info << " " << pts[fcs[faceI][vertexI]];
    //         Info << endl;
    //     }
    // Info << endl;


    label patchID(0);
    const polyPatch& pp = example.meshPtr().boundaryMesh()[patchID];
    if (isA<emptyPolyPatch>(pp))
    {
        // patch patchID is of type "empty".
        Info << "You will not see this." << endl;
    }

    // Patches may also be retrieved from the mesh using their name. This could be
    // useful if the user were to refer to a particular patch from a dictionary
    // (like when you do when calculating forces on a particular patch).
    word patchName("movingWall");
    patchID = example.meshPtr().boundaryMesh().findPatchID(patchName);
    Info << "Retrieved patch " << patchName << " at index " << patchID << " using its name only." << nl << endl;

    Info<< "End\n" << endl;

    exit(0);
}
