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
    12simpleSteadyNS.C
\*---------------------------------------------------------------------------*/

#include "unsteadyNS.H"
#include "ITHACAstream.H"
#include "ITHACAPOD.H"
#include "forces.H"
#include "IOmanip.H"
//#include "fvCFD.H" //already in SteadyNSSimple.H

#include "dynamicFvMesh.H"
//#include "singlePhaseTransportModel.H"   // already present already in SteadyNSSimple.H
#include "turbulenceModel.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H" //already present already in SteadyNSSimple.H
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
//#include "createUfIfPresent.H"

//#include "points0MotionSolver.H" // added

class tutorial22 : public unsteadyNS
{
public:
    /// Constructor
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
        #include "createMesh.H"
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
        #include "createFields.H"
        #include "createFvOptions.H"
    }

    /// Velocity field
    volVectorField& U;
    /// Pressure field
    volScalarField& p;
    ///
    surfaceScalarField& phi;
    //surfaceVectorField& Uf;

    int folderN;
    int saver;
    int middleStep;
    int middleExport;
    /// Point motion field
    //mutabe pointVectorField pointDisplacement_;



    /// Perform an Offline solve
    /// Reimplement this starting from pimplefoam
    void offlineSolve()
    {
        Vector<double> inl(0, 0, 0);
        List<scalar> mu_now(1);

        // if the offline solution is already performed read the fields
        if (offline)
        {
            ITHACAstream::read_fields(Ufield, U, "./ITHACAoutput/Offline/");
            ITHACAstream::read_fields(Pfield, p, "./ITHACAoutput/Offline/");
            mu_samples =
                ITHACAstream::readMatrix("./ITHACAoutput/Offline/mu_samples_mat.txt");
        }
        else
        {
            Vector<double> Uinl(1, 0, 0);
            label BCind = 0;

            for (label i = 0; i < mu.cols(); i++)
            {
                mu_now[0] = mu(0, i);
                change_viscosity(mu(0, i));
                assignIF(U, Uinl);
                truthSolve3(mu_now); //truthSolve2 initial
            }
        }
    }

    void truthSolve3(List<scalar> mu_now, word Folder = ".")
    {
        Time& runTime = _runTime();
        volScalarField& p = _p();
        volVectorField& U = _U();
        fvMesh& mesh = _mesh();
        surfaceScalarField& phi = _phi();
        //surfaceVectorField& Uf = _Uf();
        //simpleControl& simple = _simple();
	#include "initContinuityErrs.H"
	fv::options& fvOptions = _fvOptions();
	pimpleControl& pimple = _pimple();
	IOMRFZoneList& MRF = _MRF();
	singlePhaseTransportModel& laminarTransport = _laminarTransport();
	instantList Times = runTime.times();
	runTime.setEndTime(finalTime);
        //correctPhi = _pimple().dict().getOrDefault; // for correctPhi
        //singlePhaseTransportModel& laminarTransport = _laminarTransport();
        scalar residual = 1;
        scalar checkMeshCourantNo = 1;
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
        middleStep = para->ITHACAdict->lookupOrDefault<label>("middleStep", 20);
        middleExport = para->ITHACAdict->lookupOrDefault<bool>("middleExport", true);
        //***************************pimpleFoam algorithm******************

        turbulence->validate();
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        Info << "\nStarting time loop\n" << endl;

        while (runTime.run())
        {
            //#include "readDyMControls.H"
            #include "readTimeControls.H" // included in readDyMControls
            #include "CourantNo.H"
            #include "setDeltaT.H"

            runTime++;

            Info << "Time = " << runTime.timeName() << nl << endl;

            // --- Pressure-velocity PIMPLE corrector loop
            while (pimple.loop())
            {
                /*if (pimple.firstIter() || moveMeshOuterCorrectors)
                {
              
                    // Do any mesh changes      
                    mesh.update();

                    if (mesh.changing())
                    {
                        MRF.update();

                        if (correctPhi)
                        {
                            // Calculate absolute flux
                            // from the mapped surface velocity
                            //#include "createUf.H"
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
                }*/

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

            if (middleExport)
            {
                ITHACAstream::exportSolution(U, name(folderN + 1), Folder + name(counter));
                ITHACAstream::exportSolution(p, name(folderN + 1), Folder + name(counter));
            }
            else
            {
                ITHACAstream::exportSolution(U, name(counter), Folder);
                ITHACAstream::exportSolution(p, name(counter), Folder);
            }

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
            mu_samples.conservativeResize(mu_samples.rows() + 1, mu_now.size());

            for (label i = 0; i < mu_now.size(); i++)
            {
                mu_samples(mu_samples.rows() - 1, i) = mu_now[i];
            }

            // Resize to Unitary if not initialized by user (i.e. non-parametric problem)
            if (mu.cols() == 0)
            {
                mu.resize(1, 1);
            }

            if (mu_samples.rows() == mu.cols())
            {
                ITHACAstream::exportMatrix(mu_samples, "mu_samples", "eigen",
                                           Folder);
            }

        }
    }

};

int main(int argc, char* argv[])
{
    // Construct the tutorial object
    tutorial22 example(argc, argv);
    // Read some parameters from file
    ITHACAparameters* para = ITHACAparameters::getInstance(example._mesh(),
                             example._runTime());
    int NmodesUout = para->ITHACAdict->lookupOrDefault<int>("NmodesUout", 15);
    int NmodesPout = para->ITHACAdict->lookupOrDefault<int>("NmodesPout", 15);
    int NmodesUproj = para->ITHACAdict->lookupOrDefault<int>("NmodesUproj", 10);
    int NmodesPproj = para->ITHACAdict->lookupOrDefault<int>("NmodesPproj", 10);
    // Read the par file where the parameters are stored
    word filename("./par");
    example.mu = ITHACAstream::readMatrix(filename);
    // Set the inlet boundaries patch 0 directions x and y
    example.inletIndex.resize(1, 2);
    example.inletIndex(0, 0) = 0;
    example.inletIndex(0, 1) = 0;
    // Perform the offline solve
    example.offlineSolve();
    ITHACAstream::read_fields(example.liftfield, example.U, "./lift/");
    // Homogenize the snapshots
    example.computeLift(example.Ufield, example.liftfield, example.Uomfield);
    // Perform POD on velocity and pressure and store the first 10 modes
    ITHACAPOD::getModes(example.Uomfield, example.Umodes, example._U().name(),
                        example.podex, 0, 0,
                        NmodesUout);
    ITHACAPOD::getModes(example.Pfield, example.Pmodes, example._p().name(),
                        example.podex, 0, 0,
                        NmodesPout);
    // // Create the reduced object
    // reducedSimpleSteadyNS reduced(example);
    // PtrList<volVectorField> U_rec_list;
    // PtrList<volScalarField> P_rec_list;
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
