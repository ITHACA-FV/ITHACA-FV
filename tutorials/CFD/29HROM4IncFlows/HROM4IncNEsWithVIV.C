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
    Tutorial of flow around a moving cylinder
\*---------------------------------------------------------------------------*/
// #include <iostream>
// #include "fvCFD.H"
#include "unsteadyNS.H"
#include "ITHACAPOD.H"
#include "ITHACAstream.H"
#include "ITHACAparameters.H"
//#include "ReducedProblem.H"
#include "pointConstraints.H"
#include <chrono>
#include "Foam2Eigen.H"
#include "DEIM.H"
#include "ReducedSimpleSteadyNS.H"
#include <omp.h>

template<typename O>
class MyDEIM : public DEIM<O>
{
    public:
        using DEIM<O>::DEIM;
        // MyDEIM(PtrList<O>& SnapShots, label NModes, word FunctionName,
        //        word FieldName): DEIM(SnapShots, NModes, FunctionName, FieldName)
        // {}
        PtrList<O> fields;
        autoPtr<O> subField;
};

class tutorial27: public unsteadyNS
{
public:
        explicit tutorial27(int argc, char* argv[])
        : unsteadyNS(argc, argv), 
        U(_U()), p(_p())
    {
    //point0 = meshPtr().points();
    // DeimU = autoPtr<MyDEIM<volVectorField>> 
    //         (   new MyDEIM<volVectorField>
    //             ( Ufield, 10, U.name() + "_" + "solution", U.name()
    //             )
    //         );
     //std::cout << "################ tutorial27 ##################" << std::endl;        
    }
    /// Velocity field
    volVectorField& U;
    /// Pressure field
    volScalarField& p;
    /// pointDisplacement field
    //pointVectorField& pd;
    /// Hyper-reduced objects
    autoPtr<MyDEIM<volVectorField>> DeimU;
    autoPtr<MyDEIM<volScalarField>> DeimP;
    /// Initial coordinates of the grid points
    ///vectorField point0;
    void offlineSolve(word folder="./ITHACAoutput/Offline/")
    {
        List<scalar> mu_now(1);

        if (offline)
        {
            ITHACAstream::read_fields(Ufield, U, folder);
            ITHACAstream::read_fields(Pfield, p, folder);
            //ITHACAstream::read_fields(Dfield, pd, folder);
        }    
        else
        {
            mu_now[0] = mu(0, 0); 
            truthSolve(mu_now, folder);
        }

    }
    void HROM(int NmodesU, int NmodesP)
    {
        // std::cout << "##################################" << std::endl;
        fvMesh& mesh = _mesh();
        // Info << Ufield.size() << endl;
        // Info << Pfield.size() << endl;
        DeimU = autoPtr<MyDEIM<volVectorField>> 
        (
            new MyDEIM<volVectorField>
            (
               Ufield, NmodesU, _U().name(), _U().name()
            )
        );
        DeimP = autoPtr<MyDEIM<volScalarField>> 
        (
            new MyDEIM<volScalarField>
            (
               Pfield, NmodesP, p.name(), p.name()
            )
        );
        /// Generate the submeshes of each fields
        DeimU->generateSubmesh(1, U.mesh(), U);
        DeimP->generateSubmesh(1, p.mesh(), p);
        
    }

}; 

class HyperROMFsi: public reducedSimpleSteadyNS
{
public:
    explicit HyperROMFsi(tutorial27& FoamPb): problem(&FoamPb)
        
    {
        
            startTime =  problem->startTime;
            finalTime =  problem->finalTime;
            timeStep =   problem->timeStep;
            writeEvery = timeStep;
            nextWrite = 0.0;
            //std::cout << "################ Ctor of POD-I Fsi ##################" << std::endl;
    }

    tutorial27* problem;
    /// time control variables
    scalar startTime;
    scalar finalTime;
    scalar timeStep;
    scalar writeEvery;
    scalar nextWrite;

    PtrList<volScalarField> PredFields;
    PtrList<volVectorField> UredFields;
    label counter = problem->counter;
    Eigen::MatrixXd DeimU;
    Eigen::MatrixXd DeimP;
    double time_rom_u, time_rom_p;

    void solveOnline_Pimple(scalar mu_now,
                            int NmodesUproj, 
                            int NmodesPproj, 
                            fileName folder = "./ITHACAoutput/Online/")
    {
        Time& runTime = problem->_runTime();
        fvMesh& mesh = problem->_mesh();
        fv::options& fvOptions = problem->_fvOptions();
        pimpleControl& pimple = problem->_pimple();
        volScalarField& p = problem->_p();
        volVectorField& U = problem->_U();
        surfaceScalarField& phi = problem->_phi();
        //IOMRFZoneList& MRF = problem->_MRF();
        singlePhaseTransportModel& laminarTransport = problem->_laminarTransport();
        autoPtr<incompressible::turbulenceModel> turbulence
        (
            incompressible::turbulenceModel::New
            (
                U, 
                phi, 
                laminarTransport
            )
        );
        instantList Times = runTime.times();
        runTime.setEndTime(finalTime);
        runTime.setTime(Times[1], 1);
        runTime.setDeltaT(timeStep);
        nextWrite = startTime;
        label pRefCell = 0;
        scalar pRefValue = 0.0;
        scalar  cumulativeContErr = 0.0;
       /// SparseMatrices with corresponding vectors
       Eigen::SparseMatrix<double,Eigen::RowMajor> Su;
       Eigen::SparseMatrix<double,Eigen::RowMajor> Sp;

       Eigen::VectorXd su, sp;
       List<label> UMagicPoints = problem->DeimU->magicPoints();
       List<label> PMagicPoints = problem->DeimP->magicPoints();
       //scalarField Volumes = mesh.V();

       // Eigen::VectorXd selected_volumes_u(UMagicPoints.size());
       // Eigen::VectorXd selected_volumes_p(PMagicPoints.size());
       // /// TODO: to be compute Offline
       // for (int k = 0; k < UMagicPoints.size(); ++k) 
       // {
       //      selected_volumes_u(k) = 1.0 / mesh.V()[UMagicPoints[k]];
       // }

       // for (int k = 0; k < PMagicPoints.size(); ++k) 
       // {
       //      selected_volumes_p(k) = 1.0 / mesh.V()[PMagicPoints[k]];
       // }
       //const auto& DeimU = problem->DeimU->U.transpose();
       //const auto& DeimP = problem->DeimP->U.transpose();
        /// adjustTimeStep
        bool adjustTimeStep = problem->adjustTimeStep;
        /// maxCourant
        scalar maxCo = problem->maxCo;
        /// maxDeltaT
        scalar maxDeltaT = problem->maxDeltaT;
        // Eigen::VectorXd a0 = Foam2Eigen::field2Eigen(problem->_U());
        // Eigen::VectorXd b0 = Foam2Eigen::field2Eigen(problem->_p());
        // Eigen::VectorXd p_ref = DeimP*b0;
        // Eigen::VectorXd u_ref = DeimU*a0;

        // std::cout << "p_ref = \n" << p_ref << std::endl;
        // std::cout << "u_ref = \n" << u_ref << std::endl;
        // //std::cout << "a = \n" << a << std::endl;
        // exit(0);

        // Tikhonov regularization: to stabilize the solution: if necessary
        //B += 1e-6 * Eigen::MatrixXd::Identity(B.rows(), B.cols());
        double lambda0 = 1.0; // Initial weight
        double alpha = 1000;   // Decay rate (for exponential decay)
        double time = 0.0;    // Current simulation time
        double dt = problem->_runTime().value();      // Time step size
        /// PIMPLE algorithm starts here
Info<< "\nStarting time loop\n" << endl;
        while (runTime.run())
        {
#include "readTimeControls.H"
#include "CourantNo.H"
#include "setDeltaT.H"
            runTime.setEndTime(finalTime);
            runTime++;
            //p.storePrevIter();
            Info << "Time = " << runTime.timeName() << nl << endl;
            // --- Pressure-velocity PIMPLE corrector loop
            while (pimple.loop())
            {
// auto t1 = std::chrono::high_resolution_clock::now(); 
double lambda_t = lambda0 * std::exp(-alpha * problem->_runTime().value());               
#include "Reduced_UEqn.H"
// auto t2 = std::chrono::high_resolution_clock::now();
// auto time_span_u = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    //Eigen::VectorXd diag = Eigen::Map<Eigen::VectorXd>(UEqn.diag().data(), mesh.C().size());    

                // Pressure corrector loop
                while (pimple.correct())
                {
//auto t3 = std::chrono::high_resolution_clock::now(); 
#include "Reduced_pEqn.H"
//auto t4 = std::chrono::high_resolution_clock::now();
//auto time_span_p = std::chrono::duration_cast<std::chrono::duration<double>>(t4 - t3);
//time_rom_p +=time_span_p.count();
                }

                if (pimple.turbCorr())
                {
                    laminarTransport.correct();
                    turbulence->correct();
                }
            }
            if(checkWrite(runTime))
            {
                //ITHACAstream::exportSolution(U, name(counter), folder);
                //ITHACAstream::exportSolution(p, name(counter), folder);
                //std::ofstream of(folder + name(counter) + "/" + runTime.timeName());
                UredFields.append(U.clone());
                PredFields.append(p.clone());
                counter++;
                nextWrite += writeEvery;
            }

        }

    } 

    bool checkWrite(Time& timeObject)
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

    std::tuple<Eigen::MatrixXd, Eigen::VectorXd> HyperReducedSys(
         const Eigen::SparseMatrix<double, Eigen::RowMajor>& S,
         const Eigen::VectorXd& se,
         const List<label>& MagicPoints,
         const Eigen::MatrixXd& Modes)
        {
            const int s = MagicPoints.size();
            const int m = Modes.rows(); // number of modes
            Eigen::MatrixXd B(s,m);
            Eigen::VectorXd b(s);
            //tripletList.reserve(s*(S.nonZeros()/S.rows()) );
            //Estimate average non-zeros per row
            Eigen::VectorXd tempB(m);
            for (int k = 0; k < s; ++k) 
            {                
                int rowIdx = MagicPoints[k];
                tempB.setZero();
                // Process row directly without temporary
                for (Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator it(S,rowIdx); it; ++it)
                {
                     //tripletList.emplace_back(Trip(k, it.index(), it.value() ));
                     tempB += it.value()*Modes.col(it.index());
                }
                double volume = 1.0/problem->_mesh().V()[rowIdx];
                B.row(k) = std::sqrt(volume)*tempB;
                b(k) = std::sqrt(volume)*se(MagicPoints[k]);
            }
            return std::make_tuple(B,b);
        }
    
};
/*----------------------------------------------------------------------------------------------------------*\
                               Starting the MAIN
\*-----------------------------------------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
    // Construct the tutorial object
    tutorial27 example(argc, argv);
    tutorial27 online(argc,argv);
    //exit(0);
    std::clock_t startOff;
    double durationOff;
    // Read some parameters from file
    ITHACAparameters* para = ITHACAparameters::getInstance
                             (example._mesh(),
                              example._runTime());                        

    int NmodesUout  =  readInt(para->ITHACAdict->lookup("NmodesUout"));
    int NmodesPout  =  readInt(para->ITHACAdict->lookup("NmodesPout"));

    int NmodesUproj  = readInt(para->ITHACAdict->lookup("NmodesUproj"));
    int NmodesPproj  = readInt(para->ITHACAdict->lookup("NmodesPproj"));

    int finalTime  = readInt(para->ITHACAdict->lookup("finalTime"));
    //scalar timeStep  = readInt(para->ITHACAdict->lookup("timeStep"));
    //scalar writeEvery  = readInt(para->ITHACAdict->lookup("writeEvery"));
    // Read the par file where the parameters are stored
    word filename("./par");
    example.mu = ITHACAstream::readMatrix(filename);
    // Time parameters: We can use Ioodictionnary to access time parameters
    example.startTime  = 0;
    example.finalTime  = 5; // 2period
    example.timeStep   = 0.002; ///

    example.writeEvery = 1e-02;

    // //Perform the offline solve
    auto t1 = std::chrono::high_resolution_clock::now();
    example.offlineSolve();
    
    auto t2 = std::chrono::high_resolution_clock::now();
    auto time_span_off = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    // std::cout << "time_rom=  " << reduced.time_rom << std::endl;
    std::cout << "The Offline  phase  duration =" << time_span_off.count()<< std::endl;
    
    // if(example.podex==0 )
    // {
    //    ITHACAPOD::getModes(example.Ufield, online.Umodes, example._U().name(),
    //                 example.podex, 0, 0, NmodesUout);
    // }
    // else
    // {
    //   ITHACAstream::read_fields(online.Umodes, example._U(), "./ITHACAoutput/POD/");
    // }
    // exit(0);
    online.Ufield = example.Ufield;
    online.Pfield = example.Pfield;
    /// Construct  the DEIM object
    online.HROM(NmodesUproj,NmodesPproj);

    List<label>bCellsU,bCellsP;
    forAll(example._U().boundaryField(), patchI)
    {
        const fvPatch &pp = example._U().boundaryField()[patchI].patch();
        //if (pp.name() =="cylinder1")
        //{
            forAll(pp, faceI)
            {
                label cellI = pp.faceCells()[faceI];
                bCellsU.append(cellI);
                bCellsP.append(cellI);
            }
        //}
        
    }
    //
    for (int i = 0; i <online.DeimU->magicPoints().size() ; ++i)
    {
        bCellsU.append(online.DeimU->magicPoints()[i]);
    }
    for (int i = 0; i <online.DeimP->magicPoints().size() ; ++i)
    {
        bCellsP.append(online.DeimP->magicPoints()[i]);
    }
    /// Sorting and remove duplicates
    inplaceUniqueSort(bCellsU);
    inplaceUniqueSort(bCellsP);
    //Info << bCellsU.size() << endl;
    /// Updating magics points indices of each fields
    online.DeimU->magicPoints() = bCellsU;
    //online.DeimP->magicPoints() = bCellsP;
    //online.DeimU->uniqueMagicPoints() = bCellsU;
    List<label> xyz = online.DeimU->xyz();
    online.DeimU->setMagicPoints(bCellsU, xyz );
    //train.EliObject->uniqueMagicPoints() = ITHACAutilities::combineList(train.EliObject->totalMagicPoints());
    online.DeimU->submesh->setCellSubset(online.DeimU->uniqueMagicPoints(),true); //OKay

    volScalarField Indici
    (
        IOobject
        (
            online._U().name() + "_indices",
            online._mesh().time().timeName(),
            online._mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        online._mesh(),
        dimensionedScalar(online._U().name() + "_indices", dimensionSet(0, 0, 0, 0, 0, 0, 0),
                          0.0)
    );
    List<label> indices;
    int layers = 1;

    for (label i = 0; i < online.DeimU->magicPoints().size(); i++)
    {
        indices = ITHACAutilities::getIndices(online._mesh(), online.DeimU->magicPoints()[i], layers);
        online.DeimU->totalMagicPoints().append(indices);
    }

    online.DeimU->uniqueMagicPoints() = ITHACAutilities::combineList(online.DeimU->totalMagicPoints());
    scalar zerodot25 = 0.25;
    ITHACAutilities::assignIF(Indici, zerodot25,
                              online.DeimU->uniqueMagicPoints().List<label>::clone()());
    ITHACAutilities::assignONE(Indici, online.DeimU->magicPoints());
    ITHACAstream::exportSolution
    (Indici, "1", "./ITHACAoutput/DEIM/" + online._U().name()
    );
    //online.DeimP->magicPoints() = bCellsP;
    // online.DeimU->U = example.DeimU->U;
    // online.DeimP->U = example.DeimP->U;
    /// Restart the simulation
     //exit(0);
    /// ############### contruct the reduced the class object ###################
    HyperROMFsi hyperreduced(online);
    hyperreduced.DeimU = online.DeimU->U.transpose();
    hyperreduced.DeimP = online.DeimP->U.transpose();

    hyperreduced.startTime  =  example.startTime;
    hyperreduced.finalTime  =  example.finalTime;
    hyperreduced.timeStep   =   example.timeStep;
    hyperreduced.writeEvery = example.writeEvery;
    //Perform the online solutions
    scalar mu_now = example.mu(0, 0);
   auto t3 = std::chrono::high_resolution_clock::now();
    /// Solve the reduced problem
    hyperreduced.solveOnline_Pimple(mu_now, NmodesUproj, NmodesPproj);
auto t4 = std::chrono::high_resolution_clock::now();
auto time_span_online = std::chrono::duration_cast<std::chrono::duration<double>>(t4 - t3);
std::cout << "velocity_total_time=  " << hyperreduced.time_rom_u << std::endl;
std::cout << "pressure_total_time=  " << hyperreduced.time_rom_p << std::endl;
std::cout << "The Online  phase  duration  is  equal  to " << time_span_online.count() << std::endl;    
    for (int i = 0; i < hyperreduced.UredFields.size(); ++i)
    {
        ITHACAstream::exportSolution(hyperreduced.UredFields[i], name(i+1), "./ITHACAoutput/Online");
        ITHACAstream::exportSolution(hyperreduced.PredFields[i], name(i+1), "./ITHACAoutput/Online");
    }
    
    Eigen::MatrixXd errL2U = ITHACAutilities::errorL2Abs(example.Ufield, hyperreduced.UredFields);
    //Eigen::MatrixXd errL2U = ITHACAutilities::errorL2Rel(example.Ufield, hyperreduced.UredFields);
    std::cout << "======================= errL2U completed================================" << "\n";
    //Eigen::MatrixXd errL2P = ITHACAutilities::errorL2Rel(example.Pfield, hyperreduced.PredFields);
    Eigen::MatrixXd errL2P = ITHACAutilities::errorL2Abs(example.Pfield, hyperreduced.PredFields);


    exit(0);
}

