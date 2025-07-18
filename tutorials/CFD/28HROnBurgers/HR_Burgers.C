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
    Example of a heat transfer Reduction Problem
SourceFiles
    02thermalBlock.C
\*---------------------------------------------------------------------------*/

#include <iostream>
#include "fvCFD.H"
#include "IOmanip.H"
#include "Time.H"
#include "Burgers.H"
#include "ITHACAPOD.H"
#include "Foam2Eigen.H"
#include "DEIM.H"
#include <Eigen/Dense>
#define _USE_MATH_DEFINES
#include <cmath>
#include "simpleMatrix.H"
#include "topoSet.H"
#include "cellSet.H"
#include "fvMeshSubset.H"
//#include "cellZoneSet.H"
//#include<Eigen/SparseQR>

//#include "fvMeshSubsetter.H"  // Not fvMeshSubset (need two-step subsetting)
/// \brief Class where the tutorial number 2 is implemented.
/// \details It is a child of the laplacianProblem class and some of its
/// functions are overridden to be adapted to the specific case.

class EmpLagranInter : public DEIM<volVectorField>
{
    public:
        using DEIM::DEIM;
        PtrList<volVectorField> fields;
        autoPtr<volVectorField> subField;

};
class tutorial23: public Burgers
{
    public:
        explicit tutorial23(int argc, char* argv[])
            :
            Burgers(argc, argv),
            U(_U())
        {}

        /// Velocity field
        volVectorField& U;
        /// EIM object
        autoPtr<EmpLagranInter> EliObject;

        void offlineSolve(word folder = "./ITHACAoutput/Offline/")
        {
            if (offline)
            {
                ITHACAstream::read_fields(Ufield, U, "./ITHACAoutput/Offline/");
                //ITHACAstream::readMiddleFields(Ufield, U, "./ITHACAoutput/Offline/");

            }
            else
            {
                truthSolve(folder);
            }
        }

        void ELI(int NmodesU)
        {
            fvMesh& mesh  =  const_cast<fvMesh&>(U.mesh());
            EliObject = autoPtr<EmpLagranInter> (new EmpLagranInter(Ufield, NmodesU,"Usolution", U.name()));
            EliObject->generateSubmesh(1, U.mesh(), U);
            // EliObject->generateSubmesh(2, U.mesh(), U);
        }
        void truthSolve(word folder)
        {
            Time& runTime = _runTime();
            fvMesh& mesh = _mesh();
            volVectorField& U = _U();
            surfaceScalarField& phi = _phi();
            fv::options& fvOptions = _fvOptions();
            simpleControl& simple = _simple();
            dimensionedScalar& nu = _nu();
            counter = 1;
            //ITHACAstream::exportSolution(U, name(counter), folder);
            //counter++;
            nextWrite = startTime;
            nextWrite += writeEvery;
            while (_simple().loop())
            {
                Info << "Time = " << _runTime().timeName() << nl << endl;
                while (simple.correctNonOrthogonal())
                {
                     fvVectorMatrix  UEqn
                    (
                        fvm::ddt(U)
                        +fvm::div(phi,U)
                        - fvm::laplacian(nu,U)
                    );

                    //UEqn.relax();
                    UEqn.solve();
                }

                phi = linearInterpolate(U) & mesh.Sf();

                if (checkWrite(runTime))
                {
                    ITHACAstream::exportSolution(U, name(counter), folder );
                    counter++;
                    Ufield.append(U.clone());
                    nextWrite += writeEvery;
                }
                _runTime().printExecutionTime(Info);
            }
        }
};

class EliBurgers: public reductionProblem
{
public:

    //// Data members
    tutorial23* problem;
    scalar nextWrite;
    scalar startTime;
    scalar writeEvery;
    PtrList<volVectorField> Ufield;
    Eigen::MatrixXd Modes; // Matrix modes
    double time_rom=0.0;

    EliBurgers(){}
    explicit EliBurgers(tutorial23& FoamProblem) : problem(&FoamProblem)
    {
        nextWrite  = 0;
        startTime  = problem->startTime;
        writeEvery = problem->writeEvery;

        std::cout << "######## reduced Constructor calling ##########"<< std::endl;
    }

    ~EliBurgers(){}


    void OnlineBurgers(int NmodesUproj,word folder = "./ITHACAoutput/Online/")
    {
        Time& runTime = problem->_runTime();
        fvMesh& mesh = problem->_mesh();
        volVectorField& U = problem->_U();
        surfaceScalarField& phi = problem->_phi();
        //fv::options& fvOptions = problem->_fvOptions();
        simpleControl& simple = problem->_simple();
        //dimensionedScalar& nu = problem->_nu();
        counter = 1;
        ////
        nextWrite = startTime;
        nextWrite += writeEvery;
        
        Eigen::SparseMatrix<double,Eigen::RowMajor> S;
        Eigen::VectorXd se;

        List<label> uniqueMagicPoints = problem->EliObject->magicPoints();
        //const auto& Modes  = problem->EliObject->U.transpose();
        //Eigen::VectorXd volumes = ITHACAutilities::getMassMatrixFV(T);
        while (simple.loop() )
        {
            Info << "Time = " << problem->_runTime().timeName() << endl;

            while (simple.correctNonOrthogonal())
            {
                /// The matrix time approx = 8.08077

                fvVectorMatrix SEqn
                (
                    fvm::ddt(U)
                    +fvm::div(phi,U)
                    -fvm::laplacian(problem->_nu(),U) 
                );
                SEqn.relax();
                problem->_fvOptions().constrain(SEqn);

                /// Field Conversion time approx = 7.04161
                //Foam2Eigen::fvMatrix2Eigen(SEqn, S, se);
                Foam2Eigen::fvMat2Eigen(SEqn, S, se);
                //auto t1 = std::chrono::high_resolution_clock::now();  
                std::tuple<Eigen::MatrixXd,Eigen::VectorXd> HRSys = HyperReducedSys(S,se,
                                                            uniqueMagicPoints,Modes);
                //auto t2 = std::chrono::high_resolution_clock::now();
                //auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
                //time_rom +=time_span.count();                                                
                
                /// Solve HR problem +reconstruct + eigen2field: time approx = 3.74458
                //Eigen::VectorXd c = B.colPivHouseholderQr().solve(b);
                Eigen::VectorXd b = std::get<1>(HRSys);
                Eigen::MatrixXd B = std::get<0>(HRSys);
                
                Eigen::VectorXd c = B.partialPivLu().solve(b);  
                //auto t1 = std::chrono::high_resolution_clock::now();              
                /// Reconstruct the velocity
                Eigen::VectorXd z = problem->EliObject->U*c;
                
                //auto t2 = std::chrono::high_resolution_clock::now();
                //auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
                //time_rom +=time_span.count();
                
                /// Eigen2field
                U = Foam2Eigen::Eigen2field(problem->_U(), z );
                U.correctBoundaryConditions();
            }
            phi = linearInterpolate(U) & mesh.Sf();
            problem->_runTime().printExecutionTime(Info);

            if (checkWrite(problem->_runTime()))
            {
                //ITHACAstream::exportSolution(T, name(counter), folder);
                counter++;
                Ufield.append(U.clone());
                //OnlineCoeffs.append(a);
                nextWrite += writeEvery;
            }
        }

    }
    std::tuple<Eigen::MatrixXd, Eigen::VectorXd> HyperReducedSys(
        const Eigen::SparseMatrix<double, Eigen::RowMajor>& S,
        const Eigen::VectorXd& se, 
        const List<label>& uniqueMagicPoints, 
        const Eigen::MatrixXd& Modes)
    {

        //fvMesh& mesh = problem->_mesh();
        int s = uniqueMagicPoints.size();
        int m = Modes.rows(); // number of modes
        Eigen::MatrixXd B(s,m);
        Eigen::VectorXd b(s);
        Eigen::VectorXd tempB(m);
      
        for (int k = 0; k < uniqueMagicPoints.size(); ++k)
        {
            int magicPoint = uniqueMagicPoints[k];
            tempB.setZero();
            for ( Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(S, magicPoint); it; ++it )
            {
                tempB += it.value() * Modes.col(it.index());
            }
            // Apply scaling by volume
            double volume = 1.0/problem->_mesh().V()[magicPoint];
            B.row(k) = tempB * volume;
            b(k) = volume * se(magicPoint);
        }
       return std::make_tuple(B, b);

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
    
};

int main(int argc, char* argv[])
{
    // Create the train object of the tutorial02 type
    tutorial23 train(argc, argv);
    //tutorial23 test(argc, argv)

    // Read some parameters from file
    ITHACAparameters* para = ITHACAparameters::getInstance(train._mesh(),train._runTime());
    int NmodesUout  = readInt(para->ITHACAdict->lookup("NmodesUout"));
    int NmodesUproj  = readInt(para->ITHACAdict->lookup("NmodesUproj"));

    //startOff= std::clock();
    auto t3 = std::chrono::high_resolution_clock::now();
    train.offlineSolve();
    auto t4 = std::chrono::high_resolution_clock::now();
    auto time_span_off = std::chrono::duration_cast<std::chrono::duration<double>>(t4 - t3);

    std::cout << "The Offline  phase  duration =" << time_span_off.count()<< std::endl;
    
    /// Construct the DEIM indices.
    train.ELI(NmodesUproj);
    //List<label> test = {0,2791,10742,6295,7616,11418};//{0,10792,10824,10727,4,5,6,7,8,9}
    //Info << "test = "<< test <<endl;
    //exit(0);
    /// Restart the simulation
    train.restart();
    /// Construct the reduced problem
    EliBurgers hyperreduced(train);
    /// Transposing the matrix modes.
    hyperreduced.Modes = train.EliObject->U.transpose();
    auto t1 = std::chrono::high_resolution_clock::now();
    /// Solve the reduced problem
    hyperreduced.OnlineBurgers(NmodesUproj);
    auto t2 = std::chrono::high_resolution_clock::now();
    //auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1)
    auto time_span_online = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    std::cout << "HRSys=  " << hyperreduced.time_rom << std::endl;
    std::cout << "The Online  phase  duration =" << time_span_online.count()<< std::endl;
    //exit(0);
    /// Mapping to the whole fields

    for (int i = 0; i < hyperreduced.Ufield.size(); ++i)
    {
        ITHACAstream::exportSolution(hyperreduced.Ufield[i], name(i+1), "./ITHACAoutput/Online");
        volVectorField Ue("Ue", train.Ufield[i]-hyperreduced.Ufield[i]);
        //auto t = ITHACAutilities::L2Norm(train.Tfield[i]);
        //Te /=t;
        ITHACAstream::exportSolution(Ue, name(i+1), "./ITHACAoutput/AbsError");
    }
    
    Eigen::MatrixXd errL2U = ITHACAutilities::errorL2Abs(train.Ufield, hyperreduced.Ufield);
    std::cout << "#####################" << std::endl;
    Eigen::MatrixXd errL2RelU = ITHACAutilities::errorL2Rel(train.Ufield, hyperreduced.Ufield);
    
    exit(0);
}




//////////////////////////////////////////////////////////////////////////////////////////

// int main(int argc, char* argv[])
// {
//     // Create the train object of the tutorial02 type
//     tutorial23 train(argc, argv);
//     //tutorial23 test(argc, argv);
//     std::clock_t startOff;
//     double durationOff;
//     int dim  = 3;
//     label nCells = train._mesh().nCells();s
//     scalar zerodot25 = 0.25;

//     // Read some parameters from file
//     ITHACAparameters* para = ITHACAparameters::getInstance(train._mesh(),train._runTime());
//     int NmodesUout  = readInt(para->ITHACAdict->lookup("NmodesUout"));
//     int NmodesUproj  = readInt(para->ITHACAdict->lookup("NmodesUproj"));

//     startOff= std::clock();
//     train.offlineSolve();
//     durationOff = (std::clock() - startOff);
//     std::cout << "The Offline phase duration is equal to" << durationOff << std::endl;
//     // ITHACAPOD::getModes(train.Ufield, train.Umodes, train._U().name(),
//     //                    train.podex, 0, 0,NmodesUout);

//     Eigen::Map<Eigen::MatrixXd> test(train.Ufield[10].component(vector::X).ref().data(),1,train._U().size());
//     //Info << train.Ufield[10].component(vector::X).ref().data() << endl;
//      std::cout << test << std::endl;
//      exit(0);
//     /// Construct the DEIM indices.
//     train.ELI(NmodesUproj);
//     //exit(0);
//     /// Mask field to submesh
//     Eigen::SparseMatrix<double> field2submesh;
//     label M  = train.EliObject->uniqueMagicPoints().size();
//     field2submesh.resize(nCells*dim, M*dim);

//     for (unsigned int ith_subCell{0} ; ith_subCell < M; ith_subCell++)
//     {
//         for (unsigned int i= 0; i < dim; i++)
//         {
//           field2submesh.insert(train.EliObject->uniqueMagicPoints()[ith_subCell] + nCells*i, ith_subCell + i*M) = 1;
//         }
//     }
//     field2submesh.makeCompressed();
//     /// Important matrices
//     Eigen::MatrixXd X = train.EliObject->U; 
//     Eigen::MatrixXd P = train.EliObject->P;
//     Eigen::MatrixXd MatrixModes = train.EliObject->MatrixModes;  
//     /// Positive symmetric matrix
//     Eigen::MatrixXd phiTilde = train.EliObject->P.transpose()*train.EliObject->MatrixModes;
//     // Eigen::MatrixXd Q = phiTilde*(phiTilde.transpose()*phiTilde).inverse()*(phiTilde.transpose()*phiTilde).inverse()*phiTilde.transpose();
//     // std::cout << "phiTilde = " << phiTilde << std::endl;
//     // exit(0);
//     Eigen::MatrixXd V = (X.transpose()*field2submesh*field2submesh.transpose()*X).inverse()*X.transpose()*field2submesh; // Ok
//     /// Define the full matrix to the full field
//     Eigen::MatrixXd Y = X*V; 
//     //std::cout << "Y.rows() = " << Y.rows() << std::endl;
//     //std::cout << "Y.cols() = " << Y.cols() << std::endl;
//     train.EliObject->MatrixOnline = Y; 

//     //Eigen::MatrixXd A = X.transpose()*field2submesh*field2submesh.transpose()*X;
//     Eigen::MatrixXd A = phiTilde.transpose()*phiTilde;
// //std::cout << "A = " <<  A << std::endl;
//     Eigen::JacobiSVD<Eigen::MatrixXd> svd(A);
// double cond = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size()-1);
   
//     // PtrList<volVectorField> subfields;

//     // for (int i = 0; i < train.Ufield.size(); ++i)
//     // {

//     //   auto  subfld = train.EliObject->submesh().interpolate(train.Ufield[i] ).ref();
//     //   subfields.append(subfld);
       
//     // }
//     // ITHACAstream::exportFields(subfields, "./ITHACAoutput/SubFields", subfields[0].name() );  
    
//     //std::cout << "eigen values = " <<  svd.singularValues().transpose() << std::endl;
//     //std::cout << "condition number = " <<  cond << std::endl;

//     // forAll(train._mesh().boundary(), patch)
//     // {
//     //   //const word& patchName = train._mesh().boundary()[patch].name();  
//     //   //Info << "patchName = \t" <<  patchName << endl; 
//     // Eigen::MatrixXd B(Eigen::MatrixXd::Random(5,3)), Q;
//     // B.setRandom();
//     // Eigen::HouseholderQR<Eigen::MatrixXd> qr(B);
//     // Q = qr.householderQ();
//     // //thinQ = qr.householderQ() * thinQ; // <- here
//     // std::cout << "The complete unitary matrix Q is:\n" << Q << "\n\n";
//     //std::cout << "The thin matrix Q is:\n" << thinQ << "\n\n";
        
//     //   // // Loop over all faces of boundary patch
//     //   // forAll(train._mesh().boundary()[patch], facei)
//     //   // {
//     //   //   const label& bCell = train._mesh().boundaryMesh()[patch].faceCells()[facei];    // Boundary cell index
//     //   //   //const label& face = boundaryMesh[patch].start() + facei;        // Face index
//     //   //   //Info << "bCell = \t" <<  bCell << endl;
//     //   //   Info << "mesh.cellCells() [" << bCell << "]" << train._mesh().cellCells()[bCell] << endl;
//     //   //   // Do your calculations e.g.
//     //   //   // U.boundaryField()[patch][facei] = vector::zero;
//     //   // }
//     // }
//     // List<List<label>>bCells;
//     List<label>bCells;
//     // forAll(train._mesh().boundary(), patch)
//     // {
//     //   //const word& patchName = train._mesh().boundary()[patch].name();  
//     //   //Info << "patchName = \t" <<  patchName << endl; 

//     //   //Info << train._mesh().boundaryMesh()[patch].faceCells() << endl;
//     //   //exit(0);
//     //   //auto ll = train._mesh().boundaryMesh()[patch].faceCells();
//     //   //bCells.append(ll);
        
//     //   // Loop over all faces of boundary patch
//     //   forAll(train._mesh().boundary()[patch], facei)
//     //   {
//     //     /// Boundary cell index
//     //     const label& bCell = train._mesh().boundaryMesh()[patch].faceCells()[facei];
//     //     bCells.append(bCell);
//     //     //const label& face = boundaryMesh[patch].start() + facei;// Face index
//     //     //Info << "bCell = \t" <<  bCell << endl;
//     //     //Info <<"mesh.cellCells()["<<bCell<< "]"<< train._mesh().cellCells()[bCell] << endl;
//     //     // Do your calculations e.g.
//     //     // U.boundaryField()[patch][facei] = vector::zero;
//     //   }
//     //   // exit(0);
//     // }

//     // Assigning contribution from BC
//     forAll(train._U().boundaryField(), patchI)
//     {
//         const fvPatch &pp = train._U().boundaryField()[patchI].patch();
//         forAll(pp, faceI)
//         {
//             label cellI = pp.faceCells()[faceI];
//             bCells.append(cellI);
//         }
//     }



//     //PP.insert(306, 0) = 1;
//     for (int i = 0; i < train._U().size(); ++i)
//     {
        
//         PP.insert(i, 0) = 1;
//         Eigen::MatrixXd AA = PP.transpose()*train.EliObject->MatrixModes;
//         Eigen::MatrixXd BB = AA.transpose()*AA;

//         Eigen::JacobiSVD<Eigen::MatrixXd> svdd(BB);
//         double condd = svdd.singularValues()(0) / (1e-8 + svdd.singularValues()(svdd.singularValues().size()-1));
//         //if (condd <=5000)
//         //{
//             std::cout << "cell = " << i << std::endl;
//             std::cout << condd  << std::endl;
//         //}
        
//         //PP.conservativeResize(MatrixModes.rows(), i + 1);
//         //PP.insert(bCells[i], i) = 1;
//         //std::cout << svdd.singularValues()  << std::endl;
//     }



//     exit(0);
//     volScalarField Indici
//     (
//         IOobject
//         (
//             "indices",
//             train._mesh().time().timeName(),
//             train._mesh(),
//             IOobject::NO_READ,
//             IOobject::NO_WRITE
//         ),
//         train._mesh(),
//         dimensionedScalar("indices", dimensionSet(0, 0, 0, 0, 0, 0, 0),0.0)
//     );
    
    
//     train.EliObject->submesh->setCellSubset(bCells);
//     ITHACAutilities::assignONE(Indici, bCells );
//     ITHACAstream::exportSolution(Indici, "1", "./ITHACAoutput/BMesh/");

//     //Info << bCells << endl;


//     exit(0);
//     std::clock_t startOn;
//     double durationOn;
//     startOn = std::clock();

//     train.restart();
//     EliBurgers reduced(train);
//     /// Solve on the Hyper-reduced problem
//     reduced.OnlineBurgers(NmodesUproj);

//     durationOn = (std::clock() - startOn);
//     std::cout << "The Online  phase  duration  is  equal  to " << durationOn << std::endl;
//     //ITHACAstream::exportFields(reduced.Ufield, "./ITHACAoutput/Online", reduced.Ufield[0].name() );
//     //PtrList<volVectorField> errflds;
//     //exit(0);
// }




