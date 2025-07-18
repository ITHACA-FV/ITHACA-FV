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
    Hyper-reduction for scalar Transport Problem
SourceFiles
    Tutorial24.C
\*---------------------------------------------------------------------------*/

#include <iostream>
#include "fvCFD.H"
//#include "IOmanip.H"
#include "Time.H"
#include "ScalarTransport.H"
//#include "ITHACAPOD.H"
#include "Foam2Eigen.H"
#include "DEIM.H"
#include <Eigen/Dense>
#include <cmath>
#include <omp.h> // Include OpenMP for parallelization
#include <vector>
#include <algorithm>

 
class EmpLagranInter : public DEIM<volScalarField>
{
    public:
        using DEIM::DEIM;
        PtrList<volScalarField> fields;
        autoPtr<volScalarField> subField;

};


// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //
class tutorial24 : public ScalarTransport
{
    public:
        /// Constructor
        explicit tutorial24(int argc, char* argv[])
            :
            ScalarTransport(argc, argv),T(_T())
            {
            // curX = _mesh().points();
            // point0 = _mesh().points();
             }
            // Relevant Fields
         volScalarField& T;
         autoPtr<EmpLagranInter> EliObject;
         // volScalarField& p;
         // volScalarField& nut;
         // pointVectorField& pD;

        /// Perform an Offline solve
        void offlineSolve(word folder="./ITHACAoutput/Offline/")
        {
                if (offline)
                {
                    ITHACAstream::read_fields(Tfield, T, folder );
                     // mu_samples =
                     // ITHACAstream::readMatrix("./ITHACAoutput/Offline/mu_samples_mat.txt");
                }
                else
                {
                    truthSolve(folder);
                    //restart();
                }    
        }
        void ELI(int NmodesT)
        {
            fvMesh& mesh  =  const_cast<fvMesh&>(T.mesh());
            EliObject = autoPtr<EmpLagranInter> (new EmpLagranInter(Tfield, NmodesT,"Tsol", T.name()));
            EliObject->generateSubmesh(1, T.mesh(), T);
            // EliObject->generateSubmesh(2, U.mesh(), U);
         }

};

class EliBurgers: public reductionProblem
{
public:

    //// Data members
    tutorial24* problem;
    scalar nextWrite;
    scalar startTime;
    scalar writeEvery;
    PtrList<volScalarField> Tfield;
    Eigen::MatrixXd Tmodes;
    //volVectorModes Lagrangian
    double time_rom;
    EliBurgers(){}
    explicit EliBurgers(tutorial24& FoamProblem) : problem(&FoamProblem)
    {
        nextWrite  = 0;
        startTime  = problem->startTime;
        writeEvery = problem->writeEvery;

        std::cout << "######## reduced Constructor calling ##########"<< std::endl;
    }

    ~EliBurgers(){}


    
   void OnlineBurgers(int NmodesUproj,word folder = "./ITHACAoutput/Online/")
    {
        //Time& runTime = problem->_runTime();
        fvMesh& mesh = problem->_mesh();
        volScalarField& T = problem->_T();
        surfaceScalarField& phi = problem->_phi();
        //fv::options& fvOptions = problem->_fvOptions();
        simpleControl& simple = problem->_simple();
        //dimensionedScalar& nu = problem->_nu();
        counter = 1;
        ////
        nextWrite = startTime;
        nextWrite += writeEvery;
        
        
        //Eigen::SparseMatrix<double> S;
        Eigen::SparseMatrix<double, Eigen::RowMajor> S;
        Eigen::VectorXd se;
        //const auto& Tmodes = problem->EliObject->U.transpose();

        List<label> uniqueMagicPoints = problem->EliObject->magicPoints();
        //Eigen::VectorXd volumes = ITHACAutilities::getMassMatrixFV(T);
        while (simple.loop() )
        {
            Info << "Time = " << problem->_runTime().timeName() << endl;
            //Info << "dt = " << problem->_runTime().value() << endl;

            while (simple.correctNonOrthogonal())
            {
                /// The matrix
                fvScalarMatrix SEqn
                (
                    fvm::ddt(T)
                    +fvm::div(phi,T)
                    -fvm::laplacian(problem->_nu(),T)
                );
                SEqn.relax();

                problem->_fvOptions().constrain(SEqn);
                /// Field Conversion
                //fvMatrix2Eigen(SEqn, S, se);
                Foam2Eigen::fvMat2Eigen(SEqn, S, se);
                S.makeCompressed();
                std::tuple<Eigen::MatrixXd,Eigen::VectorXd> HRSys = HyperReducedSys(S,se,uniqueMagicPoints,Tmodes);
                /// HR system
                Eigen::VectorXd b = std::get<1>(HRSys);
                Eigen::MatrixXd B = std::get<0>(HRSys);
                //double cond = B.jacobiSvd().singularValues()(0) /
                 //B.jacobiSvd().singularValues().tail(1)(0);
           //std::cout << "cond number = " << cond << std::endl;
                Eigen::VectorXd c = B.completeOrthogonalDecomposition().solve(b);
                //Eigen::VectorXd c = B.partialPivLu().solve(b);
                //Eigen::VectorXd c = B.completeOrthogonalDecomposition().solve(b);
                //Eigen::VectorXd c = B.colPivHouseholderQr().solve(b); //.
                //Eigen::VectorXd c=B.householderQr().solve(b);
                //Eigen::VectorXd c=B.llt().solve(b);
                /// Reconstruct the temperature
                Eigen::VectorXd z = problem->EliObject->U*c;
                ///Solution 2: Copy Eigen Data Using Iterators
                std::copy(z.data(), z.data() + z.size(), T.ref().begin());
                //T = Foam2Eigen::Eigen2field(problem->_T(), z );
                T.correctBoundaryConditions();
                
            }

            problem->_runTime().printExecutionTime(Info);

            if (checkWrite(problem->_runTime()))
            {
                //ITHACAstream::exportSolution(T, name(counter), folder);
                counter++;
                Tfield.append(T.clone());
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
    tutorial24 train(argc, argv);

    // Read some parameters from file
    ITHACAparameters* para = ITHACAparameters::getInstance(train._mesh(),train._runTime());
    int NmodesTout  = readInt(para->ITHACAdict->lookup("NmodesUout"));
    int NmodesTproj  = readInt(para->ITHACAdict->lookup("NmodesUproj"));

    //startOff= std::clock();
    auto t3 = std::chrono::high_resolution_clock::now();
    train.offlineSolve();
    auto t4 = std::chrono::high_resolution_clock::now();
    auto time_span_off = std::chrono::duration_cast<std::chrono::duration<double>>(t4 - t3);
    // std::cout << "time_rom=  " << reduced.time_rom << std::endl;
    std::cout << "The Offline  phase  duration =" << time_span_off.count()<< std::endl;
    //ITHACAPOD::getModes(train.Tfield, train.Tmodes, train._T().name(),
    //                    train.podex, 0, 0, NmodesUproj);
    /// Construct the DEIM indices.
    train.ELI(NmodesTproj);

    //exit(0);
    train.restart();
    // Construct the hyper reduced object
    EliBurgers hyperreduced(train);
    hyperreduced.Tmodes = train.EliObject->U.transpose();
    
    auto t1 = std::chrono::high_resolution_clock::now();
    hyperreduced.OnlineBurgers(NmodesTproj);
    auto t2 = std::chrono::high_resolution_clock::now();
    //auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1)

    auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    //std::cout << "time_magicspoints=  " << reduced.time_rom << std::endl;
    std::cout << "The Online  phase  duration =" << time_span.count()<< std::endl;
    /// Mapping to the whole fields

    for (int i = 0; i < hyperreduced.Tfield.size(); ++i)
    {
        ITHACAstream::exportSolution(hyperreduced.Tfield[i], name(i+1), "./ITHACAoutput/Online");
        volScalarField Te("Te", train.Tfield[i]-hyperreduced.Tfield[i]);
        //auto t = ITHACAutilities::L2Norm(train.Tfield[i]);
        //Te /=t;
        ITHACAstream::exportSolution(Te, name(i+1), "./ITHACAoutput/AbsError");
    }
   
    Eigen::MatrixXd errL2T = ITHACAutilities::errorL2Abs(train.Tfield, hyperreduced.Tfield);
    std::cout << "#####################" << std::endl;
    Eigen::MatrixXd errL2RelT = ITHACAutilities::errorL2Rel(train.Tfield, hyperreduced.Tfield);
    exit(0);
}
// Eigen::Map<Eigen::VectorXd > eigenMatrix(SEqn.source().data(), T.size(), 1 );
// S = Eigen::SparseMatrix<double>::Map(T.size(), T.size(), T.size(), outer.data(), inner.data(), SEqn.diag().data());
// S.makeCompressed();
// List<label> ListupperAddr = SEqn.lduAddr().upperAddr();
// List<label> ListlowerAddr = SEqn.lduAddr().lowerAddr();

// int l = SEqn.lduAddr().upperAddr().size();
// int ll = SEqn.lduAddr().lowerAddr().size()+1;
// //ListlowerAddr.append( ListlowerAddr[l-1]);
// std::vector<int> upperAddr(ListupperAddr.begin(), ListupperAddr.end()); ///OkK
// std::vector<int> lowerAddr(ListlowerAddr.begin(), ListlowerAddr.end()); ///OkK
// // //S = Eigen::SparseMatrix<double>::Map(T.size(),T.size(),l,lowerAddr.data(), upperAddr.data(), SEqn.upper().data());///Okk
// //S = Eigen::SparseMatrix<double>::Map(T.size(),T.size(),ll,upperAddr.data(),lowerAddr.data(),SEqn.lower().data());///Okk
// // //std::cout << S << std::endl;
// Info <<ListlowerAddr.size()  << endl;
// inplaceUniqueSort(ListlowerAddr);
// Info <<ListlowerAddr.size()  << endl;
//                  exit(0);






// void OnlineBurgers(int NmodesUproj,word folder = "./ITHACAoutput/Online/")
//     {
//         autoPtr<simpleControl> _simple;
//         // /// fvOptions
//         autoPtr<fv::options> _fvOptions;
//         fvMesh& mesh = problem->_mesh();
        
//         auto submesh = problem->EliObject->submesh;
//         submesh->setCellSubset(problem->EliObject->uniqueMagicPoints() );
//         submesh->subMesh().fvSchemes::readOpt() = mesh.fvSchemes::readOpt();
//         submesh->subMesh().fvSolution::readOpt() =mesh.fvSolution::readOpt();
        
//         submesh->subMesh().fvSchemes::read();
//         submesh->subMesh().fvSolution::read();
    
//          _simple = autoPtr<simpleControl>
//               (
//                   new simpleControl
//                   (
//                       submesh().subMesh()
//                   )
//               );
//          //Info << "####################################### " << endl;
//         _fvOptions = autoPtr<fv::options>(new fv::options(submesh().subMesh()));
//         simpleControl& simple = _simple();
//         fv::options& fvOptions = _fvOptions();
//         volScalarField& T = problem->_T();
//         //// sub field
//         volScalarField S = submesh().interpolate(T).ref();
//         S.rename("S");
//         surfaceScalarField sphi = submesh().interpolate(problem->_phi()).ref();
//         sphi.rename("sphi");
//         //Info << "####################################### " << endl;
//         //List<label> uniqueMagicPoints = problem->EliObject->magicPoints();
//         Eigen::VectorXd volumes = ITHACAutilities::getMassMatrixFV(T );

//         Eigen::SparseMatrix<double> Se;
//         Eigen::VectorXd se;
//         Eigen::MatrixXd X = problem->EliObject->U; 
//         Eigen::SparseMatrix<double> P = problem->EliObject->P;

//         Eigen::VectorXd subvol = P.transpose()*volumes;
//         Eigen::MatrixXd Q = P.transpose()*X;
//         //std::cout << "Q =\n" << Q << std::endl;
        
// //Eigen::MatrixXd V = X*(P.transpose()*X).inverse(); // Ok
// Eigen::MatrixXd V = X*(X.transpose()*P*P.transpose()*X).inverse()*X.transpose()*P; // Ok

//          counter = 1;
//         ////
//         nextWrite = startTime;
//         nextWrite += writeEvery;
//         while (simple.loop() )
//         {
//             Info << "Time = " << problem->_runTime().timeName() << endl;

//             while (simple.correctNonOrthogonal())
//             {
//                 /// The matrix
//                 fvScalarMatrix SEqn
//                 (
//                     fvm::ddt(S)
//                     +fvm::div(sphi,S)
//                     -fvm::laplacian(problem->_nu(),S)
//                 );
//                 SEqn.relax();

//                 problem->_fvOptions().constrain(SEqn);
//                 /// Field Conversion
//                 Foam2Eigen::fvMatrix2Eigen(SEqn, Se, se);
//                 //// Assemble the hyper-reduced system
//                 Eigen::MatrixXd A = Q.transpose()*Se*subvol.asDiagonal()*Q;
//                 //std::cout << "A =\n" << A << std::endl;
//                 Eigen::VectorXd a = Q.transpose()*subvol.asDiagonal()*se;
//                 //std::cout << "a =\n" << a << std::endl;
//                 ///Solve the Hyper-reduced problem.
//                 Eigen::VectorXd c = A.fullPivLu().solve(a);
//                 /// Reconstruct the eigen sub solution
//                 Eigen::VectorXd x = Q*c;
//                //std::cout << "######################################" << std::endl;
//                 ///Solution 2: Copy Eigen Data Using Iterators
//                 //std::copy(x.data(), x.data() + x.size(), S.ref().begin());
//                 S = Foam2Eigen::Eigen2field(S,x);
//                 // std::cout << V.rows() << std::endl;
//                 // std::cout << V.cols() << std::endl;

//                 // std::cout << x.rows() << std::endl;
//                 // std::cout << x.cols() << std::endl;
//                 //Info <<  S << endl;
//                 Eigen::VectorXd z = X*c;
//                 T = Foam2Eigen::Eigen2field(T, z );
//                 //std::copy(z.data(), z.data() + z.size(), problem->_T().ref().begin());
//                 //S.correctBoundaryConditions();
//                 T.correctBoundaryConditions();
                
//             }

//             problem->_runTime().printExecutionTime(Info);

//             if (checkWrite(problem->_runTime()))
//             {
//                 //ITHACAstream::exportSolution(T, name(counter), folder);
//                 counter++;
//                 Tfield.append(T.clone());
//                 nextWrite += writeEvery;
//             }
//         }

//     }
