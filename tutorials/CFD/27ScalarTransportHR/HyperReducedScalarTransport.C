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
      
        Eigen::SparseMatrix<double, Eigen::RowMajor> S;
        Eigen::VectorXd se;
        const auto& Modes = problem->EliObject->U.transpose();

        List<label> uniqueMagicPoints = problem->EliObject->magicPoints();
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
                std::tuple<Eigen::MatrixXd,Eigen::VectorXd> HRSys = HyperReducedSys(S,se,uniqueMagicPoints,Modes);
                /// HR system
                Eigen::VectorXd b = std::get<1>(HRSys);
                Eigen::MatrixXd B = std::get<0>(HRSys);
                //double cond = B.jacobiSvd().singularValues()(0) /
                 //B.jacobiSvd().singularValues().tail(1)(0);
           //std::cout << "cond number = " << cond << std::endl;
                //Eigen::VectorXd c = (B.transpose()*B).inverse()*B.transpose()*b; // LS problem using minimization
                Eigen::VectorXd c = B.partialPivLu().solve(b);
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
    // std::cout << "time_rom=  " << reduced.time_rom << std::endl;
    std::cout << "The Offline  phase  duration =" << time_span_off.count()<< std::endl;
    //ITHACAPOD::getModes(train.Tfield, train.Tmodes, train._T().name(),
    //                    train.podex, 0, 0, NmodesUproj);
    /// Construct the DEIM indices.
    train.ELI(NmodesUproj);
    
    // train.EliObject->P =  field2submesh;
    //List<label> test = {0,2791,10742,6295,7616,11418};//{0,10792,10824,10727,4,5,6,7,8,9}
    //Info << "test = "<< test <<endl;
    //train.EliObject->magicPoints() = {0,163,6338, 6238};// Deim_Lebesgue
    //Info << "magicPoints() = "<< train.EliObject->magicPoints() <<endl;

    //exit(0);
    train.restart();
    EliBurgers reduced(train);
    auto t1 = std::chrono::high_resolution_clock::now();
    reduced.OnlineBurgers(NmodesUproj);
    auto t2 = std::chrono::high_resolution_clock::now();
    //auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1)

    auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    std::cout << "time_magicspoints=  " << reduced.time_rom << std::endl;
    std::cout << "The Online  phase  duration =" << time_span.count()<< std::endl;
    /// Mapping to the whole fields

    for (int i = 0; i < reduced.Tfield.size(); ++i)
    {
        ITHACAstream::exportSolution(reduced.Tfield[i], name(i+1), "./ITHACAoutput/Online");
    }
   
    Eigen::MatrixXd errL2T = ITHACAutilities::errorL2Abs(train.Tfield, reduced.Tfield);
    std::cout << "#####################" << std::endl;
    Eigen::MatrixXd errL2RelT = ITHACAutilities::errorL2Rel(train.Tfield, reduced.Tfield);
    exit(0);
}

