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
    Hyper-reduction for scalar Transport
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
    //volVectorModes Lagrangian;
    List<Eigen::VectorXd> OnlineCoeffs;
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
        
        Eigen::SparseMatrix<double> S;
        Eigen::VectorXd se;

        int n = T.size();
        int s = problem->EliObject->P.cols();
        Eigen::MatrixXd B;
        B.setZero(s,NmodesUproj);
        Eigen::VectorXd b;
        b.setZero(s);

        List<label> uniqueMagicPoints = problem->EliObject->magicPoints();
        //Eigen::VectorXd volumes = ITHACAutilities::getMassMatrixFV(T);
        while (simple.loop() )
        {
            Info << "Time = " << problem->_runTime().timeName() << endl;

            while (simple.correctNonOrthogonal())
            {
                /// The matrix
                fvScalarMatrix SEqn
                (
                    fvm::ddt(T)
                    +fvm::div(phi,T)
                    -fvm::laplacian(problem->_nu(),T)
                );
                /// Field Conversion
                Foam2Eigen::fvMatrix2Eigen(SEqn, S, se);
                Eigen::VectorXd tempB(NmodesUproj);
                /// #pragma omp parallel for
                auto t1 = std::chrono::high_resolution_clock::now();
                //#pragma omp parallel for
                for (int k = 0; k < uniqueMagicPoints.size(); ++k)
                {
                    int magicPoint = uniqueMagicPoints[k];
                    Eigen::SparseVector<double> vec = S.row(magicPoint);

                    // Temporary storage for B(k, i)
                    tempB.setZero();
                    //for (Eigen::SparseVector<double>::InnerIterator it(vec); it; ++it)
                    for (Eigen::SparseMatrix<double>::InnerIterator it(S,magicPoint); it; ++it)
                    {
                        tempB += it.value() * problem->EliObject->U.row(it.index());
                    }

                    // Apply scaling by volume
                    double volume = T.mesh().V()[magicPoint];
                    B.row(k) = tempB * volume;
                    b(k) = volume * se(magicPoint);
                }
                // for (int k=0; k < uniqueMagicPoints.size(); ++k)
                // {
                //     Eigen::SparseVector<double> vec = S.row(uniqueMagicPoints[k]);
                //     for (Eigen::SparseVector<double>::InnerIterator it(vec); it; ++it)
                //     {
                //         Eigen::RowVectorXd tempB = Eigen::RowVectorXd::Zero(NmodesUproj);
                //         for (Eigen::SparseVector<double>::InnerIterator it(vec); it; ++it)
                //         {
                //             tempB += it.value() * problem->EliObject->U.row(it.index());
                //         }
                //         // for (int i = 0; i < NmodesUproj; ++i)
                //         // {
                //         //   /// Use OpenMP atomic if B is shared
                //         //   #pragma omp atomic
                //         //   B(k, i) += it.value()*problem->EliObject->U(it.index(), i);
                //         // }
                //     }
                //     B.row(k) *=T.mesh().V()[uniqueMagicPoints[k]];
                //     //B.row(k) *=volumes(uniqueMagicPoints[k]);
                //     //b(k) = volumes(uniqueMagicPoints[k]) * se(uniqueMagicPoints[k]);
                //     b(k) = T.mesh().V()[uniqueMagicPoints[k]]*se(uniqueMagicPoints[k]);
                // } 
                /// Projections: much costly
                //Eigen::VectorXd a = problem->EliObject->P*se;
                //Eigen::MatrixXd A = problem->EliObject->P*S*problem->EliObject->U;
                auto t2 = std::chrono::high_resolution_clock::now();
                auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
                time_rom +=time_span.count();
                //auto t1 = std::chrono::high_resolution_clock::now();
                /// solve the reduced problem
                //Eigen::VectorXd c = B.colPivHouseholderQr().solve(b);
                Eigen::VectorXd c = B.householderQr().solve(b);
                //Eigen::VectorXd c = B.ldlt().solve(b);
                // auto t2 = std::chrono::high_resolution_clock::now();
                Eigen::VectorXd z = problem->EliObject->U*c;
                /// Field Conversion
                ///Solution 2: Copy Eigen Data Using Iterators
                //auto t1 = std::chrono::high_resolution_clock::now();
                std::copy(z.data(), z.data() + z.size(), T.ref().begin());
                //T = Foam2Eigen::Eigen2field(problem->_T(), z );
                T.correctBoundaryConditions();
            }

            if (checkWrite(problem->_runTime()))
            {
                //ITHACAstream::exportSolution(T, name(counter), folder);
                counter++;
                Tfield.append(T.clone());
                //OnlineCoeffs.append(a);
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
    
    
};

int main(int argc, char* argv[])
{
    // Create the train object of the tutorial02 type
    tutorial24 train(argc, argv);
    //tutorial23 test(argc, argv);
    //std::clock_t startOff;
    //double durationOff;

    // Read some parameters from file
    ITHACAparameters* para = ITHACAparameters::getInstance(train._mesh(),train._runTime());
    int NmodesUout  = readInt(para->ITHACAdict->lookup("NmodesUout"));
    int NmodesUproj  = readInt(para->ITHACAdict->lookup("NmodesUproj"));

    //startOff= std::clock();
    auto t3 = std::chrono::high_resolution_clock::now();
    train.offlineSolve();
    auto t4 = std::chrono::high_resolution_clock::now();
    //durationOff = (std::clock() - startOff);
    auto time_span_off = std::chrono::duration_cast<std::chrono::duration<double>>(t4 - t3);
    // std::cout << "time_rom=  " << reduced.time_rom << std::endl;
    std::cout << "The Offline  phase  duration =" << time_span_off.count()<< std::endl;
    /// Construct the DEIM indices.
    train.ELI(NmodesUproj);

    List<label>bCells;
    forAll(train._U().boundaryField(), patchI)
    {
        const fvPatch &pp = train._U().boundaryField()[patchI].patch();
        forAll(pp, faceI)
        {
            label cellI = pp.faceCells()[faceI];
            bCells.append(cellI);
        }
    }
    /// Sorting and remove duplicates
    inplaceUniqueSort(bCells);

    //Info << bCells << endl;
    for (int i = 0; i <train.EliObject->magicPoints().size() ; ++i)
    {
        //std::cout << train.EliObject->magicPoints()[i] << std::endl;
        bCells.append(train.EliObject->magicPoints()[i]);
    }
    /// Sorting and remove duplicates
    inplaceUniqueSort(bCells);
    //train.EliObject->uniqueMagicPoints() = bCells;

    //Info << bCells << endl;
    //exit(0);

    // int dim  = 1;
    // label nCells = train._mesh().nCells();
    // /// Mask field to submesh
    // Eigen::SparseMatrix<double> field2submesh;
    // label M  = train.EliObject->uniqueMagicPoints().size();
    // field2submesh.resize(nCells*dim, M*dim);

    // for (unsigned int ith_subCell{0}; ith_subCell < M; ith_subCell++)
    // {
    //     for (unsigned int i= 0; i < dim; i++)
    //     {
    //       field2submesh.insert(train.EliObject->uniqueMagicPoints()[ith_subCell] + nCells*i, ith_subCell + i*M) = 1;
    //     }
    // }
    // field2submesh.makeCompressed();
    // //Eigen::SparseMatrix<double>  field2submesh = train.EliObject->P;
    // /// Get the volumes of cells.
    // Eigen::VectorXd volumes = ITHACAutilities::getMassMatrixFV(train._T() );
    //     //Eigen::MatrixXd diag =  volumes.asDiagonal();
    // Eigen::SparseMatrix<double> diag(train._T().size(),train._T().size());
   
    // for (int i = 0; i < train._T().size(); ++i)
    // {
    //     diag.insert(i,i) = volumes(i); 
    //     //diag.coeffRef(i,i) = volumes(i); 
    // }
    // diag.makeCompressed();  
    // //train.EliObject->MatrixOnline  = train.EliObject->U*((field2submesh.transpose()*train.EliObject->U).inverse());
    // //auto P = problem->EliObject->P.transpose()*diag;
    // //train.EliObject->P = field2submesh.transpose()*diag;
    // train.EliObject->P = field2submesh.transpose()*diag;
    /// Important matrices
    // Eigen::MatrixXd X = train.EliObject->U; 
    // Eigen::MatrixXd V = (X.transpose()*field2submesh*field2submesh.transpose()*X).inverse()*X.transpose()*field2submesh; // Ok
    // /// Define the full matrix to the full field
    // Eigen::MatrixXd Y = X*V; 
    // train.EliObject->MatrixOnline = Y; 
 
    //exit(0);
    // Eigen::VectorXd x(8);  x << 1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 7.7, 8.8;
    // Eigen::Array4i ind_vec(0,2,4,5);
    // Eigen::Array4d result = ind_vec.unaryExpr(x);

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
    //ITHACAstream::exportFields(reduced.Tfield[i], "./ITHACAoutput/Online", train._T().name());
    exit(0);
}

//SEqn.relax();
//problem->_fvOptions().constrain(SEqn);
//C.diagonal()(= SEqn.diag().data());
/// Eigen Mapping
//Eigen::Map<Eigen::MatrixXd> test(SEqn.source().data(), SEqn.source().size(), 1);//OKKK
//Eigen::Map<Eigen::MatrixXd> output(field.ref().data(), field.size(), 1);
//std::move(output);

// int p = 2;
// int q = 5;
// int r = 4;
// Eigen::VectorXf v2(5);
// Eigen::VectorXf v(3);
// Eigen::VectorXf v1(3);
// v2 << 100, 0.2, 30, 41, 55;
// v << 1, 4, 7;
// v1 << 2, 5, 8;

// //Eigen::VectorXf vNew = Eigen::VectorXf::Zero(30);
// Eigen::VectorXf vNew = Eigen::VectorXf::Zero(9);
// int l = 0;
// int j = 2;
// int step = 3;
// //vNew(Eigen::seq(p-1, p+q*r-1, q)) = v2;
// vNew(Eigen::seq(l, l+1+j*step-1, step)) = v;
// vNew(Eigen::seq(l+1, l+2+j*step-1, step)) = v1;

// std::cout << vNew.transpose() << std::endl;
// exit(0);
/// Matrix Conversion:
// Assuming foamMatrix.data() gives access to row/column indices and values
//Eigen::Map<Eigen::SparseMatrix<double>> eigenMatrix(SEqn.ref().data(), T.size(), T.size());
