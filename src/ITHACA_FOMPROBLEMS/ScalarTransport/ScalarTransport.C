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

//#include "laplacianProblem.H"
#include "ScalarTransport.H"

//#include "fvMeshSubsetter.H"  // Not fvMeshSubset (need two-step subsetting)
/// \brief Class where the tutorial number 2 is implemented.
/// \details It is a child of the laplacianProblem class and some of its
/// functions are overridden to be adapted to the specific case.


// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //
// Constructor
ScalarTransport::ScalarTransport() {}

ScalarTransport::ScalarTransport(int argc, char* argv[])
    :
    UnsteadyProblem()
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
    _simple = autoPtr<simpleControl>
              (
                  new simpleControl
                  (
                      mesh
                  )
              );
    simpleControl& simple = _simple();
#include "createFields.H"
#include "createFvOptions.H"
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
    para = ITHACAparameters::getInstance(mesh, runTime);
    offline = ITHACAutilities::check_off();
    podex = ITHACAutilities::check_pod();
    setTimes(runTime);
}

// void offlineSolve(word folder = "./ITHACAoutput/Offline/")
// {
//     if (offline)
//     {
//         ITHACAstream::read_fields(Tfield, T, "./ITHACAoutput/Offline/");
//         //ITHACAstream::readMiddleFields(Ufield, U, "./ITHACAoutput/Offline/");

//     }
//     else
//     {
//         truthSolve(folder);
//     }
// }

// void ELI(int NmodesT)
// {
//     fvMesh& mesh  =  const_cast<fvMesh&>(T.mesh());
//     EliObject = autoPtr<EmpLagranInter> (new EmpLagranInter(Tfield, NmodesT,"Tsol", T.name()));
//     EliObject->generateSubmesh(1, T.mesh(), T);
//     // EliObject->generateSubmesh(2, U.mesh(), U);
// }
void ScalarTransport::truthSolve(word folder)
{
    Time& runTime = _runTime();
    fvMesh& mesh = _mesh();
    volScalarField& T = _T();
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
            fvScalarMatrix TEqn
            (
                fvm::ddt(T)
                + fvm::div(phi,T)
                - fvm::laplacian(nu,T)
            );

            TEqn.relax();
            fvOptions.constrain(TEqn);
            TEqn.solve();
            fvOptions.correct(T);
        }

        //phi = linearInterpolate(U) & mesh.Sf();

        if (checkWrite(runTime))
        {
            ITHACAstream::exportSolution(T, name(counter), folder );
            counter++;
            Tfield.append(T.clone());
            nextWrite += writeEvery;
        }
    }
}

void ScalarTransport::restart()
{
    _T.clear();
    _U.clear();
    _phi.clear();
    _fvOptions.clear();
    _nu.clear();
    _transportProperties.clear();
    argList& args = _args();
    Time& runTime = _runTime();
    runTime.setTime(0, 1);
    Foam::fvMesh& mesh = _mesh();
#include "createFields.H"
#include "createFvOptions.H"
}

// class EliBurgers: public reductionProblem
// {
// public:

//     //// Data members
//     tutorial23* problem;
//     scalar nextWrite;
//     scalar startTime;
//     scalar writeEvery;
//     PtrList<volVectorField> Ufield;
//     //volVectorModes Lagrangian;
//     List<Eigen::VectorXd> OnlineCoeffs;
//     EliBurgers(){}
//     explicit EliBurgers(tutorial23& FoamProblem) : problem(&FoamProblem)
//     {
//         nextWrite  = 0;
//         startTime  = problem->startTime;
//         writeEvery = problem->writeEvery;

//         std::cout << "######## reduced Constructor calling ##########"<< std::endl;
//     }

//     ~EliBurgers(){}


//     void OnlineBurgers(int NmodesUproj,word folder = "./ITHACAoutput/Online/")
//     {
//         // Time& runTime = problem->_runTime();
//         // //fvMesh& mesh = problem->_mesh();
//         //volVectorField& U = problem->_U();
//         // surfaceScalarField& phi = problem->_phi();
//         // //simpleControl& simple = problem->_simple();
//         //dimensionedScalar& nu = problem->_nu();
//         /// simpleControl
//         autoPtr<simpleControl> _simple;
//         // /// fvOptions
//         autoPtr<fv::options> _fvOptions;
        
//         auto submesh = problem->EliObject->submesh;
//         submesh->setCellSubset(problem->EliObject->uniqueMagicPoints() );
//         submesh->subMesh().fvSchemes::readOpt() = problem->_mesh().fvSchemes::readOpt();
//         submesh->subMesh().fvSolution::readOpt() = problem->_mesh().fvSolution::readOpt();
        
//         submesh->subMesh().fvSchemes::read();
//         submesh->subMesh().fvSolution::read();
//         //submesh().subMesh().write();
//          _simple = autoPtr<simpleControl>
//               (
//                   new simpleControl
//                   (
//                       submesh().subMesh()
//                   )
//               );
        
//         _fvOptions = autoPtr<fv::options>(new fv::options(submesh().subMesh()));
//         simpleControl& simple = _simple();
//         fv::options& fvOptions = _fvOptions();
    
//         autoPtr<volVectorField> _S;
//         _S = autoPtr<volVectorField>
//          (
//              new volVectorField
//              (
//                  IOobject
//                  (
//                      "S",
//                      problem->_runTime().timeName(),
//                      submesh().subMesh(),
//                      IOobject::NO_READ,
//                      IOobject::NO_WRITE
//                  ),
//                  submesh().subMesh()
//              )
//          );
//          volVectorField& S = _S();
//   //Info << "############################################" << endl;
//         autoPtr<surfaceScalarField> _sphi;
//         _sphi = autoPtr<surfaceScalarField>
//        (
//            new surfaceScalarField
//            (
//                IOobject
//                (
//                    "sphi",
//                    problem->_runTime().timeName(),
//                    submesh().subMesh(),
//                    IOobject::READ_IF_PRESENT,
//                    IOobject::NO_WRITE
//                ),
//                linearInterpolate(S) & submesh().subMesh().Sf()
//            )
//          );
//          surfaceScalarField& sphi = _sphi();
//         //Info << submesh->subMesh().schemesDict().size() << endl;
//         // Info << submesh->subMesh().solutionDict() << endl;
//         // exit(0);
//         counter = 1;
//        //  int dim  = 3;
//        //  label nCells = problem->_mesh().nCells();
//        //  //exit(0);
//        //  /// Mask field to submesh
//        //  Eigen::SparseMatrix<double> field2submesh;
//        //  label M  = problem->EliObject->uniqueMagicPoints().size();
//        //  field2submesh.resize(nCells*dim, M*dim);
 
//        //  for (unsigned int ith_subCell{0} ; ith_subCell < M; ith_subCell++)
//        //  {
//        //      for (unsigned int i= 0; i < dim; i++)
//        //      {
//        //        field2submesh.insert(problem->EliObject->uniqueMagicPoints()[ith_subCell] + nCells*i, ith_subCell + i*M) = 1;
//        //      }
//        //  }
//        //  field2submesh.makeCompressed();
//        // /// Important matrices
//        //  Eigen::MatrixXd X = problem->EliObject->U; 
//        //  Eigen::MatrixXd V = (X.transpose()*field2submesh*field2submesh.transpose()*X).inverse()*X.transpose()*field2submesh; // Ok
//        //  /// Define the full matrix to the full field
//        //  Eigen::MatrixXd Y = X*V; 
//         ////
//         nextWrite = startTime;
//         nextWrite += writeEvery;
        
//         Eigen::SparseMatrix<double> Se;

//         Eigen::VectorXd se;
//         Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > solver; //0k
//         // performs a Cholesky factorization of A
//         //Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> chol(Se);    
//         //Eigen::SparseQR<Eigen::SparseMatrix<double> > solver;
//         //Eigen::SimplicialLDLT <Eigen::SparseMatrix<double> > solver;
//        //Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
//         // Eigen::VectorXd a;
//         //int r = problem->EliObject->MatrixOnline.cols();
//         //a.resize(r,1);
//         while (simple.loop() )
//         {
//             Info << "Time = " << problem->_runTime().timeName() << endl;

//             while (simple.correctNonOrthogonal())
//             {

//                 fvVectorMatrix SEqn
//                 (
//                     fvm::ddt(S) 
//                     + fvm::div(sphi,S)
//                     - fvm::laplacian(problem->_nu(),S)
//                 );
//                 //SEqn.relax();
//                 Foam2Eigen::fvMatrix2Eigen(SEqn, Se, se);
//                 /// Hyper-reduced system
//                 /// TODO: make sure to include the values of the volumes.
//                 Eigen::VectorXd a = solver.compute(Se).solve(se);
//                 //std::cout << "a.rows() = " << a.rows() << std::endl;
//                 //std::cout << "a.cols() = " << a.cols() << std::endl;
//                 //Eigen::VectorXd z = Y*a; 
//                 //exit(0);
//                 /// Solution in eigen format
//                 S = Foam2Eigen::Eigen2field(_S(), a);
//             }
//             sphi = linearInterpolate(S) & submesh->subMesh().Sf();

//             if (checkWrite(problem->_runTime()))
//             {
//                 //std::cout << "########## line 190" << std::endl;
//                 //ITHACAstream::exportSolution(S, name(counter), folder);
//                 counter++;
//                 Ufield.append(S.clone());
//                 //OnlineCoeffs.append(a);
//                 nextWrite += writeEvery;
//             }
//         }

//     }
   
//     bool checkWrite(Time& timeObject)
//     {
//         scalar diffnow = mag(nextWrite - atof(timeObject.timeName().c_str()));
//         scalar diffnext = mag(nextWrite - atof(timeObject.timeName().c_str()) -
//                               timeObject.deltaTValue());

//         if ( diffnow < diffnext)
//         {
//             return true;
//         }
//         else
//         {
//             return false;
//         }
//     }
    
// };




// int main(int argc, char* argv[])
// {
//     // Create the train object of the tutorial02 type
//     tutorial24 train(argc, argv);
//     //tutorial23 test(argc, argv);
//     std::clock_t startOff;
//     double durationOff;

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

//     // /// Construct the DEIM indices.
//     // train.ELI(NmodesUproj);

//     // int dim  = 3;
//     // label nCells = train._mesh().nCells();
//     // //exit(0);
//     // /// Mask field to submesh
//     // Eigen::SparseMatrix<double> field2submesh;
//     // label M  = train.EliObject->uniqueMagicPoints().size();
//     // field2submesh.resize(nCells*dim, M*dim);

//     // for (unsigned int ith_subCell{0} ; ith_subCell < M; ith_subCell++)
//     // {
//     //     for (unsigned int i= 0; i < dim; i++)
//     //     {
//     //       field2submesh.insert(train.EliObject->uniqueMagicPoints()[ith_subCell] + nCells*i, ith_subCell + i*M) = 1;
//     //     }
//     // }
//     // field2submesh.makeCompressed();
    
//     // /// Important matrices
//     // Eigen::MatrixXd X = train.EliObject->U; 
//     // Eigen::MatrixXd V = (X.transpose()*field2submesh*field2submesh.transpose()*X).inverse()*X.transpose()*field2submesh; // Ok
//     // /// Define the full matrix to the full field
//     // Eigen::MatrixXd Y = X*V; 
//     // //std::cout << "Y.rows() = " << Y.rows() << std::endl;
//     // //std::cout << "Y.cols() = " << Y.cols() << std::endl;
//     // train.EliObject->MatrixOnline = Y; 
   
//     // // PtrList<volVectorField> subfields;

//     // // for (int i = 0; i < train.Ufield.size(); ++i)
//     // // {

//     // //   auto  subfld = train.EliObject->submesh().interpolate(train.Ufield[i] ).ref();
//     // //   subfields.append(subfld);
       
//     // // }
//     // // ITHACAstream::exportFields(subfields, "./ITHACAoutput/SubFields", subfields[0].name() );  

//     // //exit(0);
//     // std::clock_t startOn;
//     // double durationOn;
//     // startOn = std::clock();

//     // train.restart();
//     // EliBurgers reduced(train);

//     // reduced.OnlineBurgers(NmodesUproj);

//     // durationOn = (std::clock() - startOn);
//     // std::cout << "The Online  phase  duration  is  equal  to " << durationOn << std::endl;
//     // //ITHACAstream::exportFields(reduced.Ufield, "./ITHACAoutput/Online", reduced.Ufield[0].name() );
//     // //PtrList<volVectorField> errflds;
//     // //exit(0);
//     // /// Mapping to the whole fields
//     // for (int i = 0; i < reduced.Ufield.size(); ++i)
//     // {
//     //     Eigen::VectorXd c = Foam2Eigen::field2Eigen(reduced.Ufield[i]);
//     //     Eigen::VectorXd z = Y*c; //reduced.OnlineCoeffs[i];
//     //     auto fld = Foam2Eigen::Eigen2field(train._U(), z);
//     //     fld.correctBoundaryConditions();
//     //     ITHACAstream::exportSolution(fld, name(i+1), "./ITHACAoutput/Online");
//     // }
//     //ITHACAstream::exportFields(reduced.Ufield, "./ITHACAoutput/Online", "S");
//     exit(0);
// }

// // ########################################################
    // label n = train.EliObject->uniqueMagicPoints().size();
    // // int size = train.EliObject->magicPoints().size();
    // Info << train.EliObject->magicPoints() << endl;
    // Info << train.EliObject->uniqueMagicPoints() << endl;
    // // //Info << train.EliObject->uniqueMagicPoints()[5] << endl;

    // // for (int i = 0; i <n ; ++i)
    // // {
    // //     //Info << "mesh.cells() = " << train._mesh().cells()[train.EliObject->uniqueMagicPoints()[5]] << endl;

    // //     auto cell  = train._mesh().cells()[train.EliObject->uniqueMagicPoints()[i]];
    // //     Info << cell << endl;
    // //     forAll (cell, faceI)
    // //     {
    // //         Info << cell[faceI] << endl;
    // //     }
    
    // // }
    // // // ##############################################################



//for(label i = 0; i < M; i++)
// {

//      CoeffMat[i][i] = UEqn.diag()[problem->EliObject->magicPoints()[i]];
//      CoeffMat.source()[i] = UEqn.source()[problem->EliObject->magicPoints()[i]];
// }
// Info << "Time2 = " << problem->_runTime().timeName() << nl << endl;
// // Assigning off-diagonal coefficients

// for (label k = 0; k < M; k++)
// {
//     //Info << "mesh.cells() = " << train._mesh().cells()[train.EliObject->uniqueMagicPoints()[5]] << endl;
//     auto cell  = mesh.cells()[problem->EliObject->magicPoints()[k]];
//     //Info << cell << endl;
//     forAll (cell, faceI)
//     {
//         label l = UEqn.lduAddr().lowerAddr()[faceI];
//         label u = UEqn.lduAddr().upperAddr()[faceI];
//         CoeffMat[l][u] = UEqn.upper()[faceI];
//         CoeffMat[u][l] = UEqn.lower()[faceI];
//         //Info << cell[faceI] << endl;
//     }

// }

// auto x  = CoeffMat.LUsolve();
// Info << x << endl;
// //volVectorField x(xx);

// std::cout << problem->EliObject->MatrixOnline.cols() << std::endl;
// std::cout << problem->EliObject->MatrixOnline.rows() << std::endl;



// Eigen::VectorXd y = Foam2Eigen::field2Eigen(x);

// std::cout << y.cols() << std::endl;
// std::cout << y.rows() << std::endl;

//  exit(0);

 
// Eigen::VectorXd z  = problem->EliObject->MatrixOnline*y;

//U = Foam2Eigen::Eigen2field(problem->_U(),z, true);
//Info<< CoeffMat.solve() << endl;
//Info<< x << endl;

// Assigning contribution from BC
// forAll(U.boundaryField(), patchI)
// {
//     const fvPatch &pp = U.boundaryField()[patchI].patch();
//     forAll(pp, faceI)
//     {
//         label cellI = pp.faceCells()[faceI];
//         //CoeffMat[cellI][cellI]   += UEqn.internalCoeffs()[patchI][faceI];
//         //CoeffMat.source()[cellI] +=UEqn.boundaryCoeffs()[patchI][faceI];
//     }
// }
// Info << CoeffMat << endl;
// Foam2Eigen::fvMatrix2Eigen(UEqn, Ae, be);
// auto Ar = P.transpose() * Ae * P;
// auto br = P.transpose() * be;
// auto Ap = Q.transpose()* Ar *Q;
// auto bp = Q.transpose() * br;

//std::cout <<  "P = " << P << std::endl;

//std::cout <<  "Ar = " << Q.transpose()* Ar *Q << std::endl;
//std::cout <<  "br = " << Q.transpose() * br << std::endl;
//auto a = Ap.colPivHouseholderQr().solve(bp);