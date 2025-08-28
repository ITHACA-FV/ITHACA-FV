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

\*---------------------------------------------------------------------------*/

/// \file
/// Source file of the sequentialIHTP class.


#include "sequentialIHTP.H"

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructors
sequentialIHTP::sequentialIHTP() {}

sequentialIHTP::sequentialIHTP(int argc, char* argv[])
    :
    DT("DT", dimensionSet(0, 2, -1, 0, 0, 0, 0), 1.0)
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
#include "createFields.H"
#include "createThermocouples.H"
    thermocouplesPos = TCpos;
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
#include "createRegularization.H"
    para = ITHACAparameters::getInstance(mesh, runTime);
    offline = ITHACAutilities::check_off();
    podex = ITHACAutilities::check_pod();
    startTime = runTime.startTime().value();
    deltaTime = runTime.deltaTValue();
    endTime = runTime.endTime().value();
    Ntimes = (endTime - startTime) / deltaTime;
    timeSteps.resize( Ntimes + 1 );
    forAll(timeSteps, timeI)
    {
        timeSteps[timeI] = startTime + (timeI) * deltaTime;
    }
    //Info << "debug: timeSteps = " << timeSteps << endl;
    //Info << "debug: startTime = " << startTime << endl;
    //Info << "debug: endTime = " << endTime << endl;
    //Info << "debug: deltaTime = " << deltaTime << endl;
}

// * * * * * * * * * * * * * * Full Order Methods * * * * * * * * * * * * * * //

void sequentialIHTP::setDiffusivity(scalar _diff)
{
    diffusivity = _diff;
}

void sequentialIHTP::setSpaceBasis(word type, scalar shapeParameter, label Npod)
{
    if (!thermocouplesRead)  // Kabir: It checks if the thermocouples have been read by evaluating the thermocouplesRead variable. If they haven't been read, it calls the readThermocouples function.
    {
        readThermocouples(); // Kabir:Defining positions of thermocouples
    }
    volScalarField& T = _T();
    fvMesh& mesh = _mesh();
    
    // ################### Kabir: Read coordinates of those selected thermocouples from ITHACAdict
    Time& runTime = _runTime();
    ITHACAparameters* para1 = new ITHACAparameters(mesh, runTime);                                // Create an instance of ITHACAparameters with the required arguments
    label NbasisInSpace = para1->ITHACAdict->lookupOrDefault<label>("sizeOfTheParameter", 0);     // label sizeOfParameter = 5; 
    // NbasisInSpace = 5; //thermocouplesNum;      // RBF <= measurements

    Eigen::VectorXd thermocoupleXValues(NbasisInSpace);
    Eigen::VectorXd thermocoupleZValues(NbasisInSpace);
    Eigen::VectorXd thermocoupleYValues(NbasisInSpace);
    
    List<scalar> thermocoupleXList = para1->ITHACAdict->subDict("thermocoupleCoordinates").lookup("XValues");
    List<scalar> thermocoupleZList = para1->ITHACAdict->subDict("thermocoupleCoordinates").lookup("ZValues");
    List<scalar> thermocoupleYList = para1->ITHACAdict->subDict("thermocoupleCoordinates").lookup("YValues");

    //Kabir: Convert the lists to Eigen::VectorXd
    for (label i = 0; i < NbasisInSpace; i++) 
    {
        thermocoupleXValues(i) = thermocoupleXList[i];
        thermocoupleZValues(i) = thermocoupleZList[i];
        thermocoupleYValues(i) = thermocoupleYList[i];
    }
    // ################### Kabir: Read coordinates of those selected thermocouples from ITHACAdict

    // Kabir: For defining the location of RBFs, we just project some of those selected measurement points on the boundary hotSide.
    // Kabir: RBFs are center at thermocouple projection.
    Info << "\nRadial Basis Functions are used." << endl;
    Info << "The center of each function is at the projection " << endl;
    Info << "of some selected thermocouple on the boundary hotSide." << endl; 
    Info << "RBFs must be less than measurements or at least are equal to measurements." << endl; 
    Info << "If RBFs = measurements it means that we have one basis function in front of each thermocouples.\n\n"; 
    // Kabir: Because RBFs are defined by a shape parameter and the center of the kernel(center of a radial basis function).

    Eigen::MatrixXd radius_kb(NbasisInSpace, T.boundaryField()[hotSide_ind].size());  // Kabir:
    heatFluxSpaceBasis.resize(NbasisInSpace);

    int thermocouplesCounter = 0;   // I changed the code, we do not need this counter. 
    int rbfCenterTimeI = 0;         // I changed the code, we do not need this counter. 

    // Kabir: It calculates the maximum X and Z coordinates of the face centers on the hot side boundary. No need to do that.
    scalar maxX = 1.0; // Foam::max(mesh.boundaryMesh()[hotSide_ind].faceCentres().component(Foam::vector::X));
    scalar maxZ = 1.0;  // Foam::max(mesh.boundaryMesh()[hotSide_ind].faceCentres().component(Foam::vector::Z));

    forAll(heatFluxSpaceBasis, funcI)  //5, Kabir: this loop will execute NbasisInSpace(5) times, once for each basis function.
    {
        // Kabir: Get the x and z coordinates of the thermocouple location
        scalar thermocoupleX =thermocoupleXValues(funcI);
        scalar thermocoupleZ =thermocoupleZValues(funcI);
        scalar thermocoupleY =thermocoupleYValues(funcI);

        // Kabir: Resize the heatFluxSpaceBasis[funcI] vector to match the size of hotSide boundary faces
        heatFluxSpaceBasis[funcI].resize(T.boundaryField()[hotSide_ind].size());
        forAll (T.boundaryField()[hotSide_ind], faceI) 
        {
            scalar faceX = mesh.boundaryMesh()[hotSide_ind].faceCentres()[faceI].x();                                           // Kabir: Get the x and z coordinates of the face center
            scalar faceZ = mesh.boundaryMesh()[hotSide_ind].faceCentres()[faceI].z();

            // Kabir: R is the distance value for the kernel. Calculate the distance (radius) between those selected thermocouples and face center
            scalar radius = Foam::sqrt((faceX - thermocoupleX) * (faceX - thermocoupleX) / maxX / maxX + (faceZ - thermocoupleZ) * (faceZ - thermocoupleZ) / maxZ / maxZ);
            // Kabir: The normalization factors maxX and maxZ are used to ensure that the components (differences in x and z coordinates) are scaled properly be.fore calculating the radius. 
            // Kabir: I think the normalization factor should not be used and it must utilize just Euclidean distance itself. So we put maxX and maxZ equalt to 
            // Kabir: Store the calculated radius for later use
            radius_kb(funcI,faceI)=radius;
            
            // Kabir Calculate the value of the kernel function and assign it to heatFluxSpaceBasis[funcI][faceI]

            heatFluxSpaceBasis[funcI][faceI] = Foam::sqrt(1 + (shapeParameter * radius) * (shapeParameter * radius));            // Kabir: Multiquadric kernel 
            // heatFluxSpaceBasis[funcI][faceI] = Foam::exp(-1.0 * (shapeParameter * shapeParameter) * (radius * radius));       // Kabir: Gaussian kernel
            // heatFluxSpaceBasis[funcI][faceI] = 1.0 / (1.0 + (shapeParameter * radius) * (shapeParameter * radius));           // Kabir: Inverse Quadratic Kernel
            // heatFluxSpaceBasis[funcI][faceI] = 1.0 / Foam::sqrt(1.0 + (shapeParameter * radius) * (shapeParameter * radius)); // Kabir: Inverse Multiquadric Kernel 

        }
        // Kabir: Increment the thermocouplesCounter to move to the next thermocouple location. I changed the code, we do not need this counter. 
        thermocouplesCounter++;
    }
    /// ###################Kabir: Export radius and the location of those selected thermocouples
    cnpy::save(radius_kb, "radius_kb.npy");
    cnpy::save(thermocoupleXValues, "thermocoupleXValues.npy");
    cnpy::save(thermocoupleZValues, "thermocoupleZValues.npy");
    cnpy::save(thermocoupleYValues, "thermocoupleYValues.npy");
    /// ###################Kabir: Export radius and the location of those selected thermocouples

    /// ###################Kabir: Export the heatFluxSpaceBasis data in order to plot the reconstructed heat flux, ITHACAoutput/projection/HeatFluxSpaceRBF
    std::string heat_flux_space_basis="heat_flux_space_basis";

    std::string folderNamePrRb = "HeatFluxSpaceRBF";
    std::string folderPathPrRb = "ITHACAoutput/projection/" + folderNamePrRb;

               // Convert heatFluxSpaceBasis which is Foam::List<Foam::List<double>> object to heatFluxSpaceBasisMatrix which is Eigen::Matrix object before calling the exportMatrix function by using Eigen::Map constructor 
    Eigen::MatrixXd heatFluxSpaceBasisMatrix(heatFluxSpaceBasis.size(), heatFluxSpaceBasis[0].size());
    for (label i = 0; i < heatFluxSpaceBasis.size(); ++i)
        heatFluxSpaceBasisMatrix.row(i) = Eigen::Map<Eigen::VectorXd>(heatFluxSpaceBasis[i].begin(), heatFluxSpaceBasis[i].size());

    ITHACAstream::exportMatrix(heatFluxSpaceBasisMatrix, heat_flux_space_basis, "eigen", folderPathPrRb);

    // ###################Kabir: Export the heatFluxSpaceBasis data in order to plot the reconstructed heat flux, ITHACAoutput/projection/HeatFluxSpaceRBF


    if (type == "pod")
    {
        Info << "POD basis are not yet implemented, exiting" << endl;
        exit(10);
        //List<scalar> massVector(T.boundaryField()[hotSide_ind].size());
        //forAll (T.boundaryField()[hotSide_ind], faceI)
        //{
        //    massVector[faceI] = mesh.boundary()[hotSide_ind].magSf()[faceI];
        //}
        //List<List<scalar>> tempBasis;
        //word debugFolder = "./ITHACAoutput/debugParameterizedBasis/";
        //ITHACAPOD::getModesSVD(heatFluxSpaceBasis, tempBasis, "PODbase", false, false, false, Npod, true);

        //forAll(heatFluxSpaceBasis, baseI)
        //{
        //    volScalarField base = list2Field(heatFluxSpaceBasis[baseI], 0.0);
        //    ITHACAstream::exportSolution(base,
        //                                 std::to_string(1),
        //                                 debugFolder,
        //                                 "RBFbase" + std::to_string(baseI + 1));
        //}
        //forAll(tempBasis, baseI)
        //{
        //    volScalarField base = list2Field(tempBasis[baseI], 0.0);
        //    ITHACAstream::exportSolution(base,
        //                                 std::to_string(1),
        //                                 debugFolder,
        //                                 "PODbase" + std::to_string(baseI + 1));
        //}
        //heatFluxSpaceBasis = tempBasis;
    }
}

void sequentialIHTP::set_gParametrized(word spaceBaseFuncType,
        scalar shapeParameter_space)
{
    volScalarField& T = _T();
    setSpaceBasis(spaceBaseFuncType, shapeParameter_space);

    Info << "Using CONSTANT basis in time" << endl;
    NbasisInTime = 1;
    NsamplesWindow = 1;
    Nbasis = NbasisInTime * NbasisInSpace;
    gBaseFunctions.resize(Nbasis);
    NtimestepsInSequence = NtimeStepsBetweenSamples * NbasisInTime + 1;
    forAll(gBaseFunctions, baseI)
    {
        gBaseFunctions[baseI].resize(NtimestepsInSequence);
        for(int timeI = 0; timeI < NtimestepsInSequence; timeI++)
        {
            gBaseFunctions[baseI][timeI] = heatFluxSpaceBasis[baseI];
        }
    }

    g.resize(timeSteps.size());
    gWeights.resize(Nbasis);
    forAll (gWeights, weigthI)
    {
        gWeights[weigthI] = 0;
    }
    forAll(timeSteps, timeI)
    {
        g[timeI].resize(T.boundaryField()[hotSide_ind].size(), 0.0);
    }
}

List<List<scalar>> sequentialIHTP::interpolateWeights(List<scalar> Wold, List<scalar> Wnew)
{
    M_Assert(Wold.size() == Wnew.size(), "Input weights vectors must have the same size");

    double t0 = 0;
    double t1 = NtimeStepsBetweenSamples * deltaTime;
    List<List<scalar>> Wout;
    Wout.resize(Wold.size());
    forAll (Wold, wI)
    {
        Wout[wI].resize(NtimeStepsBetweenSamples + 1);
        for(int timeI = 0; timeI < NtimeStepsBetweenSamples + 1; timeI++)
        {
            double time = (timeI + 1) * deltaTime;
            double a = Wold[wI] - (Wnew[wI] - Wold[wI]) / (t1 - t0) * t0;
            double b = (Wnew[wI] - Wold[wI]) / (t1 - t0);
            Wout[wI][timeI] = a + b * time;
        }
    }
    return Wout;
}

void sequentialIHTP::update_gParametrized(List<scalar> weights)
{
    M_Assert(weights.size() == Nbasis,
             "weigths size different from basis functions size");
    volScalarField& T = _T();
    label firstTimeI = timeSampleI * NtimeStepsBetweenSamples;
    List<List<scalar>> interpolatedWeights;

    if(interpolationFlag && !offlineFlag)
    {
        if(timeSampleI > 0)
        {
            interpolatedWeights = interpolateWeights(gWeightsOld, weights);
        }
    }

    int lastTimeStep = firstTimeI + NtimeStepsBetweenSamples + 1;

    label shortTime = 0;
    for(int timeI = firstTimeI; timeI < lastTimeStep; timeI++)
    {
        forAll (T.boundaryField()[hotSide_ind], faceI)
        {
            g[timeI][faceI] = 0.0;
            forAll (weights, weightI)
            {
                if(interpolationFlag && !offlineFlag)
                {
                    if(timeSampleI > 0)
                    {
                        g[timeI][faceI] += interpolatedWeights[weightI][shortTime] * 
                            gBaseFunctions[weightI][shortTime][faceI];
                    }
                    else
                    {
                        g[timeI][faceI] += weights[weightI] * 
                            gBaseFunctions[weightI][shortTime][faceI];
                    }
                }
                else
                {
                    g[timeI][faceI] += weights[weightI] * 
                        gBaseFunctions[weightI][shortTime][faceI];
                }
            }
        }
	shortTime++;
    }
}

volScalarField sequentialIHTP::list2Field(List<scalar> list,
        scalar innerField)
{
    volScalarField& T = _T();
    fvMesh& mesh = _mesh();
    volScalarField field(T);
    ITHACAutilities::assignIF(field, innerField);
    //Access the mesh information for the boundary
    const polyPatch& cPatch = mesh.boundaryMesh()[hotSide_ind];
    //List of cells close to a boundary
    const labelUList& faceCells = cPatch.faceCells();
    forAll(cPatch, faceI)
    {
        //id of the owner cell having the face
        label faceOwner = faceCells[faceI] ;
        field[faceOwner] = list[faceI];
    }
    return field;
}

void sequentialIHTP::parameterizedBCoffline(bool force)
{
    fvMesh& mesh = _mesh();
    Tbasis.resize(0);
    T0field.resize(0);
    M_Assert(diffusivity > 0.0, "Call setDiffusivity to set up the diffusivity");

    offlineTimestepsSize = NtimeStepsBetweenSamples;
    offlineEndTime = NtimeStepsBetweenSamples * NbasisInTime * deltaTime;

    if (ITHACAutilities::check_file(folderOffline + "/Theta_mat.txt") && force == 0)
    {
        Info << "\nOffline already computed." << endl;
        Info << "Check that the basis used for the parameterized BC are correct (RBF, POD, etc.)"
             << endl;
        Theta = ITHACAstream::readMatrix(folderOffline + "Theta_mat.txt");
        addSol = ITHACAstream::readMatrix(folderOffline + "addSol_mat.txt");
        ITHACAstream::read_fields(Tad_time, "Tad", folderOffline);

        for (label baseI = 0; baseI < Theta.cols(); baseI++)
        {
            Ttime.resize(0);
            ITHACAstream::read_fields(Ttime, "T" + std::to_string(baseI + 1),
                                      folderOffline);
            Tbasis.append(Ttime.clone());
        }
    }
    else
    {
        Info << "\nComputing offline" << endl;
        Theta.resize(Nbasis, gWeights.size());
	offlineFlag = 1;
        Info << "Theta size = " << Theta.rows() << ", " << Theta.cols() << endl;
        solveAdditional();

        for (label baseI = 0; baseI < Theta.cols(); baseI++)
        {
            Info << "\n--------------------------------------\n" << endl;
            Info << "Base " << baseI + 1 << " of " << Theta.cols() << endl;
            Info << "\n--------------------------------------\n" << endl;
            restart();
            Ttime.resize(0);
            gWeights = Foam::zero();
            gWeights[baseI] =  1;
	    timeSampleI = 0;
            update_gParametrized(gWeights);
            solveDirect();

            for(int timeI = 0; timeI < offlineTimestepsSize; timeI++)
            {
                volScalarField& T = Ttime[timeI];
                /// Saving basis
                volScalarField gParametrizedField = list2Field(g[timeI]);
                ITHACAstream::exportSolution(gParametrizedField,
                                             std::to_string(timeSteps[timeI + 1]),
                                             folderOffline,
                                             "g" + std::to_string(baseI + 1));
                ITHACAstream::exportSolution(T, std::to_string(timeSteps[timeI + 1]),
                                             folderOffline,
                                             "T" + std::to_string(baseI + 1));
            }
            Tbasis.append(Ttime.clone());
            Tcomp = fieldValueAtThermocouples(Ttime);
            M_Assert(Tcomp.size() == addSol.size(), "Something wrong in reading values at the observations points");
            for(int i = 0; i < Tcomp.size(); i++)
            {
                Theta(i, baseI) = Tcomp(i) + addSol(i);
            }
        }

        ITHACAstream::exportMatrix(Theta, "Theta", "eigen", folderOffline);
        ITHACAstream::exportMatrix(addSol, "addSol", "eigen", folderOffline);
    }

    Eigen::MatrixXd A = Theta;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A,
                                          Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::MatrixXd singularValues = svd.singularValues();
    double conditionNumber = singularValues.maxCoeff() / singularValues.minCoeff();
    Info << "Theta Condition number = " << conditionNumber << endl;
    ITHACAstream::exportMatrix(singularValues, "ThetaSingularValues", "eigen",
                               folderOffline);
    Eigen::MatrixXd ThetaTI = Theta * Eigen::MatrixXd::Identity(Theta.rows(), Theta.cols());
    ITHACAstream::exportMatrix(ThetaTI, "ThetaTI", "eigen",
            folderOffline);
    offlineFlag = 0;
    Info << "\nOffline ENDED" << endl;
}

void sequentialIHTP::reconstrucT(word outputFolder)
{
    Info << "Reconstructing field T" << endl;
    Ttime.resize(0);
    Info << "\nExporting solution in the time domain (" << timeSteps[NtimeStepsBetweenSamples * timeSampleI] << ", " << timeSteps[NtimeStepsBetweenSamples + NtimeStepsBetweenSamples * timeSampleI] << "]\n" << endl;
    restart();
    if(timeSampleI == 0)
    {
        ITHACAstream::exportSolution(T0_time[0], std::to_string(timeSteps[0]),
                                     outputFolder,
                                     "Treconstructed");
        volScalarField gParametrizedField = list2Field(g[0]);
        ITHACAstream::exportSolution(gParametrizedField,
                                     std::to_string(timeSteps[0]),
                                     outputFolder,
                                     "gReconstructed");
    }
    for(int timeI = 0; timeI < NtimeStepsBetweenSamples; timeI++)
    {
        volScalarField T=_T.clone();
        ITHACAutilities::assignIF(T, homogeneousBC);
        forAll(Tbasis, baseI)
        {
            T += gWeights[baseI] * (Tbasis[baseI][timeI] + Tad_time[timeI]);
        }
        T += - Tad_time[timeI] + T0_time[timeI];

        label realTimeStep = timeI + NtimeStepsBetweenSamples * timeSampleI + 1;
        ITHACAstream::exportSolution(T, std::to_string(timeSteps[realTimeStep]),
                                     outputFolder,
                                     "Treconstructed");
        volScalarField gParametrizedField = list2Field(g[realTimeStep - 1]);
        ITHACAstream::exportSolution(gParametrizedField,
                                     std::to_string(timeSteps[realTimeStep]),
                                     outputFolder,
                                     "gReconstructed");
        Ttime.append(T.clone());
    }
}

Eigen::VectorXd sequentialIHTP::reconstrucT(Eigen::VectorXi cells)
{
    Eigen::VectorXd Tout(cells.size());
    int timeI = NtimeStepsBetweenSamples - 1;
    for(int cellI = 0; cellI < cells.size(); cellI++)
    {
        forAll(Tbasis, baseI)
        {
            Tout(cellI) += gWeights[baseI] * (Tbasis[baseI][timeI].internalField()[cellI] 
                    + Tad_time[timeI].internalField()[cellI]);
        }
        Tout(cellI) += - Tad_time[timeI].internalField()[cellI] + T0_time[timeI].internalField()[cellI];
    }

    return Tout;
}

void sequentialIHTP::parameterizedBC(word outputFolder, volScalarField initialField)
{
    Info << endl << "Using quasilinearity of direct problem ::" << endl;
    Info << "Using " << linSys_solver << " to solve the linear system" << endl;
    timeSampleI = 0;

    List<Eigen::MatrixXd> linSys;
    linSys.resize(2);
    linSys[0] = Theta.transpose() * Theta;
    while(timeSampleI < timeSamplesNum)
    {
        Info << "\nTime sample " << timeSampleI + 1 << endl;
        auto t_start = std::chrono::high_resolution_clock::now();


        if(timeSampleI > 0)
        {
             reconstrucT("./ITHACAoutput/debugReconstrucT/");
             initialField = Ttime[NtimeStepsBetweenSamples - 1];
        }
	solveT0(initialField);

	TmeasShort = Tmeas.segment(thermocouplesNum * timeSampleI, thermocouplesNum * NsamplesWindow); // Kabir: thermocouplesNum comes from createThermocouples.H inside /u/k/kbakhsha/ITHACA-FV-KF/src/ITHACA_FOMPROBLEMS/sequentialIHTP
        linSys[1] = Theta.transpose() * (TmeasShort + addSol - T0_vector);
        Eigen::VectorXd weigths;

        if (linSys_solver == "fullPivLU")
        {
            weigths = linSys[0].fullPivLu().solve(linSys[1]);
        }
        else if (linSys_solver == "jacobiSvd")
        {
            Eigen::JacobiSVD<Eigen::MatrixXd> svd(linSys[0],
                                                  Eigen::ComputeThinU | Eigen::ComputeThinV);
            weigths = svd.solve(linSys[1]);
        }
        else if (linSys_solver == "householderQr")
        {
            weigths = linSys[0].householderQr().solve(linSys[1]);
        }
        else if (linSys_solver == "ldlt")
        {
            weigths = linSys[0].ldlt().solve(linSys[1]);
        }
        else if (linSys_solver == "inverse")
        {
            weigths = linSys[0].inverse() * linSys[1];
        }
        else if (linSys_solver == "TSVD")
        {
            Info << "Using TSVD" << endl;
            weigths = ITHACAregularization::TSVD(Theta, TmeasShort + addSol - T0_vector, TSVD_filter);
        }
        else if (linSys_solver == "Tikhonov")
        {
            linSys[0] = Theta;
            linSys[1] = (TmeasShort + addSol - T0_vector);
            Eigen::JacobiSVD<Eigen::MatrixXd> svd(linSys[0],
                                                Eigen::ComputeThinU | Eigen::ComputeThinV);

            weigths = ITHACAregularization::Tikhonov(Theta, linSys[1], Tikhonov_filter);
        }
        else
        {
            Info << "Select a linear system solver in this list:" << endl
                 << "fullPivLU, jacobiSvd, householderQr, ldlt, TSVD, Tikhonov, conjugateGradient" << endl;
            exit(1);
        }
        gWeightsOld = gWeights;
        gWeights.resize(weigths.size());
        forAll(gWeights, weightI)
        {
            gWeights[weightI] = weigths(weightI);
        }
        Info << "Weights = \n" << gWeights << endl;
        update_gParametrized(gWeights);
        label verbose = 0;
        parameterizedBC_postProcess(linSys, weigths, outputFolder, verbose);
	timeSampleI++;
        auto t_end = std::chrono::high_resolution_clock::now();
        double elapsed_time_ms = std::chrono::duration<double, std::milli>(t_end-t_start).count();
        Info << "CPU time = " << elapsed_time_ms << " milliseconds" << endl << endl;
    }
    ITHACAstream::exportMatrix(Jlist, "costFunction", "eigen", outputFolder);
    Info << "End" << endl;
    Info << endl;
}

void sequentialIHTP::set_valueFraction()
{
    fvMesh& mesh = _mesh();
    valueFraction.resize(mesh.boundaryMesh()["coldSide"].size());
    homogeneousBCcoldSide.resize(mesh.boundaryMesh()["coldSide"].size());
    Eigen::VectorXd faceCellDist =
        ITHACAutilities::boudaryFaceToCellDistance(mesh, coldSide_ind);
    forAll (valueFraction, faceI)
    {
        valueFraction[faceI] = 1 / (1 + (thermalCond / HTC / faceCellDist(faceI)));
        homogeneousBCcoldSide[faceI] =  0;
    }
    refGrad = homogeneousBCcoldSide;
}

void sequentialIHTP::assignDirectBC(label timeI)
{
    fvMesh& mesh = _mesh();
    volScalarField& T = _T();
    set_valueFraction();
    forAll(mesh.boundaryMesh(), patchI)
    {
        if (patchI == mesh.boundaryMesh().findPatchID("coldSide"))
        {
            ITHACAutilities::assignMixedBC(T, patchI, Tf, refGrad, valueFraction);
        }
        else if (patchI == mesh.boundaryMesh().findPatchID("hotSide"))
        {
            ITHACAutilities::assignBC(T, patchI, - g[timeI] / thermalCond);
        }
        else
        {
            ITHACAutilities::assignBC(T, patchI, homogeneousBC);
        }
    }
}

void sequentialIHTP::solveT0(volScalarField initialField)
{
    Info << "\nSolving FULL T0 problem" << endl;
    restartOffline();
    fvMesh& mesh = _mesh();
    simpleControl& simple = _simple();
    fv::options& fvOptions(_fvOptions());
    volScalarField T0= _T.clone();
    Foam::Time& runTime = _runTime();
    set_valueFraction();
    List<scalar> RobinBC = Tf * 0.0;
    word outputFolder = "./ITHACAoutput/debugT0/";

    T0 = initialField;

    T0field.append(T0.clone());
    T0_time.resize(0);
    label timeI = 0;
    forAll(mesh.boundaryMesh(), patchI)
    {
        if (patchI == mesh.boundaryMesh().findPatchID("coldSide"))
        {
            ITHACAutilities::assignMixedBC(T0, patchI, RobinBC, refGrad,
                                           valueFraction);
        }
        else
        {
            ITHACAutilities::assignBC(T0, patchI, homogeneousBC);
        }
    }

    while (runTime.loop())
    {
        Info << "Time = " << runTime.timeName() << nl << endl;
        timeI++;

        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix TEqn
            (
                fvm::ddt(T0) - fvm::laplacian(DT * diffusivity, T0)
            );
            fvOptions.constrain(TEqn);
            TEqn.solve();
            fvOptions.correct(T0);
        }

        T0_time.append(T0.clone());
        T0field.append(T0.clone());
        ITHACAstream::exportSolution(T0, std::to_string(
                    timeSteps[samplingSteps[timeSampleI] - NtimeStepsBetweenSamples + 
                    timeI]), outputFolder, "T0");
        runTime.printExecutionTime(Info);
        runTime.write();
    }

    T0_vector = fieldValueAtThermocouples(T0_time);
    Info << "SolveT0 ENDED\n" << endl;
}

void sequentialIHTP::getT0modes()
{
    word outputFolder = "./ITHACAoutput/modes/";
    Info << "Computing " << NmodesT0 << " T0 modes" << endl;

    ITHACAPOD::getModes(T0field, T0modes, "T0",
                        0, 0, 0,
                        NmodesT0);
    PtrList<volScalarField> modes = T0modes.toPtrList();
    forAll(modes, modeI)
    {
        ITHACAstream::exportSolution(modes[modeI], std::to_string(modeI + 1),
                outputFolder, "T0modes");
    }
    Info << "T0 modes COMPUTED\n" << endl;
}

void sequentialIHTP::projectT0()
{
    Info << "\n*****************************************************" << endl;
    Info << "Computing projection matrices" << endl;
    fvMesh& mesh = _mesh();
    simpleControl& simple = _simple();
    fv::options& fvOptions(_fvOptions());
    volScalarField T0=_T.clone();
    Foam::Time& runTime = _runTime();
    set_valueFraction();
    List<scalar> RobinBC = Tf * 0.0;

    forAll(mesh.boundaryMesh(), patchI)
    {
        if (patchI == mesh.boundaryMesh().findPatchID("coldSide"))
        {
            ITHACAutilities::assignMixedBC(T0, patchI, RobinBC, refGrad,
                                           valueFraction);
        }
        else
        {
            ITHACAutilities::assignBC(T0, patchI, homogeneousBC);
        }
    }

    Eigen::SparseMatrix<double> T0implicitMatrix;
    Eigen::SparseMatrix<double> T0explicitMatrix;

    ITHACAutilities::assignIF(T0, 1.0);
    fvScalarMatrix Teq(fvm::ddt(T0) - fvm::laplacian(DT * diffusivity, T0));
    Eigen::VectorXd b;
    Foam2Eigen::fvMatrix2Eigen(Teq, T0implicitMatrix, b);

    T0modes.toEigen();
    T0implicitMatrix_red = T0modes.EigenModes[0].transpose() 
        * T0implicitMatrix * T0modes.EigenModes[0];
    T0explicitMatrix_red = T0modes.EigenModes[0].transpose() 
        * b.asDiagonal() * T0modes.EigenModes[0];
    Info << "Projection matrices COMPUTED\n" << endl;
}

void sequentialIHTP::projectDirectOntoT0()
{
    /// Creation of the matrices to project direct solution at the last timestep
    /// onto the T0 reduced space
    Info << "Computing direct problem projection matrices" << endl;
    int internalFieldSize = Tbasis[0][0].internalField().size();
    Eigen::MatrixXd Tbasis_Eigen(internalFieldSize, Nbasis);
    Eigen::VectorXd Tad_Eigen = 
        Foam2Eigen::field2Eigen(Tad_time[NtimeStepsBetweenSamples - 1]);

    PtrList<volScalarField> Tbasis_lastTime;
    M_Assert(Tbasis[0].size() == NtimeStepsBetweenSamples, 
            "The basis for the direct problem have wrong dimention in time");

    forAll(Tbasis, baseI)
    {
        volScalarField temp = Tbasis[baseI][NtimeStepsBetweenSamples - 1] 
            + Tad_time[NtimeStepsBetweenSamples - 1];
        Tbasis_lastTime.append(temp.clone());
    }

    Tbasis_projectionMat = T0modes.project(Tbasis_lastTime); 
    Tad_projected = T0modes.project(Tad_time[NtimeStepsBetweenSamples - 1]);
}

void sequentialIHTP::pointProjectionOffline()
{
    int Ncells = magicPoints.size();
    M_Assert(Ncells > 0, "Set the number of magic points");

    int lastTimestepID = NtimeStepsBetweenSamples - 1;

    pointsReconstructMatrix.resize(Ncells, NmodesT0);
    for(int cellI = 0; cellI < Ncells; cellI++)
    {
        pointsReconstructMatrix.row(cellI) = 
            T0modes.EigenModes[0].row(magicPoints[cellI]);
    }
    pointTbasis_reconstructionMat = pointsReconstructMatrix * Tbasis_projectionMat;
    pointTad_reconstructed = pointsReconstructMatrix * Tad_projected;
}

void sequentialIHTP::projectionErrorOffline()
{
    /// I compute the L2 norm of the Tbasis and Tad perpendicular to the projection
    Info << "Computing the offline part for the projection error" << endl;

    //TODO You can consider all modes when computing the error and then 
    //choose the right ammount of modes based on the behaviour of the error
    int lastTimestepID = NtimeStepsBetweenSamples - 1;
    projectionErrorTbasis.resize(0);
    word outputFolder = "./ITHACAoutput/projectionError";


    forAll(Tbasis, baseI)
    {
        volScalarField base = Tbasis[baseI][lastTimestepID];
        volScalarField baseProj(base);
        baseProj = T0modes.projectSnapshot(base, NmodesT0);//, "L2");
        volScalarField temp = base - baseProj;
        projectionErrorTbasis.append(temp.clone());
        ITHACAstream::exportSolution(projectionErrorTbasis[baseI], std::to_string(baseI + 1),
                outputFolder, "projectionErrorTbasis");
    }
    projectionErrorTad.resize(0);
    volScalarField TadProj(Tbasis[0][lastTimestepID]);
    TadProj = T0modes.projectSnapshot(Tad_time[lastTimestepID], NmodesT0);//, "L2");
    volScalarField temp(Tad_time[lastTimestepID] - TadProj);
    projectionErrorTad.append(temp.clone());
    ITHACAstream::exportSolution(projectionErrorTad[0], std::to_string(1),
            outputFolder, "projectionErrorTad");

    //forAll(Tbasis, baseI)
    //{
    //    volScalarField base = Tbasis[baseI][lastTimestepID];
    //    volScalarField baseProj(base);
    //    T0modes.projectSnapshot(base, baseProj, NmodesT0, "L2");
    //    volScalarField basePerp = base - baseProj;
    //    projectionErrorTbasis(baseI) = ITHACAutilities::L2Norm(basePerp);
    //    Info << "Projection error Tbase["<< baseI << "] = " << projectionErrorTbasis(baseI) << endl;
    //}
    //volScalarField TadProj(Tad_time[lastTimestepID]);
    //T0modes.projectSnapshot(Tad_time[lastTimestepID], TadProj, NmodesT0, "L2");
    //volScalarField TadPerp = Tad_time[lastTimestepID] - TadProj;
    //projectionErrorTad = ITHACAutilities::L2Norm(TadPerp);
    //Info << "Projection error Tad = " << projectionErrorTad << endl;
}

void sequentialIHTP::T0offline(int NmagicPoints)
{
    getT0modes();
    findMagicPoints(NmagicPoints);
    projectT0();
    projectDirectOntoT0();
    projectionErrorOffline();
    pointProjectionOffline();
}

void sequentialIHTP::solveAdditional()
{
    Info << "Solving additional problem" << endl;
    restartOffline();
    fvMesh& mesh = _mesh();
    simpleControl& simple = _simple();
    fv::options& fvOptions(_fvOptions());
    volScalarField Tad=_T.clone();
    Foam::Time& runTime = _runTime();
    set_valueFraction();
    Tad_time.resize(0);
    ITHACAutilities::assignIF(Tad, homogeneousBC);
    label timeI = 0;
    assignDirectBC(timeI);
    List<scalar> RobinBC = - Tf;
    forAll(mesh.boundaryMesh(), patchI)
    {
        if (patchI == mesh.boundaryMesh().findPatchID("coldSide"))
        {
            ITHACAutilities::assignMixedBC(Tad, patchI, RobinBC, refGrad, valueFraction);
        }
        else
        {
            ITHACAutilities::assignBC(Tad, patchI, homogeneousBC);
        }
    }

    while (runTime.loop())
    {
        Info << "Time = " << runTime.timeName() << nl << endl;
        timeI++;
        assignDirectBC(timeI);
        RobinBC = - Tf;
        forAll(mesh.boundaryMesh(), patchI)
        {
            if (patchI == mesh.boundaryMesh().findPatchID("coldSide"))
            {
                ITHACAutilities::assignMixedBC(Tad, patchI, RobinBC, refGrad, valueFraction);
            }
            else
            {
                ITHACAutilities::assignBC(Tad, patchI, homogeneousBC);
            }
        }

        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix TEqn
            (
                fvm::ddt(Tad) - fvm::laplacian(DT * diffusivity, Tad)
            );
            fvOptions.constrain(TEqn);
            TEqn.solve();
            fvOptions.correct(Tad);
        }
        Tad_time.append(Tad.clone());
        ITHACAstream::exportSolution(Tad, std::to_string(timeSteps[timeI]),
                                     folderOffline,
                                     "Tad");
        runTime.printExecutionTime(Info);
        runTime.write();
    }

    Info << "Tad_time.size = " << Tad_time.size() << endl;
    Info << "Ntime = " << Ntimes << endl;
    addSol = fieldValueAtThermocouples(Tad_time);
    Info << "END \n" << endl;
}

void sequentialIHTP::solveDirect()
{
    if(offlineFlag)
    {
        restartOffline();
    }
    else
    {
        restart();
    }
    M_Assert(diffusivity>1e-36, "Set the diffusivity value");
    volScalarField& T = _T();
    assignDirectBC(0);
    if(offlineFlag)
    {
        ITHACAutilities::assignIF(T, homogeneousBC);
    }
    simpleControl& simple = _simple();
    Foam::Time& runTime = _runTime();
    fv::options& fvOptions(_fvOptions());
    label timeI = 0;
    Ttime.resize(0);

    while (runTime.loop())
    {
        Info << "Time = " << runTime.timeName() << nl << endl;
        timeI++;
        assignDirectBC(timeI);

        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix TEqn
            (
                fvm::ddt(T) - fvm::laplacian(DT * diffusivity, T)
            );
            fvOptions.constrain(TEqn);
            TEqn.solve();
            fvOptions.correct(T);
        }
        Ttime.append(T.clone());

        runTime.printExecutionTime(Info);
        runTime.write();
    }
    Info << "Direct computation ENDED" << endl;
    
}

void sequentialIHTP::readThermocouples()
{
    Info << "Defining positions of thermocouples" << endl;

    if (!thermocouplesRead)
    {
        word fileName = "./thermocouplesCellsID";

        if (ITHACAutilities::check_file(fileName + "_mat.txt"))
        {
            Info << "Reading thermocouples cells from file" << endl;
            Eigen::MatrixXi TCmatrix = ITHACAstream::readMatrix(fileName + "_mat.txt").cast
                                       <int> ();
            thermocouplesCellID = Foam2Eigen::EigenMatrix2List(TCmatrix);
        }
        else
        {
            Info << "Defining positions of thermocouples" << endl;
            fvMesh& mesh = _mesh();
            volScalarField& T = _T();
            thermocouplesCellID.resize(thermocouplesPos.size());
            forAll(thermocouplesPos, tcI)
            {
                thermocouplesCellID[tcI] = mesh.findCell(thermocouplesPos[tcI]);
            }
            volScalarField thermocouplesField(T);
            ITHACAutilities::assignIF(thermocouplesField, homogeneousBC);
            forAll(thermocouplesCellID, tcI)
            {
                thermocouplesField.ref()[thermocouplesCellID[tcI]] = 1;
            }
            ITHACAstream::exportSolution(thermocouplesField, "1", "./ITHACAoutput/thermocouplesField/",
                                         "thermocouplesField,");
            Eigen::MatrixXi thermocouplesCellID_eigen = Foam2Eigen::List2EigenMatrix(
                        thermocouplesCellID);
            ITHACAstream::exportMatrix(thermocouplesCellID_eigen, fileName,
                                       "eigen", "./");
        }

        thermocouplesRead = 1;
        samplingTime.resize(timeSamplesNum);
        forAll(samplingTime, timeI)
        {
            samplingTime[timeI] = timeSamplesT0 + timeI * timeSamplesDeltaT;
        }
        sampling2symulationTime();
	NtimeStepsBetweenSamples = timeSamplesDeltaT / deltaTime;
        residual.resize(thermocouplesNum * timeSamplesNum);
    }
    else
    {
        WarningInFunction << "readThermocouples function called twice." << endl;
        WarningInFunction << "I am not doing the second reading." << endl;
    }
}

Eigen::VectorXd sequentialIHTP::fieldValueAtThermocouples(
    volScalarField& field)
{
    if (!thermocouplesRead)
    {
        readThermocouples();
    }

    fvMesh& mesh = _mesh();
    dictionary interpolationDict =
        mesh.solutionDict().subDict("interpolationSchemes");
    autoPtr<Foam::interpolation<scalar>> fieldInterp =
                                          Foam::interpolation<scalar>::New(interpolationDict, field);
    Eigen::VectorXd fieldInt;
    fieldInt.resize(thermocouplesPos.size());
    forAll(thermocouplesPos, tcI)
    {
        fieldInt(tcI) = fieldInterp->interpolate(thermocouplesPos[tcI],
                        thermocouplesCellID[tcI]);
    }
    return fieldInt;
}

Eigen::VectorXd sequentialIHTP::fieldValueAtThermocouples(
    PtrList<volScalarField> fieldList, label fieldI)
{
    Eigen::VectorXd fieldInt = fieldValueAtThermocouples(fieldList[fieldI]);
    return fieldInt;
}

Eigen::VectorXd sequentialIHTP::fieldValueAtThermocouples(
    PtrList<volScalarField> fieldList)
{
    Eigen::VectorXd fieldInt;
    if ( fieldList.size() == Ntimes + 1 )
    {
        Info << "\n Sampling for ALL sampling times \n\n" << endl;
        fieldInt.resize(timeSamplesNum * thermocouplesNum);
        forAll(samplingSteps, sampleTimeI)
        {
            fieldInt.segment(sampleTimeI * thermocouplesNum, thermocouplesNum) =
                fieldValueAtThermocouples(fieldList, samplingSteps[sampleTimeI]);
        }
    }
    else if ( fieldList.size() == NtimeStepsBetweenSamples )
    {
        Info << "\nField size = " << fieldList.size() << ".\nSampling ONLY the last timestep\n\n" << endl;
        fieldInt =
            fieldValueAtThermocouples(fieldList, NtimeStepsBetweenSamples - 1);
    }
    else
    {
        Info << "The input fieldList of sequentialIHTP::fieldValueAtThermocouples can have size Ntimes + 1 (=" << Ntimes + 1 << ") or\n";
        Info << " NtimeStepsBetweenSamples  (=" <<  NtimeStepsBetweenSamples << ") but has size " << fieldList.size() << endl;
        Info << "Exiting." << endl;
        exit(23);
    }
    Info << "\nSampling done \n" << endl;
    return fieldInt;
}


void sequentialIHTP::restart()
{
    Time& runTime = _runTime();
    instantList Times = runTime.times();
    runTime.setTime(Times[1], 0);
    _simple.clear();
   _T.clear();

    Foam::fvMesh& mesh = _mesh();
    _simple = autoPtr<simpleControl>
              (
                  new simpleControl
                  (
                      mesh
                  )
              );

    _T = autoPtr<volScalarField>
         (
             new volScalarField
             (
                 IOobject
                 (
                     "T",
                     runTime.timeName(),
                     mesh,
                     IOobject::MUST_READ,
                     IOobject::AUTO_WRITE
                 ),
                 mesh
             )
         );
    

    Info << "Ready for new computation" << endl;
}

void sequentialIHTP::restartOffline()
{
    Info << "Setting endTime to offlineEndTime" << endl;
    restart();
    Time& runTime = _runTime();
    instantList Times = runTime.times();
    runTime.setTime(0.0, 0);
    runTime.setEndTime(offlineEndTime);
    Info << "Ready for new offline computation" << endl;
}

void sequentialIHTP::restartT0()
{
    Info << "Setting endTime to offlineEndTime" << endl;
    restart();
    Time& runTime = _runTime();
    instantList Times = runTime.times();
    runTime.setTime(0.0, 0);
    runTime.setEndTime(timeSamplesDeltaT);
    Info << "Ready for new T0 computation" << endl;
}

void sequentialIHTP::sampling2symulationTime()
{
    scalar EPSILON = 2e-16;
    label deltaTimeQuotient = std::floor(timeSamplesDeltaT / deltaTime);
    //Info << "debug: deltaTimeQuotient= " <<deltaTimeQuotient<< endl;
    //Info << "debug: timeSamplesDeltaT= " <<timeSamplesDeltaT<< endl;
    //Info << "debug: deltaTime= " <<deltaTime<< endl;
    M_Assert(std::fabs(timeSamplesDeltaT / deltaTime - std::trunc(
                           timeSamplesDeltaT / deltaTime)) < EPSILON,
             "timeSamplesDeltaT should be a multiple of deltaTime");
    label n0 = (timeSamplesT0 - startTime) / deltaTime;
    M_Assert(n0 > 0, "First sampling step cannot be 0");
    //Info << "debug: n0 = " << n0 << endl;
    //Info << "debug: (timeSamplesT0 - startTime) / deltaTime = " << (timeSamplesT0 - startTime) / deltaTime << endl;
    M_Assert(std::fabs(n0 * deltaTime - timeSamplesT0) < EPSILON,
             "The first sampling time must coincide with a simulation timestep");
    scalar samplingEndTime = timeSamplesDeltaT * (timeSamplesNum - 1) + timeSamplesT0;
    //Info << "debug: samplingEndTime = " << samplingEndTime << endl;
    //Info << "debug: EndTime = " << endTime << endl;
    M_Assert(!(endTime + EPSILON < samplingEndTime
               && std::fabs(endTime - samplingEndTime) > EPSILON),
             "The samplingEndTime cannot be later than the symulation endTime");
    samplingSteps.resize(timeSamplesNum);
    forAll(samplingTime, sampleI)
    {
        samplingSteps[sampleI] = n0 + sampleI * deltaTimeQuotient;
    }
    //Info << "debug: samplingSteps = " << samplingSteps << endl;
}

void sequentialIHTP::parameterizedBC_postProcess(
    List<Eigen::MatrixXd> linSys, Eigen::VectorXd weigths, word outputFolder,
    label verbose)
{
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(Theta,
                                          Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::MatrixXd singVal = svd.singularValues();
    ITHACAstream::exportMatrix(singVal, "singularValues", "eigen", outputFolder);
    if (verbose)
    {
        // Printing outputs at screen
        std::cout << "Singular values of Theta.transpose() * Theta are " << std::endl;
        std::cout << svd.singularValues() << std::endl;
        std::cout << "weigths = " << std::endl;
        std::cout << weigths << std::endl;
        std::cout << "linSys[1] = " << std::endl;
        std::cout << linSys[1] << std::endl;
        std::cout << "Theta = " << std::endl;
        std::cout << Theta << std::endl;
        residual =  linSys[0] * weigths - linSys[1];
        //std::cout << "Residual  = " << std::endl;
        //std::cout << residual << std::endl;
        std::cout << "Residual 2-norm = " << std::endl;
        std::cout << residual.squaredNorm() << std::endl;
        std::cout << "\n addSol = " << std::endl;
        std::cout << addSol << std::endl;
        std::cout << "T0_vector = " << std::endl;
        std::cout << T0_vector << std::endl;
        std::cout << "Tmeas = " << std::endl;
        std::cout << Tmeas << std::endl;
    }

    reconstrucT(outputFolder);
    Tcomp = fieldValueAtThermocouples(Ttime);
    std::cout << "Tcomp = \n" << Tcomp.transpose() << std::endl;
    std::cout << "TmeasShort = \n" << TmeasShort.transpose() << std::endl;
    J = 0.5 * Foam::sqrt((Tcomp - TmeasShort).dot(Tcomp - TmeasShort));
    Info << "J = " << J << endl;
    Jlist.conservativeResize(Jlist.size() + 1);
    Jlist(Jlist.size() - 1) = J;
}

void sequentialIHTP::findMagicPoints(int NmagicPoints)
{
    M_Assert(NmagicPoints > 0, "Set number of magic points");
    magicPoints.clear();
    M_Assert(NmagicPoints <= NmodesT0, 
            "Number of magic points bigger than number of modes");
    Eigen::MatrixXd A;
    Eigen::VectorXd b;
    Eigen::VectorXd c;
    Eigen::VectorXd r;
    Eigen::VectorXd rho(1);
    Eigen::MatrixXd MatrixModes = T0modes.toEigen()[0];
    label ind_max, c1;
    double max = MatrixModes.cwiseAbs().col(0).maxCoeff(&ind_max, &c1);
    rho(0) = max;
    magicPoints.append(ind_max);
    Eigen::MatrixXd U = MatrixModes.col(0);
    Eigen::SparseMatrix<double> P;
    P.resize(MatrixModes.rows(), 1);
    P.insert(ind_max, 0) = 1;

    for (label i = 1; i < NmagicPoints; i++)
    {
        A = P.transpose() * U;
        b = P.transpose() * MatrixModes.col(i);
        c = A.fullPivLu().solve(b);
        r = MatrixModes.col(i) - U * c;
        max = r.cwiseAbs().maxCoeff(&ind_max, &c1);
        P.conservativeResize(MatrixModes.rows(), i + 1);
        P.insert(ind_max, i) = 1;
        U.conservativeResize(MatrixModes.rows(), i + 1);
        U.col(i) =  MatrixModes.col(i);
        rho.conservativeResize(i + 1);
        rho(i) = max;
        magicPoints.append(ind_max);
    }

    Info << "magicPoints:\n" << magicPoints << endl;
}
