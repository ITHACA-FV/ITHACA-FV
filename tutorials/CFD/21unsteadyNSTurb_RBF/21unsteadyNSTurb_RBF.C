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
    Example of an unsteady NS Reduction Problem
SourceFiles
    04unsteadyNS.C
\*---------------------------------------------------------------------------*/

#include "unsteadyNS.H"
#include "UnsteadyNSTurb.H"
#include "ITHACAPOD.H"
#include "ReducedUnsteadyNS.H"
#include "ReducedUnsteadyNSTurb.H"
#include "ITHACAstream.H"
#include <chrono>
#include<math.h>
#include<iomanip>

#include "fvCFD.H"
#include "fvOptions.H"
#include "IOmanip.H"
#include "ITHACAutilities.H"
#include "Foam2Eigen.H"
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <Eigen/SparseLU>
#include "EigenFunctions.H"
#include "reductionProblem.H"
#include "steadyNS.H"

class tutorial21 : public UnsteadyNSTurb
{
public:
    explicit tutorial21(int argc, char* argv[])
        : UnsteadyNSTurb(argc, argv),
          U(_U()),
          p(_p()),
          nut(_nut())
    {}

    // Fields
    volVectorField& U;
    volScalarField& p;
    volScalarField& nut;

    void offlineSolve(std::string offlinepath);
    void computeSnapshots(std::string inputPath, std::string avgOutputPath, std::string fluctOutputPath);
    void computeModes(
    std::string avgInputPath, 
    std::string fluctInputPath,
    std::string avgModesOutputPath, 
    std::string fluctModesOutputPath,
    int NmodesU_avg, 
    int NmodesNut_avg, 
    int NmodesU_fluct, 
    int NmodesNut_fluct
);
    void computeCoefficients(std::string avgSnapshotsInputPath, std::string fluctSnapshotsInputPath,
                             std::string avgModesInputPath, std::string fluctModesInputPath,
                             std::string coeffsOutputPath);

};

void tutorial21::offlineSolve(std::string offlinepath)
{
    Vector<double> inl(1, 0, 0);
    List<scalar> mu_now(1);
    label BCind = 0;

    if ((offline) && (ITHACAutilities::check_folder(offlinepath) == true))
    {
        ITHACAstream::read_fields(Ufield, U, offlinepath);
        ITHACAstream::read_fields(Pfield, p, offlinepath);
        ITHACAstream::read_fields(nutFields, nut, offlinepath);
    }
    else
    {
        for (label i = 0; i < mu.cols(); i++)
        {
            inl[0] = mu(0, i);
            mu_now[0] = mu(0, i);
            assignBC(U, BCind, inl);
            truthSolve(mu_now, offlinepath);
        }
    }
}
std::vector<std::pair<int, int>> velocity_ranges = {
    {0, 199},
    {200, 399},
    {400, 599},
    {600, 799},
    {800, 999},
    {1000, 1199},
    {1200, 1399},
    {1400, 1599},
    {1600, 1799},
    {1800, 1999},
    {2000, 2199},
    {2200, 2399}
};
void tutorial21::computeSnapshots(std::string inputPath, std::string avgOutputPath, std::string fluctOutputPath) {
    PtrList<volVectorField> avgUList(velocity_ranges.size());
    PtrList<volVectorField> fluctUList;
    PtrList<volScalarField> avgNutList(velocity_ranges.size());
    PtrList<volScalarField> fluctNutList;

    for (size_t group = 0; group < velocity_ranges.size(); ++group) {
        int start = velocity_ranges[group].first;
        int end = velocity_ranges[group].second;
        int numSnapshots = end - start + 1;

        // Load the U snapshots
        PtrList<volVectorField> Ufield;
        ITHACAstream::read_fields(Ufield, U, inputPath, start, numSnapshots);

        // Load the nut snapshots
        PtrList<volScalarField> nutField;
        ITHACAstream::read_fields(nutField, nut, inputPath, start, numSnapshots);

        if (Ufield.empty() || nutField.empty()) {
            std::cerr << " No snapshots found for group " << group << "\n";
            continue;
        }

        // Compute U average
        volVectorField avgU = Ufield[0];
        for (int i = 1; i < numSnapshots; ++i)
            avgU.internalFieldRef() += Ufield[i].internalField();
        avgU.internalFieldRef() /= numSnapshots;
        avgUList.set(group, new volVectorField(avgU));

        // Compute nut average
        volScalarField avgNut = nutField[0];
        for (int i = 1; i < numSnapshots; ++i)
            avgNut.internalFieldRef() += nutField[i].internalField();
        avgNut.internalFieldRef() /= numSnapshots;
        avgNutList.set(group, new volScalarField(avgNut));

        // Store U fluct and nut fluct
        for (int i = 0; i < numSnapshots; ++i) {
            volVectorField fluctU = Ufield[i];
            fluctU.internalFieldRef() -= avgU.internalField();
            fluctUList.append(new volVectorField(fluctU));

            volScalarField fluctNut = nutField[i];
            fluctNut.internalFieldRef() -= avgNut.internalField();
            fluctNutList.append(new volScalarField(fluctNut));
        }
    }

    // Export U and nut averages
    ITHACAstream::exportFields(avgUList, avgOutputPath, "U_avg");
    ITHACAstream::exportFields(avgNutList, avgOutputPath, "Nut_avg");

    // Export fluctuations
    ITHACAstream::exportFields(fluctUList, fluctOutputPath, "U_fluct");
    ITHACAstream::exportFields(fluctNutList, fluctOutputPath, "Nut_fluct");

    std::cout << " U and nut: average and fluctuation snapshots written to disk.\n";
}



void tutorial21::computeModes(
    std::string avgInputPath, 
    std::string fluctInputPath,
    std::string avgModesOutputPath, 
    std::string fluctModesOutputPath,
    int NmodesU_avg, 
    int NmodesNut_avg, 
    int NmodesU_fluct, 
    int NmodesNut_fluct)
{
    // Load averaged snapshots explicitly
    PtrList<volVectorField> avgUList;
    PtrList<volScalarField> avgNutList;

    ITHACAstream::read_fields(avgUList, "U_avg", avgInputPath, 0, 12);
    ITHACAstream::read_fields(avgNutList, "Nut_avg", avgInputPath, 0, 12);

    PtrList<volVectorField> UavgModes;
    PtrList<volScalarField> NutavgModes;

    NmodesU_avg = std::min(NmodesU_avg, static_cast<int>(avgUList.size()));
    NmodesNut_avg = std::min(NmodesNut_avg, static_cast<int>(avgNutList.size()));

    ITHACAPOD::getModes(avgUList, UavgModes, U.name(), podex, supex, 0, NmodesU_avg);
    ITHACAPOD::getModes(avgNutList, NutavgModes, nut.name(), podex, supex, 0, NmodesNut_avg);

    ITHACAstream::exportFields(UavgModes, avgModesOutputPath, "U_avg_modes");
    ITHACAstream::exportFields(NutavgModes, avgModesOutputPath, "Nut_avg_modes");

    // Clearly handle fluctuating snapshots explicitly
    PtrList<volVectorField> fluctUList;
    PtrList<volScalarField> fluctNutList;

    ITHACAstream::read_fields(fluctUList, "U_fluct", fluctInputPath, 0, 2400);
    ITHACAstream::read_fields(fluctNutList, "Nut_fluct", fluctInputPath, 0, 2400);

    PtrList<volVectorField> UfluctModes;
    PtrList<volScalarField> NutfluctModes;

    NmodesU_fluct = std::min(NmodesU_fluct, static_cast<int>(fluctUList.size()));
    NmodesNut_fluct = std::min(NmodesNut_fluct, static_cast<int>(fluctNutList.size()));

    ITHACAPOD::getModes(fluctUList, UfluctModes, U.name(), podex, supex, 0, NmodesU_fluct);
    ITHACAPOD::getModes(fluctNutList, NutfluctModes, nut.name(), podex, supex, 0, NmodesNut_fluct);

    ITHACAstream::exportFields(UfluctModes, fluctModesOutputPath, "U_fluct_modes");
    ITHACAstream::exportFields(NutfluctModes, fluctModesOutputPath, "Nut_fluct_modes");


}
void tutorial21::computeCoefficients(
    std::string avgSnapshotsInputPath, 
    std::string fluctSnapshotsInputPath,
    std::string avgModesInputPath, 
    std::string fluctModesInputPath,
    std::string coeffsOutputPath)
{
    // --- Loading snapshots explicitly ---
    PtrList<volVectorField> UavgSnapshots, UfluctSnapshots;
    PtrList<volScalarField> NutavgSnapshots, NutfluctSnapshots;

    ITHACAstream::read_fields(UavgSnapshots, "U_avg", avgSnapshotsInputPath, 0, 12);
    ITHACAstream::read_fields(NutavgSnapshots, "Nut_avg", avgSnapshotsInputPath, 0, 12);
    ITHACAstream::read_fields(UfluctSnapshots, "U_fluct", fluctSnapshotsInputPath, 0, 2400);
    ITHACAstream::read_fields(NutfluctSnapshots, "Nut_fluct", fluctSnapshotsInputPath, 0, 2400);

    std::cout << "Loaded snapshots explicitly:\n"
              << "U_avg: " << UavgSnapshots.size() << "\n"
              << "Nut_avg: " << NutavgSnapshots.size() << "\n"
              << "U_fluct: " << UfluctSnapshots.size() << "\n"
              << "Nut_fluct: " << NutfluctSnapshots.size() << std::endl;

    // --- Loading modes explicitly ---
    PtrList<volVectorField> UavgModes, UfluctModes;
    PtrList<volScalarField> NutavgModes, NutfluctModes;

    ITHACAstream::read_fields(UavgModes, "U_avg_modes", avgModesInputPath, 0, 12);
    ITHACAstream::read_fields(NutavgModes, "Nut_avg_modes", avgModesInputPath, 0, 12);
    ITHACAstream::read_fields(UfluctModes, "U_fluct_modes", fluctModesInputPath, 0, 15);
    ITHACAstream::read_fields(NutfluctModes, "Nut_fluct_modes", fluctModesInputPath, 0,15);
    std::cout << "Loaded modes explicitly:\n"
              << "U_avg_modes: " << UavgModes.size() << "\n"
              << "Nut_avg_modes: " << NutavgModes.size() << "\n"
              << "U_fluct_modes: " << UfluctModes.size() << "\n"
              << "Nut_fluct_modes: " << NutfluctModes.size() << std::endl;

    // --- Ensure modes are explicitly loaded ---
    if (UavgModes.empty() || NutavgModes.empty() || UfluctModes.empty() || NutfluctModes.empty()) {
        FatalErrorIn("computeCoefficients")
            << "Modes are missing or not loaded correctly."
            << exit(FatalError);
    }

    // --- Compute coefficients explicitly ---
    Eigen::MatrixXd UavgCoefficients = ITHACAutilities::getCoeffs(UavgSnapshots, UavgModes);
    Eigen::MatrixXd NutavgCoefficients = ITHACAutilities::getCoeffs(NutavgSnapshots, NutavgModes);
    Eigen::MatrixXd UfluctCoefficients = ITHACAutilities::getCoeffs(UfluctSnapshots, UfluctModes);
    Eigen::MatrixXd NutfluctCoefficients = ITHACAutilities::getCoeffs(NutfluctSnapshots, NutfluctModes);

    // --- Debugging explicitly coefficient sizes ---
    std::cout << "Computed coefficient sizes explicitly:\n"
              << "U_avg: " << UavgCoefficients.rows() << "x" << UavgCoefficients.cols() << "\n"
              << "Nut_avg: " << NutavgCoefficients.rows() << "x" << NutavgCoefficients.cols() << "\n"
              << "U_fluct: " << UfluctCoefficients.rows() << "x" << UfluctCoefficients.cols() << "\n"
              << "Nut_fluct: " << NutfluctCoefficients.rows() << "x" << NutfluctCoefficients.cols() << std::endl;

    // --- Export coefficients explicitly to both Python and MATLAB ---
    ITHACAstream::exportMatrix(UavgCoefficients, "U_avg_coeffs", "python", coeffsOutputPath);
    ITHACAstream::exportMatrix(NutavgCoefficients, "Nut_avg_coeffs", "python", coeffsOutputPath);
    ITHACAstream::exportMatrix(UfluctCoefficients, "U_fluct_coeffs", "python", coeffsOutputPath);
    ITHACAstream::exportMatrix(NutfluctCoefficients, "Nut_fluct_coeffs", "python", coeffsOutputPath);

    ITHACAstream::exportMatrix(UavgCoefficients, "U_avg_coeffs", "matlab", coeffsOutputPath);
    ITHACAstream::exportMatrix(NutavgCoefficients, "Nut_avg_coeffs", "matlab", coeffsOutputPath);
    ITHACAstream::exportMatrix(UfluctCoefficients, "U_fluct_coeffs", "matlab", coeffsOutputPath);
    ITHACAstream::exportMatrix(NutfluctCoefficients, "Nut_fluct_coeffs", "matlab", coeffsOutputPath);
    
    ITHACAstream::exportMatrix(UavgCoefficients, "U_avg_coeffs", "eigen", coeffsOutputPath);
    ITHACAstream::exportMatrix(NutavgCoefficients, "Nut_avg_coeffs", "eigen", coeffsOutputPath);
    ITHACAstream::exportMatrix(UfluctCoefficients, "U_fluct_coeffs", "eigen", coeffsOutputPath);
    ITHACAstream::exportMatrix(NutfluctCoefficients, "Nut_fluct_coeffs", "eigen", coeffsOutputPath);

    std::cout << "Exported all coefficients explicitly and successfully to: " << coeffsOutputPath << std::endl;
};

/*---------------------------------------------------------------------------*\
                               Starting the MAIN
\*---------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
    // Construct the tutorial21 object
    tutorial21 example(argc, argv);

    // Read parameters from ITHACAdict file
    ITHACAparameters* para =
        ITHACAparameters::getInstance(example._mesh(), example._runTime());

    const int NmodesU       = para->ITHACAdict->lookupOrDefault<int>("NmodesU", 15);
    const int NmodesP       = para->ITHACAdict->lookupOrDefault<int>("NmodesP", 15);
    const int NmodesSUP     = para->ITHACAdict->lookupOrDefault<int>("NmodesSUP", 15);
    const int NmodesNUT     = para->ITHACAdict->lookupOrDefault<int>("NmodesNUT", 15);
    const int NmodesProject = para->ITHACAdict->lookupOrDefault<int>("NmodesProject", 15);

    // --- Fixed settings (Book Table 8.1 extended to 12.5 and 13.0) ---
    // Keep writeEvery = 20 * dt for integer alignment.
    const std::vector<scalar> inletVel = {
        7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0
    };

    const std::vector<scalar> timeSteps = {
        0.0004, 0.0004, 0.00035, 0.0003, 0.0003, 0.0003,
        0.0003, 0.0003, 0.00025, 0.00025, 0.00025, 0.00025
    };

    const std::vector<scalar> writeEverys = {
        0.008, 0.008, 0.007, 0.006, 0.006, 0.006,
        0.006, 0.006, 0.005, 0.005, 0.005, 0.005
    };

    const scalar startTimeAll = 10.0;  // start saving after warm-up (periodic regime)
    const label  snapsPerU    = 200;   // target snapshots per velocity

    // Prepare parameters scaffold (single parameter)
    example.Pnumber = 1;
    example.Tnumber = 1;
    example.setParameters();

    // Inlet patch index
    example.inletIndex.resize(1, 2);
    example.inletIndex(0, 0) = 0;
    example.inletIndex(0, 1) = 0;

    // Accumulate all mu_samples across all U into one matrix [time, mu]
    Eigen::MatrixXd allMuSamples(0, 2);

    for (label i = 0; i < static_cast<label>(inletVel.size()); ++i)
    {
        // Set the single-parameter sample to this U
        example.mu_range.resize(1, 2);
        example.mu_range(0, 0) = inletVel[i];
        example.mu_range(0, 1) = inletVel[i];
        example.genEquiPar();

        // Per-U time controls
        example.startTime  = startTimeAll;
        example.timeStep   = timeSteps[i];
        example.writeEvery = writeEverys[i];

        // Exactly 200 snapshots:
        //  - first at startTime + writeEvery
        //  - last  at startTime + 200*writeEvery
        //
        // Use a tiny guard (-SMALL) to avoid an extra snapshot on an exact final-time write
        // while not accidentally skipping the last intended write due to float/time-control logic.
        example.finalTime = example.startTime + snapsPerU * example.writeEvery - SMALL;

        Info<< "\n[OFFLINE] U_in=" << inletVel[i]
            << "   dt=" << example.timeStep
            << "   writeEvery=" << example.writeEvery
            << "   ratio=" << (example.writeEvery / example.timeStep)  // should be 20
            << "   start=" << example.startTime
            << "   final=" << example.finalTime
            << "   (snapshots planned = " << snapsPerU << ")\n"
            << endl;

        // Reset mu_samples (time, mu) for this U before running
        example.mu_samples.resize(0, 2);

        // Run and collect
        example.offlineSolve("./ITHACAoutput/Offline/");

        // Append to global mu_samples
        const label oldRows = allMuSamples.rows();
        allMuSamples.conservativeResize(oldRows + example.mu_samples.rows(), 2);
        allMuSamples.block(oldRows, 0, example.mu_samples.rows(), 2) = example.mu_samples;
    }

    // Export combined mu_samples
    ITHACAstream::exportMatrix(allMuSamples, "mu_samples", "eigen",  "./ITHACAoutput/Offline/");
    ITHACAstream::exportMatrix(allMuSamples, "mu_samples", "python", "./ITHACAoutput/Offline/");
    ITHACAstream::exportMatrix(allMuSamples, "mu_samples", "matlab", "./ITHACAoutput/Offline/");

    Info<< "\nAll snapshots and full mu_samples exported for all 12 inlet velocities.\n"
        << endl;


   example.computeSnapshots("./ITHACAoutput/Offline/",
                             "./ITHACAoutput/Offline_avg/",
                             "./ITHACAoutput/Offline_fluct/");

    // Call explicitly using existing variables:
example.computeModes(
        "./ITHACAoutput/Offline_avg/",
        "./ITHACAoutput/Offline_fluct/",
        "./ITHACAoutput/Modes_avg/",
        "./ITHACAoutput/Modes_fluct/",
        NmodesU, NmodesNUT, NmodesU, NmodesNUT);
     std::string avgSnapshotsInputPath = "./ITHACAoutput/Offline_avg/";
    std::string fluctSnapshotsInputPath = "./ITHACAoutput/Offline_fluct/";
    std::string avgModesInputPath = "./ITHACAoutput/Modes_avg/";
    std::string fluctModesInputPath = "./ITHACAoutput/Modes_fluct/";
    std::string coeffsOutputPath = "./ITHACAoutput/Coefficients/";

    example.computeCoefficients(avgSnapshotsInputPath, fluctSnapshotsInputPath,
                                avgModesInputPath, fluctModesInputPath,
                                coeffsOutputPath);

    // Define velRBF
    auto mu_mat =
        ITHACAstream::readMatrix("./ITHACAoutput/Offline/mu_samples_mat.txt");
    example.velRBF = mu_mat.col(1);
    // Solve the supremizer problem
    example.solvesupremizer();
    // Search the lift function
    example.liftSolve();
    // Normalize the lifting function
    ITHACAutilities::normalizeFields(example.liftfield);
    // Create homogeneous basis functions for velocity
    example.computeLift(example.Ufield, example.liftfield, example.Uomfield);
    // Export the homogeneous velocity snapshots
    ITHACAstream::exportFields(example.Uomfield, "./ITHACAoutput/Offline",
                               "Uofield");
    // Perform a POD decomposition for the velocity, the pressure and the eddy viscosity
    ITHACAPOD::getModes(example.Uomfield, example.Umodes, example._U().name(),
                        example.podex,
                        example.supex, 0, NmodesProject);
    // -----------------------------------------------------------------------------
// Compute projection coefficients of Uomfield onto U modes (20, 3000)
// This will be used for eddy viscosity interpolation (Section 8.4.2 of Rozza)
Eigen::MatrixXd velCoeffs = ITHACAutilities::getCoeffs(example.Uomfield, example.Umodes);

// Debugging info
std::cout << " vel.txt shape: (" << velCoeffs.rows() << ", " << velCoeffs.cols() << ")" << std::endl;
// Export to Python and MATLAB formats
ITHACAstream::exportMatrix(velCoeffs, "vel", "python", "./ITHACAoutput/Coefficients/");
ITHACAstream::exportMatrix(velCoeffs, "vel", "matlab", "./ITHACAoutput/Coefficients/");
ITHACAstream::exportMatrix(velCoeffs, "vel", "eigen", "./ITHACAoutput/Coefficients/");

    ITHACAPOD::getModes(example.Pfield, example.Pmodes, example.p().name(),
                        example.podex,
                        example.supex, 0, NmodesProject);
    ITHACAPOD::getModes(example.nutFields, example.nutModes, example._nut().name(),
                        example.podex,
                        example.supex, 0, NmodesProject);
    // Solve the supremizer problem based on the pressure modes
    example.solvesupremizer("modes");
    
// ----------------- Load nut avg/fluct fields and modes (IMPORTANT!) ----------------- //
    PtrList<volScalarField> nutAvgFields;
    ITHACAstream::read_fields(nutAvgFields, "Nut_avg", "./ITHACAoutput/Offline_avg/", 0, 12);

    PtrList<volScalarField> nutAvgModes;
    ITHACAstream::read_fields(nutAvgModes, "Nut_avg_modes", "./ITHACAoutput/Modes_avg/", 0, 12);

    PtrList<volScalarField> nutFluctFields;
    ITHACAstream::read_fields(nutFluctFields, "Nut_fluct", "./ITHACAoutput/Offline_fluct/", 0, 2400);

    PtrList<volScalarField> nutFluctModes;
    ITHACAstream::read_fields(nutFluctModes, "Nut_fluct_modes", "./ITHACAoutput/Modes_fluct/", 0, 15);
std::cout << "[DEBUG] nutFluctModes.size() = " << nutFluctModes.size() << std::endl;
example.nutAve = nutAvgModes;             // Average nut fields
example.nutFluctModes = nutFluctModes;   // Fluctuation nut fields

// (Optional: modes if used elsewhere)
example.nutAvgModes = nutAvgModes;
example.nutFluctModes = nutFluctModes;
example.nNutAvgModes = nutAvgModes.size();
example.nNutFluctModes = nutFluctModes.size();


    // Debug print
    std::cout << "[DEBUG] nutAve.size()    = " << example.nutAve.size()    << std::endl;
    std::cout << "[DEBUG] nutModes.size()  = " << example.nutModes.size()  << std::endl;
    std::cout << "[DEBUG] nutAvgModes.size() = " << example.nutAvgModes.size() << std::endl;

    // ---------------- Project (build tensors & matrices) ---------------- //
  example.projectSUP("./Matrices", NmodesU, NmodesP, NmodesSUP,
                       NmodesNUT, true);
  std::cout << "[MAIN] rbfSplinesNutAvg in FOM: " << example.rbfSplinesNutAvg.size() << std::endl;
std::cout << "[MAIN] rbfSplinesNutFluct in FOM: " << example.rbfSplinesNutFluct.size() << std::endl;


// Create an object of the turbulent class
ReducedUnsteadyNSTurb pod_rbf(example);

// Set value of the reduced viscosity and the penalty factor
pod_rbf.nu = 1e-04;
pod_rbf.tauU.resize(1, 1);
pod_rbf.tstart = 0;
pod_rbf.finalTime = 1.194;
pod_rbf.dt = 0.006;
pod_rbf.storeEvery = 0.006;
pod_rbf.exportEvery = 0.006;

// For single online velocity (k = 0)
Eigen::MatrixXd velNow(1, 1);
velNow(0, 0) = 7.75;
pod_rbf.tauU(0, 0) = 0;

// If you want to set the initial condition explicitly
//pod_rbf.initCond = Eigen::MatrixXd::Zero(pod_rbf.Nphi_u + pod_rbf.Nphi_p, 1);

// === Perform the online solve ===
pod_rbf.solveOnlineSUPAve(velNow);

// === Diagnostics: print shapes before copying ===
std::cout << "[CHECK] rbfCoeffMat size: "
          << pod_rbf.rbfCoeffMat.rows() << " x " << pod_rbf.rbfCoeffMat.cols() << std::endl;
std::cout << "[CHECK] online_solution.size(): " << pod_rbf.online_solution.size() << std::endl;

// === SAFELY copy all interpolated eddy viscosity coefficients ===
// The full time series is usually needed for post-processing and reconstruction
Eigen::MatrixXd rbfCoeff = pod_rbf.rbfCoeffMat;

// === Export matrices ===
ITHACAstream::exportMatrix(rbfCoeff, "rbfCoeff", "python", "./ITHACAoutput/Matrices/");
ITHACAstream::exportMatrix(rbfCoeff, "rbfCoeff", "matlab", "./ITHACAoutput/Matrices/");
ITHACAstream::exportMatrix(pod_rbf.online_solution, "red_coeff", "python", "./ITHACAoutput/red_coeff");
ITHACAstream::exportMatrix(pod_rbf.online_solution, "red_coeff", "matlab", "./ITHACAoutput/red_coeff");
ITHACAstream::exportMatrix(pod_rbf.online_solution, "red_coeff", "eigen", "./ITHACAoutput/red_coeff");

// === Call reconstruction ===
pod_rbf.reconstruct(true, "./ITHACAoutput/Reconstruction/");


std::cout << "[INFO] Main completed successfully." << std::endl;

    // Solve the full order problem for the online velocity values for the purpose of comparison
    tutorial21 example2(argc, argv);
    /// Set the number of parameters
    example2.Pnumber = 1;
    /// Set samples
    example2.Tnumber = 1;
    /// Set the parameters infos
    example2.setParameters();
    // Set the parameter ranges
    example2.mu_range(0, 0) = 7.75;
    example2.mu_range(0, 1) = 7.75;
    // Generate equispaced samples inside the parameter range
    example2.genEquiPar();
    // Set the inlet boundaries where we have non homogeneous boundary conditions
    example2.inletIndex.resize(1, 2);
    example2.inletIndex(0, 0) = 0;
    example2.inletIndex(0, 1) = 0;
    // Time parameters
    example2.startTime = 9.1;
    example2.finalTime = 10.30;
    example2.timeStep = 0.0004;
    example2.writeEvery = 0.006;
    example2.offlineSolve("./ITHACAoutput/Offline_check/");
// -------------------- Standard relative errors (as before) --------------------
// -------------------- Standard relative errors (as before) --------------------
Eigen::MatrixXd errFrobU  = ITHACAutilities::errorFrobRel(example2.Ufield,  pod_rbf.uRecFields);
Eigen::MatrixXd errFrobP  = ITHACAutilities::errorFrobRel(example2.Pfield,  pod_rbf.pRecFields);
Eigen::MatrixXd errFrobNu = ITHACAutilities::errorFrobRel(example2.nutFields, pod_rbf.nutRecFields);

ITHACAstream::exportMatrix(errFrobU,  "errFrobU",  "matlab", "./ITHACAoutput/ErrorsFrob/");
ITHACAstream::exportMatrix(errFrobP,  "errFrobP",  "matlab", "./ITHACAoutput/ErrorsFrob/");
ITHACAstream::exportMatrix(errFrobNu, "errFrobNut","matlab", "./ITHACAoutput/ErrorsFrob/");

Eigen::MatrixXd errL2U  = ITHACAutilities::errorL2Rel(example2.Ufield,  pod_rbf.uRecFields);
Eigen::MatrixXd errL2P  = ITHACAutilities::errorL2Rel(example2.Pfield,  pod_rbf.pRecFields);
Eigen::MatrixXd errL2Nu = ITHACAutilities::errorL2Rel(example2.nutFields, pod_rbf.nutRecFields);

ITHACAstream::exportMatrix(errL2U,  "errL2U",  "matlab", "./ITHACAoutput/ErrorsL2/");
ITHACAstream::exportMatrix(errL2P,  "errL2P",  "matlab", "./ITHACAoutput/ErrorsL2/");
ITHACAstream::exportMatrix(errL2Nu, "errL2Nut","matlab", "./ITHACAoutput/ErrorsL2/");

// --------------------------- Helpers (lambdas) ---------------------------
// Domain volume |Ω| (dimensioned) and its numeric value
auto VtotDim = [&](const fvMesh& mesh) -> dimensionedScalar { return sum(mesh.V()); };

// Proper dimensioned handling for L2 norms (no .value() inside Foam::sqrt!)
auto L2s = [&](const volScalarField& f) -> scalar
{
    const dimensionedScalar I  = fvc::domainIntegrate(sqr(f)); // [f]^2·[V]
    const dimensionedScalar L2 = Foam::sqrt(I);                // [f]·[V]^{1/2}
    return L2.value();
};
auto L2v = [&](const volVectorField& f) -> scalar
{
    const dimensionedScalar I  = fvc::domainIntegrate(magSqr(f)); // [U]^2·[V]
    const dimensionedScalar L2 = Foam::sqrt(I);                   // [U]·[V]^{1/2}
    return L2.value();
};

// ---------------- Absolute L2 errors (unnormalized) for U and p -------------
auto absL2SeriesVector = [&](const PtrList<volVectorField>& FOM,
                             const PtrList<volVectorField>& ROM) -> Eigen::MatrixXd
{
    const label n = Foam::min(FOM.size(), ROM.size());
    if (FOM.size() != ROM.size())
        Info<< "[WARN] absL2(U): snapshot mismatch FOM=" << FOM.size()
            << " ROM=" << ROM.size() << " -> trunc " << n << nl;

    Eigen::MatrixXd out(n, 1);
    for (label i = 0; i < n; ++i)
    {
        tmp<volVectorField> tDiff = FOM[i] - ROM[i];
        out(i, 0) = L2v(tDiff());
    }
    return out;
};

auto absL2SeriesScalar = [&](const PtrList<volScalarField>& FOM,
                             const PtrList<volScalarField>& ROM) -> Eigen::MatrixXd
{
    const label n = Foam::min(FOM.size(), ROM.size());
    if (FOM.size() != ROM.size())
        Info<< "[WARN] absL2(p): snapshot mismatch FOM=" << FOM.size()
            << " ROM=" << ROM.size() << " -> trunc " << n << nl;

    Eigen::MatrixXd out(n, 1);
    for (label i = 0; i < n; ++i)
    {
        tmp<volScalarField> tDiff = FOM[i] - ROM[i];
        out(i, 0) = L2s(tDiff());
    }
    return out;
};

// Compute and export absolute L2 time series
Eigen::MatrixXd errL2U_abs = absL2SeriesVector(example2.Ufield, pod_rbf.uRecFields);
Eigen::MatrixXd errL2P_abs = absL2SeriesScalar(example2.Pfield, pod_rbf.pRecFields);

ITHACAstream::exportMatrix(errL2U_abs, "errL2U_abs", "matlab", "./ITHACAoutput/ErrorsL2/");
ITHACAstream::exportMatrix(errL2P_abs, "errL2P_abs", "matlab", "./ITHACAoutput/ErrorsL2/");

// ---------------------- RMS (divide by sqrt(|Ω|)) -------------------------
const fvMesh& meshCheck = example2.Ufield[0].mesh();
const scalar sqrtV = Foam::sqrt(VtotDim(meshCheck)).value();

Eigen::MatrixXd errL2U_rms(errL2U_abs.rows(), 1);
Eigen::MatrixXd errL2P_rms(errL2P_abs.rows(), 1);
for (label i = 0; i < errL2U_abs.rows(); ++i) errL2U_rms(i,0) = errL2U_abs(i,0) / (sqrtV + SMALL);
for (label i = 0; i < errL2P_abs.rows(); ++i) errL2P_rms(i,0) = errL2P_abs(i,0) / (sqrtV + SMALL);

ITHACAstream::exportMatrix(errL2U_rms, "errL2U_rms", "matlab", "./ITHACAoutput/ErrorsL2/");
ITHACAstream::exportMatrix(errL2P_rms, "errL2P_rms", "matlab", "./ITHACAoutput/ErrorsL2/");

// ---------------------- Option A: scaled RMS (u*, p*) ----------------------
const scalar rho  = 1.0;   // set appropriately if not 1
const scalar Vinf = 7.75;  // your online inlet speed (or read from input)
const scalar qInf = 0.5 * rho * Vinf * Vinf;

Eigen::MatrixXd errU_rms_star(errL2U_rms.rows(), 1);
Eigen::MatrixXd errP_rms_star(errL2P_rms.rows(), 1);

for (label i = 0; i < errL2U_rms.rows(); ++i) errU_rms_star(i,0) = errL2U_rms(i,0) / (Vinf + SMALL);
for (label i = 0; i < errL2P_rms.rows(); ++i) errP_rms_star(i,0) = errL2P_rms(i,0) / (qInf + SMALL);

ITHACAstream::exportMatrix(errU_rms_star, "errU_rms_star", "matlab", "./ITHACAoutput/ErrorsL2/");
ITHACAstream::exportMatrix(errP_rms_star, "errP_rms_star", "matlab", "./ITHACAoutput/ErrorsL2/");

// -------------------------- Console summaries --------------------------------
Info<< "[ABS-L2] U: n=" << errL2U_abs.rows()
    << " min=" << errL2U_abs.minCoeff()
    << " max=" << errL2U_abs.maxCoeff()
    << " mean="<< errL2U_abs.mean() << nl;

Info<< "[ABS-L2] p: n=" << errL2P_abs.rows()
    << " min=" << errL2P_abs.minCoeff()
    << " max=" << errL2P_abs.maxCoeff()
    << " mean="<< errL2P_abs.mean() << nl;

Info<< "[RMS] U: mean=" << errL2U_rms.mean()
    << " | p: mean=" << errL2P_rms.mean() << nl;

Info<< "[Scaled RMS] u*: mean=" << errU_rms_star.mean()
    << " | p*: mean=" << errP_rms_star.mean() << nl;



}
