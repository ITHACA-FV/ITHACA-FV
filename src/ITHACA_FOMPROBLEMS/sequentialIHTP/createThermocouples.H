IOdictionary thermocouplesPosition  // Kabir: "Right now, thermocouplesPosition has as many thermocouples as thermocouplesDict has. It means that thermocouplesNum is equal to the number of thermocouples in thermocouplesDict."
(    
    IOobject
    (
        "thermocouplesDict",        // inside 06enKFwDF_3dIHTP.H the function named int getObservationSize() reads measurementsDict
        runTime.constant(),         // inside 06enKFwDF_3dIHTP.H the function named Eigen::VectorXd observe(volScalarField field) reads measurementsDict
        mesh,                       // inside 06enKFwDF_3dIHTP.H the function named Eigen::VectorXd observe(Eigen::VectorXd field) reads measurementsDict
        IOobject::MUST_READ,        // inside 06enKFwDF_3dIHTP.H the function named void readObservationTimes() reads observationDeltaTimesteps and observationStartTimestep
        IOobject::NO_WRITE
    )
);

List<vector> TCpos(thermocouplesPosition.lookup("positions"));
thermocouplesNum = TCpos.size();
thermocouplesPosition.lookup("timeSamplesT0") >> timeSamplesT0;                        // timeSamplesT0     1.0;
thermocouplesPosition.lookup("timeSamplesDeltaT") >> timeSamplesDeltaT;                // timeSamplesDeltaT 1.0;
timeSamplesNum = (runTime.endTime().value() - timeSamplesT0) / timeSamplesDeltaT + 1;  // 20, This is wrong. (It should be runTime.endTime().value()/Ntimes)/observationDeltaTimesteps
Info << "First sampling time = " << timeSamplesT0 << endl;                             // 1
Info << "Delta time between samples = " << timeSamplesDeltaT << endl;                  // 1
Info << "Number of samplings = " << timeSamplesNum << endl;                            // 20, This is wrong. (It should be runTime.endTime().value()/Ntimes)/observationDeltaTimesteps
Info << "Number of samplings should be calculated by (endtime / Ntimes) / observationDeltaTimesteps, so, the above calculated one is inncorect " << endl;



// timeSamplesT0     1.0;      First sampling time
// timeSamplesDeltaT 1.0;      Delta time between samples


// observationDeltaTimesteps 2;       inside 06enKFwDF_3dIHTP.H
// observationStartTimestep  1;       inside 06enKFwDF_3dIHTP.H

// Number of samples in time = 50     it is equal to equal measurements.cols() inside 06enKFwDF_3dIHTP.H,   Ntimes/observationDeltaTimesteps
// Number of sample in space = 100    it is equal to equal measurements.rows() inside 06enKFwDF_3dIHTP.H,   number of thermocouples in thermocouplesDict."

