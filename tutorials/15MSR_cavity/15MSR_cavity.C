#include "usmsrProblem.H"
#include "ReducedUnsteadyMSR.H"
#include "ITHACAstream.H"
#include "LRSensitivity.H"
#include "Tm.H"
#include "Ptot.H"
#include <Eigen/Dense>
#include <math.h>
#include <iomanip>
#include <chrono>
#include <ctime>
#include "Tm_time.H"
#include "Ptot_time.H"

class msr : public usmsrProblem
{
    public:
        explicit msr(int argc, char* argv[])
            :
            usmsrProblem(argc, argv),
            U(_U()),
            p(_p()),
            flux(_flux()),
            prec1(_prec1()),
            prec2(_prec2()),
            prec3(_prec3()),
            prec4(_prec4()),
            prec5(_prec5()),
            prec6(_prec6()),
            prec7(_prec7()),
            prec8(_prec8()),
            T(_T()),
            dec1(_dec1()),
            dec2(_dec2()),
            dec3(_dec3()),
            v(_v()),
            D(_D()),
            NSF(_NSF()),
            A(_A()),
            SP(_SP()),
            TXS(_TXS())
        {}

        volVectorField& U;
        volScalarField& p;
        volScalarField& flux;
        volScalarField& prec1;
        volScalarField& prec2;
        volScalarField& prec3;
        volScalarField& prec4;
        volScalarField& prec5;
        volScalarField& prec6;
        volScalarField& prec7;
        volScalarField& prec8;
        volScalarField& T;
        volScalarField& dec1;
        volScalarField& dec2;
        volScalarField& dec3;
        volScalarField& v;
        volScalarField& D;
        volScalarField& NSF;
        volScalarField& A;
        volScalarField& SP;
        volScalarField& TXS;

        ITHACAparameters para0;
        int NUproj = para0.ITHACAdict->lookupOrDefault<int>("NUproj", 2);
        int NPproj = para0.ITHACAdict->lookupOrDefault<int>("NPproj", 2);
        int NFluxproj = para0.ITHACAdict->lookupOrDefault<int>("NFluxproj", 2);
        int NPrecproj1 = para0.ITHACAdict->lookupOrDefault<int>("NPrecproj1", 2);
        int NPrecproj2 = para0.ITHACAdict->lookupOrDefault<int>("NPrecproj2", 2);
        int NPrecproj3 = para0.ITHACAdict->lookupOrDefault<int>("NPrecproj3", 2);
        int NPrecproj4 = para0.ITHACAdict->lookupOrDefault<int>("NPrecproj4", 2);
        int NPrecproj5 = para0.ITHACAdict->lookupOrDefault<int>("NPrecproj5", 2);
        int NPrecproj6 = para0.ITHACAdict->lookupOrDefault<int>("NPrecproj6", 2);
        int NPrecproj7 = para0.ITHACAdict->lookupOrDefault<int>("NPrecproj7", 2);
        int NPrecproj8 = para0.ITHACAdict->lookupOrDefault<int>("NPrecproj8", 2);
        int NTproj = para0.ITHACAdict->lookupOrDefault<int>("NTproj", 2);
        int NDecproj1 = para0.ITHACAdict->lookupOrDefault<int>("NDecproj1", 2);
        int NDecproj2 = para0.ITHACAdict->lookupOrDefault<int>("NDecproj2", 2);
        int NDecproj3 = para0.ITHACAdict->lookupOrDefault<int>("NDecproj3", 2);
        int NCproj = para0.ITHACAdict->lookupOrDefault<int>("NCproj", 2);

        double r1 = _beta1().value() / _betaTot().value();
        double r2 = _beta2().value() / _betaTot().value();
        double r3 = _beta3().value() / _betaTot().value();
        double r4 = _beta4().value() / _betaTot().value();
        double r5 = _beta5().value() / _betaTot().value();
        double r6 = _beta6().value() / _betaTot().value();
        double r7 = _beta7().value() / _betaTot().value();
        double r8 = _beta8().value() / _betaTot().value();

        void offlineSolve()
        {
            if (offline == false)
            {
                List<scalar> mu_now(3);

                for (int i = 0; i < mu.cols(); i++)
                {
                    _nu().value() = mu(0, i);
                    change_viscosity(mu(0, i));
                    _betaTot().value() = mu(1, i);
                    _beta1().value() = r1 * mu(1, i);
                    _beta2().value() = r2 * mu(1, i);
                    _beta3().value() = r3 * mu(1, i);
                    _beta4().value() = r4 * mu(1, i);
                    _beta5().value() = r5 * mu(1, i);
                    _beta6().value() = r6 * mu(1, i);
                    _beta7().value() = r7 * mu(1, i);
                    _beta8().value() = r8 * mu(1, i);
                    _decLam3().value() = mu(2, i);
                    mu_now[0] = mu(0, i);
                    mu_now[1] = mu(1, i);
                    mu_now[2] = mu(2, i);
                    truthSolve(mu_now);
                    restart();
                }
            }
            else
            {
                readMSRfields();
            }
        }

};

class Tmlocal : public Tm
{
    public:
        explicit Tmlocal(int argc, char* argv[], int Nsampled)
            :
            Tm(argc, argv, Nsampled),
            T(_T())
        {}

        volScalarField& T;
};

class Ploct: public Ptot_time
{
    public:
        explicit Ploct(int argc, char* argv[], int Nsampled)
            :
            Ptot_time(argc, argv, Nsampled),
            powerDens(_powerDens())
        {}

        volScalarField& powerDens;
};
class Plocal: public Ptot
{
    public:
        explicit Plocal(int argc, char* argv[], int Nsampled)
            :
            Ptot(argc, argv, Nsampled),
            powerDens(_powerDens())
        {}

        volScalarField& powerDens;
};
class Tmloct: public Tm_time
{
    public:
        explicit Tmloct(int argc, char* argv[], int Nsampled)
            :
            Tm_time(argc, argv, Nsampled),
            T(_T())
        {}

        volScalarField& T;
};

int main(int argc, char* argv[])
{
    //object initialization
    msr prova(argc, argv);
    //inletIndex set
    prova.inletIndex.resize(1, 2);
    prova.inletIndex(0, 0) = 0;
    prova.inletIndex(0, 1) = 0;
    prova.inletIndexT.resize(4, 1);
    prova.inletIndexT(0, 0) = 0;
    prova.inletIndexT(1, 0) = 1;
    prova.inletIndexT(2, 0) = 2;
    prova.inletIndexT(3, 0) = 3;
    //set the number of parameters and values of each
    prova.Pnumber = 3;
    prova.Tnumber = 31;
    prova.setParameters();
    int central = static_cast<int>(prova.Tnumber / 2);
    double nu0 = 2.46E-06;
    double betatot0 = 321.8E-05;
    double dlam30 = 3.58e-04;
    //+-10% from default values
    prova.mu_range(0, 0) = 1.1 * nu0; //nu range
    prova.mu_range(0, 1) = 0.9 * nu0;
    prova.mu_range(1, 0) = 1.1 * betatot0; //betaTot range
    prova.mu_range(1, 1) = 0.9 * betatot0;
    prova.mu_range(2, 0) = 1.1 * dlam30; //decLam3 range
    prova.mu_range(2, 1) = 0.9 * dlam30;
    prova.genRandPar();
    //override central column of mu with reference values of the parameters
    prova.mu(0, central) = nu0;
    prova.mu(1, central) = betatot0;
    prova.mu(2, central) = dlam30;
    double tstart = std::time(0);
    Eigen::MatrixXd arange(prova.Pnumber, 2);
    Eigen::MatrixXd mu_f = prova.mu;

    if (prova.offline == false)
    {
        for (int i = 0; i < prova.Pnumber; i++)
        {
            arange(i, 0) = prova.mu.row(i).minCoeff();
            arange(i, 1) = prova.mu.row(i).maxCoeff();
        }

        ITHACAstream::exportMatrix(arange, "range", "eigen", "./ITHACAoutput");
        ITHACAstream::exportMatrix(prova.mu, "mu_off", "eigen", "./ITHACAoutput");
    }
    else
    {
        mu_f = ITHACAstream::readMatrix("./ITHACAoutput/mu_off_mat.txt");
        arange = ITHACAstream::readMatrix("./ITHACAoutput/range_mat.txt");
    }

    //perform the offline step
    prova.offlineSolve();
    prova._nu().value() = nu0;
    prova.change_viscosity(nu0);
    double tend = std::time(0);
    double deltat_offline;
    deltat_offline = difftime(tend, tstart);

    //save the time needed to perform offline stage
    if (prova.offline == false)
    {
        std::ofstream diff_t("./ITHACAoutput/t_offline=" + name(deltat_offline));
    }

    std::cout << "Offline stage time= " << deltat_offline << std::endl;
    //compute the liftfunction for the velocity
    prova.liftSolve();
    //compute the liftfuncion for the temperature
    prova.liftSolveT();
    //homogenize the velocity field
    prova.homogenizeU();
    //homogenize the temperature field
    prova.homogenizeT();
    //compute the spatial modes
    prova.msrgetModesEVD();
    // perform the projection
    Eigen::VectorXi NPrec(8);
    NPrec(0) = prova.NPrecproj1;
    NPrec(1) = prova.NPrecproj2;
    NPrec(2) = prova.NPrecproj3;
    NPrec(3) = prova.NPrecproj4;
    NPrec(4) = prova.NPrecproj5;
    NPrec(5) = prova.NPrecproj6;
    NPrec(6) = prova.NPrecproj7;
    NPrec(7) = prova.NPrecproj8;
    Eigen::VectorXi NDec(3);
    NDec(0) = prova.NDecproj1;
    NDec(1) = prova.NDecproj2;
    NDec(2) = prova.NDecproj3;
    prova.projectPPE("./Matrices", prova.NUproj, prova.NPproj, prova.NFluxproj,
                     NPrec, prova.NTproj, NDec, prova.NCproj);
    //define the moving wall bc condition for U
    double umedio = 2.46E-03;
    //define the temperature at the moving wall
    Eigen::MatrixXd Twall(4, 2);
    Twall.setZero();
    Twall(0, 0) = 850;
    Twall(1, 0) = 950;
    Twall(2, 0) = 950;
    Twall(3, 0) = 950;
    //define the ROM object and its time settings
    reducedusMSR ridotto(prova);
    ridotto.tstart = 0;
    ridotto.finalTime = 10;
    ridotto.dt = 0.1;
    //vel_now for the ROM object
    Eigen::MatrixXd vel_now(1, 1);
    vel_now(0, 0) = umedio;
    //temp_now for the ROM object
    Eigen::MatrixXd temp_now = Twall;
    // set the online value of mu equal to mu(:,central)
    Eigen::VectorXd mu_on(3);
    ridotto.nu = mu_f(0, central);
    ridotto.btot = mu_f(1, central);
    ridotto.b1 = prova.r1 * ridotto.btot;
    ridotto.b2 = prova.r2 * ridotto.btot;
    ridotto.b3 = prova.r3 * ridotto.btot;
    ridotto.b4 = prova.r4 * ridotto.btot;
    ridotto.b5 = prova.r5 * ridotto.btot;
    ridotto.b6 = prova.r6 * ridotto.btot;
    ridotto.b7 = prova.r7 * ridotto.btot;
    ridotto.b8 = prova.r8 * ridotto.btot;
    ridotto.dl3 = mu_f(2, central);
    // index corresponding at central solution
    int ref_start = 465;
    mu_on(0) = ridotto.nu;
    mu_on(1) = ridotto.btot;
    mu_on(2) = ridotto.dl3;
    ridotto.solveOnline(vel_now, temp_now, mu_on, ref_start);
    ridotto.reconstructAP("./ITHACAoutput/Reconstruction/", 10);
    //ridotto.reconstructAP("./ITHACAoutput/Reconstruction/",50);
    //final snapshots for ROM and FOM respectively
    int tf = 21;
    int tf_fom = 31;
    //Some post process calculations
    int c = 0;
    Eigen::MatrixXd err(10, tf);
    Eigen::MatrixXd TmFOM(prova.Tnumber, tf);
    Eigen::MatrixXd PtotFOM(prova.Tnumber, tf);
    Eigen::MatrixXd TmROM(prova.Tnumber, tf);
    Eigen::MatrixXd PtotROM(prova.Tnumber, tf);
    std::cout << "Computing errors on reference case..." << std::endl;

    for (int i = ref_start; i < ref_start + tf; i++)
    {
        err(0, c) = ITHACAutilities::error_fields(prova.Ufield[i], ridotto.UREC[c]);
        err(1, c) = ITHACAutilities::error_fields(prova.Fluxfield[i],
                    ridotto.FLUXREC[c]);
        err(2, c) = ITHACAutilities::error_fields(prova.Prec1field[i],
                    ridotto.PREC1REC[c]);
        err(3, c) = ITHACAutilities::error_fields(prova.Prec4field[i],
                    ridotto.PREC4REC[c]);
        err(4, c) = ITHACAutilities::error_fields(prova.Prec8field[i],
                    ridotto.PREC8REC[c]);
        err(5, c) = ITHACAutilities::error_fields(prova.Tfield[i], ridotto.TREC[c]);
        err(6, c) = ITHACAutilities::error_fields(prova.Dec1field[i],
                    ridotto.DEC1REC[c]);
        err(7, c) = ITHACAutilities::error_fields(prova.Dec3field[i],
                    ridotto.DEC3REC[c]);
        err(8, c) = ITHACAutilities::error_fields(prova.PowerDensfield[i],
                    ridotto.POWERDENSREC[c]);
        err(9, c) = ITHACAutilities::error_fields(prova.TXSFields[i],
                    ridotto.TXSREC[c]);
        c++;
    }

    std::cout << "End" << std::endl;
    ITHACAstream::exportMatrix(err, "err", "eigen",
                               "./ITHACAoutput");
    std::cout << "computing FOM output..." << std::endl;
    double tmpt;
    double tmpp;
    int rig = 0;
    int col = 0;

    for (int i = 0; i < prova.Tfield.size(); i++)
    {
        tmpt = prova.Tfield[i].weightedAverage(prova._mesh().V()).value();
        tmpp = fvc::domainIntegrate(prova.PowerDensfield[i]).value();

        if (col < tf)
        {
            TmFOM(rig, col) = tmpt;
            PtotFOM(rig, col) = tmpp;
            col++;
        }
        else if (col == tf_fom)
        {
            col = 0;
            rig++;
        }
        else
        {
            col++;
        }
    }

    std::cout << "End" << std::endl;
    ITHACAstream::exportMatrix(TmFOM, "fom_resT", "eigen",
                               "./ITHACAoutput");
    ITHACAstream::exportMatrix(PtotFOM, "fom_resP", "eigen",
                               "./ITHACAoutput");
    // perform online solve for all the values of mu_offline
    std::string folder = {"./ITHACAoutput/ROMOutput/"};

    for (int i = 0; i < prova.Tnumber; i++)
    {
        ridotto.nu = mu_f(0, i);
        ridotto.btot = mu_f(1, i);
        ridotto.b1 = prova.r1 * ridotto.btot;
        ridotto.b2 = prova.r2 * ridotto.btot;
        ridotto.b3 = prova.r3 * ridotto.btot;
        ridotto.b4 = prova.r4 * ridotto.btot;
        ridotto.b5 = prova.r5 * ridotto.btot;
        ridotto.b6 = prova.r6 * ridotto.btot;
        ridotto.b7 = prova.r7 * ridotto.btot;
        ridotto.b8 = prova.r8 * ridotto.btot;
        ridotto.dl3 = mu_f(2, i);
        mu_on(0) = ridotto.nu;
        mu_on(1) = ridotto.btot;
        mu_on(2) = ridotto.dl3;
        folder.append(std::to_string(i));
        ridotto.solveOnline(vel_now, temp_now, mu_on, ref_start);
        ridotto.reconstructAP(folder, 50);
        folder = {"./ITHACAoutput/ROMOutput/"};
        ridotto.clearFields();
    }

    //compute the figures of merit for every "snapshot" time
    Tmloct Tmedia(argc, argv, prova.Tnumber);
    Ploct Ptotale(argc, argv, prova.Tnumber);
    c = 0;

    for (int i = 0; i < tf; i++)
    {
        Tmedia.buildMO(folder, i);
        Ptotale.buildMO(folder, i);
        TmROM.col(c) = Tmedia.modelOutput;
        PtotROM.col(c) = Ptotale.modelOutput;
        c++;
    }

    ITHACAstream::exportMatrix(TmROM, "rom_resT", "eigen",
                               "./ITHACAoutput");
    ITHACAstream::exportMatrix(PtotROM, "rom_resP", "eigen",
                               "./ITHACAoutput");
    //Sensitivity Analysis, two figures of merit considered:
    // average temperature and total power at the last time instant, i.e rom.finalTime
    //define the folder where save the output
    folder = {"./ITHACAoutput/Sensitivity/ModelOutput/"};
    int Nparameters = prova.Pnumber;
    std::vector<std::string> pdflist = {"normal", "normal", "normal"};
    int samplingPoints = 1000;
    LRSensitivity analisi(Nparameters, samplingPoints);
    //define the training range, values sampled outside are rejected
    analisi.trainingRange = arange;
    //set the parameters of each distributions
    Eigen::MatrixXd Mp(Nparameters, 2);
    Mp(0, 0) = nu0;
    Mp(0, 1) = nu0 * 0.1 / 3;
    Mp(1, 0) = betatot0;
    Mp(1, 1) = betatot0 * 0.1 / 3;
    Mp(2, 0) = dlam30;
    Mp(2, 1) = dlam30 * 0.1 / 3;
    //build and export the sampling sets
    analisi.buildSamplingSet(pdflist, Mp);
    ITHACAstream::exportMatrix(analisi.MatX, "sampled", "eigen", "./ITHACAoutput");
    //compute parameters statistics
    analisi.getXstats();
    //perform the simulations with the reduced object
    int counter = 0;
    tstart = std::time(0);
    ridotto.clearFields();

    for (int j = 0; j < samplingPoints; j++)
    {
        ridotto.nu = analisi.MatX(j, 0);
        ridotto.btot = analisi.MatX(j, 1);
        ridotto.b1 = prova.r1 * ridotto.btot;
        ridotto.b2 = prova.r2 * ridotto.btot;
        ridotto.b3 = prova.r3 * ridotto.btot;
        ridotto.b4 = prova.r4 * ridotto.btot;
        ridotto.b5 = prova.r5 * ridotto.btot;
        ridotto.b6 = prova.r6 * ridotto.btot;
        ridotto.b7 = prova.r7 * ridotto.btot;
        ridotto.b8 = prova.r8 * ridotto.btot;
        ridotto.dl3 = analisi.MatX(j, 2);
        mu_on(0) = ridotto.nu;
        mu_on(1) = ridotto.btot;
        mu_on(2) = ridotto.dl3;
        folder.append(std::to_string(counter));
        ridotto.solveOnline(vel_now, temp_now, mu_on, ref_start);
        ridotto.reconstructAP(folder, 1000);
        folder = {"./ITHACAoutput/Sensitivity/ModelOutput/"};
        ridotto.clearFields();
        counter++;
    }

    tend = std::time(0);
    int Ntot = samplingPoints;
    double deltat_online = difftime(tend, tstart);
    // save SA online time
    std::ofstream diff_ton("./ITHACAoutput/t_onlineSA_" + name(Ntot) + "=" + name(
                               deltat_online));
    //initialize the figure of merit objects
    Tmlocal TmediaSA(argc, argv, Ntot);
    Plocal PtotSA(argc, argv, Ntot);
    //build the Model Output
    TmediaSA.buildMO(folder);
    PtotSA.buildMO(folder);
    //assign FofM to LRSensitivity object, first figure of merit considered is the average temperature
    analisi.M = TmediaSA;
    //load the model output in the LRSensitivity object
    analisi.load_output();
    //compute output statistics
    analisi.getYstat();
    std::cout << "Mean value and variance of the output:" << std::endl;
    std::cout << analisi.Ey << "\t" << analisi.Vy << std::endl;
    std::cout << "-------------" << std::endl;
    // compute linear regression coefficients
    analisi.getBetas();
    std::cout << "analisi regression coefficients:" << std::endl;
    std::cout << analisi.betas << std::endl;
    std::cout << "-------------" << std::endl;
    // assess the quality of the linear regression, i.e. R^2
    analisi.assessQuality();
    std::cout << "quality: " << analisi.QI << std::endl;
    // save the ROM output, linear output, statistics
    Eigen::MatrixXd saveyout(samplingPoints, 1);
    Eigen::MatrixXd saveymodel(samplingPoints, 1);
    saveyout.col(0) = analisi.y;
    saveymodel.col(0) = analisi.ylin;
    ITHACAstream::exportMatrix(saveyout, "T_out", "eigen", "./ITHACAoutput");
    ITHACAstream::exportMatrix(saveymodel, "T_model", "eigen", "./ITHACAoutput");
    Eigen::MatrixXd coeffsT(4, 3);
    coeffsT.setZero();
    coeffsT(0, 0) = analisi.Ey;
    coeffsT(0, 1) = analisi.Vy;
    coeffsT(0, 2) = analisi.QI;
    coeffsT.row(1) = analisi.EX;
    coeffsT.row(2) = analisi.VX;
    coeffsT.row(3) = analisi.betas;
    ITHACAstream::exportMatrix(coeffsT, "coeffsT", "eigen", "./ITHACAoutput");
    //assign FofM object to LRSensitivity object
    analisi.M = PtotSA;
    //load the model output in the LRSensitivity object, repeat all the previous step
    // for the total power
    analisi.load_output();
    analisi.getYstat();
    std::cout << "Mean value and variance of the output:" << std::endl;
    std::cout << analisi.Ey << "\t" << analisi.Vy << std::endl;
    std::cout << "-------------" << std::endl;
    analisi.getBetas();
    std::cout << "analisi regression coefficients:" << std::endl;
    std::cout << analisi.betas << std::endl;
    std::cout << "-------------" << std::endl;
    analisi.assessQuality();
    std::cout << "quality: " << analisi.QI << std::endl;
    Eigen::MatrixXd saveyoutP(samplingPoints, 1);
    Eigen::MatrixXd saveymodelP(samplingPoints, 1);
    saveyoutP.col(0) = analisi.y;
    saveymodelP.col(0) = analisi.ylin;
    ITHACAstream::exportMatrix(saveyoutP, "P_out", "eigen", "./ITHACAoutput");
    ITHACAstream::exportMatrix(saveymodelP, "P_model", "eigen", "./ITHACAoutput");
    Eigen::MatrixXd coeffsP(4, 3);
    coeffsP.setZero();
    coeffsP(0, 0) = analisi.Ey;
    coeffsP(0, 1) = analisi.Vy;
    coeffsP(0, 2) = analisi.QI;
    coeffsP.row(1) = analisi.EX;
    coeffsP.row(2) = analisi.VX;
    coeffsP.row(3) = analisi.betas;
    ITHACAstream::exportMatrix(coeffsP, "coeffsP", "eigen", "./ITHACAoutput");
    return 0;
}

