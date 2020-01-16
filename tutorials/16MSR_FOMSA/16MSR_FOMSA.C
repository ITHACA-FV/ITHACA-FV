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

        void offlineSolve(std::string dir)
        {
            if (offline == false)
            {
                List<scalar> mu_now(3);
                std::string folder = dir;

                for (int i = 0; i < mu.cols(); i++)
                {
                    folder.append(std::to_string(i));
                    folder.append("/");
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
                    truthSolve(mu_now, folder);
                    folder = dir;
                    restart();
                }
            }
            else
            {
                readMSRfields(dir);
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

int main(int argc, char* argv[])
{
    //object initialization
    msr prova(argc, argv);
    //inletIndex set
    prova.inletIndex.resize(1, 2);
    prova.inletIndex(0, 0) = 0;
    prova.inletIndex(0, 1) = 0;
    prova.Pnumber = 3;
    prova.Tnumber = 1000;
    prova.setParameters();
    // define the distribution central and std. deviation values
    double nu0 = 2.46E-06;
    double betatot0 = 321.8E-05;
    double dlam30 = 3.58e-04;
    double signu = 0.1 / 3 * nu0;
    double sigbeta = 0.1 / 3 * betatot0;
    double sigdlam3 = 0.1 / 3 * dlam30;
    // read range to adopt as training range
    prova.mu_range = ITHACAstream::readMatrix("./range_mat.txt");
    std::string dist = {"normal"};
    // sample the values using ITHACA::sampling
    prova.mu.row(0) = ITHACAsampling::samplingMC(dist, prova.mu_range(0, 0),
                      prova.mu_range(0, 1), nu0, signu, prova.Tnumber);
    prova.mu.row(1) = ITHACAsampling::samplingMC(dist, prova.mu_range(1, 0),
                      prova.mu_range(1, 1), betatot0, sigbeta, prova.Tnumber);
    prova.mu.row(2) = ITHACAsampling::samplingMC(dist, prova.mu_range(2, 0),
                      prova.mu_range(2, 1), dlam30, sigdlam3, prova.Tnumber);
    double tstart = std::time(0);
    Eigen::MatrixXd mu_f = prova.mu;

    if (prova.offline == false)
    {
        ITHACAstream::exportMatrix(prova.mu, "mu_off", "eigen", "./ITHACAoutput");
    }
    else
    {
        mu_f = ITHACAstream::readMatrix("./ITHACAoutput/mu_off_mat.txt");
    }

    std::string folder = {"./ITHACAoutput/FOMoutput/"};
    //perform the offline step
    prova.offlineSolve(folder);
    double tend = std::time(0);
    double deltat_offline;
    deltat_offline = difftime(tend, tstart);

    if (prova.offline == false)
    {
        std::ofstream diff_t("./ITHACAoutput/t_SA=" + name(deltat_offline));
    }

    std::cout << "Offline stage time= " << deltat_offline << std::endl;
    //sensitivity analysis
    int Nparameters = prova.Pnumber;
    int samplingPoints = prova.Tnumber;
    LRSensitivity analisi(Nparameters, samplingPoints);
    //load the sampling set
    analisi.MatX.col(0) = mu_f.row(0);
    analisi.MatX.col(1) = mu_f.row(1);
    analisi.MatX.col(2) = mu_f.row(2);
    analisi.getXstats();
    //initialize Tm object
    Tmlocal TmediaSA(argc, argv, samplingPoints);
    Plocal PtotSA(argc, argv, samplingPoints);
    //build the Model Output
    TmediaSA.buildMO(folder);
    PtotSA.buildMO(folder);
    //assign it to VDSensitivity object
    analisi.M = TmediaSA;
    //load the model output in the VDSensitivity object
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
    //assign it to VDSensitivity object
    analisi.M = PtotSA;
    //load the model output in the VDSensitivity object
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

