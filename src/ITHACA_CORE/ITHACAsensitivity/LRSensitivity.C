#include "LRSensitivity.H"


using std::vector;
using std::string;

// Constructors
LRSensitivity::LRSensitivity() {};

LRSensitivity::LRSensitivity(label Npara, label Np)
{
    No_parameters = Npara;
    Npoints = Np;
    setAll();
}


void LRSensitivity::buildSamplingSet(std::vector<std::string>& pdflist,
                                     Eigen::MatrixXd plist)
{
    for (label i = 0; i < No_parameters; i++)
    {
        MatX.col(i) = ITHACAsampling::samplingMC(pdflist[i], trainingRange(i, 0),
                      trainingRange(i, 1), plist(i, 0), plist(i, 1), Npoints);
    }
}


void LRSensitivity::setAll()
{
    MatX.resize(Npoints, No_parameters);
    MatX.setZero();
    MatXn.resize(Npoints, No_parameters);
    MatXn.setZero();
    trainingRange.resize(No_parameters, 2);
    trainingRange.setZero();
    y.resize(Npoints);
    y.setZero();
    yn.resize(Npoints);
    yn.setZero();
    ylin.resize(Npoints);
    ylin.setZero();
    EX.resize(No_parameters);
    EX.setZero();
    VX.resize(No_parameters);
    VX.setZero();
    betas.resize(No_parameters);
    betas.setZero();
}

void LRSensitivity::load_output()
{
    if (M.MObuilt == true)
    {
        y = M.modelOutput;
    }
    else
    {
        std::cout << "Model output not computed or assigned, program aborted" <<
                  std::endl;
        exit(0);
    }
}

void LRSensitivity::getYstat()
{
    Ey = 0;
    Vy = 0;
    Ey = y.sum() / Npoints;

    for (label i = 0; i < Npoints; i++)
    {
        Vy += (y(i) - Ey) * (y(i) - Ey);
    }

    Vy = Vy / (Npoints - 1);
    Ydone = true;
}

void LRSensitivity::getXstats()
{
    for (label i = 0; i < No_parameters; i++)
    {
        EX(i) = MatX.col(i).sum() / Npoints;
        VX(i) = 0;

        for (label j = 0; j < Npoints; j++)
        {
            VX(i) += (MatX(j, i) - EX(i)) * (MatX(j, i) - EX(i));
        }

        VX(i) = VX(i) / (Npoints - 1);
    }

    Xdone = true;
}


void LRSensitivity::getBetas()
{
    if (Ydone == true && Xdone == true)
    {
        //Normalize independent variables and output
        for (label j = 0; j < Npoints; j++)
        {
            yn(j) = (y(j) - Ey) / std::sqrt(Vy);

            for (label i = 0; i < No_parameters; i++)
            {
                MatXn(j, i) = (MatX(j, i) - EX(i)) / std::sqrt(VX(i));
            }
        }

        Eigen::MatrixXd A = MatXn.transpose() * MatXn;
        Eigen::VectorXd sol = MatXn.transpose() * yn;
        betas = A.colPivHouseholderQr().solve(sol);
        bdone = true;
    }
    else
    {
        std::cout <<
                  "Statistics about inputs or output are not computed yet, nothing to do ..." <<
                  std::endl;
    }
}

void LRSensitivity::assessQuality()
{
    if (bdone == true)
    {
        QI = 0;
        double  N = 0;
        double  D = 0;

        for (label j = 0; j < Npoints; j++)
        {
            ylin(j) = Ey;

            for (label i = 0; i < No_parameters; i++)
            {
                ylin(j) += betas(i) * MatXn(j, i) * std::sqrt(Vy);
            }

            N += (ylin(j) - Ey) * (ylin(j) - Ey);
            D += (y(j) - Ey) * (y(j) - Ey);
        }

        QI = N / D;
    }
    else
    {
        std::cout <<
                  "Linear regression coefficients are not computed yet, nothing to do ..." <<
                  std::endl;
    }
}



