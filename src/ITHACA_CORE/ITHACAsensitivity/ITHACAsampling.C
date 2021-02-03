#include "ITHACAsampling.H"

std::vector<std::string> ITHACAsampling::distributions = {"UNIFORM", "NORMAL", "POISSON", "EXPONENTIAL"};

Eigen::VectorXd ITHACAsampling::samplingMC(std::string pdftype, double& lowerE,
        double& upperE, double& distpara1, double& distpara2, label& Npoints)
{
    std::random_device rd;
    std::mt19937 generator(rd());

    //to make it non-case sensitive
    for (label i = 0; i < pdftype.size(); i++)
    {
        pdftype[i] = toupper(pdftype[i]);
    }

    label pos = -1;
    bool pdf_found = false;
    Eigen::VectorXd samplingVector;
    samplingVector.resize(Npoints);
    std::uniform_real_distribution<> dist0(distpara1, distpara2);
    std::normal_distribution<> dist1(distpara1, distpara2);
    std::poisson_distribution<> dist2(distpara1);
    std::exponential_distribution<> dist3(distpara1);

    for (label j = 0; j < distributions.size(); j++)
    {
        if (pdftype == distributions[j])
        {
            pdf_found = true;
            pos = j;
        }
    }

    double random;
    label p_counter = 0;

    if (pdf_found == true)
    {
        while (p_counter < samplingVector.size())
        {
            switch (pos)
            {
                case 0:
                    random = dist0(generator);
                    break;

                case 1:
                    random = dist1(generator);
                    break;

                case 2:
                    random = dist2(generator);
                    break;

                case 3:
                    random = dist3(generator);
                    break;
            }

            if (random >= lowerE && random <= upperE)
            {
                samplingVector(p_counter) = random;
                p_counter++;
            }
        }
    }
    else
    {
        std::cout << "pdf '" << pdftype << "' not implemented, programm aborted" <<
                  std::endl;
        exit(0);
    }

    return samplingVector;
}
