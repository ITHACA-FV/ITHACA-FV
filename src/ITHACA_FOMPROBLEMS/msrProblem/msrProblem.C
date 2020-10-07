#include "msrProblem.H"


// Constructor
msrProblem::msrProblem() {}

msrProblem::msrProblem(int argc, char* argv[])
{
#include "setRootCase.H"
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
#include "createFields_Neutronics.H"
#include "createFields_Thermal.H"
#include "createConstants.H"
#include "createFvOptions.H"
    turbulence->validate();
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
    tolerance = ITHACAdict->lookupOrDefault<scalar>("tolerance", 1e-5);
    maxIter = ITHACAdict->lookupOrDefault<scalar>("maxIter", 1000);
    offline = ITHACAutilities::check_off();
    podex = ITHACAutilities::check_pod();
}

void msrProblem::truthSolve(List<scalar> mu_now)
{
    Time& runTime = _runTime();
    fvMesh& mesh = _mesh();
    volScalarField& p = _p();
    volVectorField& U = _U();
    surfaceScalarField& phi = _phi();
    fv::options& fvOptions = _fvOptions();
    simpleControl& simple = _simple();
    IOMRFZoneList& MRF = _MRF();
    singlePhaseTransportModel& laminarTransport = _laminarTransport();
    dimensionedScalar& IV1 = _IV1();
    dimensionedScalar& Keff = _Keff();
    volScalarField& flux = _flux();
    volScalarField flux_old = _flux();
    dimensionedScalar& betaTot = _betaTot();
    volScalarField& SP = _SP();
    dimensionedScalar& SP1_0 = _SP1_0();
    dimensionedScalar& alfa_SP1 = _alfa_SP1();
    volScalarField& D = _D();
    dimensionedScalar& D1_0 = _D1_0();
    dimensionedScalar& alfa_D1 = _alfa_D1();
    volScalarField& NSF = _NSF();
    dimensionedScalar& NSF1_0 = _NSF1_0();
    dimensionedScalar& alfa_NSF1 = _alfa_NSF1();
    volScalarField& A = _A();
    dimensionedScalar& A1_0 = _A1_0();
    dimensionedScalar& alfa_A1 = _alfa_A1();
    volScalarField& prec1 = _prec1();
    dimensionedScalar& Sc = _Sc();
    dimensionedScalar& Sct = _Sct();
    dimensionedScalar& lam1 = _lam1();
    dimensionedScalar& beta1 = _beta1();
    volScalarField& prec2 = _prec2();
    dimensionedScalar& lam2 = _lam2();
    dimensionedScalar& beta2 = _beta2();
    volScalarField& prec3 = _prec3();
    dimensionedScalar& lam3 = _lam3();
    dimensionedScalar& beta3 = _beta3();
    volScalarField& prec4 = _prec4();
    dimensionedScalar& lam4 = _lam4();
    dimensionedScalar& beta4 = _beta4();
    volScalarField& prec5 = _prec5();
    dimensionedScalar& lam5 = _lam5();
    dimensionedScalar& beta5 = _beta5();
    volScalarField& prec6 = _prec6();
    dimensionedScalar& lam6 = _lam6();
    dimensionedScalar& beta6 = _beta6();
    volScalarField& prec7 = _prec7();
    dimensionedScalar& lam7 = _lam7();
    dimensionedScalar& beta7 = _beta7();
    volScalarField& prec8 = _prec8();
    dimensionedScalar& lam8 = _lam8();
    dimensionedScalar& beta8 = _beta8();
    volScalarField& T = _T();
    dimensionedScalar& Pr = _Pr();
    dimensionedScalar& Prt = _Prt();
    volScalarField& dec1 = _dec1();
    dimensionedScalar& decLam1 = _decLam1();
    dimensionedScalar& decBeta1 = _decBeta1();
    volScalarField& dec2 = _dec2();
    dimensionedScalar& decLam2 = _decLam2();
    dimensionedScalar& decBeta2 = _decBeta2();
    volScalarField& dec3 = _dec3();
    dimensionedScalar& decLam3 = _decLam3();
    dimensionedScalar& decBeta3 = _decBeta3();
    dimensionedScalar& decbetaTot = _decbetaTot();
    dimensionedScalar& rhoRef = _rhoRef();
    dimensionedScalar& CpRef = _CpRef();
    volScalarField v = _v();
    volScalarField TXS = _TXS();
    dimensionedScalar& nu = _nu();
    dimensionedScalar& betaTE = _betaTE();
    dimensionedScalar& Tref = _Tref();
    dimensionedScalar& TrefXS = _TrefXS();
    volScalarField& logT = _logT();
    volScalarField& alphat = _alphat();
    volScalarField& difft = _difft();
    volScalarField powerDens = (1 - decbetaTot) * flux * SP +
                               (decLam1 * dec1 + decLam2 * dec2 + decLam3 * dec3);
    powerDens.rename("powerDens");
#include "NLsolve.H"
    counter++;
    writeMu(mu_now);
    // --- Fill in the mu_samples with parameters (mu) to be used for the PODI sample points
    mu_samples.conservativeResize(mu_samples.rows() + 1, mu_now.size());

    for (int i = 0; i < mu_now.size(); i++)
    {
        mu_samples(mu_samples.rows() - 1, i) = mu_now[i];
    }

    // Resize to Unitary if not initialized by user (i.e. non-parametric problem)
    if (mu.cols() == 0)
    {
        mu.resize(1, 1);
    }

    if (mu_samples.rows() == mu.cols())
    {
        ITHACAstream::exportMatrix(mu_samples, "mu_samples", "eigen",
                                   "./ITHACAoutput/Offline");
    }
}


// Get the modes function SVD
void msrProblem::msrgetModesSVD()
{
    int NU = para->ITHACAdict->lookupOrDefault<int>("NUout", 10);
    int NP = para->ITHACAdict->lookupOrDefault<int>("NPout", 10);
    int NF = para->ITHACAdict->lookupOrDefault<int>("NFluxout", 10);
    int NPrec1 = para->ITHACAdict->lookupOrDefault<int>("NPrecout1", 10);
    int NPrec2 = para->ITHACAdict->lookupOrDefault<int>("NPrecout2", 10);
    int NPrec3 = para->ITHACAdict->lookupOrDefault<int>("NPrecout3", 10);
    int NPrec4 = para->ITHACAdict->lookupOrDefault<int>("NPrecout4", 10);
    int NPrec5 = para->ITHACAdict->lookupOrDefault<int>("NPrecout5", 10);
    int NPrec6 = para->ITHACAdict->lookupOrDefault<int>("NPrecout6", 10);
    int NPrec7 = para->ITHACAdict->lookupOrDefault<int>("NPrecout7", 10);
    int NPrec8 = para->ITHACAdict->lookupOrDefault<int>("NPrecout8", 10);
    int NT = para->ITHACAdict->lookupOrDefault<int>("NTout", 10);
    int NDec1 = para->ITHACAdict->lookupOrDefault<int>("NDecout1", 10);
    int NDec2 = para->ITHACAdict->lookupOrDefault<int>("NDecout2", 10);
    int NDec3 = para->ITHACAdict->lookupOrDefault<int>("NDecout3", 10);
    int NC = para->ITHACAdict->lookupOrDefault<int>("NCout", 10);

    if (homboolU == true)
    {
        ITHACAPOD::getModesSVD(Uomfield, Umodes, _U().name(), podex, 0, 0, NU);
    }
    else
    {
        ITHACAPOD::getModesSVD(Ufield, Umodes, _U().name(), podex, 0, 0, NU);
    }

    ITHACAPOD::getModesSVD(Pfield, Pmodes, _p().name(), podex, 0, 0, NP);
    ITHACAPOD::getModesSVD(Fluxfield, Fluxmodes, _flux().name(), podex, 0, 0, NF);
    ITHACAPOD::getModesSVD(Prec1field, Prec1modes, _prec1().name(), podex, 0, 0,
                           NPrec1);
    ITHACAPOD::getModesSVD(Prec2field, Prec2modes, _prec2().name(), podex, 0, 0,
                           NPrec2);
    ITHACAPOD::getModesSVD(Prec3field, Prec3modes, _prec3().name(), podex, 0, 0,
                           NPrec3);
    ITHACAPOD::getModesSVD(Prec4field, Prec4modes, _prec4().name(), podex, 0, 0,
                           NPrec4);
    ITHACAPOD::getModesSVD(Prec5field, Prec5modes, _prec5().name(), podex, 0, 0,
                           NPrec5);
    ITHACAPOD::getModesSVD(Prec6field, Prec6modes, _prec6().name(), podex, 0, 0,
                           NPrec6);
    ITHACAPOD::getModesSVD(Prec7field, Prec7modes, _prec7().name(), podex, 0, 0,
                           NPrec7);
    ITHACAPOD::getModesSVD(Prec8field, Prec8modes, _prec8().name(), podex, 0, 0,
                           NPrec8);

    if (homboolT == true)
    {
        ITHACAPOD::getModesSVD(Tomfield, Tmodes, _T().name(), podex, 0, 0, NT);
    }
    else
    {
        ITHACAPOD::getModesSVD(Tfield, Tmodes, _T().name(), podex, 0, 0, NT);
    }

    ITHACAPOD::getModesSVD(Dec1field, Dec1modes, _dec1().name(), podex, 0, 0,
                           NDec1);
    ITHACAPOD::getModesSVD(Dec2field, Dec2modes, _dec2().name(), podex, 0, 0,
                           NDec2);
    ITHACAPOD::getModesSVD(Dec3field, Dec3modes, _dec3().name(), podex, 0, 0,
                           NDec3);
    ITHACAPOD::getModesSVD(vFields, vmodes, _v().name(), podex, 0, 0, NC);
    ITHACAPOD::getModesSVD(DFields, Dmodes, _D().name(), podex, 0, 0, NC);
    ITHACAPOD::getModesSVD(NSFFields, NSFmodes, _NSF().name(), podex, 0, 0, NC);
    ITHACAPOD::getModesSVD(AFields, Amodes, _A().name(), podex, 0, 0, NC);
    ITHACAPOD::getModesSVD(SPFields, SPmodes, _SP().name(), podex, 0, 0, NC);
    ITHACAPOD::getModesSVD(TXSFields, TXSmodes, _TXS().name(), podex, 0, 0, NC);
    Info << "End\n" << endl;
}

// Get the modes function Compliance matrix method
void msrProblem::msrgetModesEVD()
{
    int NU = para->ITHACAdict->lookupOrDefault<int>("NUout", 10);
    int NP = para->ITHACAdict->lookupOrDefault<int>("NPout", 10);
    int NF = para->ITHACAdict->lookupOrDefault<int>("NFluxout", 10);
    int NPrec1 = para->ITHACAdict->lookupOrDefault<int>("NPrecout1", 10);
    int NPrec2 = para->ITHACAdict->lookupOrDefault<int>("NPrecout2", 10);
    int NPrec3 = para->ITHACAdict->lookupOrDefault<int>("NPrecout3", 10);
    int NPrec4 = para->ITHACAdict->lookupOrDefault<int>("NPrecout4", 10);
    int NPrec5 = para->ITHACAdict->lookupOrDefault<int>("NPrecout5", 10);
    int NPrec6 = para->ITHACAdict->lookupOrDefault<int>("NPrecout6", 10);
    int NPrec7 = para->ITHACAdict->lookupOrDefault<int>("NPrecout7", 10);
    int NPrec8 = para->ITHACAdict->lookupOrDefault<int>("NPrecout8", 10);
    int NT = para->ITHACAdict->lookupOrDefault<int>("NTout", 10);
    int NDec1 = para->ITHACAdict->lookupOrDefault<int>("NDecout1", 10);
    int NDec2 = para->ITHACAdict->lookupOrDefault<int>("NDecout2", 10);
    int NDec3 = para->ITHACAdict->lookupOrDefault<int>("NDecout3", 10);
    int NC = para->ITHACAdict->lookupOrDefault<int>("NCout", 10);

    if (homboolU == true)
    {
        ITHACAPOD::getModes(Uomfield, Umodes, _U().name(), podex, 0, 0, NU);
    }
    else
    {
        ITHACAPOD::getModes(Ufield, Umodes, _U().name(), podex, 0, 0, NU);
    }

    ITHACAPOD::getModes(Pfield, Pmodes, _U().name(), podex, 0, 0, NP);
    ITHACAPOD::getModes(Fluxfield, Fluxmodes, _U().name(), podex, 0, 0, NF);
    ITHACAPOD::getModes(Prec1field, Prec1modes, _prec1().name(), podex, 0, 0,
                        NPrec1);
    ITHACAPOD::getModes(Prec2field, Prec2modes, _prec2().name(), podex, 0, 0,
                        NPrec2);
    ITHACAPOD::getModes(Prec3field, Prec3modes, _prec3().name(), podex, 0, 0,
                        NPrec3);
    ITHACAPOD::getModes(Prec4field, Prec4modes, _prec4().name(), podex, 0, 0,
                        NPrec4);
    ITHACAPOD::getModes(Prec5field, Prec5modes, _prec5().name(), podex, 0, 0,
                        NPrec5);
    ITHACAPOD::getModes(Prec6field, Prec6modes, _prec6().name(), podex, 0, 0,
                        NPrec6);
    ITHACAPOD::getModes(Prec7field, Prec7modes, _prec7().name(), podex, 0, 0,
                        NPrec7);
    ITHACAPOD::getModes(Prec8field, Prec8modes, _prec8().name(), podex, 0, 0,
                        NPrec8);

    if (homboolT == true)
    {
        ITHACAPOD::getModes(Tomfield, Tmodes, _T().name(), podex, 0, 0, NT);
    }
    else
    {
        ITHACAPOD::getModes(Tfield, Tmodes, _T().name(), podex, 0, 0, NT);
    }

    ITHACAPOD::getModes(Dec1field, Dec1modes, _dec1().name(), podex, 0, 0, NDec1);
    ITHACAPOD::getModes(Dec2field, Dec2modes, _dec2().name(), podex, 0, 0, NDec2);
    ITHACAPOD::getModes(Dec3field, Dec3modes, _dec3().name(), podex, 0, 0, NDec3);
    ITHACAPOD::getModes(vFields, vmodes, _v().name(), podex, 0, 0, NC);
    ITHACAPOD::getModes(DFields, Dmodes, _D().name(), podex, 0, 0, NC);
    ITHACAPOD::getModes(NSFFields, NSFmodes, _NSF().name(), podex, 0, 0, NC);
    ITHACAPOD::getModes(AFields, Amodes, _A().name(), podex, 0, 0, NC);
    ITHACAPOD::getModes(SPFields, SPmodes, _SP().name(), podex, 0, 0, NC);
    ITHACAPOD::getModes(TXSFields, TXSmodes, _TXS().name(), podex, 0, 0, NC);
    Info << "End\n" << endl;
}

void msrProblem::liftSolve()
{
    for (label k = 0; k < inletIndex.rows(); k++)
    {
        Time& runTime = _runTime();
        surfaceScalarField& phi = _phi();
        fvMesh& mesh = _mesh();
        volScalarField p = _p();
        volVectorField U = _U();
        IOMRFZoneList& MRF = _MRF();
        label BCind = inletIndex(k, 0);
        volVectorField Ulift("Ulift" + name(k), U);
        instantList Times = runTime.times();
        runTime.setTime(Times[1], 1);
        pisoControl potentialFlow(mesh, "potentialFlow");
        Info << "Solving a lifting Problem" << endl;
        Vector<double> v1(0, 0, 0);
        v1[inletIndex(k, 1)] = 1;
        Vector<double> v0(0, 0, 0);

        for (label j = 0; j < U.boundaryField().size(); j++)
        {
            if (j == BCind)
            {
                assignBC(Ulift, j, v1);
            }
            else if (U.boundaryField()[BCind].type() == "fixedValue")
            {
                assignBC(Ulift, j, v0);
            }
            else
            {
            }

            assignIF(Ulift, v0);
            phi = linearInterpolate(Ulift) & mesh.Sf();
        }

        Info << "Constructing velocity potential field Phi\n" << endl;
        volScalarField Phi
        (
            IOobject
            (
                "Phi",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("Phi", dimLength * dimVelocity, 0),
            p.boundaryField().types()
        );
        label PhiRefCell = 0;
        scalar PhiRefValue = 0;
        setRefCell
        (
            Phi,
            potentialFlow.dict(),
            PhiRefCell,
            PhiRefValue
        );
        mesh.setFluxRequired(Phi.name());
        runTime.functionObjects().start();
        MRF.makeRelative(phi);
        adjustPhi(phi, Ulift, p);

        while (potentialFlow.correctNonOrthogonal())
        {
            fvScalarMatrix PhiEqn
            (
                fvm::laplacian(dimensionedScalar("1", dimless, 1), Phi)
                ==
                fvc::div(phi)
            );
            PhiEqn.setReference(PhiRefCell, PhiRefValue);
            PhiEqn.solve();

            if (potentialFlow.finalNonOrthogonalIter())
            {
                phi -= PhiEqn.flux();
            }
        }

        MRF.makeAbsolute(phi);
        Info << "Continuity error = "
             << mag(fvc::div(phi))().weightedAverage(mesh.V()).value()
             << endl;
        Ulift = fvc::reconstruct(phi);
        Ulift.correctBoundaryConditions();
        Info << "Interpolated velocity error = "
             << (sqrt(sum(sqr((fvc::interpolate(U) & mesh.Sf()) - phi)))
                 / sum(mesh.magSf())).value()
             << endl;
        Ulift.write();
        liftfield.append(Ulift);
    }
}




// * * * * * * * * * * * * * * Projection Methods * * * * * * * * * * * * * * //

void msrProblem::projectPPE(fileName folder, label NU, label NP, label NF,
                            Eigen::VectorXi NPrec, label NT, Eigen::VectorXi NDec, label NC)
{
    if (NPrec.size() != 8 || NDec.size() != 3)
    {
        std::cout <<
                  "The model assumes 8 groups of precursors and 3 of decay heat, check NDrec and NDec dimensions..."
                  << std::endl;
        exit(0);
    }

    NUmodes = NU;
    NPmodes = NP;
    NFluxmodes = NF;
    NPrecmodes = NPrec;
    int NPrec1 = NPrecmodes(0);
    int NPrec2 = NPrecmodes(1);
    int NPrec3 = NPrecmodes(2);
    int NPrec4 = NPrecmodes(3);
    int NPrec5 = NPrecmodes(4);
    int NPrec6 = NPrecmodes(5);
    int NPrec7 = NPrecmodes(6);
    int NPrec8 = NPrecmodes(7);
    NTmodes = NT;
    NDecmodes = NDec;
    int NDec1 = NDecmodes(0);
    int NDec2 = NDecmodes(1);
    int NDec3 = NDecmodes(2);
    NCmodes = NC;
    Info << "\n Computing fluid-dynamics matrices\n" << endl;
    B_matrix = diffusive_term(NUmodes, NPmodes);
    C_matrix = convective_term(NUmodes, NPmodes);
    M_matrix = mass_term(NUmodes, NPmodes);
    K_matrix = pressure_gradient_term(NUmodes, NPmodes);
    D_matrix = laplacian_pressure(NPmodes);
    G_matrix = div_momentum(NUmodes, NPmodes);
    BC1_matrix = pressure_BC1(NUmodes, NPmodes);
    BC2_matrix = pressure_BC2(NUmodes, NPmodes);
    BC3_matrix = pressure_BC3(NUmodes, NPmodes);
    Info << "\n End \n" << endl;
    Info << "\n Computing neutronics matrices\n" << endl;
    LF_matrix = laplacian_flux(NFluxmodes, NCmodes);
    MF_matrix = mass_flux(NFluxmodes);
    PF_matrix = prod_flux(NFluxmodes, NCmodes);
    AF_matrix = abs_flux(NFluxmodes, NCmodes);
    PS1_matrix = prec_source(NFluxmodes, NPrec1, 1);
    PS2_matrix = prec_source(NFluxmodes, NPrec2, 2);
    PS3_matrix = prec_source(NFluxmodes, NPrec3, 3);
    PS4_matrix = prec_source(NFluxmodes, NPrec4, 4);
    PS5_matrix = prec_source(NFluxmodes, NPrec5, 5);
    PS6_matrix = prec_source(NFluxmodes, NPrec6, 6);
    PS7_matrix = prec_source(NFluxmodes, NPrec7, 7);
    PS8_matrix = prec_source(NFluxmodes, NPrec8, 8);
    ST1_matrix = stream_term(NUmodes, NPrec1, 1);
    ST2_matrix = stream_term(NUmodes, NPrec2, 2);
    ST3_matrix = stream_term(NUmodes, NPrec3, 3);
    ST4_matrix = stream_term(NUmodes, NPrec4, 4);
    ST5_matrix = stream_term(NUmodes, NPrec5, 5);
    ST6_matrix = stream_term(NUmodes, NPrec6, 6);
    ST7_matrix = stream_term(NUmodes, NPrec7, 7);
    ST8_matrix = stream_term(NUmodes, NPrec8, 8);
    MP1_matrix = prec_mass(NPrec1, 1);
    MP2_matrix = prec_mass(NPrec2, 2);
    MP3_matrix = prec_mass(NPrec3, 3);
    MP4_matrix = prec_mass(NPrec4, 4);
    MP5_matrix = prec_mass(NPrec5, 5);
    MP6_matrix = prec_mass(NPrec6, 6);
    MP7_matrix = prec_mass(NPrec7, 7);
    MP8_matrix = prec_mass(NPrec8, 8);
    LP1_matrix = laplacian_prec(NPrec1, 1);
    LP2_matrix = laplacian_prec(NPrec2, 2);
    LP3_matrix = laplacian_prec(NPrec3, 3);
    LP4_matrix = laplacian_prec(NPrec4, 4);
    LP5_matrix = laplacian_prec(NPrec5, 5);
    LP6_matrix = laplacian_prec(NPrec6, 6);
    LP7_matrix = laplacian_prec(NPrec7, 7);
    LP8_matrix = laplacian_prec(NPrec8, 8);
    FS1_matrix = flux_source(NFluxmodes, NPrec1, NCmodes, 1);
    FS2_matrix = flux_source(NFluxmodes, NPrec2, NCmodes, 2);
    FS3_matrix = flux_source(NFluxmodes, NPrec3, NCmodes, 3);
    FS4_matrix = flux_source(NFluxmodes, NPrec4, NCmodes, 4);
    FS5_matrix = flux_source(NFluxmodes, NPrec5, NCmodes, 5);
    FS6_matrix = flux_source(NFluxmodes, NPrec6, NCmodes, 6);
    FS7_matrix = flux_source(NFluxmodes, NPrec7, NCmodes, 7);
    FS8_matrix = flux_source(NFluxmodes, NPrec8, NCmodes, 8);
    Info << "\n End \n" << endl;
    Info << "\n Computing thermal matrices\n" << endl;
    SD1_matrix = stream_dec(NUmodes, NDec1, 1);
    SD2_matrix = stream_dec(NUmodes, NDec2, 2);
    SD3_matrix = stream_dec(NUmodes, NDec3, 3);
    MD1_matrix = dec_mass(NDec1, 1);
    MD2_matrix = dec_mass(NDec2, 2);
    MD3_matrix = dec_mass(NDec3, 3);
    LD1_matrix = laplacian_dec(NDec1, 1);
    LD2_matrix = laplacian_dec(NDec2, 2);
    LD3_matrix = laplacian_dec(NDec3, 3);
    DFS1_matrix = dec_fluxsource(NFluxmodes, NDec1, NCmodes, 1);
    DFS2_matrix = dec_fluxsource(NFluxmodes, NDec2, NCmodes, 2);
    DFS3_matrix = dec_fluxsource(NFluxmodes, NDec3, NCmodes, 3);
    TM_matrix = mass_temp(NTmodes);
    TS_matrix = temp_stream(NUmodes, NTmodes);
    LT_matrix = laplacian_temp(NTmodes);
    TXS_matrix = temp_XSfluxsource(NTmodes, NFluxmodes, NCmodes);
    THS1_matrix = temp_heatsource(NTmodes, NDec1, NCmodes, 1);
    THS2_matrix = temp_heatsource(NTmodes, NDec2, NCmodes, 2);
    THS3_matrix = temp_heatsource(NTmodes, NDec3, NCmodes, 3);
    Info << "\n End \n" << endl;
    msrcoeff(NCmodes);
}

// * * * * * * * * * * * * * * Momentum Eq. Methods * * * * * * * * * * * * * //

Eigen::MatrixXd msrProblem::diffusive_term(label NUmodes, label NPmodes)
{
    label Bsize = NUmodes + liftfield.size();
    Eigen::MatrixXd B_matrix;
    B_matrix.resize(Bsize, Bsize);
    PtrList<volVectorField> Together(0);

    if (liftfield.size() != 0)
    {
        for (label k = 0; k < liftfield.size(); k++)
        {
            Together.append(liftfield[k]);
        }
    }

    if (NUmodes != 0)
    {
        for (label k = 0; k < NUmodes; k++)
        {
            Together.append(Umodes[k]);
        }
    }

    // Project everything
    for (label i = 0; i < Bsize; i++)
    {
        for (label j = 0; j < Bsize; j++)
        {
            B_matrix(i, j) = fvc::domainIntegrate(Together[i] & fvc::laplacian(
                    dimensionedScalar("1", dimless, 1), Together[j])).value();
        }
    }

    // Export the matrix
    ITHACAstream::exportMatrix(B_matrix, "B", "matlab",
                               "./ITHACAoutput/Matrices/fluid_dynamics/");
    return B_matrix;
}

Eigen::MatrixXd msrProblem::pressure_gradient_term(label NUmodes, label NPmodes)
{
    label K1size = NUmodes + liftfield.size();
    label K2size = NPmodes;
    Eigen::MatrixXd K_matrix(K1size, K2size);
    PtrList<volVectorField> Together(0);

    if (liftfield.size() != 0)
    {
        for (label k = 0; k < liftfield.size(); k++)
        {
            Together.append(liftfield[k]);
        }
    }

    if (NUmodes != 0)
    {
        for (label k = 0; k < NUmodes; k++)
        {
            Together.append(Umodes[k]);
        }
    }

    // Project everything
    for (label i = 0; i < K1size; i++)
    {
        for (label j = 0; j < K2size; j++)
        {
            K_matrix(i, j) = fvc::domainIntegrate(Together[i] & fvc::grad(
                    Pmodes[j])).value();
        }
    }

    // Export the matrix
    ITHACAstream::exportMatrix(K_matrix, "K", "matlab",
                               "./ITHACAoutput/Matrices/fluid_dynamics/");
    return K_matrix;
}

List <Eigen::MatrixXd> msrProblem::convective_term(label NUmodes,
        label NPmodes)
{
    label Csize = NUmodes + liftfield.size();
    List <Eigen::MatrixXd> C_matrix;
    C_matrix.setSize(Csize);

    for (label j = 0; j < Csize; j++)
    {
        C_matrix[j].resize(Csize, Csize);
    }

    PtrList<volVectorField> Together(0);

    if (liftfield.size() != 0)
    {
        for (label k = 0; k < liftfield.size(); k++)
        {
            Together.append(liftfield[k]);
        }
    }

    if (NUmodes != 0)
    {
        for (label k = 0; k < NUmodes; k++)
        {
            Together.append(Umodes[k]);
        }
    }

    for (label i = 0; i < Csize; i++)
    {
        for (label j = 0; j < Csize; j++)
        {
            for (label k = 0; k < Csize; k++)
            {
                C_matrix[i](j, k) = fvc::domainIntegrate(Together[i] & fvc::div(
                                        linearInterpolate(Together[j]) & Together[j].mesh().Sf(), Together[k])).value();
            }
        }
    }

    // Export the matrix
    ITHACAstream::exportMatrix(C_matrix, "C", "matlab",
                               "./ITHACAoutput/Matrices/fluid_dynamics/");
    return C_matrix;
}

Eigen::MatrixXd msrProblem::mass_term(label NUmodes, label NPmodes)
{
    label Msize = NUmodes + liftfield.size();
    Eigen::MatrixXd M_matrix(Msize, Msize);
    PtrList<volVectorField> Together(0);

    if (liftfield.size() != 0)
    {
        for (label k = 0; k < liftfield.size(); k++)
        {
            Together.append(liftfield[k]);
        }
    }

    if (NUmodes != 0)
    {
        for (label k = 0; k < NUmodes; k++)
        {
            Together.append(Umodes[k]);
        }
    }

    // Project everything
    for (label i = 0; i < Msize; i++)
    {
        for (label j = 0; j < Msize; j++)
        {
            M_matrix(i, j) = fvc::domainIntegrate(Together[i] & Together[j]).value();
        }
    }

    // Export the matrix
    ITHACAstream::exportMatrix(M_matrix, "M", "matlab",
                               "./ITHACAoutput/Matrices/fluid_dynamics/");
    return M_matrix;
}

// * * * * * * * * * * * * * * Continuity Eq. Methods * * * * * * * * * * * * * //

Eigen::MatrixXd msrProblem::divergence_term(label NUmodes, label NPmodes)
{
    label P1size = NPmodes;
    label P2size = NUmodes + liftfield.size();
    Eigen::MatrixXd P_matrix(P1size, P2size);
    PtrList<volVectorField> Together(0);

    if (liftfield.size() != 0)
    {
        for (label k = 0; k < liftfield.size(); k++)
        {
            Together.append(liftfield[k]);
        }
    }

    if (NUmodes != 0)
    {
        for (label k = 0; k < NUmodes; k++)
        {
            Together.append(Umodes[k]);
        }
    }

    // Project everything
    for (label i = 0; i < P1size; i++)
    {
        for (label j = 0; j < P2size; j++)
        {
            P_matrix(i, j) = fvc::domainIntegrate(Pmodes[i] * fvc::div (
                    Together[j])).value();
        }
    }

    //Export the matrix
    ITHACAstream::exportMatrix(P_matrix, "P", "matlab",
                               "./ITHACAoutput/Matrices/fluid_dynamics/");
    return P_matrix;
}


List <Eigen::MatrixXd> msrProblem::div_momentum(label NUmodes, label NPmodes)
{
    label G1size = NPmodes;
    label G2size = NUmodes + liftfield.size();
    List <Eigen::MatrixXd> G_matrix;
    G_matrix.setSize(G1size);

    for (label j = 0; j < G1size; j++)
    {
        G_matrix[j].resize(G2size, G2size);
    }

    PtrList<volVectorField> Together(0);

    if (liftfield.size() != 0)
    {
        for (label k = 0; k < liftfield.size(); k++)
        {
            Together.append(liftfield[k]);
        }
    }

    if (NUmodes != 0)
    {
        for (label k = 0; k < NUmodes; k++)
        {
            Together.append(Umodes[k]);
        }
    }

    for (label i = 0; i < G1size; i++)
    {
        for (label j = 0; j < G2size; j++)
        {
            for (label k = 0; k < G2size; k++)
            {
                G_matrix[i](j, k) = fvc::domainIntegrate(fvc::grad(Pmodes[i]) & (fvc::div(
                                        fvc::interpolate(Together[j]) & Together[j].mesh().Sf(), Together[k]))).value();
            }
        }
    }

    // Export the matrix
    ITHACAstream::exportMatrix(G_matrix, "G", "matlab",
                               "./ITHACAoutput/Matrices/fluid_dynamics/");
    return G_matrix;
}

Eigen::MatrixXd msrProblem::laplacian_pressure(label NPmodes)
{
    label Dsize = NPmodes;
    Eigen::MatrixXd D_matrix(Dsize, Dsize);

    // Project everything
    for (label i = 0; i < Dsize; i++)
    {
        for (label j = 0; j < Dsize; j++)
        {
            D_matrix(i, j) = fvc::domainIntegrate(fvc::grad(Pmodes[i])&fvc::grad(
                    Pmodes[j])).value();
        }
    }

    //Export the matrix
    ITHACAstream::exportMatrix(D_matrix, "D", "matlab",
                               "./ITHACAoutput/Matrices/fluid_dynamics/");
    return D_matrix;
}

Eigen::MatrixXd msrProblem::pressure_BC1(label NUmodes, label NPmodes)
{
    label P_BC1size = NPmodes;
    label P_BC2size = NUmodes + liftfield.size();
    Eigen::MatrixXd BC1_matrix(P_BC1size, P_BC2size);
    fvMesh& mesh = _mesh();
    PtrList<volVectorField> Together(0);

    if (liftfield.size() != 0)
    {
        for (label k = 0; k < liftfield.size(); k++)
        {
            Together.append(liftfield[k]);
        }
    }

    if (NUmodes != 0)
    {
        for (label k = 0; k < NUmodes; k++)
        {
            Together.append(Umodes[k]);
        }
    }

    for (label i = 0; i < P_BC1size; i++)
    {
        for (label j = 0; j < P_BC2size; j++)
        {
            surfaceScalarField lpl((fvc::interpolate(fvc::laplacian(
                                        Together[j]))&mesh.Sf())*fvc::interpolate(Pmodes[i]));
            double s = 0;

            for (label k = 0; k < lpl.boundaryField().size(); k++)
            {
                s += gSum(lpl.boundaryField()[k]);
            }

            BC1_matrix(i, j) = s;
        }
    }

    ITHACAstream::exportMatrix(BC1_matrix, "BC1", "matlab",
                               "./ITHACAoutput/Matrices/fluid_dynamics/");
    return BC1_matrix;
}


List <Eigen::MatrixXd> msrProblem::pressure_BC2(label NUmodes, label NPmodes)
{
    label P2_BC1size = NPmodes;
    label P2_BC2size = NUmodes + liftfield.size();
    List <Eigen::MatrixXd> BC2_matrix;
    fvMesh& mesh = _mesh();
    BC2_matrix.setSize(P2_BC1size);

    for (label j = 0; j < P2_BC1size; j++)
    {
        BC2_matrix[j].resize(P2_BC2size, P2_BC2size);
    }

    PtrList<volVectorField> Together(0);

    if (liftfield.size() != 0)
    {
        for (label k = 0; k < liftfield.size(); k++)
        {
            Together.append(liftfield[k]);
        }
    }

    if (NUmodes != 0)
    {
        for (label k = 0; k < NUmodes; k++)
        {
            Together.append(Umodes[k]);
        }
    }

    for (label i = 0; i < P2_BC1size; i++)
    {
        for (label j = 0; j < P2_BC2size; j++)
        {
            for (label k = 0; k < P2_BC2size; k++)
            {
                surfaceScalarField div_m(fvc::interpolate(fvc::div(fvc::interpolate(
                                             Together[j]) & mesh.Sf(), Together[k]))&mesh.Sf()*fvc::interpolate(Pmodes[i]));
                double s = 0;

                for (label k = 0; k < div_m.boundaryField().size(); k++)
                {
                    s += gSum(div_m.boundaryField()[k]);
                }

                BC2_matrix[i](j, k) = s;
            }
        }
    }

    // Export the matrix
    ITHACAstream::exportMatrix(BC2_matrix, "BC2", "matlab",
                               "./ITHACAoutput/Matrices/fluid_dynamics/");
    return BC2_matrix;
}

Eigen::MatrixXd msrProblem::pressure_BC3(label NUmodes, label NPmodes)
{
    label P3_BC1size = NPmodes;
    label P3_BC2size = NUmodes + liftfield.size();
    Eigen::MatrixXd BC3_matrix(P3_BC1size, P3_BC2size);
    fvMesh& mesh = _mesh();
    PtrList<volVectorField> Together(0);

    if (liftfield.size() != 0)
    {
        for (label k = 0; k < liftfield.size(); k++)
        {
            Together.append(liftfield[k]);
        }
    }

    if (NUmodes != 0)
    {
        for (label k = 0; k < NUmodes; k++)
        {
            Together.append(Umodes[k]);
        }
    }

    surfaceVectorField n(mesh.Sf() / mesh.magSf());

    for (label i = 0; i < P3_BC1size; i++)
    {
        for (label j = 0; j < P3_BC2size; j++)
        {
            surfaceVectorField BC3 = fvc::interpolate(fvc::curl(Together[j]));
            surfaceVectorField BC4 = n ^ fvc::interpolate(fvc::grad(Pmodes[i]));
            surfaceScalarField BC5 = (BC3 & BC4) * mesh.magSf();
            double s = 0;

            for (label k = 0; k < BC5.boundaryField().size(); k++)
            {
                s += gSum(BC5.boundaryField()[k]);
            }

            BC3_matrix(i, j) = s;
        }
    }

    ITHACAstream::exportMatrix(BC3_matrix, "BC3", "matlab",
                               "./ITHACAoutput/Matrices/fluid_dynamics/");
    return BC3_matrix;
}


// * * * * * * * * * * * * * * Diffusion Eq. Methods * * * * * * * * * * * * * //

List<Eigen::MatrixXd> msrProblem::laplacian_flux(label NFluxmodes,
        label NCmodes)
{
    label LFsize = NFluxmodes;
    List<Eigen::MatrixXd> LF_matrix;
    LF_matrix.setSize(LFsize);

    for (label j = 0; j < LFsize; j++)
    {
        LF_matrix[j].resize(NCmodes, LFsize);
        LF_matrix[j].setZero();
    }

    // Project everything
    for (label i = 0; i < LFsize; i++)
    {
        for (label j = 0; j < NCmodes; j++)
        {
            for (label k = 0; k < LFsize; k++)
            {
                LF_matrix[i](j, k) = fvc::domainIntegrate(Fluxmodes[i] * fvc::laplacian(
                                         Dmodes[j], Fluxmodes[k])).value();
            }
        }
    }

    // Export the matrix
    ITHACAstream::exportMatrix(LF_matrix, "LF", "matlab",
                               "./ITHACAoutput/Matrices/neutronics/");
    return LF_matrix;
}

Eigen::MatrixXd msrProblem::mass_flux(label NFluxmodes)
{
    label MFsize = NFluxmodes;
    Eigen::MatrixXd MF_matrix;
    MF_matrix.resize(MFsize, MFsize);

    // Project everything
    for (label i = 0; i < MFsize; i++)
    {
        for (label j = 0; j < MFsize; j++)
        {
            MF_matrix(i, j) = fvc::domainIntegrate(Fluxmodes[i] * Fluxmodes[j]).value();
        }
    }

    // Export the matrix
    ITHACAstream::exportMatrix(MF_matrix, "MF", "matlab",
                               "./ITHACAoutput/Matrices/neutronics/");
    return MF_matrix;
}

List<Eigen::MatrixXd> msrProblem::prod_flux(label NFluxmodes, label NCmodes)
{
    label PFsize = NFluxmodes;
    List<Eigen::MatrixXd> PF_matrix;
    PF_matrix.setSize(PFsize);

    for (label j = 0; j < PFsize; j++)
    {
        PF_matrix[j].resize(NCmodes, PFsize);
        PF_matrix[j].setZero();
    }

    //Project everything
    for (label i = 0; i < PFsize; i++)
    {
        for (label j = 0; j < NCmodes; j++)
        {
            for (label k = 0; k < PFsize; k++)
            {
                PF_matrix[i](j, k) = fvc::domainIntegrate(Fluxmodes[i] * NSFmodes[j] *
                                     Fluxmodes[k]).value();
            }
        }
    }

    // Export the matrix
    ITHACAstream::exportMatrix(PF_matrix, "PF", "matlab",
                               "./ITHACAoutput/Matrices/neutronics/");
    return PF_matrix;
}

List<Eigen::MatrixXd> msrProblem::abs_flux(label NFluxmodes, label NCmodes)
{
    label AFsize = NFluxmodes;
    List<Eigen::MatrixXd> AF_matrix;
    AF_matrix.setSize(AFsize);

    for (label j = 0; j < AFsize; j++)
    {
        AF_matrix[j].resize(NCmodes, AFsize);
        AF_matrix[j].setZero();
    }

    //Project everything
    for (label i = 0; i < AFsize; i++)
    {
        for (label j = 0; j < NCmodes; j++)
        {
            for (label k = 0; k < AFsize; k++)
            {
                AF_matrix[i](j, k) = fvc::domainIntegrate(Fluxmodes[i] * Amodes[j] *
                                     Fluxmodes[k]).value();
            }
        }
    }

    // Export the matrix
    ITHACAstream::exportMatrix(AF_matrix, "AF", "matlab",
                               "./ITHACAoutput/Matrices/neutronics/");
    return AF_matrix;
}

Eigen::MatrixXd msrProblem::prec_source(label NFluxmodes, label NPrecmodes,
                                        label family)
{
    label p = family;
    PtrList<volScalarField> Precmodes = choose_group("prec", p);
    label PS1size = NFluxmodes;
    label PS2size = NPrecmodes;
    Eigen::MatrixXd PS_matrix;
    PS_matrix.resize(PS1size, PS2size);

    // Project everything
    for (label i = 0; i < PS1size; i++)
    {
        for (label j = 0; j < PS2size; j++)
        {
            PS_matrix(i, j) = fvc::domainIntegrate(Fluxmodes[i] * Precmodes[j]).value();
        }
    }

    savegroupMatrix("PS", p, "./ITHACAoutput/Matrices/neutronics/", PS_matrix);
    return PS_matrix;
}

// * * * * * * * * * * * * * * Precursor Eq. Methods * * * * * * * * * * * * * //

List<Eigen::MatrixXd> msrProblem::stream_term(label NUmodes, label NPrecmodes,
        label family)
{
    label p = family;
    PtrList<volScalarField> Precmodes = choose_group("prec", p);
    label ST1size = NPrecmodes; //Qsizet
    label ST2size = NUmodes + liftfield.size();
    List<Eigen::MatrixXd> ST_matrix;
    ST_matrix.setSize(ST1size);

    for (label j = 0; j < ST1size; j++)
    {
        ST_matrix[j].resize(ST2size, ST1size);
        ST_matrix[j].setZero();
    }

    PtrList<volVectorField> Together(0);

    if (liftfield.size() != 0)
    {
        for (label k = 0; k < liftfield.size(); k++)
        {
            Together.append(liftfield[k]);
        }
    }

    if (NUmodes != 0)
    {
        for (label k = 0; k < NUmodes; k++)
        {
            Together.append(Umodes[k]);
        }
    }

    // Project everything
    for (label i = 0; i < ST1size; i++)
    {
        for (label j = 0; j < ST2size; j++)
        {
            for (label k = 0; k < ST1size; k++)
            {
                ST_matrix[i](j, k) = fvc::domainIntegrate(Precmodes[i] * fvc::div(
                                         fvc::interpolate(Together[j]) & Together[j].mesh().Sf(), Precmodes[k])).value();
            }
        }
    }

    savegroupMatrix("ST", p, "./ITHACAoutput/Matrices/neutronics/", ST_matrix);
    return ST_matrix;
}

Eigen::MatrixXd msrProblem::prec_mass(label NPrecmodes, label family)
{
    label p = family;
    PtrList<volScalarField> Precmodes = choose_group("prec", p);
    label MPsize = NPrecmodes;
    Eigen::MatrixXd MP_matrix;
    MP_matrix.resize(MPsize, MPsize);

    // Project everything
    for (label i = 0; i < MPsize; i++)
    {
        for (label j = 0; j < MPsize; j++)
        {
            MP_matrix(i, j) = fvc::domainIntegrate(Precmodes[i] * Precmodes[j]).value();
        }
    }

    savegroupMatrix("MP", p, "./ITHACAoutput/Matrices/neutronics/", MP_matrix);
    return MP_matrix;
}

Eigen::MatrixXd msrProblem::laplacian_prec(label NPrecmodes, label family)
{
    label p = family;
    PtrList<volScalarField> Precmodes = choose_group("prec", p);
    label LPsize = NPrecmodes;
    Eigen::MatrixXd LP_matrix;
    LP_matrix.resize(LPsize, LPsize);

    for (label i = 0; i < LPsize; i++)
    {
        for (label j = 0; j < LPsize; j++)
        {
            LP_matrix(i, j) = fvc::domainIntegrate(Precmodes[i] * fvc::laplacian(
                    dimensionedScalar("1", dimless, 1), Precmodes[j])).value();
        }
    }

    savegroupMatrix("LP", p, "./ITHACAoutput/Matrices/neutronics/", LP_matrix);
    return LP_matrix;
}

List<Eigen::MatrixXd> msrProblem::flux_source(label NFluxmodes,
        label NPrecmodes, label NCmodes, label family)
{
    label p = family;
    PtrList<volScalarField> Precmodes = choose_group("prec", p);
    label FSsize = NPrecmodes;
    List<Eigen::MatrixXd> FS_matrix;
    FS_matrix.setSize(FSsize);

    for (label j = 0; j < FSsize; j++)
    {
        FS_matrix[j].resize(NCmodes, NFluxmodes);
        FS_matrix[j].setZero();
    }

    // Project everything
    for (label i = 0; i < FSsize; i++)
    {
        for (label j = 0; j < NCmodes; j++)
        {
            for (label k = 0; k < NFluxmodes; k++)
            {
                FS_matrix[i](j, k) = fvc::domainIntegrate(Precmodes[i] * NSFmodes[j] *
                                     Fluxmodes[k]).value();
            }
        }
    }

    savegroupMatrix("FS", p, "./ITHACAoutput/Matrices/neutronics/", FS_matrix);
    return FS_matrix;
}



// * * * * * * * * * *  Decay heat Eq. Methods * * * * * * * * * * * * * //



List<Eigen::MatrixXd> msrProblem::stream_dec(label NUmodes, label NDecmodes,
        label decgroup)
{
    label g = decgroup;
    PtrList<volScalarField> Decmodes = choose_group("dec", g);
    label SD1size = NDecmodes; //Qsizet
    label SD2size = NUmodes + liftfield.size();
    List<Eigen::MatrixXd> SD_matrix;
    SD_matrix.setSize(SD1size);

    for (label j = 0; j < SD1size; j++)
    {
        SD_matrix[j].resize(SD2size, SD1size);
        SD_matrix[j].setZero();
    }

    PtrList<volVectorField> Together(0);

    if (liftfield.size() != 0)
    {
        for (label k = 0; k < liftfield.size(); k++)
        {
            Together.append(liftfield[k]);
        }
    }

    if (NUmodes != 0)
    {
        for (label k = 0; k < NUmodes; k++)
        {
            Together.append(Umodes[k]);
        }
    }

    // Project everything
    for (label i = 0; i < SD1size; i++)
    {
        for (label j = 0; j < SD2size; j++)
        {
            for (label k = 0; k < SD1size; k++)
            {
                SD_matrix[i](j, k) = fvc::domainIntegrate(Decmodes[i] * fvc::div(
                                         fvc::interpolate(Together[j]) & Together[j].mesh().Sf(), Decmodes[k])).value();
            }
        }
    }

    savegroupMatrix("SD", g, "./ITHACAoutput/Matrices/thermal/", SD_matrix);
    return SD_matrix;
}

Eigen::MatrixXd msrProblem::dec_mass(label NDecmodes, label decgroup)
{
    label g = decgroup;
    PtrList<volScalarField> Decmodes = choose_group("dec", g);
    label MDsize = NDecmodes;
    Eigen::MatrixXd MD_matrix;
    MD_matrix.resize(MDsize, MDsize);

    // Project everything
    for (label i = 0; i < MDsize; i++)
    {
        for (label j = 0; j < MDsize; j++)
        {
            MD_matrix(i, j) = fvc::domainIntegrate(Decmodes[i] * Decmodes[j]).value();
        }
    }

    savegroupMatrix("MD", g, "./ITHACAoutput/Matrices/thermal/", MD_matrix);
    return MD_matrix;
}


Eigen::MatrixXd msrProblem::laplacian_dec(label NDecmodes, label decgroup)
{
    label g = decgroup;
    PtrList<volScalarField> Decmodes = choose_group("dec", g);
    label LDsize = NDecmodes;
    Eigen::MatrixXd LD_matrix;
    LD_matrix.resize(LDsize, LDsize);

    // Project everything
    for (label i = 0; i < LDsize; i++)
    {
        for (label j = 0; j < LDsize; j++)
        {
            LD_matrix(i, j) = fvc::domainIntegrate(Decmodes[i] * fvc::laplacian(
                    dimensionedScalar("1", dimless, 1), Decmodes[j])).value();
        }
    }

    savegroupMatrix("LD", g, "./ITHACAoutput/Matrices/thermal/", LD_matrix);
    return LD_matrix;
}

List<Eigen::MatrixXd> msrProblem::dec_fluxsource(label NFluxmodes,
        label NDecmodes, label NCmodes, label decgroup)
{
    label g = decgroup;
    PtrList<volScalarField> Decmodes = choose_group("dec", g);
    label DFSsize = NDecmodes;
    List<Eigen::MatrixXd> DFS_matrix;
    DFS_matrix.setSize(DFSsize);

    for (label j = 0; j < DFSsize; j++)
    {
        DFS_matrix[j].resize(NCmodes, NFluxmodes);
        DFS_matrix[j].setZero();
    }

    // Project everything
    for (label i = 0; i < DFSsize; i++)
    {
        for (label j = 0; j < NCmodes; j++)
        {
            for (label k = 0; k < NFluxmodes; k++)
            {
                DFS_matrix[i](j, k) = fvc::domainIntegrate(Decmodes[i] * SPmodes[j] *
                                      Fluxmodes[k]).value();
            }
        }
    }

    savegroupMatrix("DFS", g, "./ITHACAoutput/Matrices/thermal/", DFS_matrix);
    return DFS_matrix;
}



// * * * * * * * * * * * Temperature Eq. Methods * * * * * * * * * * * * * //

Eigen::MatrixXd msrProblem::mass_temp(label NTmodes)
{
    label TMsize = NTmodes + liftfieldT.size();
    Eigen::MatrixXd TM_matrix;
    TM_matrix.resize(TMsize, TMsize);
    PtrList<volScalarField> TogetherT(0);

    if (liftfieldT.size() != 0)
    {
        for (label k = 0; k < liftfieldT.size(); k++)
        {
            TogetherT.append(liftfieldT[k]);
        }
    }

    if (NTmodes != 0)
    {
        for (label k = 0; k < NTmodes; k++)
        {
            TogetherT.append(Tmodes[k]);
        }
    }

    // Project everything
    for (label i = 0; i < TMsize; i++)
    {
        for (label j = 0; j < TMsize; j++)
        {
            TM_matrix(i, j) = fvc::domainIntegrate(TogetherT[i] * TogetherT[j]).value();
        }
    }

    // Export the matrix
    ITHACAstream::exportMatrix(TM_matrix, "TM", "matlab",
                               "./ITHACAoutput/Matrices/thermal/");
    return TM_matrix;
}

List<Eigen::MatrixXd> msrProblem::temp_stream(label NUmodes, label NTmodes)
{
    label TS1size = NTmodes + liftfieldT.size(); //Qsizet
    label TS2size = NUmodes + liftfield.size();
    List<Eigen::MatrixXd> TS_matrix;
    TS_matrix.setSize(TS1size);

    for (label j = 0; j < TS1size; j++)
    {
        TS_matrix[j].resize(TS2size, TS1size);
        TS_matrix[j].setZero();
    }

    PtrList<volVectorField> Together(0);

    if (liftfield.size() != 0)
    {
        for (label k = 0; k < liftfield.size(); k++)
        {
            Together.append(liftfield[k]);
        }
    }

    if (NUmodes != 0)
    {
        for (label k = 0; k < NUmodes; k++)
        {
            Together.append(Umodes[k]);
        }
    }

    PtrList<volScalarField> TogetherT(0);

    if (liftfieldT.size() != 0)
    {
        for (label k = 0; k < liftfieldT.size(); k++)
        {
            TogetherT.append(liftfieldT[k]);
        }
    }

    if (NTmodes != 0)
    {
        for (label k = 0; k < NTmodes; k++)
        {
            TogetherT.append(Tmodes[k]);
        }
    }

    // Project everything
    for (label i = 0; i < TS1size; i++)
    {
        for (label j = 0; j < TS2size; j++)
        {
            for (label k = 0; k < TS1size; k++)
            {
                TS_matrix[i](j, k) = fvc::domainIntegrate(TogetherT[i] * fvc::div(
                                         fvc::interpolate(Together[j]) & Together[j].mesh().Sf(), TogetherT[k])).value();
            }
        }
    }

    // Export the matrix
    ITHACAstream::exportMatrix(TS_matrix, "TS", "matlab",
                               "./ITHACAoutput/Matrices/thermal/");
    return TS_matrix;
}


Eigen::MatrixXd msrProblem::laplacian_temp(label NTmodes)
{
    label LTsize = NTmodes + liftfieldT.size();
    Eigen::MatrixXd LT_matrix;
    LT_matrix.resize(LTsize, LTsize);
    PtrList<volScalarField> TogetherT(0);

    if (liftfieldT.size() != 0)
    {
        for (label k = 0; k < liftfieldT.size(); k++)
        {
            TogetherT.append(liftfieldT[k]);
        }
    }

    if (NTmodes != 0)
    {
        for (label k = 0; k < NTmodes; k++)
        {
            TogetherT.append(Tmodes[k]);
        }
    }

    // Project everything
    for (label i = 0; i < LTsize; i++)
    {
        for (label j = 0; j < LTsize; j++)
        {
            LT_matrix(i, j) = fvc::domainIntegrate(TogetherT[i] * fvc::laplacian(
                    dimensionedScalar("1", dimless, 1), TogetherT[j])).value();
        }
    }

    // Export the matrix
    ITHACAstream::exportMatrix(LT_matrix, "LT", "matlab",
                               "./ITHACAoutput/Matrices/thermal/");
    return LT_matrix;
}

List<Eigen::MatrixXd> msrProblem::temp_XSfluxsource(label NTmodes,
        label NFluxmodes, label NCmodes)
{
    label TXSsize = NTmodes + liftfieldT.size();
    List<Eigen::MatrixXd> TXS_matrix;
    TXS_matrix.setSize(TXSsize);

    for (label j = 0; j < TXSsize; j++)
    {
        TXS_matrix[j].resize(NCmodes, NFluxmodes);
        TXS_matrix[j].setZero();
    }

    PtrList<volScalarField> TogetherT(0);

    if (liftfieldT.size() != 0)
    {
        for (label k = 0; k < liftfieldT.size(); k++)
        {
            TogetherT.append(liftfieldT[k]);
        }
    }

    if (NTmodes != 0)
    {
        for (label k = 0; k < NTmodes; k++)
        {
            TogetherT.append(Tmodes[k]);
        }
    }

    // Project everything
    for (label i = 0; i < TXSsize; i++)
    {
        for (label j = 0; j < NCmodes; j++)
        {
            for (label k = 0; k < NFluxmodes; k++)
            {
                TXS_matrix[i](j, k) = fvc::domainIntegrate(TogetherT[i] * TXSmodes[j] *
                                      Fluxmodes[k]).value();
            }
        }
    }

    ITHACAstream::exportMatrix(TXS_matrix, "TXS", "matlab",
                               "./ITHACAoutput/Matrices/thermal/");
    return TXS_matrix;
}

List<Eigen::MatrixXd> msrProblem::temp_heatsource(label NTmodes,
        label NDecmodes, label NCmodes, label decgroup)
{
    label  g = decgroup;
    PtrList<volScalarField> Decmodes = choose_group("dec", g);
    label THSsize = NTmodes + liftfieldT.size();
    List<Eigen::MatrixXd> THS_matrix;
    THS_matrix.setSize(THSsize);

    for (label j = 0; j < THSsize; j++)
    {
        THS_matrix[j].resize(NCmodes, NDecmodes);
        THS_matrix[j].setZero();
    }

    PtrList<volScalarField> TogetherT(0);

    if (liftfieldT.size() != 0)
    {
        for (label k = 0; k < liftfieldT.size(); k++)
        {
            TogetherT.append(liftfieldT[k]);
        }
    }

    if (NTmodes != 0)
    {
        for (label k = 0; k < NTmodes; k++)
        {
            TogetherT.append(Tmodes[k]);
        }
    }

    // Project everything
    for (label i = 0; i < THSsize; i++)
    {
        for (label j = 0; j < NCmodes; j++)
        {
            for (label k = 0; k < NDecmodes; k++)
            {
                THS_matrix[i](j, k) = fvc::domainIntegrate(TogetherT[i] * vmodes[j] *
                                      Decmodes[k]).value();
            }
        }
    }

    savegroupMatrix("THS", g, "./ITHACAoutput/Matrices/thermal/", THS_matrix);
    return THS_matrix;
}

PtrList<volScalarField> msrProblem::choose_group(string field, int ith)
{
    if (field == "prec")
    {
        switch (ith)
        {
            case 1:
                return Prec1modes;
                break;

            case 2:
                return Prec2modes;
                break;

            case 3:
                return Prec3modes;
                break;

            case 4:
                return Prec4modes;
                break;

            case 5:
                return Prec5modes;
                break;

            case 6:
                return Prec6modes;
                break;

            case 7:
                return Prec7modes;
                break;

            case 8:
                return Prec8modes;
                break;
        }
    }

    if (field == "dec")
    {
        switch (ith)
        {
            case 1:
                return Dec1modes;
                break;

            case 2:
                return Dec2modes;
                break;

            case 3:
                return Dec3modes;
                break;
        }
    }
}


template<typename M>
void msrProblem::savegroupMatrix(string nome, int n, word folder, M matrice)
{
    nome.append(std::to_string(n));
    word name = nome;
    ITHACAstream::exportMatrix(matrice, name, "matlab", folder);
}

void msrProblem::homogenizeU()
{
    computeLift(Ufield, liftfield, Uomfield);
    homboolU = true;
}

void msrProblem::homogenizeT()
{
    computeLiftT(Tfield, liftfieldT, Tomfield);
    homboolT = true;
}

void msrProblem::readMSRfields()
{
    volVectorField& U = _U();
    volScalarField& p = _p();
    volScalarField& flux = _flux();
    volScalarField& prec1 = _prec1();
    volScalarField& prec2 = _prec2();
    volScalarField& prec3 = _prec3();
    volScalarField& prec4 = _prec4();
    volScalarField& prec5 = _prec5();
    volScalarField& prec6 = _prec6();
    volScalarField& prec7 = _prec7();
    volScalarField& prec8 = _prec8();
    volScalarField& T = _T();
    volScalarField& dec1 = _dec1();
    volScalarField& dec2 = _dec2();
    volScalarField& dec3 = _dec3();
    volScalarField& v = _v();
    volScalarField& D = _D();
    volScalarField& NSF = _NSF();
    volScalarField& A = _A();
    volScalarField& SP = _SP();
    volScalarField& TXS = _TXS();
    volScalarField powerDens = flux * (1 - _decbetaTot()) * _SP() + _decLam1() *
                               dec1 + _decLam2() * dec2 + _decLam3() * dec3;
    powerDens.rename("powerDens");
    ITHACAstream::read_fields(Ufield, U, "./ITHACAoutput/Offline/");
    ITHACAstream::read_fields(Pfield, p, "./ITHACAoutput/Offline/");
    ITHACAstream::read_fields(Fluxfield, flux, "./ITHACAoutput/Offline/");
    ITHACAstream::read_fields(Tfield, T, "./ITHACAoutput/Offline/");
    ITHACAstream::read_fields(Prec1field, prec1, "./ITHACAoutput/Offline/");
    ITHACAstream::read_fields(Prec2field, prec2, "./ITHACAoutput/Offline/");
    ITHACAstream::read_fields(Prec3field, prec3, "./ITHACAoutput/Offline/");
    ITHACAstream::read_fields(Prec4field, prec4, "./ITHACAoutput/Offline/");
    ITHACAstream::read_fields(Prec5field, prec5, "./ITHACAoutput/Offline/");
    ITHACAstream::read_fields(Prec6field, prec6, "./ITHACAoutput/Offline/");
    ITHACAstream::read_fields(Prec7field, prec7, "./ITHACAoutput/Offline/");
    ITHACAstream::read_fields(Prec8field, prec8, "./ITHACAoutput/Offline/");
    ITHACAstream::read_fields(Dec1field, dec1, "./ITHACAoutput/Offline/");
    ITHACAstream::read_fields(Dec2field, dec2, "./ITHACAoutput/Offline/");
    ITHACAstream::read_fields(Dec3field, dec3, "./ITHACAoutput/Offline/");
    ITHACAstream::read_fields(PowerDensfield, powerDens, "./ITHACAoutput/Offline/");
    ITHACAstream::read_fields(vFields, v, "./ITHACAoutput/Offline/");
    ITHACAstream::read_fields(DFields, D, "./ITHACAoutput/Offline/");
    ITHACAstream::read_fields(NSFFields, NSF, "./ITHACAoutput/Offline/");
    ITHACAstream::read_fields(AFields, A, "./ITHACAoutput/Offline/");
    ITHACAstream::read_fields(SPFields, SP, "./ITHACAoutput/Offline/");
    ITHACAstream::read_fields(TXSFields, TXS, "./ITHACAoutput/Offline/");
    mu_samples =
        ITHACAstream::readMatrix("./ITHACAoutput/Offline/mu_samples_mat.txt");
}

void msrProblem::readMSRfields(std::string& dir)
{
    volVectorField& U = _U();
    volScalarField& p = _p();
    volScalarField& flux = _flux();
    volScalarField& prec1 = _prec1();
    volScalarField& prec2 = _prec2();
    volScalarField& prec3 = _prec3();
    volScalarField& prec4 = _prec4();
    volScalarField& prec5 = _prec5();
    volScalarField& prec6 = _prec6();
    volScalarField& prec7 = _prec7();
    volScalarField& prec8 = _prec8();
    volScalarField& T = _T();
    volScalarField& dec1 = _dec1();
    volScalarField& dec2 = _dec2();
    volScalarField& dec3 = _dec3();
    volScalarField& v = _v();
    volScalarField& D = _D();
    volScalarField& NSF = _NSF();
    volScalarField& A = _A();
    volScalarField& SP = _SP();
    volScalarField& TXS = _TXS();
    volScalarField powerDens = flux * (1 - _decbetaTot()) * _SP() + _decLam1() *
                               dec1 + _decLam2() * dec2 + _decLam3() * dec3;
    powerDens.rename("powerDens");
    std::string folder = dir;

    if (ITHACAutilities::check_folder(folder) == true)
    {
        for (int i = 0; i < mu.cols(); i++)
        {
            folder.append(std::to_string(i));
            folder.append("/");
            ITHACAstream::read_fields(Ufield, U, folder);
            ITHACAstream::read_fields(Pfield, p, folder);
            ITHACAstream::read_fields(Fluxfield, flux, folder);
            ITHACAstream::read_fields(Tfield, T, folder);
            ITHACAstream::read_fields(Prec1field, prec1, folder);
            ITHACAstream::read_fields(Prec2field, prec2, folder);
            ITHACAstream::read_fields(Prec3field, prec3, folder);
            ITHACAstream::read_fields(Prec4field, prec4, folder);
            ITHACAstream::read_fields(Prec5field, prec5, folder);
            ITHACAstream::read_fields(Prec6field, prec6, folder);
            ITHACAstream::read_fields(Prec7field, prec7, folder);
            ITHACAstream::read_fields(Prec8field, prec8, folder);
            ITHACAstream::read_fields(Dec1field, dec1, folder);
            ITHACAstream::read_fields(Dec2field, dec2, folder);
            ITHACAstream::read_fields(Dec3field, dec3, folder);
            ITHACAstream::read_fields(PowerDensfield, powerDens, folder);
            ITHACAstream::read_fields(vFields, v, folder);
            ITHACAstream::read_fields(DFields, D, folder);
            ITHACAstream::read_fields(NSFFields, NSF, folder);
            ITHACAstream::read_fields(AFields, A, folder);
            ITHACAstream::read_fields(SPFields, SP, folder);
            ITHACAstream::read_fields(TXSFields, TXS, folder);
            folder = dir;
        }
    }
    else
    {
        std::cout << "Error" << std::endl;
    }
}


void msrProblem::change_viscosity(double mu)
{
    const volScalarField& nu =  _laminarTransport().nu();
    volScalarField& ciao = const_cast<volScalarField&>(nu);
    this->assignIF(ciao, mu);

    for (int i = 0; i < ciao.boundaryFieldRef().size(); i++)
    {
        this->assignBC(ciao, i, mu);
    }
}

void msrProblem::msrcoeff(label& NC)
{
    NCmodes = NC;
    Eigen::MatrixXd Ncoeff_v = ITHACAutilities::getCoeffs(vFields, vmodes);
    ITHACAstream::exportMatrix(Ncoeff_v, "Ncoeff_v", "matlab",
                               "./ITHACAoutput/Matrices/");
    int Ncol = Ncoeff_v.cols();
    SAMPLES_v.resize(NCmodes);
    rbfsplines_v.resize(NCmodes);

    for (int i = 0; i < NCmodes; i++)
    {
        std::cout << "Constructing v RadialBasisFunction for mode " << i + 1 <<
                  std::endl;
        SAMPLES_v[i] = new SPLINTER::DataTable(1, 1);

        for (int j = 0; j < Ncol; j++)
        {
            SAMPLES_v[i]->addSample(mu_samples.row(j), Ncoeff_v(i, j));
        }

        rbfsplines_v[i] = new SPLINTER::RBFSpline(*SAMPLES_v[i],
                SPLINTER::RadialBasisFunctionType::GAUSSIAN);
    }

    Eigen::MatrixXd Ncoeff_D = ITHACAutilities::getCoeffs(DFields, Dmodes);
    ITHACAstream::exportMatrix(Ncoeff_D, "Ncoeff_D", "matlab",
                               "./ITHACAoutput/Matrices/");
    SAMPLES_D.resize(NCmodes);
    rbfsplines_D.resize(NCmodes);

    for (int i = 0; i < NCmodes; i++)
    {
        std::cout << "Constructing D RadialBasisFunction for mode " << i + 1 <<
                  std::endl;
        SAMPLES_D[i] = new SPLINTER::DataTable(1, 1);

        for (int j = 0; j < Ncol; j++)
        {
            SAMPLES_D[i]->addSample(mu_samples.row(j), Ncoeff_D(i, j));
        }

        rbfsplines_D[i] = new SPLINTER::RBFSpline(*SAMPLES_D[i],
                SPLINTER::RadialBasisFunctionType::GAUSSIAN);
    }

    Eigen::MatrixXd Ncoeff_NSF = ITHACAutilities::getCoeffs(NSFFields,
                                 NSFmodes);
    ITHACAstream::exportMatrix(Ncoeff_NSF, "Ncoeff_NSF", "matlab",
                               "./ITHACAoutput/Matrices/");
    SAMPLES_NSF.resize(NCmodes);
    rbfsplines_NSF.resize(NCmodes);

    for (int i = 0; i < NCmodes; i++)
    {
        std::cout << "Constructing NSF RadialBasisFunction for mode " << i + 1 <<
                  std::endl;
        SAMPLES_NSF[i] = new SPLINTER::DataTable(1, 1);

        for (int j = 0; j < Ncol; j++)
        {
            SAMPLES_NSF[i]->addSample(mu_samples.row(j), Ncoeff_NSF(i, j));
        }

        rbfsplines_NSF[i] = new SPLINTER::RBFSpline(*SAMPLES_NSF[i],
                SPLINTER::RadialBasisFunctionType::GAUSSIAN);
    }

    Eigen::MatrixXd Ncoeff_A = ITHACAutilities::getCoeffs(AFields, Amodes);
    ITHACAstream::exportMatrix(Ncoeff_A, "Ncoeff_A", "matlab",
                               "./ITHACAoutput/Matrices/");
    SAMPLES_A.resize(NCmodes);
    rbfsplines_A.resize(NCmodes);

    for (int i = 0; i < NCmodes; i++)
    {
        std::cout << "Constructing A RadialBasisFunction for mode " << i + 1 <<
                  std::endl;
        SAMPLES_A[i] = new SPLINTER::DataTable(1, 1);

        for (int j = 0; j < Ncol; j++)
        {
            SAMPLES_A[i]->addSample(mu_samples.row(j), Ncoeff_A(i, j));
        }

        rbfsplines_A[i] = new SPLINTER::RBFSpline(*SAMPLES_A[i],
                SPLINTER::RadialBasisFunctionType::GAUSSIAN);
    }

    Eigen::MatrixXd Ncoeff_SP = ITHACAutilities::getCoeffs(SPFields,
                                SPmodes);
    ITHACAstream::exportMatrix(Ncoeff_SP, "Ncoeff_SP", "matlab",
                               "./ITHACAoutput/Matrices/");
    SAMPLES_SP.resize(NCmodes);
    rbfsplines_SP.resize(NCmodes);

    for (int i = 0; i < NCmodes; i++)
    {
        std::cout << "Constructing  SP RadialBasisFunction for mode " << i + 1 <<
                  std::endl;
        SAMPLES_SP[i] = new SPLINTER::DataTable(1, 1);

        for (int j = 0; j < Ncol; j++)
        {
            SAMPLES_SP[i]->addSample(mu_samples.row(j), Ncoeff_SP(i, j));
        }

        rbfsplines_SP[i] = new SPLINTER::RBFSpline(*SAMPLES_SP[i],
                SPLINTER::RadialBasisFunctionType::GAUSSIAN);
    }

    Eigen::MatrixXd Ncoeff_TXS = ITHACAutilities::getCoeffs(TXSFields,
                                 TXSmodes);
    ITHACAstream::exportMatrix(Ncoeff_TXS, "Ncoeff_TXS", "matlab",
                               "./ITHACAoutput/Matrices/");
    SAMPLES_TXS.resize(NCmodes);
    rbfsplines_TXS.resize(NCmodes);

    for (int i = 0; i < NCmodes; i++)
    {
        std::cout << "Constructing  TXS RadialBasisFunction for mode " << i + 1 <<
                  std::endl;
        SAMPLES_TXS[i] = new SPLINTER::DataTable(1, 1);

        for (int j = 0; j < Ncol; j++)
        {
            SAMPLES_TXS[i]->addSample(mu_samples.row(j), Ncoeff_TXS(i, j));
        }

        rbfsplines_TXS[i] = new SPLINTER::RBFSpline(*SAMPLES_TXS[i],
                SPLINTER::RadialBasisFunctionType::GAUSSIAN);
    }
}


void msrProblem::liftSolveT()
{
    for (label k = 0; k < inletIndexT.rows(); k++)
    {
        Time& runTime = _runTime();
        fvMesh& mesh = _mesh();
        volScalarField& T = _T();
        volVectorField& U = _U();
        surfaceScalarField& phi = _phi();
        phi = linearInterpolate(U) & mesh.Sf();
        //pisoControl potentialFlow(mesh,"potentialFLow");
        simpleControl simple(mesh);
        IOMRFZoneList& MRF = _MRF();
        singlePhaseTransportModel& laminarTransport = _laminarTransport();
        turbulence = autoPtr<incompressible::turbulenceModel>
                     (
                         incompressible::turbulenceModel::New(U, phi, laminarTransport)
                     );
        dimensionedScalar& nu = _nu();
        dimensionedScalar& Pr = _Pr();
        dimensionedScalar& Prt = _Prt();
        volScalarField& alphat = _alphat();
        volScalarField& v = _v();
        dimensionedScalar& cp = _CpRef();
        dimensionedScalar& rho = _rhoRef();
        dimensionedScalar& sp = _SP1_0();
        dimensionedScalar& dbtot = _decbetaTot();
        label BCind = inletIndexT(k, 0);
        volScalarField Tlift("Tlift" + name(k), T);
        instantList Times = runTime.times();
        runTime.setTime(Times[1], 1);
        Info << "Solving a lifting Problem" << endl;
        scalar t1 = 1;
        scalar t0 = 0;
        alphat = turbulence->nut() / Prt;
        alphat.correctBoundaryConditions();
        volScalarField source
        (
            IOobject
            (
                "source",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("source", dimensionSet(0, -2, -1, 0, 0, 0, 0), 1)
        );

        for (label j = 0; j < T.boundaryField().size(); j++)
        {
            if (j == BCind)
            {
                assignBC(Tlift, j, t1);
                assignIF(Tlift, t0);
            }
            else if (T.boundaryField()[BCind].type() == "fixedValue")
            {
                assignBC(Tlift, j, t0);
                assignIF(Tlift, t0);
            }
            else
            {
            }
        }

        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix TEqn
            (
                fvm::div(phi, Tlift) == fvm::laplacian(turbulence->nu() / Pr + alphat, Tlift)
                + sp * (1 - dbtot) / (cp * rho)*source
            );
            TEqn.solve();
            Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                 << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                 << nl << endl;
        }

        Tlift.write();
        liftfieldT.append(Tlift);
    }
}


void msrProblem::restart()
{
    volScalarField& p = _p();
    volScalarField& p0 = _p0();
    volVectorField& U = _U();
    volVectorField& U0 = _U0();
    surfaceScalarField& phi = _phi();
    surfaceScalarField& phi0 = _phi0();
    volScalarField& flux = _flux();
    volScalarField& flux0 = _flux0();
    volScalarField& prec1 = _prec1();
    volScalarField& prec10 = _prec10();
    volScalarField& prec2 = _prec2();
    volScalarField& prec20 = _prec20();
    volScalarField& prec3 = _prec3();
    volScalarField& prec30 = _prec30();
    volScalarField& prec4 = _prec4();
    volScalarField& prec40 = _prec40();
    volScalarField& prec5 = _prec5();
    volScalarField& prec50 = _prec50();
    volScalarField& prec6 = _prec6();
    volScalarField& prec60 = _prec60();
    volScalarField& prec7 = _prec7();
    volScalarField& prec70 = _prec70();
    volScalarField& prec8 = _prec8();
    volScalarField& prec80 = _prec80();
    volScalarField& T = _T();
    volScalarField& T0 = _T0();
    volScalarField& dec1 = _dec1();
    volScalarField& dec10 = _dec10();
    volScalarField& dec2 = _dec2();
    volScalarField& dec20 = _dec20();
    volScalarField& dec3 = _dec3();
    volScalarField& dec30 = _dec30();
    dimensionedScalar& Keff = _Keff();
    dimensionedScalar& K0 = _K0();
    volScalarField& v = _v();
    volScalarField& v0 = _v0();
    volScalarField& NSF = _NSF();
    volScalarField& NSF0 = _NSF0();
    volScalarField& A = _A();
    volScalarField& A0 = _A0();
    volScalarField& D = _D();
    volScalarField& D0 = _D0();
    volScalarField& SP = _SP();
    volScalarField& SP0 = _SP0();
    volScalarField& TXS = _TXS();
    volScalarField& TXS0 = _TXS0();
    fvMesh& mesh = _mesh();
    p = p0;
    U = U0;
    phi = phi0;
    turbulence.reset(
        (incompressible::turbulenceModel::New(U, phi, _laminarTransport())).ptr()
    );
    flux = flux0;
    prec1 = prec10;
    prec2 = prec20;
    prec3 = prec30;
    prec4 = prec40;
    prec5 = prec50;
    prec6 = prec60;
    prec7 = prec70;
    prec8 = prec80;

    if (precInBool == true)
    {
        prec1.boundaryFieldRef().set(precinIndex,
                                     fvPatchField<scalar>::New(prec10.boundaryField()[precinIndex].type(),
                                             mesh.boundary()[precinIndex], prec1));
        prec2.boundaryFieldRef().set(precinIndex,
                                     fvPatchField<scalar>::New(prec20.boundaryField()[precinIndex].type(),
                                             mesh.boundary()[precinIndex], prec2));
        prec3.boundaryFieldRef().set(precinIndex,
                                     fvPatchField<scalar>::New(prec30.boundaryField()[precinIndex].type(),
                                             mesh.boundary()[precinIndex], prec3));
        prec4.boundaryFieldRef().set(precinIndex,
                                     fvPatchField<scalar>::New(prec40.boundaryField()[precinIndex].type(),
                                             mesh.boundary()[precinIndex], prec4));
        prec5.boundaryFieldRef().set(precinIndex,
                                     fvPatchField<scalar>::New(prec50.boundaryField()[precinIndex].type(),
                                             mesh.boundary()[precinIndex], prec5));
        prec6.boundaryFieldRef().set(precinIndex,
                                     fvPatchField<scalar>::New(prec60.boundaryField()[precinIndex].type(),
                                             mesh.boundary()[precinIndex], prec6));
        prec7.boundaryFieldRef().set(precinIndex,
                                     fvPatchField<scalar>::New(prec70.boundaryField()[precinIndex].type(),
                                             mesh.boundary()[precinIndex], prec7));
        prec8.boundaryFieldRef().set(precinIndex,
                                     fvPatchField<scalar>::New(prec80.boundaryField()[precinIndex].type(),
                                             mesh.boundary()[precinIndex], prec8));
    }

    T = T0;
    dec1 = dec10;
    dec2 = dec20;
    dec3 = dec30;
    Keff = K0;
    v = v0;
    NSF = NSF0;
    A = A0;
    D = D0;
    SP = SP0;
    TXS = TXS0;
}







