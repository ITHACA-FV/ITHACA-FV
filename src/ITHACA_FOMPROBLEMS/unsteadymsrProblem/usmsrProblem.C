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
#include "usmsrProblem.H"
// Construct Null
usmsrProblem::usmsrProblem() {}

usmsrProblem::usmsrProblem(int argc, char* argv[])
{
#include "setRootCase.H"
#include "createTime.H"
#include "createMesh.H"
    _pimple = autoPtr<pimpleControl>
              (
                  new pimpleControl
                  (
                      mesh
                  )
              );
    _npimple = autoPtr<pimpleControl>
               (
                   new pimpleControl
                   (
                       mesh,
                       "NPIMPLE"
                   )
               );
#include "createFields.H"
#include "createFields_Neutronics.H"
#include "createFields_Thermal.H"
#include "createConstants.H"
#include "createFvOptions.H"
    offline = ITHACAutilities::check_off();
    podex = ITHACAutilities::check_pod();
}

void usmsrProblem::truthSolve(List<scalar> mu_now)
{
    label cc_p = 0;
    Time& runTime = _runTime();
    fvMesh& mesh = _mesh();
#include "initContinuityErrs.H"
    volScalarField& p = _p();
    volVectorField& U = _U();
    surfaceScalarField& phi = _phi();
    fv::options& fvOptions = _fvOptions();
    pimpleControl& pimple = _pimple();
    pimpleControl& npimple = _npimple();
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
    dimensionedScalar& tau = _tau();
    volScalarField powerDens = (1 - decbetaTot) * flux * SP +
                               (decLam1 * dec1 + decLam2 * dec2 + decLam3 * dec3);
    powerDens.rename("powerDens");
    para = ITHACAparameters::getInstance(mesh, runTime);
    startTime = para->ITHACAdict->lookupOrDefault("startTime", 0);
    finalTime = para->ITHACAdict->lookupOrDefault("finalTime", 1);
    timeStep = para->ITHACAdict->lookupOrDefault("timeStep", 0.1);
    writeEvery = para->ITHACAdict->lookupOrDefault("writeEvery", 0.1);
    instantList Times = runTime.times();
    runTime.setEndTime(finalTime);
    runTime.setTime(Times[1], 1);
    runTime.setDeltaT(timeStep);
    nextWrite = startTime;
    label Ntau = static_cast<int> (tau.value() / timeStep);
    label Ntot = static_cast<int> (finalTime / timeStep);
    bc_prec.resize(8, Ntot + 1);
    bool flagBC = false;
    label nsnapshots = 0;

    while (runTime.run())
    {
#include "readTimeControls.H"
#include "CourantNo.H"
#include "setDeltaT.H"
        runTime.setEndTime(finalTime + timeStep);
        Info << "Time = " << runTime.timeName() << nl << endl;

        while (pimple.loop())
        {
#include "UEqn.H"
#include "TEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
#include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }

        flux_old = flux;

        while (npimple.loop())
        {
#include "updateConsts.H"
#include "DiffEqn.H"
#include "precEqns.H"
#include "TEqn.H"
#include "decEqns.H"
        }

#include "updateK.H"
        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;
        powerDens = (1 - decbetaTot) * flux * SP + (decLam1 * dec1 + decLam2 * dec2 +
                    decLam3 * dec3);

        if (checkWrite(runTime))
        {
            nsnapshots += 1;
            ITHACAstream::exportSolution(U, name(counter), "./ITHACAoutput/Offline/");
            ITHACAstream::exportSolution(p, name(counter), "./ITHACAoutput/Offline/");
            ITHACAstream::exportSolution(flux, name(counter), "./ITHACAoutput/Offline/");
            ITHACAstream::exportSolution(prec1, name(counter), "./ITHACAoutput/Offline/");
            ITHACAstream::exportSolution(prec2, name(counter), "./ITHACAoutput/Offline/");
            ITHACAstream::exportSolution(prec3, name(counter), "./ITHACAoutput/Offline/");
            ITHACAstream::exportSolution(prec4, name(counter), "./ITHACAoutput/Offline/");
            ITHACAstream::exportSolution(prec5, name(counter), "./ITHACAoutput/Offline/");
            ITHACAstream::exportSolution(prec6, name(counter), "./ITHACAoutput/Offline/");
            ITHACAstream::exportSolution(prec7, name(counter), "./ITHACAoutput/Offline/");
            ITHACAstream::exportSolution(prec8, name(counter), "./ITHACAoutput/Offline/");
            ITHACAstream::exportSolution(T, name(counter), "./ITHACAoutput/Offline/");
            ITHACAstream::exportSolution(dec1, name(counter), "./ITHACAoutput/Offline/");
            ITHACAstream::exportSolution(dec2, name(counter), "./ITHACAoutput/Offline/");
            ITHACAstream::exportSolution(dec3, name(counter), "./ITHACAoutput/Offline/");
            ITHACAstream::exportSolution(powerDens, name(counter),
                                         "./ITHACAoutput/Offline/");
            ITHACAstream::exportSolution(v, name(counter), "./ITHACAoutput/Offline/");
            ITHACAstream::exportSolution(D, name(counter), "./ITHACAoutput/Offline/");
            ITHACAstream::exportSolution(NSF, name(counter), "./ITHACAoutput/Offline/");
            ITHACAstream::exportSolution(A, name(counter), "./ITHACAoutput/Offline/");
            ITHACAstream::exportSolution(SP, name(counter), "./ITHACAoutput/Offline/");
            ITHACAstream::exportSolution(TXS, name(counter), "./ITHACAoutput/Offline/");
            std::ofstream of("./ITHACAoutput/Offline/" + name(counter) + "/" +
                             runTime.timeName());
            std::ofstream ofk("./ITHACAoutput/Offline/" + name(counter) + "/" + name(
                                  Keff.value()));
            Ufield.append(tmp<volVectorField>(U));
            Pfield.append(tmp<volScalarField>(p));
            Fluxfield.append(tmp<volScalarField>(flux));
            Prec1field.append(tmp<volScalarField>(prec1));
            Prec2field.append(tmp<volScalarField>(prec2));
            Prec3field.append(tmp<volScalarField>(prec3));
            Prec4field.append(tmp<volScalarField>(prec4));
            Prec5field.append(tmp<volScalarField>(prec5));
            Prec6field.append(tmp<volScalarField>(prec6));
            Prec7field.append(tmp<volScalarField>(prec7));
            Prec8field.append(tmp<volScalarField>(prec8));
            Tfield.append(tmp<volScalarField>(T));
            Dec1field.append(tmp<volScalarField>(dec1));
            Dec2field.append(tmp<volScalarField>(dec2));
            Dec3field.append(tmp<volScalarField>(dec3));
            PowerDensfield.append(tmp<volScalarField>(powerDens));
            vFields.append(tmp<volScalarField>(v));
            DFields.append(tmp<volScalarField>(D));
            NSFFields.append(tmp<volScalarField>(NSF));
            AFields.append(tmp<volScalarField>(A));
            SPFields.append(tmp<volScalarField>(SP));
            TXSFields.append(tmp<volScalarField>(TXS));
            counter++;
            nextWrite += writeEvery;
            writeMu(mu_now);
            // --- Fill in the mu_samples with parameters (time, mu) to be used for the PODI sample points
            mu_samples.conservativeResize(mu_samples.rows() + 1, mu_now.size() + 1);
            mu_samples(mu_samples.rows() - 1, 0) = atof(runTime.timeName().c_str());

            for (label i = 0; i < mu_now.size(); i++)
            {
                mu_samples(mu_samples.rows() - 1, i + 1) = mu_now[i];
            }
        }

        if (precInBool == true)
        {
            computePrecsBC(cc_p);
        }

        if (runTime.value() >= tau.value() && precInBool == true)
        {
            if (flagBC == false)
            {
                flagBC = true;
                changePrecsBC();
            }

            assignPrecsBC(cc_p, Ntau);
        }

        cc_p++;
        runTime++;
    }

    // Resize to Unitary if not initialized by user (i.e. non-parametric problem)
    if (mu.cols() == 0)
    {
        mu.resize(1, 1);
    }

    if (mu_samples.rows() == nsnapshots * mu.cols())
    {
        ITHACAstream::exportMatrix(mu_samples, "mu_samples", "eigen",
                                   "./ITHACAoutput/Offline");
    }
}


bool usmsrProblem::checkWrite(Time& timeObject)
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



void usmsrProblem::truthSolve(List<scalar> mu_now, std::string folder)
{
    mkDir(folder);
    ITHACAutilities::createSymLink(folder);
    counter = 1;
#include "initContinuityErrs.H"
    label cc_p = 0;
    Time& runTime = _runTime();
    fvMesh& mesh = _mesh();
    volScalarField& p = _p();
    volVectorField& U = _U();
    surfaceScalarField& phi = _phi();
    fv::options& fvOptions = _fvOptions();
    pimpleControl& pimple = _pimple();
    pimpleControl& npimple = _npimple();
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
    dimensionedScalar& tau = _tau();
    para = ITHACAparameters::getInstance(mesh, runTime);
    startTime = para->ITHACAdict->lookupOrDefault("startTime", 0);
    finalTime = para->ITHACAdict->lookupOrDefault("finalTime", 1);
    timeStep = para->ITHACAdict->lookupOrDefault("timeStep", 0.1);
    writeEvery = para->ITHACAdict->lookupOrDefault("writeEvery", 0.1);
    instantList Times = runTime.times();
    runTime.setEndTime(finalTime);
    runTime.setTime(Times[1], 1);
    runTime.setDeltaT(timeStep);
    nextWrite = startTime;
    label Ntau = static_cast<int> (tau.value() / timeStep);
    label Ntot = static_cast<int> (finalTime / timeStep);
    bc_prec.resize(8, Ntot + 1);
    bool flagBC = false;
    label nsnapshots = 0;

    while (runTime.run())
    {
#include "readTimeControls.H"
#include "CourantNo.H"
#include "setDeltaT.H"
        runTime.setEndTime(finalTime + timeStep);
        Info << "Time = " << runTime.timeName() << nl << endl;

        while (pimple.loop())
        {
#include "UEqn.H"
#include "TEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
#include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }

        flux_old = flux;

        while (npimple.loop())
        {
#include "updateConsts.H"
#include "DiffEqn.H"
#include "precEqns.H"
#include "TEqn.H"
#include "decEqns.H"
        }

#include "updateK.H"
        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;
        powerDens = (1 - decbetaTot) * flux * SP + (decLam1 * dec1 + decLam2 * dec2 +
                    decLam3 * dec3);

        if (checkWrite(runTime))
        {
            nsnapshots += 1;
            ITHACAstream::exportSolution(U, name(counter), folder);
            ITHACAstream::exportSolution(p, name(counter), folder);
            ITHACAstream::exportSolution(flux, name(counter), folder);
            ITHACAstream::exportSolution(prec1, name(counter), folder);
            ITHACAstream::exportSolution(prec2, name(counter), folder);
            ITHACAstream::exportSolution(prec3, name(counter), folder);
            ITHACAstream::exportSolution(prec4, name(counter), folder);
            ITHACAstream::exportSolution(prec5, name(counter), folder);
            ITHACAstream::exportSolution(prec6, name(counter), folder);
            ITHACAstream::exportSolution(prec7, name(counter), folder);
            ITHACAstream::exportSolution(prec8, name(counter), folder);
            ITHACAstream::exportSolution(T, name(counter), folder);
            ITHACAstream::exportSolution(dec1, name(counter), folder);
            ITHACAstream::exportSolution(dec2, name(counter), folder);
            ITHACAstream::exportSolution(dec3, name(counter), folder);
            ITHACAstream::exportSolution(powerDens, name(counter), folder);
            ITHACAstream::exportSolution(v, name(counter), folder);
            ITHACAstream::exportSolution(D, name(counter), folder);
            ITHACAstream::exportSolution(NSF, name(counter), folder);
            ITHACAstream::exportSolution(A, name(counter), folder);
            ITHACAstream::exportSolution(SP, name(counter), folder);
            ITHACAstream::exportSolution(TXS, name(counter), folder);
            std::ofstream of(folder + "/" + name(counter) + "/" + runTime.timeName());
            std::ofstream ofk(folder + "/" + name(counter) + "/" + name(Keff.value()));
            Ufield.append(tmp<volVectorField>(U));
            Pfield.append(tmp<volScalarField>(p));
            Fluxfield.append(tmp<volScalarField>(flux));
            Prec1field.append(tmp<volScalarField>(prec1));
            Prec2field.append(tmp<volScalarField>(prec2));
            Prec3field.append(tmp<volScalarField>(prec3));
            Prec4field.append(tmp<volScalarField>(prec4));
            Prec5field.append(tmp<volScalarField>(prec5));
            Prec6field.append(tmp<volScalarField>(prec6));
            Prec7field.append(tmp<volScalarField>(prec7));
            Prec8field.append(tmp<volScalarField>(prec8));
            Tfield.append(tmp<volScalarField>(T));
            Dec1field.append(tmp<volScalarField>(dec1));
            Dec2field.append(tmp<volScalarField>(dec2));
            Dec3field.append(tmp<volScalarField>(dec3));
            PowerDensfield.append(tmp<volScalarField>(powerDens));
            vFields.append(tmp<volScalarField>(v));
            DFields.append(tmp<volScalarField>(D));
            NSFFields.append(tmp<volScalarField>(NSF));
            AFields.append(tmp<volScalarField>(A));
            SPFields.append(tmp<volScalarField>(SP));
            TXSFields.append(tmp<volScalarField>(TXS));
            counter++;
            nextWrite += writeEvery;
            writeMu(mu_now);
            // --- Fill in the mu_samples with parameters (time, mu) to be used for the PODI sample points
            mu_samples.conservativeResize(mu_samples.rows() + 1, mu_now.size() + 1);
            mu_samples(mu_samples.rows() - 1, 0) = atof(runTime.timeName().c_str());

            for (label i = 0; i < mu_now.size(); i++)
            {
                mu_samples(mu_samples.rows() - 1, i + 1) = mu_now[i];
            }
        }

        if (precInBool == true)
        {
            computePrecsBC(cc_p);
        }

        if (runTime.value() >= tau.value() && precInBool == true)
        {
            if (flagBC == false)
            {
                flagBC = true;
                changePrecsBC();
            }

            assignPrecsBC(cc_p, Ntau);
        }

        cc_p++;
        runTime++;
    }
}


void usmsrProblem::changePrecsBC()
{
    fvMesh& mesh = _mesh();
    volScalarField& prec1 = _prec1();
    volScalarField& prec2 = _prec2();
    volScalarField& prec3 = _prec3();
    volScalarField& prec4 = _prec4();
    volScalarField& prec5 = _prec5();
    volScalarField& prec6 = _prec6();
    volScalarField& prec7 = _prec7();
    volScalarField& prec8 = _prec8();
    prec1.boundaryFieldRef().set(precinIndex,
                                 fvPatchField<scalar>::New("fixedValue", mesh.boundary()[precinIndex], prec1));
    prec2.boundaryFieldRef().set(precinIndex,
                                 fvPatchField<scalar>::New("fixedValue", mesh.boundary()[precinIndex], prec2));
    prec3.boundaryFieldRef().set(precinIndex,
                                 fvPatchField<scalar>::New("fixedValue", mesh.boundary()[precinIndex], prec3));
    prec4.boundaryFieldRef().set(precinIndex,
                                 fvPatchField<scalar>::New("fixedValue", mesh.boundary()[precinIndex], prec4));
    prec5.boundaryFieldRef().set(precinIndex,
                                 fvPatchField<scalar>::New("fixedValue", mesh.boundary()[precinIndex], prec5));
    prec6.boundaryFieldRef().set(precinIndex,
                                 fvPatchField<scalar>::New("fixedValue", mesh.boundary()[precinIndex], prec6));
    prec7.boundaryFieldRef().set(precinIndex,
                                 fvPatchField<scalar>::New("fixedValue", mesh.boundary()[precinIndex], prec7));
    prec8.boundaryFieldRef().set(precinIndex,
                                 fvPatchField<scalar>::New("fixedValue", mesh.boundary()[precinIndex], prec8));
}

void usmsrProblem::computePrecsBC(label call)
{
    fvMesh& mesh = _mesh();
    volScalarField& prec1 = _prec1();
    volScalarField& prec2 = _prec2();
    volScalarField& prec3 = _prec3();
    volScalarField& prec4 = _prec4();
    volScalarField& prec5 = _prec5();
    volScalarField& prec6 = _prec6();
    volScalarField& prec7 = _prec7();
    volScalarField& prec8 = _prec8();
    dimensionedScalar& lam1 = _lam1();
    dimensionedScalar& lam2 = _lam2();
    dimensionedScalar& lam3 = _lam3();
    dimensionedScalar& lam4 = _lam4();
    dimensionedScalar& lam5 = _lam5();
    dimensionedScalar& lam6 = _lam6();
    dimensionedScalar& lam7 = _lam7();
    dimensionedScalar& lam8 = _lam8();
    dimensionedScalar& tau = _tau();
    bc_prec(0, call) = gSum(prec1.boundaryField()[precoutIndex] *
                            mesh.magSf().boundaryField()[precoutIndex]) / gSum(
                           mesh.magSf().boundaryField()[precoutIndex]) * std::exp(-lam1.value() *
                                   tau.value());
    bc_prec(1, call) = gSum(prec2.boundaryField()[precoutIndex] *
                            mesh.magSf().boundaryField()[precoutIndex]) / gSum(
                           mesh.magSf().boundaryField()[precoutIndex]) * std::exp(-lam2.value() *
                                   tau.value());
    bc_prec(2, call) = gSum(prec3.boundaryField()[precoutIndex] *
                            mesh.magSf().boundaryField()[precoutIndex]) / gSum(
                           mesh.magSf().boundaryField()[precoutIndex]) * std::exp(-lam3.value() *
                                   tau.value());
    bc_prec(3, call) = gSum(prec4.boundaryField()[precoutIndex] *
                            mesh.magSf().boundaryField()[precoutIndex]) / gSum(
                           mesh.magSf().boundaryField()[precoutIndex]) * std::exp(-lam4.value() *
                                   tau.value());
    bc_prec(4, call) = gSum(prec5.boundaryField()[precoutIndex] *
                            mesh.magSf().boundaryField()[precoutIndex]) / gSum(
                           mesh.magSf().boundaryField()[precoutIndex]) * std::exp(-lam5.value() *
                                   tau.value());
    bc_prec(5, call) = gSum(prec6.boundaryField()[precoutIndex] *
                            mesh.magSf().boundaryField()[precoutIndex]) / gSum(
                           mesh.magSf().boundaryField()[precoutIndex]) * std::exp(-lam6.value() *
                                   tau.value());
    bc_prec(6, call) = gSum(prec7.boundaryField()[precoutIndex] *
                            mesh.magSf().boundaryField()[precoutIndex]) / gSum(
                           mesh.magSf().boundaryField()[precoutIndex]) * std::exp(-lam7.value() *
                                   tau.value());
    bc_prec(7, call) = gSum(prec8.boundaryField()[precoutIndex] *
                            mesh.magSf().boundaryField()[precoutIndex]) / gSum(
                           mesh.magSf().boundaryField()[precoutIndex]) * std::exp(-lam8.value() *
                                   tau.value());
}

void usmsrProblem::assignPrecsBC(label call, label Ntau)
{
    volScalarField& prec1 = _prec1();
    volScalarField& prec2 = _prec2();
    volScalarField& prec3 = _prec3();
    volScalarField& prec4 = _prec4();
    volScalarField& prec5 = _prec5();
    volScalarField& prec6 = _prec6();
    volScalarField& prec7 = _prec7();
    volScalarField& prec8 = _prec8();
    assignBC(prec1, precinIndex, bc_prec(0, call - Ntau));
    assignBC(prec2, precinIndex, bc_prec(1, call - Ntau));
    assignBC(prec3, precinIndex, bc_prec(2, call - Ntau));
    assignBC(prec4, precinIndex, bc_prec(3, call - Ntau));
    assignBC(prec5, precinIndex, bc_prec(4, call - Ntau));
    assignBC(prec6, precinIndex, bc_prec(5, call - Ntau));
    assignBC(prec7, precinIndex, bc_prec(6, call - Ntau));
    assignBC(prec8, precinIndex, bc_prec(7, call - Ntau));
}
