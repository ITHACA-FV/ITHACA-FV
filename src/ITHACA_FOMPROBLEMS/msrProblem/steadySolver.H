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

while (simple.loop(runTime))
{
    Info << "Time = " << runTime.timeName() << nl << endl;
    // --- Pressure-velocity SIMPLE corrector
    {
#include "UEqn.H"
#include "TEqn.H"
#include "pEqn.H"
        rho = rhoRef * (1 - betaTE * (T - Tref));
        logT = log(max(0.5, (T / TrefXS)));
        NSF1 = (1 / Keff) * (NSF1_0 + alfa_NSF1 * logT);
        M11 = NSF1 * (1 - betaTot) - (A1_0 + alfa_A1 * logT);
        D1 = (D1_0 + alfa_D1 * logT);
        SP1 = (1 / Keff) * (SP1_0 + alfa_SP1 * logT);
        alphaEff = turbulence->nu() / Pr + turbulence->nut() / Prt;
        diffEff = turbulence->nu() / Sc + turbulence->nut() / Sct;
#include "DiffEqn.H"
#include "precEqns.H"
#include "decEqns.H"
#include "TEqn.H"
    }
    laminarTransport.correct();
    turbulence->correct();
    runTime.write();
    Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
         << "  ClockTime = " << runTime.elapsedClockTime() << " s"
         << nl << endl;
    powerDens = (1 - decbetaTot) * flux * SP1 + (decLam1 * dec1 + decLam2 * dec2 +
                decLam3 * dec3);
    exportSolution(U, name(counter), "./ITHACAoutput/Offline/");
    exportSolution(p, name(counter), "./ITHACAoutput/Offline/");
    exportSolution(flux, name(counter), "./ITHACAoutput/Offline/");
    exportSolution(prec1, name(counter), "./ITHACAoutput/Offline/");
    exportSolution(prec2, name(counter), "./ITHACAoutput/Offline/");
    exportSolution(prec3, name(counter), "./ITHACAoutput/Offline/");
    exportSolution(prec4, name(counter), "./ITHACAoutput/Offline/");
    exportSolution(prec5, name(counter), "./ITHACAoutput/Offline/");
    exportSolution(prec6, name(counter), "./ITHACAoutput/Offline/");
    exportSolution(prec7, name(counter), "./ITHACAoutput/Offline/");
    exportSolution(prec8, name(counter), "./ITHACAoutput/Offline/");
    exportSolution(T, name(counter), "./ITHACAoutput/Offline/");
    exportSolution(dec1, name(counter), "./ITHACAoutput/Offline/");
    exportSolution(dec2, name(counter), "./ITHACAoutput/Offline/");
    exportSolution(dec3, name(counter), "./ITHACAoutput/Offline/");
    exportSolution(powerDens, name(counter), "./ITHACAoutput/Offline/");
    exportSolution(rho, name(counter), "./ITHACAoutput/Offline/");
    exportSolution(NSF1, name(counter), "./ITHACAoutput/Offline/");
    exportSolution(M11, name(counter), "./ITHACAoutput/Offline/");
    exportSolution(D1, name(counter), "./ITHACAoutput/Offline/");
    exportSolution(SP1, name(counter), "./ITHACAoutput/Offline/");
    exportSolution(alphaEff, name(counter), "./ITHACAoutput/Offline/");
    exportSolution(diffEff, name(counter), "./ITHACAoutput/Offline/");
    Ufield.append(U);
    Pfield.append(p);
    Fluxfield.append(flux);
    Prec1field.append(prec1);
    Prec2field.append(prec2);
    Prec3field.append(prec3);
    Prec4field.append(prec4);
    Prec5field.append(prec5);
    Prec6field.append(prec6);
    Prec7field.append(prec7);
    Prec8field.append(prec8);
    Tfield.append(T);
    Dec1field.append(dec1);
    Dec2field.append(dec2);
    Dec3field.append(dec3);
    rhofield.append(rho);
    NSF1field.append(NSF1);
    M11field.append(M11);
    D1field.append(D1);
    SP1field.append(SP1);
    aEfield.append(alphaEff);
    dEfield.append(diffEff);
    counter++;
}



