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
/// Source file of the steadyNS class.

#include "compressibleSteadyNS.H"
#include "viscosityModel.H"

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructor
compressibleSteadyNS::compressibleSteadyNS() {}
compressibleSteadyNS::compressibleSteadyNS(int argc, char* argv[])
{
    #include "setRootCaseLists.H"
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
//#include "createFvOptions.H"
    //supex = ITHACAutilities::check_sup();
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
    
    para = new ITHACAparameters;
}


// * * * * * * * * * * * * * * Full Order Methods * * * * * * * * * * * * * * //

// Method to perform a truthSolve
void compressibleSteadyNS::truthSolve(List<scalar> mu_now)
{
    Time& runTime = _runTime();
    fvMesh& mesh = _mesh();
    volScalarField& E = _E();
    volScalarField& rho = _rho();
    volScalarField& psi = _psi();    
    volVectorField& U = _U();
    pressureControl& pressureControl=_pressureControl(); 
    dimensionedScalar& initialMass = _initialMass();
    surfaceScalarField& phi = _phi();
    fv::options& fvOptions = _fvOptions();
    simpleControl& simple = _simple();
    IOMRFZoneList& MRF = _MRF();
    fluidThermo& thermo = pThermo();
    volScalarField& p = pThermo().p();
    volScalarField _nut(turbulence->nut());
    //Info << thermo.mu() << endl;
    //Info << thermo.nu() << endl;
#include "NLsolve.H"
    ITHACAstream::exportSolution(U, name(counter), "./ITHACAoutput/Offline/");
    ITHACAstream::exportSolution(p, name(counter), "./ITHACAoutput/Offline/");
    ITHACAstream::exportSolution(E, name(counter), "./ITHACAoutput/Offline/");  
    ITHACAstream::exportSolution(_nut, name(counter), "./ITHACAoutput/Offline/");
    Ufield.append(U);
    Pfield.append(p);
    Efield.append(E);
    nutFields.append(_nut);
    counter++;
    writeMu(mu_now); 
    mu_samples.conservativeResize(mu_samples.rows() + 1, mu_now.size());
}

void compressibleSteadyNS::change_viscosity(double mu_new)
{
    const volScalarField& mu =  pThermo().mu();
    volScalarField& mu_field = const_cast<volScalarField&>(mu);
    this->assignIF(mu_field, mu_new);

    for (int i = 0; i < mu_field.boundaryFieldRef().size(); i++)
    {
        this->assignBC(mu_field, i, mu_new);
    }
}

fvVectorMatrix compressibleSteadyNS::get_Umatrix(volVectorField& U, volScalarField& p, Vector<double>& uresidual_v)
  {
    IOMRFZoneList& MRF = _MRF();
    surfaceScalarField& phi = _phi();
    volScalarField & rho = _rho();
    fv::options& fvOptions = _fvOptions();
    simpleControl& simple = _simple();

    MRF.correctBoundaryVelocity(U);

    fvVectorMatrix UEqn
    (
        fvm::div(phi, U)
      + MRF.DDt(rho, U)
      + turbulence->divDevRhoReff(U)
     ==
        fvOptions(rho, U)
    );

    UEqn.relax();
    fvOptions.constrain(UEqn);

    if (simple.momentumPredictor())
    {
        uresidual_v = solve(UEqn == -fvc::grad(p)).initialResidual();
        fvOptions.correct(U);
    }
    Ueqn_global = &UEqn;

    return UEqn;
  }

fvScalarMatrix compressibleSteadyNS::get_Ematrix(volVectorField& U, volScalarField& he){}

fvScalarMatrix compressibleSteadyNS::get_Pmatrix(volVectorField& U, volScalarField& p){}

fvScalarMatrix compressibleSteadyNS::get_Pcmatrix(volVectorField& U, volScalarField& p){}