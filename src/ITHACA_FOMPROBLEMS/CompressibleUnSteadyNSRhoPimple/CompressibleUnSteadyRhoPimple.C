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
/// Source file of the Compressible UnSteady RhoPimple with dynamic mesh.

#include "CompressibleUnSteadyRhoPimple.H"
//#include "viscosityModel.H"

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //
// Constructor
CompressibleUnSteadyRhoPimple::CompressibleUnSteadyRhoPimple() {}
CompressibleUnSteadyRhoPimple::CompressibleUnSteadyRhoPimple(int argc,
        char* argv[])
    : UnsteadyProblem()
{
    _args = autoPtr<argList>
            (
                new argList(argc, argv)
            );

    if (!_args->checkRootCase())
    {
        Foam::FatalError.exit();
    }

    argList& args = _args();
#include "createTime.H"
    Info << "Create a dynamic mesh for time = "
         << runTime.timeName() << nl << endl;
    meshPtr = autoPtr<dynamicFvMesh> (dynamicFvMesh::New(args, runTime));
    dynamicFvMesh& mesh = meshPtr();
    _pimple = autoPtr<pimpleControl>
              (
                  new pimpleControl
                  (
                      mesh
                  )
              );
    pimpleControl& pimple = _pimple();
#include "createFields.H"
    //#include "createFvOptions.H"
#include "initContinuityErrs.H"
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
    ITHACAparameters* para = ITHACAparameters::getInstance(mesh, _runTime());
    offline = ITHACAutilities::check_off();
    podex = ITHACAutilities::check_pod();
    //middleExport = para->ITHACAdict->lookupOrDefault<bool>("middleExport", true);
    //middleStep = para->ITHACAdict->lookupOrDefault<label>("middleStep", 20);
}


// * * * * * * * * * * * * * * Full Order Methods * * * * * * * * * * * * * * //

/// Method to perform a truthSolve
void CompressibleUnSteadyRhoPimple::truthSolve(word folder)
{
    Time& runTime = _runTime();
    volVectorField& U = _U();
    volScalarField& p = pThermo().p();
    volScalarField& E = pThermo().he();
    pimpleControl& pimple = _pimple();
    volScalarField _nut(turbulence->nut());
    instantList Times = runTime.times();
    runTime.setEndTime(finalTime);
    runTime.setTime(Times[1], 1);
    runTime.setDeltaT(timeStep);
    nextWrite = startTime; // timeStep initialization
    int counter = 1;
#include "NLsolve.H"
}

void CompressibleUnSteadyRhoPimple::changeViscosity(double mu_new)
{
    const volScalarField& mu =  pThermo().mu();
    volScalarField& mu_field = const_cast<volScalarField&>(mu);
    this->assignIF(mu_field, mu_new);

    for (int i = 0; i < mu_field.boundaryFieldRef().size(); i++)
    {
        this->assignBC(mu_field, i, mu_new);
    }
}

fvVectorMatrix CompressibleUnSteadyRhoPimple::getNLTerm(volVectorField& U)
{
    surfaceScalarField& phi = _phi();
    fvVectorMatrix NLTerm = fvm::div(phi, U);
    return NLTerm;
}

fvVectorMatrix CompressibleUnSteadyRhoPimple::getViscTerm(volVectorField& U)
{
    volScalarField& rho = _rho();
    volScalarField nuEff = turbulence->nuEff();
    //fvVectorMatrix viscTerm = turbulence->divDevRhoReff(U);
    fvVectorMatrix viscTerm =  - fvc::div((rho * nuEff) * dev2(T(fvc::grad(
            U)))) - fvm::laplacian(rho * nuEff, U);
    return viscTerm;
}

volVectorField CompressibleUnSteadyRhoPimple::getGradP(volScalarField& p)
{
    volVectorField gradP = fvc::grad(p);
    return gradP;
}


fvVectorMatrix CompressibleUnSteadyRhoPimple::getUmatrix(volVectorField& U)
{
    volScalarField& rho = _rho();
    fv::options& fvOptions = _fvOptions();
    Ueqn_global.reset(
        new fvVectorMatrix
        ( fvm::ddt(rho, U)
          +
          getNLTerm(U)
          +
          getViscTerm(U)
          ==
          fvOptions(rho, U)
        )
    );
    Ueqn_global().relax();
    fvOptions.constrain(Ueqn_global());
    return Ueqn_global();
}

fvScalarMatrix CompressibleUnSteadyRhoPimple::getFluxTerm()
{
    volScalarField& he = pThermo().he();
    surfaceScalarField& phi = _phi();
    fvScalarMatrix fluxTerm = fvm::div(phi, he);
    return fluxTerm;
}

volScalarField CompressibleUnSteadyRhoPimple::getKinEnTerm(volVectorField& U,
        volScalarField& p)
{
    surfaceScalarField& phi = _phi();
    volScalarField& rho = _rho();
    volScalarField kinEn = fvc::div(phi, volScalarField("Ekp",
                                    0.5 * magSqr(U) + p / rho));
    return kinEn;
}

fvScalarMatrix CompressibleUnSteadyRhoPimple::getDiffTerm()
{
    volScalarField& he = pThermo().he();
    fvScalarMatrix diffTerm = fvm::laplacian(turbulence->alphaEff(), he);
    return diffTerm;
}

fvScalarMatrix CompressibleUnSteadyRhoPimple::getEmatrix(volVectorField& U,
        volScalarField& p)
{
    volScalarField& he = pThermo().he();
    volScalarField& rho = _rho();
    fv::options& fvOptions = _fvOptions();
    Eeqn_global.reset(new fvScalarMatrix(
                          getFluxTerm() + getKinEnTerm(U, p) - getDiffTerm()
                          ==
                          fvOptions(rho, he)
                      ));
    Eeqn_global().relax();
    fvOptions.constrain(Eeqn_global());
    return Eeqn_global();
}

surfaceScalarField CompressibleUnSteadyRhoPimple::getPhiHbyA(
    fvVectorMatrix& Ueqn,
    volVectorField& U, volScalarField& p)
{
    volScalarField& rho = _rho();
    volScalarField rAU(1.0 /
                       Ueqn.A()); // Inverse of the diagonal part of the U equation matrix
    HbyA.reset(new volVectorField(constrainHbyA(rAU * Ueqn.H(), U,
                                  p))); // H is the extra diagonal part summed to the r.h.s. of the U equation
    phiHbyA.reset(new surfaceScalarField("phiHbyA",
                                         fvc::interpolate(rho)*fvc::flux(HbyA.ref())));
    return phiHbyA;
}

volScalarField CompressibleUnSteadyRhoPimple::getDivPhiHbyA(
    fvVectorMatrix& Ueqn,
    volVectorField& U, volScalarField& p)
{
    volScalarField divPhiHbyA = fvc::div(getPhiHbyA(Ueqn, U, p));
    return divPhiHbyA;
}

surfaceScalarField CompressibleUnSteadyRhoPimple::getRhorAUf(
    fvVectorMatrix& Ueqn)
{
    volScalarField& rho = _rho();
    volScalarField rAU(1.0 /
                       Ueqn.A()); // Inverse of the diagonal part of the U equation matrix
    rhorAUf.reset(new surfaceScalarField("rhorAUf", fvc::interpolate(rho * rAU)));
    return rhorAUf;
}

fvScalarMatrix CompressibleUnSteadyRhoPimple::getPoissonTerm(
    fvVectorMatrix& Ueqn, volScalarField& p)
{
    fvScalarMatrix poissonTerm = fvm::laplacian(getRhorAUf(Ueqn), p);
    return poissonTerm;
}

fvScalarMatrix CompressibleUnSteadyRhoPimple::getPmatrix(fvVectorMatrix& Ueqn,
        volVectorField& U, volScalarField& p)
{
    volScalarField& rho = _rho();
    fv::options& fvOptions = _fvOptions();
    volScalarField& psi = _psi();
    pressureControl& pressureControl = _pressureControl();
    Peqn_global.reset(new fvScalarMatrix(
                          getDivPhiHbyA(Ueqn, U, p)
                          - getPoissonTerm(Ueqn, p)
                          ==
                          fvOptions(psi, p, rho.name())
                      ));
    Peqn_global().setReference
    (
        pressureControl.refCell(),
        pressureControl.refValue()
    );
    return Peqn_global();
}

void CompressibleUnSteadyRhoPimple::restart()
{
    _runTime().objectRegistry::clear();
    _mesh().objectRegistry::clear();
    // _mesh.clear();
    // _runTime.clear();
    _simple.clear();
    pThermo.clear();
    _p.clear();
    _rho.clear();
    _E.clear();
    _U.clear();
    _phi.clear();
    turbulence.clear();
    _initialMass.clear();
    _pressureControl.clear();
    _psi.clear();
    _fvOptions.clear();
    argList& args = _args();
    Info << "ReCreate time\n" << Foam::endl;
    Time& runTime = _runTime();
    runTime.setTime(0, 1);
    Foam::fvMesh& mesh = _mesh();
    _simple = autoPtr<simpleControl>
              (
                  new simpleControl
                  (
                      mesh
                  )
              );
    simpleControl& simple = _simple();
    pThermo = autoPtr<fluidThermo>
              (
                  fluidThermo::New(mesh)
              );
    fluidThermo& thermo = pThermo();
    thermo.validate(_args().executable(), "h", "e");
    _p = autoPtr<volScalarField>
         (
             new volScalarField(thermo.p())
         );
    volScalarField& p = thermo.p();
    _rho = autoPtr<volScalarField>
           (
               new volScalarField
               (
                   IOobject
                   (
                       "rho",
                       runTime.timeName(),
                       mesh,
                       IOobject::READ_IF_PRESENT,
                       IOobject::AUTO_WRITE
                   ),
                   thermo.rho()
               )
           );
    volScalarField& rho = _rho();
    _E = autoPtr<volScalarField>
         (
             new volScalarField(thermo.he())
         );
    Info << "ReReading field U\n" << endl;
    _U = autoPtr<volVectorField>
         (
             new volVectorField
             (
                 IOobject
                 (
                     "U",
                     runTime.timeName(),
                     mesh,
                     IOobject::MUST_READ,
                     IOobject::AUTO_WRITE
                 ),
                 mesh
             )
         );
    volVectorField& U = _U();
    Info << "ReReading/calculating face flux field phi\n" << endl;
    _phi = autoPtr<surfaceScalarField>
           (
               new surfaceScalarField
               (
                   IOobject
                   (
                       "phi",
                       runTime.timeName(),
                       mesh,
                       IOobject::READ_IF_PRESENT,
                       IOobject::AUTO_WRITE
                   ),
                   linearInterpolate(rho * U) & mesh.Sf()
               )
           );
    surfaceScalarField& phi = _phi();
    _pressureControl = autoPtr<pressureControl>
                       (
                           new pressureControl(p, rho, _simple().dict())
                       );
    mesh.setFluxRequired(p.name());
    Info << "ReCreating turbulence model\n" << endl;
    turbulence = autoPtr<compressible::turbulenceModel>
                 (
                     compressible::turbulenceModel::New
                     (
                         rho,
                         U,
                         phi,
                         thermo
                     )
                 );
    _initialMass = autoPtr<dimensionedScalar>
                   (
                       new dimensionedScalar(fvc::domainIntegrate(rho))
                   );
    _psi = autoPtr<volScalarField>
           (
               new volScalarField(thermo.psi())
           );
    _fvOptions = autoPtr<fv::options>(new fv::options(mesh));
#include "initContinuityErrs.H"
    //turbulence->validate();
}
