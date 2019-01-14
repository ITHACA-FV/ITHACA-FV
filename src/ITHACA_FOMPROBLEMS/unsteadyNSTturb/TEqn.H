// Solve temperature equation
{
    alphat = turbulence->nut()/Prt;
    alphat.correctBoundaryConditions();
    
//    volScalarField kappaEff("alphaEff", turbulence->nu()/Pr + alphat);
    volScalarField alphaEff("alphaEff",  nu/Pr + alphat);

    fvScalarMatrix TEqn
    (
        fvm::ddt(T)
      + fvm::div(phi, T)
      - fvm::laplacian(alphaEff, T)
    );

    TEqn.relax();
   
    TEqn.solve();



  //  h=thermo.Cp()*T;     //enthalpy

  //  thermo.correct();

  //  rho=thermo.rho();
  //  nu=thermo.mu()/rho;

}
