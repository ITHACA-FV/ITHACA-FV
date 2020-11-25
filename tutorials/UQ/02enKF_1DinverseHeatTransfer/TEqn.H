fv::options& fvOptions(_fvOptions());

fvScalarMatrix TEqn
(
    fvm::ddt(T) - fvm::laplacian(diffusivity, T)
);
TEqn.relax();

fvOptions.constrain(TEqn);

TEqn.solve();

fvOptions.correct(T);

