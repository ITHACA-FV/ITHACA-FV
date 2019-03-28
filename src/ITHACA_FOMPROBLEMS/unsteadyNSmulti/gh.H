Info<< "Calculating field g.h\n" << endl;
dimensionedScalar ghRef(- mag(g)*hRef);

//volScalarField gh("gh", (g & mesh.C()) - ghRef);
//surfaceScalarField ghf("ghf", (g & mesh.Cf()) - ghRef);

_gh = autoPtr<volScalarField>
(
	new volScalarField("gh", (g & mesh.C()) - ghRef)
);
volScalarField& gh = _gh();

_ghf = autoPtr<surfaceScalarField>
(
	new surfaceScalarField("ghf", (g & mesh.Cf()) - ghRef)
);
surfaceScalarField& ghf = _ghf();

