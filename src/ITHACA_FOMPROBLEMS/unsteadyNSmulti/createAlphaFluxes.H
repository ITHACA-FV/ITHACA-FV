IOobject alphaPhi10Header
(
    IOobject::groupName("alphaPhi0", alpha1.group()),
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
);

alphaRestart = alphaPhi10Header.typeHeaderOk<surfaceScalarField>();

if (alphaRestart)
{
    Info << "Restarting alpha" << endl;
}

// MULES flux from previous time-step
_alphaPhi10 = autoPtr<surfaceScalarField>
( 
new surfaceScalarField
(
    alphaPhi10Header,
    phi*fvc::interpolate(alpha1)
)
);

surfaceScalarField& alphaPhi10 = _alphaPhi10();

// MULES Correction
_talphaPhi1Corr0 = autoPtr<tmp<surfaceScalarField>>
(
	new tmp<surfaceScalarField> 
	);
tmp<surfaceScalarField>& talphaPhi1Corr0 = _talphaPhi1Corr0();