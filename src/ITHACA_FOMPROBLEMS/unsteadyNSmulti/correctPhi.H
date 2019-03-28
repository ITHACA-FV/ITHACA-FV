CorrectPhi
(
    U,
    phi,
    p_rgh,
    surfaceScalarField("rAUf", fvc::interpolate(_rAU())),
    geometricZeroField(),
    pimple
);

#include "continuityErrs.H"
