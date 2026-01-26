//- Store the motion state at the beginning of the time-step
const Time& t = mesh.time();
bool firstIter = false;

if (curTimeIndex_ != mesh.time().timeIndex())
{
    sDRBM.newTime();
    curTimeIndex_ = mesh.time().timeIndex();
    firstIter = true;
}

Foam::dimensionedVector g("g", dimAcceleration, Zero);

if (mesh.time().foundObject<uniformDimensionedVectorField>("g"))
{
    g = mesh.time().lookupObject<uniformDimensionedVectorField>("g");
}
else
{
    dictCoeffs.readIfPresent("g", g);
}

const scalar ramp = 1.0;
dictionary forcesDict;
forcesDict.add("type", functionObjects::forces::typeName);
forcesDict.add("patches", dictCoeffs.get<wordRes>("patches"));
forcesDict.add("rhoInf", 1.0);
forcesDict.add("rho", dictCoeffs.getOrDefault<word>("rho", "rho"));
forcesDict.add("CofR", sDRBM.centreOfRotation());

//Foam::functionObjects::forces romforces("romforces", mesh, forcesDict);
romforces.execute(); //calcForcesMoment(); //calcForcesMoments()
//romforcex.append(romforces.forceEff().x());
//romforcey.append(romforces.forceEff().y());
// Solving the sixRigidMotion problem
sDRBM.update
(
    firstIter,
    ramp * (romforces.forceEff() + sDRBM.mass()*g.value()),
    ramp
    * (
        romforces.momentEff()
        + sDRBM.mass() * (sDRBM.momentArm() ^ g.value())
    ),
    t.deltaTValue(),
    t.deltaT0Value()
);

muEval(0, 0) =  sDRBM.centreOfMass().y();

for (int i = 0; i < NmodesDproj; i++)
{
    //New value of the parameter
    //pdCoeff(i, 0) = problem->rbfSplines[i]->eval(muEval);
    pdCoeff(i, 0) = problem->rbfSplines[i]->eval(muEval);
}

// Reconstruction of the pointdisplacement
Dmodes.reconstruct(pointDisplacement, pdCoeff, "pointDisplacement");

/// Correct again the boundary fields
//pointDisplacement.correctBoundaryConditions();
/// update the new pointDisplacement
problem->sDRBMS().pointDisplacement().primitiveFieldRef()
    = pointDisplacement.primitiveFieldRef();
/// Displacement has changed. Update boundary conditions //
pointConstraints::New
(
    pointDisplacement.mesh()

).constrainDisplacement(problem->sDRBMS().pointDisplacement());

mesh.movePoints(problem->sDRBMS().curPoints());
