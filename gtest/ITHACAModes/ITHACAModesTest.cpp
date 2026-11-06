// Test file relative to the Modes class
#include "fvCFD.H"
#include "Modes.H"
#include "ITHACAPOD.H"
#include "ITHACAutilities.H"
#include "processorFvPatch.H"
#include "../OpenFOAMTest.hpp"
#include <gtest/gtest.h>
#include <type_traits>

/*
The objective here is to generate some fields, not necessary to solve a physical problem.

Vol fields (p_rgh, p): the snapshots keep the original boundary conditions of fixedFluxPressure, so that the
POD computation (and reconstruction) checks the correct treatment of the boundary conditions. If called in the
standard way, modes.reconstruct(...) calls correctBoundaryConditions(), which raises a fatal error in the case
of non-initialized fixedFluxPressure. This can be an issue in the first iterations (reconstructing before the
pressure equation has been solved) or when reconstructing from coefficients.

Surface fields (phi): the faces on processor patches are internal faces of the undecomposed mesh, so the
reconstruction has to assign them from the modes. Old versions of ITHACA did not do it. This is only testable
in a decomposed case.
*/

namespace
{

constexpr label nSnapshots = 20;
constexpr label nPODmodes = 4; // The snapshots can be exactly represented by 4 modes
constexpr scalar tolerance = 1e-10;

class Case
{
public:
    Case(int argc, char* argv[])
    {
        argList::noBanner(); // Little trick to avoid printing the OpenFOAM banner in the test output
        _args = autoPtr<argList>(new argList(argc, argv, true, true, false));
        if (!_args->checkRootCase())
        {
            Foam::FatalError.exit();
        }
        _runTime = autoPtr<Time>(new Time(Time::controlDictName, _args()));
        _mesh = autoPtr<fvMesh>(
            new fvMesh(
                IOobject(
                    fvMesh::defaultRegion,
                    _runTime->timeName(),
                    _runTime->time(),
                    IOobject::MUST_READ)));
        ITHACAparameters::getInstance(*_mesh, *_runTime);
    }

    Time& runTime() { return *_runTime; }
    fvMesh& mesh() { return *_mesh; }

private:
    autoPtr<argList> _args;
    autoPtr<Time> _runTime;
    autoPtr<fvMesh> _mesh;
};

std::unique_ptr<Case> foamCase;

template <class F>
constexpr bool isVol = std::is_same_v<F, volScalarField> || std::is_same_v<F, volVectorField>;
template <class F>
constexpr bool isSurface = std::is_same_v<F, surfaceScalarField>;

// Modes type matching a field type
template <class F>
struct ModesOf;
template <class Type, template <class> class PatchField, class GeoMesh>
struct ModesOf<GeometricField<Type, PatchField, GeoMesh>> {
    using type = Modes<Type, PatchField, GeoMesh>;
};

// Type of the values of a field (scalar or vector)
template <class F>
struct ValueOf;
template <class Type, template <class> class PatchField, class GeoMesh>
struct ValueOf<GeometricField<Type, PatchField, GeoMesh>> {
    using type = Type;
};

// Name of the field (read from 0/) used for each field type
template <class F>
const word fieldNameOf = isSurface<F> ? "phi" : (std::is_same_v<F, volVectorField> ? "U" : "p");

// Value of a field of type Type built from the scalar s: a vector has two different non-zero components
template <class Type>
Type toValue(const scalar s)
{
    if constexpr (std::is_same_v<Type, vector>)
    {
        return vector(s, -0.5 * s, 0.0);
    }
    else
    {
        return s;
    }
}


// Field with the boundary conditions defined in 0/<fieldName> and zero internal values
template <class F>
F readTemplate(fvMesh& mesh, const word& fieldName)
{
    F field(IOobject(fieldName, "0", mesh, IOobject::MUST_READ, IOobject::NO_WRITE), mesh);
    field.primitiveFieldRef() = Zero;
    return field;
}

void requireFixedFluxPressure(const volScalarField& field)
{
    bool found = false;
    forAll(field.boundaryField(), patchI)
    {
        found = found || isA<fixedFluxPressureFvPatchScalarField>(field.boundaryField()[patchI]);
    }
    if (!found)
    {
        FatalErrorInFunction << "The test case must have at least one fixedFluxPressure boundary condition" << exit(FatalError);
    }
}

// Value of the i-th snapshot at the position c: a combination of 4 spatial/temporal modes
scalar snapshotValue(const vector& c, const label i)
{
    constexpr label K = 4;
    constexpr scalar phases[K] = { 0.5, 1.0, 1.5, 2.0 };
    constexpr scalar amplitudes[K] = { 1.0, 0.6, 0.3, 0.15 };
    const scalar pi = constant::mathematical::pi;

    scalar value = 0.0;
    for (label k = 0; k < K; k++)
    {
        const scalar spatialMode = 2.0 * Foam::sin(phases[k] * pi * c.x()) * Foam::sin(phases[k] * pi * c.y());
        const scalar temporalMode = Foam::sqrt(2.0) * Foam::cos(2.0 * pi * scalar(k + 1) * i / scalar(nSnapshots));
        value += amplitudes[k] * spatialMode * temporalMode;
    }
    return value;
}

// Sets the values of a surface field on every patch from the same function as the internal faces.
void setBoundaryValues(fvMesh& mesh, surfaceScalarField& field, const label i)
{
    forAll(field.boundaryFieldRef(), patchI)
    {
        fvsPatchScalarField& patchField = field.boundaryFieldRef()[patchI];
        const vectorField& faceCentres = mesh.Cf().boundaryField()[patchI];

        scalar sign = 1.0;
        if (auto* procPatch = dynamic_cast<const processorFvPatch*>(&patchField.patch()))
        {
            sign = (procPatch->myProcNo() < procPatch->neighbProcNo()) ? 1.0 : -1.0;
        }
        forAll(patchField, faceI)
        {
            patchField[faceI] = sign * snapshotValue(faceCentres[faceI], i);
        }
    }
}

// nSnapshots snapshots of a vol (scalar or vector) or surface scalar field, evaluated at the cell or face centres.
// Surface fields also get their values on every patch (processor ones included).
template <class F>
PtrList<F> makeSnapshots(Case& c, const word& fieldName, const bool checkBC = true)
{
    static_assert(isVol<F> || isSurface<F>, "Only vol scalar/vector and surface scalar fields are supported");
    using Type = typename ValueOf<F>::type;
    fvMesh& mesh = c.mesh();

    F boundaryTemplate = readTemplate<F>(mesh, fieldName);
    if constexpr (std::is_same_v<F, volScalarField>)
    {
        if (checkBC)
        {
            requireFixedFluxPressure(boundaryTemplate);
        }
    }

    // Where the internal field is evaluated: cell centres for vol fields, internal face centres for surface fields
    const vectorField& centres = isVol<F> ? mesh.C().primitiveField() : mesh.Cf().primitiveField();

    PtrList<F> snapshots(nSnapshots);
    c.runTime().setTime(0, 0);
    for (label i = 0; i < nSnapshots; i++)
    {
        c.runTime()++;
        snapshots.set(
            i,
            new F(
                IOobject(fieldName, c.runTime().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
                boundaryTemplate));

        forAll(centres, j)
        {
            snapshots[i].primitiveFieldRef()[j] = toValue<Type>(snapshotValue(centres[j], i));
        }
        if constexpr (isSurface<F>)
        {
            setBoundaryValues(mesh, snapshots[i], i);
        }
    }
    return snapshots;
}

volScalarModes computePOD(PtrList<volScalarField>& snapshots, const word& fieldName, const bool correctBC)
{
    volScalarModes modes;
    ITHACAPOD::getModes(snapshots, modes, fieldName, false, 0, 0, nPODmodes, correctBC);
    return modes;
}

// Reconstructs one field per column of coeffs, starting from a zeroed copy of the field read from 0/
template <class F>
PtrList<F> reconstructAll(
    Case& c, typename ModesOf<F>::type& modes, const Eigen::MatrixXd& coeffs,
    const word& fieldName, const bool treatBC)
{
    F inputTemplate = readTemplate<F>(c.mesh(), fieldName);

    PtrList<F> reconstructed(coeffs.cols());
    for (label i = 0; i < coeffs.cols(); i++)
    {
        reconstructed.set(
            i,
            modes.reconstruct(inputTemplate, coeffs.col(i), "reconstructed" + fieldName, treatBC).clone());
    }
    return reconstructed;
}

// Uses the snapshots themselves as modes, with the unit vectors as coefficients: the reconstruction has to give
// back the snapshots (ITHACAPOD::getModes weighs them with the cell volumes, which is not meaningful for them)
template <class F>
PtrList<F> reproduceSnapshots(Case& c, PtrList<F>& snapshots, const word& fieldName, const bool treatBC)
{
    typename ModesOf<F>::type modes;
    modes.resize(snapshots.size());
    forAll(snapshots, i)
    {
        modes.set(i, snapshots[i].clone());
    }
    return reconstructAll<F>(c, modes, Eigen::MatrixXd::Identity(snapshots.size(), snapshots.size()), fieldName, treatBC);
}

// Worst relative L2 error over the snapshots. It includes the internal field and every patch, or only the
// processor patches if onlyProcessorPatches is true. The sums are reduced exactly once per snapshot, so that
// all the ranks do the same number of collective calls whatever their number of processor patches.
template <class F>
scalar relativeError(PtrList<F>& reference, PtrList<F>& approximation, const bool onlyProcessorPatches = false)
{
    using Type = typename ValueOf<F>::type;
    scalar worst = 0.0;
    for (label i = 0; i < reference.size(); i++)
    {
        scalar diff2 = 0.0, ref2 = 0.0;
        auto accumulate = [&](const Field<Type>& ref, const Field<Type>& approx) {
            diff2 += sum(magSqr(ref - approx));
            ref2 += sum(magSqr(ref));
        };

        if (!onlyProcessorPatches)
        {
            accumulate(reference[i].primitiveField(), approximation[i].primitiveField());
        }
        forAll(reference[i].boundaryField(), patchI)
        {
            const auto& refPatch = reference[i].boundaryField()[patchI];
            if (!onlyProcessorPatches || isA<processorFvPatch>(refPatch.patch()))
            {
                accumulate(refPatch, approximation[i].boundaryField()[patchI]);
            }
        }
        reduce(diff2, sumOp<scalar>());
        reduce(ref2, sumOp<scalar>());
        worst = max(worst, Foam::sqrt(diff2 / max(ref2, VSMALL)));
    }
    return worst;
}

// Total number of faces on processor patches (0 in a serial run)
label nProcessorFaces(const fvMesh& mesh)
{
    label n = 0;
    for (const fvPatch& patch : mesh.boundary())
    {
        if (isA<processorFvPatch>(patch))
        {
            n += patch.size();
        }
    }
    reduce(n, sumOp<label>());
    return n;
}

class ModesTest : public testing::Test
{
protected:
    void SetUp() override
    {
        Foam::FatalError.throwExceptions(); // So we don't stop the test on a fatal error
    }

    void TearDown() override
    {
        Foam::FatalError.dontThrowExceptions();
        // In parallel runs all the ranks write in the same folder: wait for them before removing it
        if (Pstream::parRun())
        {
            UPstream::barrier(UPstream::worldComm);
        }
        if (Pstream::master())
        {
            system("rm -r ./ITHACAoutput");
        }
    }

    Case& testCase() { return *foamCase; }
};

} // namespace


TEST_F(ModesTest, PODGivesTheRequestedNumberOfModes)
{
    PtrList<volScalarField> prgh = makeSnapshots<volScalarField>(testCase(), "p_rgh");
    PtrList<volScalarField> p = makeSnapshots<volScalarField>(testCase(), "p", false);

    ASSERT_EQ(computePOD(prgh, "p_rgh", false).size(), nPODmodes);
    ASSERT_EQ(computePOD(p, "p", true).size(), nPODmodes);
}

TEST_F(ModesTest, CoefficientsHaveOneColumnPerSnapshot)
{
    PtrList<volScalarField> prgh = makeSnapshots<volScalarField>(testCase(), "p_rgh");
    volScalarModes modes = computePOD(prgh, "p_rgh", false);

    Eigen::MatrixXd coeffs = ITHACAutilities::getCoeffs(prgh, modes, nPODmodes);
    ASSERT_EQ(coeffs.rows(), nPODmodes);
    ASSERT_EQ(coeffs.cols(), nSnapshots);
}

// Treating the BC of p_rgh is expected to fail: no gradient info for fixedFluxPressure. The error can be
// raised in any step of the chain (POD, projection, reconstruction), so all of them are in the checked block
TEST_F(ModesTest, TreatingBCFailsForFixedFluxPressure)
{
    PtrList<volScalarField> prgh = makeSnapshots<volScalarField>(testCase(), "p_rgh");

    ASSERT_THROW(
        {
            volScalarModes modes = computePOD(prgh, "p_rgh", true);
            Eigen::MatrixXd coeffs = ITHACAutilities::getCoeffs(prgh, modes, nPODmodes);
            reconstructAll<volScalarField>(testCase(), modes, coeffs, "p_rgh", true);
        },
        std::exception);
}

TEST_F(ModesTest, ReconstructionIsExactWithBC)
{
    PtrList<volScalarField> p = makeSnapshots<volScalarField>(testCase(), "p", false);
    volScalarModes modes = computePOD(p, "p", true);
    Eigen::MatrixXd coeffs = ITHACAutilities::getCoeffs(p, modes, nPODmodes);

    PtrList<volScalarField> reconstructed = reconstructAll<volScalarField>(testCase(), modes, coeffs, "p", true);
    EXPECT_LT(ITHACAutilities::errorL2Rel(p, reconstructed).norm(), tolerance);
}

TEST_F(ModesTest, ReconstructionIsExactWithoutBC)
{
    PtrList<volScalarField> prgh = makeSnapshots<volScalarField>(testCase(), "p_rgh");
    volScalarModes modes = computePOD(prgh, "p_rgh", false);
    Eigen::MatrixXd coeffs = ITHACAutilities::getCoeffs(prgh, modes, nPODmodes);

    PtrList<volScalarField> reconstructed = reconstructAll<volScalarField>(testCase(), modes, coeffs, "p_rgh", false);
    EXPECT_LT(ITHACAutilities::errorL2Rel(prgh, reconstructed).norm(), tolerance);
}

// Trick to run the same test for every field type
template <class F>
class RoundTripTest : public ModesTest
{
};

typedef testing::Types<volScalarField, volVectorField, surfaceScalarField> FieldTypes;
TYPED_TEST_SUITE(RoundTripTest, FieldTypes);

TYPED_TEST(RoundTripTest, ReconstructionAssignsEveryPatch)
{
    const word fieldName = fieldNameOf<TypeParam>;
    if (Pstream::parRun())
    {
        ASSERT_GT(nProcessorFaces(this->testCase().mesh()), 0) << "The decomposed mesh has no processor faces";
    }

    PtrList<TypeParam> snapshots = makeSnapshots<TypeParam>(this->testCase(), fieldName, false);
    std::vector<bool> treatBCOptions = { false };
    if constexpr (isVol<TypeParam>)
    {
        forAll(snapshots, i)
        {
            snapshots[i].correctBoundaryConditions();
        }
        treatBCOptions.push_back(true);
    }

    for (const bool treatBC : treatBCOptions)
    {
        SCOPED_TRACE(treatBC ? "with BC treatment" : "without BC treatment");
        PtrList<TypeParam> reconstructed = reproduceSnapshots(this->testCase(), snapshots, fieldName, treatBC);

        // assignBC leaves the patches whose values are derived (e.g. zeroGradient) to the BC evaluation. The
        // whole vol field is then only comparable if the reconstruction evaluates the BCs. The processor
        // patches are always compared.
        if (!isVol<TypeParam> || treatBC)
        {
            EXPECT_LT(relativeError(snapshots, reconstructed), tolerance)
                << "Reconstruction error on the whole field is too high";
        }
        EXPECT_LT(relativeError(snapshots, reconstructed, true), tolerance)
            << "Reconstruction error on the processor patches is too high";
    }
}

int main(int argc, char** argv)
{
    OpenFOAMTest::setArguments(argc, argv);
    testing::InitGoogleTest(&argc, argv);

    foamCase = std::make_unique<Case>(argc, argv);
    const int result = RUN_ALL_TESTS();
    foamCase.reset(); // Finalizes MPI in a parallel run, after all the tests

    return result;
}
