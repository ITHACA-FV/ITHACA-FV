// Test file relative to the Modes class
#include "fvCFD.H"
#include "Modes.H"
#include "ITHACAPOD.H"
#include "ITHACAutilities.H"
#include "../OpenFOAMTest.hpp"
#include <gtest/gtest.h>

/*
The objective here is to generate some fields, not necessary to solve a physical problem.
It is important that the snapshot keep the original boundary conditions of fixedFluxPreesure
so that the POD computation (and reconstruction) checks the correct treatment of the boundary conditions.
If called in the standard way, the modes.reconstruct(...) will call correctBoundaryConditions(), which
should raise a fatal error in the case of non-initialized fixedFluxPressure.
This can be an issue in first iterations (when reconstructing before the pressure equation has been solved)
or reconstruction from coefficients.
*/

class TestCase
{
public:
    TestCase(int argc, char* argv[])
    {
        // Standard OpenFOAM boilerplate to set up the Time and fvMesh objects
        argList::noBanner(); // Little trick to avoid printing the OpenFOAM banner in the test output
        autoPtr<argList> _args = autoPtr<argList>(new argList(argc, argv, true, true, false));
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
        _para = ITHACAparameters::getInstance(*_mesh, *_runTime);
    }

    label nSnapshots;
    label nPODmodes;

    PtrList<volScalarField> generateSnapshots(const word fieldName, const bool checkBC = true)
    {
        const label K = 4;
        const scalarList phi = { 0.5, 1.0, 1.5, 2.0 }; // fixed phases
        const scalarList A = { 1.0, 0.6, 0.3, 0.15 }; // fixed amplitudes
        if (nSnapshots < 2)
        {
            FatalErrorInFunction << "Number of snapshots must be at least 2" << exit(FatalError);
        }

        volScalarField boundaryTemplate(
            IOobject(
                fieldName,
                "0",
                *_mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE),
            *_mesh);
        boundaryTemplate.primitiveFieldRef() = 0.0;

        if (checkBC)
        {
            bool fixedfluxpressure_bc = false;
            forAll(boundaryTemplate.boundaryField(), patchI)
            {
                fixedfluxpressure_bc = fixedfluxpressure_bc || isA<fixedFluxPressureFvPatchScalarField>(boundaryTemplate.boundaryField()[patchI]);
            }
            if (!fixedfluxpressure_bc)
            {
                FatalErrorInFunction << "The test case must have at least one fixedFluxPressure boundary condition" << exit(FatalError);
            }
        }

        PtrList<volScalarField> snapshots(nSnapshots);
        const scalar pi = constant::mathematical::pi;
        runTime().setTime(0, 0);
        for (label i = 0; i < nSnapshots; i++)
        {
            runTime()++;
            const scalar mu = 2.0 * pi * scalar(i) / scalar(nSnapshots);
            snapshots.set(
                i,
                new volScalarField(
                    IOobject(
                        fieldName,
                        _runTime->timeName(),
                        *_mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE),
                    boundaryTemplate));

            scalarField& values = snapshots[i].primitiveFieldRef();
            forAll(values, cellI)
            {
                const scalar x = _mesh->C()[cellI].x();
                const scalar y = _mesh->C()[cellI].y();

                for (label k = 0; k < K; k++)
                {
                    const scalar spatialMode = 2.0 * Foam::sin(phi[k] * pi * x) * Foam::sin(phi[k] * pi * y);
                    const scalar temporalMode = Foam::sqrt(2.0) * Foam::cos(2.0 * pi * scalar(k + 1) * i / scalar(nSnapshots));

                    values[cellI] += A[k] * spatialMode * temporalMode;
                }
            }
            // snapshots[i].write();
        }
        return snapshots;
    }

    volScalarModes computePOD(PtrList<volScalarField>& snapshots, const word fieldName, const bool correctBC)
    {
        volScalarModes pod_modes;
        ITHACAPOD::getModes(snapshots, pod_modes, fieldName, false, 0, 0, nPODmodes, correctBC);
        return pod_modes;
    }

    Eigen::MatrixXd getCoefficients(
        PtrList<volScalarField>& snapshots,
        volScalarModes& pod_modes)
    {
        Eigen::MatrixXd coeffs = ITHACAutilities::getCoeffs(
            snapshots, pod_modes, nPODmodes);
        return coeffs;
    }

    PtrList<volScalarField> reconstructSnapshots(
        volScalarModes& pod_modes,
        const Eigen::MatrixXd& coeffs,
        const bool treatBC, const word inputFieldName = "p_rgh")
    {
        autoPtr<volScalarField> field = autoPtr<volScalarField>(
            new volScalarField(
                IOobject(
                    inputFieldName,
                    "0",
                    *_mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE),
                *_mesh));
        field->primitiveFieldRef() = 0.0;
        PtrList<volScalarField> reconstructedSnapshots(nSnapshots);
        volScalarField& inputTemplate = *field;

        const word reconstructedFieldName = "reconstructed" + inputFieldName;
        for (label i = 0; i < nSnapshots; i++)
        {
            reconstructedSnapshots.set(
                i,
                pod_modes.reconstruct(inputTemplate, coeffs.col(i), reconstructedFieldName, treatBC).clone());
            // word timeName = Foam::name(i+1);
            // ITHACAstream::exportSolution(reconstructedSnapshots[i], timeName, "./", reconstructedFieldName);
        }
        return reconstructedSnapshots;
    }

    Time& runTime() { return *_runTime; }
    fvMesh& mesh() { return *_mesh; }

private:
    ITHACAparameters* _para;
    autoPtr<Time> _runTime;
    autoPtr<fvMesh> _mesh;
};


class ModeTestFixture : public testing::Test
{
protected:
    void SetUp() override
    {
        testCase = std::make_unique<TestCase>(OpenFOAMTest::argc, OpenFOAMTest::argv);
        testCase->nSnapshots = 20;
        testCase->nPODmodes = 4; // The test function can be exactly represented by 4 modes
        snapshots = testCase->generateSnapshots("p_rgh");
        snapshotsP = testCase->generateSnapshots("p", false);
        Foam::FatalError.throwExceptions(); // So we don't stop the test on a fatal error
    }

    ~ModeTestFixture() override
    {
        Foam::FatalError.dontThrowExceptions();
        system("rm -r ./ITHACAoutput");
    }

    std::unique_ptr<TestCase> testCase;
    PtrList<volScalarField> snapshots;
    PtrList<volScalarField> snapshotsP;
};

TEST_F(ModeTestFixture, ModesComputation)
{
    volScalarModes pod_modes = testCase->computePOD(snapshots, "p_rgh", false);
    volScalarModes pod_modesP = testCase->computePOD(snapshotsP, "p", true);
    ASSERT_EQ(pod_modes.size(), testCase->nPODmodes);
}

TEST_F(ModeTestFixture, CoefficientsComputation)
{
    volScalarModes pod_modes = testCase->computePOD(snapshots, "p_rgh", false);
    Eigen::MatrixXd coeffs = testCase->getCoefficients(snapshots, pod_modes);
    ASSERT_EQ(coeffs.rows(), testCase->nPODmodes);
    ASSERT_EQ(coeffs.cols(), testCase->nSnapshots);
}

// Here we check both the reconstruction and its accuracy.
TEST_F(ModeTestFixture, ReconstructFromModesWithBC)
{
    // Expect a fatal error due to lack of gradient info for fixedFluxPressure with p_rgh
    ASSERT_THROW({
        volScalarModes pod_modes = testCase->computePOD(snapshots, "p_rgh", true);
        Eigen::MatrixXd coeffs = testCase->getCoefficients(snapshots, pod_modes);
        testCase->reconstructSnapshots(pod_modes, coeffs, true);
    }, std::exception) << "Expected a fatal error to be thrown during reconstruction with BC for 'p_rgh'.";


    ASSERT_NO_THROW({
        volScalarModes pod_modes = testCase->computePOD(snapshotsP, "p", true);
        Eigen::MatrixXd coeffs = testCase->getCoefficients(snapshotsP, pod_modes);
        PtrList<volScalarField> reconstructedSnapshots = testCase->reconstructSnapshots(pod_modes, coeffs, true, "p");

        Eigen::MatrixXd relative_error = ITHACAutilities::errorL2Rel(snapshotsP, reconstructedSnapshots);
        double error_norm = relative_error.norm();
        ASSERT_LT(error_norm, 1e-10) << "Reconstruction error norm is too high: " << error_norm;
    });
}

TEST_F(ModeTestFixture, ReconstructFromModesWithoutBC)
{
    ASSERT_NO_THROW({
        volScalarModes pod_modes = testCase->computePOD(snapshots, "p_rgh", false);
        Eigen::MatrixXd coeffs = testCase->getCoefficients(snapshots, pod_modes);
        PtrList<volScalarField> reconstructedSnapshots = testCase->reconstructSnapshots(pod_modes, coeffs, false, "p_rgh");
        Eigen::MatrixXd relative_error = ITHACAutilities::errorL2Rel(snapshots, reconstructedSnapshots);
        double error_norm = relative_error.norm();
        ASSERT_LT(error_norm, 1e-10) << "Reconstruction error norm is too high: " << error_norm;
    });
}


int main(int argc, char** argv)
{
    OpenFOAMTest::setArguments(argc, argv);
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
