// Tests of the Foam2Eigen conversions: field2Eigen, field2EigenMap, field2EigenBC, field2EigenMapBC
// and PtrList2EigenBC. The expected values are read element by element from the OpenFOAM fields,
// with the components of each element contiguous
#include "fvCFD.H"
#include "pointFields.H"
#include "Foam2Eigen.H"
#include "../OpenFOAMTest.hpp"
#include <gtest/gtest.h>
#include <chrono>
#include <functional>
#include <iomanip>
#include <type_traits>
#include <utility>

namespace
{
// Distinct, smooth values for every component
template<class Type>
void fillValues(Field<Type>& f, const vectorField& x, const label seed)
{
    forAll(f, i)
    {
        for (direction d = 0; d < pTraits<Type>::nComponents; d++)
        {
            // Function where every value is different
            setComponent(f[i], d) =
                Foam::sin((d + 1) * x[i].x() + seed) * Foam::cos((d + 2) * x[i].y()) + 0.1 * d;
        }
    }
}

template<class Type>
Eigen::VectorXd reference(const UList<Type>& values) // UList is a common base for List, Field and fvPatchField
{
    // Get reference values without using Field2Eigen stuff
    const label nCmpts = pTraits<Type>::nComponents;
    Eigen::VectorXd out(values.size() * nCmpts);
    forAll(values, i)
    {
        for (direction d = 0; d < nCmpts; d++)
        {
            out(i * nCmpts + d) = component(values[i], d);
        }
    }
    return out;
}

// Point patches store no values: their values are the internal values at the patch points
template<class Type, template<class> class PatchField, class GeoMesh>
Eigen::VectorXd patchReference(const GeometricField<Type, PatchField, GeoMesh>& f, const label patchi)
{
    if constexpr (std::is_same<GeoMesh, pointMesh>::value)
    {
        const labelList& meshPoints = f.boundaryField()[patchi].patch().meshPoints();
        Field<Type> values(meshPoints.size());
        forAll(meshPoints, j)
        {
            values[j] = f.primitiveField()[meshPoints[j]];
        }
        return reference(values);
    }
    else
    {
        return reference(f.boundaryField()[patchi]);
    }
}

// This is used to test that field2EigenMap(temporary) is rejected. A bit obscure
template<class FieldType, class = void>
struct isMappable : std::false_type {};

template<class FieldType>
struct isMappable<FieldType,
    std::void_t<decltype(Foam2Eigen::field2EigenMap(std::declval<FieldType>()))>> : std::true_type {};
}
// std::declval<FieldType> pretends to produce a value OF TYPE FieldType, without actually constructing it.
// decltype() gives the type of the expression
// std::void_t<> produces void if the expression is valid
// If the call is invalid, the specialization is not used and the default isMappable is used, which is false.


class TestCase
{
public:
    TestCase(int argc, char* argv[])
    {
        argList::noBanner();
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
    }

    Time& runTime() { return *_runTime; }

    IOobject io(const word& name) const
    {
        return IOobject(name, _runTime->timeName(), *_mesh, IOobject::NO_READ, IOobject::NO_WRITE, false);
    }

    template<class Type>
    GeometricField<Type, fvPatchField, volMesh> volField(const word& name, const label seed) const
    {
        GeometricField<Type, fvPatchField, volMesh> f(
            io(name), *_mesh, dimensioned<Type>("zero", dimless, pTraits<Type>::zero));
        const volVectorField& C = _mesh->C();
        fillValues(f.primitiveFieldRef(), C.primitiveField(), seed);
        forAll(f.boundaryField(), patchi)
        {
            Field<Type> values(f.boundaryField()[patchi].size());
            fillValues(values, C.boundaryField()[patchi], seed + 1);
            f.boundaryFieldRef()[patchi] == values;
        }
        return f;
    }

    surfaceScalarField surfaceField(const word& name, const label seed) const
    {
        surfaceScalarField f(io(name), *_mesh, dimensionedScalar(dimless, 0));
        const surfaceVectorField& Cf = _mesh->Cf();
        fillValues(f.primitiveFieldRef(), Cf.primitiveField(), seed);
        forAll(f.boundaryField(), patchi)
        {
            scalarField values(f.boundaryField()[patchi].size());
            fillValues(values, Cf.boundaryField()[patchi], seed + 1);
            f.boundaryFieldRef()[patchi] == values;
        }
        return f;
    }

    template<class Type>
    GeometricField<Type, pointPatchField, pointMesh> pointField(const word& name, const label seed) const
    {
        GeometricField<Type, pointPatchField, pointMesh> f(
            io(name), pointMesh::New(*_mesh), dimensioned<Type>("zero", dimless, pTraits<Type>::zero));
        fillValues(f.primitiveFieldRef(), _mesh->points(), seed);
        return f;
    }

private:
    autoPtr<argList> _args;
    autoPtr<Time> _runTime;
    autoPtr<fvMesh> _mesh;
};


class Foam2EigenFixture : public testing::Test
{
protected:
    // One argList for the whole suite: destroying it finalizes MPI, which cannot be initialized again
    static void SetUpTestSuite()
    {
        testCase = std::make_unique<TestCase>(OpenFOAMTest::argc, OpenFOAMTest::argv);
    }

    static void TearDownTestSuite()
    {
        testCase.reset();
    }

    void SetUp() override
    {
        std::srand(42);
    }

    template<class Type, template<class> class PatchField, class GeoMesh>
    void checkCopy(GeometricField<Type, PatchField, GeoMesh>& f)
    {
        const Eigen::VectorXd expected = reference(f.primitiveField());
        Eigen::VectorXd copy = Foam2Eigen::field2Eigen(f);

        ASSERT_EQ(copy.size(), f.size() * pTraits<Type>::nComponents);
        EXPECT_TRUE(copy == expected);
        EXPECT_TRUE(Foam2Eigen::field2Eigen(f.primitiveField()) == expected);
        EXPECT_NE(copy.data(), reinterpret_cast<const scalar*>(f.primitiveField().cdata()));

        // The copy and the field are independent
        copy.setConstant(1e6);
        EXPECT_TRUE(reference(f.primitiveField()) == expected);
        const Eigen::VectorXd before = Foam2Eigen::field2Eigen(f);
        Foam2Eigen::field2EigenMap(f) *= 2;
        EXPECT_TRUE(before == expected);
    }

    template<class Type, template<class> class PatchField, class GeoMesh>
    void checkMap(GeometricField<Type, PatchField, GeoMesh>& f)
    {
        const Eigen::VectorXd expected = reference(f.primitiveField());
        Eigen::Map<Eigen::VectorXd> map = Foam2Eigen::field2EigenMap(f);
        Eigen::Map<const Eigen::VectorXd> cmap = Foam2Eigen::field2EigenMap(std::as_const(f));

        ASSERT_EQ(map.size(), f.size() * pTraits<Type>::nComponents);
        ASSERT_EQ(cmap.size(), map.size());
        EXPECT_TRUE(map == expected);
        EXPECT_TRUE(cmap == expected);

        const scalar* storage = reinterpret_cast<const scalar*>(f.primitiveField().cdata());
        EXPECT_EQ(map.data(), storage);
        EXPECT_EQ(cmap.data(), storage);

        const Eigen::VectorXd values = Eigen::VectorXd::Random(map.size());
        map = values;
        EXPECT_TRUE(reference(f.primitiveField()) == values);
        EXPECT_TRUE(cmap == values);
    }

    template<class Type, template<class> class PatchField, class GeoMesh>
    void checkBoundary(GeometricField<Type, PatchField, GeoMesh>& f)
    {
        const List<Eigen::VectorXd> copies = Foam2Eigen::field2EigenBC(f);
        ASSERT_EQ(copies.size(), f.boundaryField().size());

        forAll(f.boundaryField(), patchi)
        {
            const Eigen::VectorXd expected = patchReference(f, patchi);
            ASSERT_EQ(copies[patchi].size(), f.boundaryField()[patchi].size() * pTraits<Type>::nComponents);
            EXPECT_TRUE(copies[patchi] == expected) << "patch " << patchi;

            if constexpr (!std::is_same<GeoMesh, pointMesh>::value)
            {
                const Field<Type>& patch = f.boundaryField()[patchi];
                Eigen::Map<Eigen::VectorXd> map = Foam2Eigen::field2EigenMapBC(f, patchi);
                Eigen::Map<const Eigen::VectorXd> cmap = Foam2Eigen::field2EigenMapBC(std::as_const(f), patchi);

                EXPECT_TRUE(map == expected);
                EXPECT_TRUE(cmap == expected);
                EXPECT_TRUE(Foam2Eigen::field2Eigen(patch) == expected);
                EXPECT_EQ(cmap.data(), reinterpret_cast<const scalar*>(patch.cdata()));

                const Eigen::VectorXd values = Eigen::VectorXd::Random(map.size());
                map = values;
                EXPECT_TRUE(reference(patch) == values);
            }
        }
    }

    template<class FieldType>
    void checkPtrList(const PtrList<FieldType>& fields)
    {
        typedef typename FieldType::value_type Type;
        for (const label nFields : {label(-1), fields.size() - 1})
        {
            const label Nf = nFields == -1 ? fields.size() : nFields;
            const List<Eigen::MatrixXd> matrices = Foam2Eigen::PtrList2EigenBC(fields, nFields);
            ASSERT_EQ(matrices.size(), fields[0].boundaryField().size());

            forAll(matrices, patchi)
            {
                ASSERT_EQ(matrices[patchi].rows(),
                          fields[0].boundaryField()[patchi].size() * pTraits<Type>::nComponents);
                ASSERT_EQ(matrices[patchi].cols(), Nf);
                for (label k = 0; k < Nf; k++)
                {
                    EXPECT_TRUE(matrices[patchi].col(k) == patchReference(fields[k], patchi))
                        << "patch " << patchi << ", field " << k;
                }
            }
        }
    }

    template<class FieldType>
    PtrList<FieldType> makeList(const std::function<FieldType(label)>& make)
    {
        PtrList<FieldType> fields(3);
        forAll(fields, k)
        {
            fields.set(k, new FieldType(make(k)));
        }
        return fields;
    }

    static std::unique_ptr<TestCase> testCase;
};

std::unique_ptr<TestCase> Foam2EigenFixture::testCase;


TEST(Foam2EigenOverloads, ConstAndRvalueResolution)
{
    // Here we check that trying to use a Eigen::Map to a temporary is rejected by the compiler, since it's quite risky
    static_assert(std::is_same<decltype(Foam2Eigen::field2EigenMap(std::declval<volVectorField&>())),
                  Eigen::Map<Eigen::VectorXd>>::value, "A non-const field gives a writable map");
    static_assert(std::is_same<decltype(Foam2Eigen::field2EigenMap(std::declval<const volVectorField&>())),
                  Eigen::Map<const Eigen::VectorXd>>::value, "A const field gives a read-only map");
    static_assert(isMappable<volVectorField&>::value && isMappable<const volVectorField&>::value,
                  "Lvalues can be mapped");
    static_assert(!isMappable<volVectorField>::value && !isMappable<const volVectorField>::value,
                  "Temporaries cannot be mapped");
    static_assert(std::is_same<decltype(Foam2Eigen::field2Eigen(std::declval<volVectorField>())),
                  Eigen::VectorXd>::value, "Temporaries can be copied");
    static_assert(std::is_same<decltype(Foam2Eigen::field2Eigen(std::declval<const volVectorField&>())),
                  Eigen::VectorXd>::value, "Const fields can be copied");
    SUCCEED();
}

TEST_F(Foam2EigenFixture, Field2Eigen)
{
    volScalarField s(testCase->volField<scalar>("s", 0));
    volVectorField U(testCase->volField<vector>("U", 1));
    volTensorField T(testCase->volField<tensor>("T", 2));
    surfaceScalarField phi(testCase->surfaceField("phi", 3));
    pointVectorField d(testCase->pointField<vector>("d", 4));
    checkCopy(s);
    checkCopy(U);
    checkCopy(T);
    checkCopy(phi);
    checkCopy(d);
}

TEST_F(Foam2EigenFixture, Field2EigenMap)
{
    volScalarField s(testCase->volField<scalar>("s", 0));
    volVectorField U(testCase->volField<vector>("U", 1));
    volTensorField T(testCase->volField<tensor>("T", 2));
    surfaceScalarField phi(testCase->surfaceField("phi", 3));
    pointVectorField d(testCase->pointField<vector>("d", 4));
    checkMap(s);
    checkMap(U);
    checkMap(T);
    checkMap(phi);
    checkMap(d);
}

TEST_F(Foam2EigenFixture, BoundaryField)
{
    volScalarField s(testCase->volField<scalar>("s", 0));
    volVectorField U(testCase->volField<vector>("U", 1));
    volTensorField T(testCase->volField<tensor>("T", 2));
    surfaceScalarField phi(testCase->surfaceField("phi", 3));
    pointVectorField d(testCase->pointField<vector>("d", 4));
    pointScalarField ps(testCase->pointField<scalar>("ps", 5));
    checkBoundary(s);
    checkBoundary(U);
    checkBoundary(T);
    checkBoundary(phi);
    checkBoundary(d);
    checkBoundary(ps);
}

TEST_F(Foam2EigenFixture, PtrList2EigenBC)
{
    checkPtrList(makeList<volScalarField>([this](label k) { return testCase->volField<scalar>("s", k); }));
    checkPtrList(makeList<volVectorField>([this](label k) { return testCase->volField<vector>("U", k); }));
    checkPtrList(makeList<volTensorField>([this](label k) { return testCase->volField<tensor>("T", k); }));
    checkPtrList(makeList<surfaceScalarField>([this](label k) { return testCase->surfaceField("phi", k); }));
    checkPtrList(makeList<pointVectorField>([this](label k) { return testCase->pointField<vector>("d", k); }));
    checkPtrList(makeList<pointScalarField>([this](label k) { return testCase->pointField<scalar>("ps", k); }));
}

TEST_F(Foam2EigenFixture, ReadsDoNotModifyField)
{
    // At a new time step, only a write access stores the old time and updates the time index
    volVectorField U(testCase->volField<vector>("U", 1));
    U.oldTime();
    testCase->runTime()++;
    const label timeIndex = U.timeIndex();
    ASSERT_NE(timeIndex, testCase->runTime().timeIndex());

    Foam2Eigen::field2Eigen(U);
    Foam2Eigen::field2EigenMap(std::as_const(U));
    Foam2Eigen::field2EigenBC(U);
    Foam2Eigen::field2EigenMapBC(std::as_const(U), 0);
    EXPECT_EQ(U.timeIndex(), timeIndex);

    Foam2Eigen::field2EigenMap(U);
    EXPECT_EQ(U.timeIndex(), testCase->runTime().timeIndex());
}


// *** Timings *** //

namespace
{
struct Timing
{
    std::string type;
    label size;
    double field2EigenMs;
    double mapReadMs;
    double eigen2fieldMs;
    double mapWriteMs;
};

// Mean time per call. Reads sum the values so that both versions use them. The field has two
// old-time levels, as in a transient solver, which Eigen2field copies with the field
template<class Type, template<class> class PatchField, class GeoMesh>
Timing timeConversions(
    const std::string& type,
    GeometricField<Type, PatchField, GeoMesh>& f,
    const label nRepetitions)
{
    typedef std::chrono::steady_clock clock;
    f.oldTime().oldTime();
    Eigen::VectorXd values = Eigen::VectorXd::Random(f.size() * pTraits<Type>::nComponents);
    scalar sink = 0;

    // Warm up
    sink += Foam2Eigen::field2Eigen(f).sum() + Foam2Eigen::field2EigenMap(std::as_const(f)).sum();

    const auto t0 = clock::now();
    for (label r = 0; r < nRepetitions; r++)
    {
        const Eigen::VectorXd copy = Foam2Eigen::field2Eigen(f);
        sink += copy.sum();
    }
    const auto t1 = clock::now();
    for (label r = 0; r < nRepetitions; r++)
    {
        sink += Foam2Eigen::field2EigenMap(std::as_const(f)).sum();
    }
    const auto t2 = clock::now();
    for (label r = 0; r < nRepetitions; r++)
    {
        f = Foam2Eigen::Eigen2field(f, values, false);
    }
    const auto t3 = clock::now();
    for (label r = 0; r < nRepetitions; r++)
    {
        Foam2Eigen::field2EigenMap(f) = values;
    }
    const auto t4 = clock::now();

    EXPECT_TRUE(std::isfinite(sink));
    const auto ms = [nRepetitions](const clock::duration& d)
    {
        return std::chrono::duration_cast<std::chrono::microseconds>(d).count() * 1e-3 / nRepetitions;
    };
    return Timing{type, f.size(), ms(t1 - t0), ms(t2 - t1), ms(t3 - t2), ms(t4 - t3)};
}
}

TEST_F(Foam2EigenFixture, Timing)
{
    const label nRepetitions = 200;
    volScalarField s(testCase->volField<scalar>("s", 0));
    volVectorField U(testCase->volField<vector>("U", 1));
    volTensorField T(testCase->volField<tensor>("T", 2));
    surfaceScalarField phi(testCase->surfaceField("phi", 3));
    pointVectorField d(testCase->pointField<vector>("d", 4));

    std::vector<Timing> timings;
    timings.push_back(timeConversions("volScalarField", s, nRepetitions));
    timings.push_back(timeConversions("volVectorField", U, nRepetitions));
    timings.push_back(timeConversions("volTensorField", T, nRepetitions));
    timings.push_back(timeConversions("surfaceScalarField", phi, nRepetitions));
    timings.push_back(timeConversions("pointVectorField", d, nRepetitions));

    for (const Timing& t : timings)
    {
        RecordProperty(t.type + "_field2Eigen_ms", std::to_string(t.field2EigenMs));
        RecordProperty(t.type + "_mapRead_ms", std::to_string(t.mapReadMs));
        RecordProperty(t.type + "_Eigen2field_ms", std::to_string(t.eigen2fieldMs));
        RecordProperty(t.type + "_mapWrite_ms", std::to_string(t.mapWriteMs));
    }
    if (!Pstream::master())
    {
        return;
    }

    std::cout << "\nMean time per call [ms] on " << (Pstream::parRun() ? "the master rank" : "a serial run")
              << " (" << nRepetitions << " repetitions, 2 old-time levels)\n"
              << std::left << std::setw(20) << "field" << std::right << std::setw(8) << "size"
              << std::setw(13) << "field2Eigen" << std::setw(10) << "mapRead" << std::setw(9) << "speedup"
              << std::setw(13) << "Eigen2field" << std::setw(10) << "mapWrite" << std::setw(9) << "speedup" << "\n";
    for (const Timing& t : timings)
    {
        std::cout << std::left << std::setw(20) << t.type << std::right << std::setw(8) << t.size
                  << std::fixed << std::setprecision(4)
                  << std::setw(13) << t.field2EigenMs << std::setw(10) << t.mapReadMs
                  << std::setprecision(1) << std::setw(8) << t.field2EigenMs / t.mapReadMs << "x"
                  << std::setprecision(4)
                  << std::setw(13) << t.eigen2fieldMs << std::setw(10) << t.mapWriteMs
                  << std::setprecision(1) << std::setw(8) << t.eigen2fieldMs / t.mapWriteMs << "x\n";
    }
    std::cout << std::endl;
}


int main(int argc, char** argv)
{
    // InitGoogleTest removes the gtest flags from argv, OpenFOAM gets the remaining arguments
    testing::InitGoogleTest(&argc, argv);
    OpenFOAMTest::setArguments(argc, argv);
    return RUN_ALL_TESTS();
}
