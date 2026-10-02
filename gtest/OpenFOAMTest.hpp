#ifndef ITHACA_OPENFOAM_TEST_HPP
#define ITHACA_OPENFOAM_TEST_HPP

#include <gtest/gtest.h>
#include "fvCFD.H"

namespace OpenFOAMTest
{
    inline int argc = 0;
    inline char** argv = nullptr;

    inline void setArguments(int testArgc, char** testArgv)
    {
        argc = testArgc;
        argv = testArgv;
    }

    template<class CaseType>
    class Fixture : public testing::Test
    {
    protected:
        std::unique_ptr<CaseType> testCase;

        void SetUp() override
        {
            ASSERT_NE(OpenFOAMTest::argv, nullptr);
            ASSERT_GT(OpenFOAMTest::argc, 0);

            testCase = std::make_unique<CaseType>(
                OpenFOAMTest::argc,
                OpenFOAMTest::argv
            );
        }
    };
}

#endif
