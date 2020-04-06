//
// Created by Logan Morrison on 3/19/20.
//

//
// Created by Logan Morrison on 2019-05-24.
//

#include "lanre/cosmology/standard_model.hpp"
#include "gtest/gtest.h"
#include <vector>
#include <ctime>
#include <cmath>

using namespace lanre::cosmology;

const double log_T_min = SM_LOG_TEMP_MIN;
const double log_T_max = SM_LOG_TEMP_MAX;
const size_t num_Ts = SM_DATA_SQRT_GSTAR.size();


class StandardModelTest : public ::testing::Test {
protected:
    void SetUp() override {
        sm_Ts.reserve(num_Ts);
        double step_log_T = (log_T_max - log_T_min) / double(num_Ts - 1);
        for (size_t n = 0; n < num_Ts; n++) {
            double log_T = log_T_min + n * step_log_T;
            sm_Ts.push_back(pow(10.0, log_T));
        }
    }

    // void TearDown() override {}
    std::vector<double> sm_Ts;
};

TEST_F(StandardModelTest, TestSqrtGStar) {
    for (size_t i = 0; i < num_Ts; i++) {
        double approx = sm_sqrt_gstar(sm_Ts[i]);
        double exact = SM_DATA_SQRT_GSTAR[i];
        double frac_diff = fabs(approx - exact) / exact;
        ASSERT_LE(frac_diff, 5e-2);
    }
}

TEST_F(StandardModelTest, TestGeff) {
    for (size_t i = 0; i < num_Ts; i++) {
        double approx = sm_geff(sm_Ts[i]);
        double exact = SM_DATA_GEFF[i];
        double frac_diff = fabs(approx - exact) / exact;
        ASSERT_LE(frac_diff, 5e-2);
    }
}

TEST_F(StandardModelTest, TestHeff) {
    for (size_t i = 0; i < num_Ts; i++) {
        double approx = sm_heff(sm_Ts[i]);
        double exact = SM_DATA_HEFF[i];
        double frac_diff = fabs(approx - exact) / exact;
        ASSERT_LE(frac_diff, 5e-2);
    }
}


int main(int argc, char *argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
