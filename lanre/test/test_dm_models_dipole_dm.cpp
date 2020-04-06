//
// Created by Logan Morrison on 3/29/20.
//

#include "lanre/dm_models/dipole_dm.hpp"
#include "gtest/gtest.h"
#include <iostream>

using namespace lanre;
using namespace lanre::dm_models;

TEST(TestDipoleDM, TestGamma) {
    DipoleDM model{100.0, 200.0, 1e-3, 1e-3, 1e5};
    std::cout << model.get_width2() << std::endl;
    std::cout << model.gamma_integrand(10.0, 10.0) << std::endl;
    std::cout << model.gamma(10.0) << std::endl;
}

int main(int argc, char *argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
