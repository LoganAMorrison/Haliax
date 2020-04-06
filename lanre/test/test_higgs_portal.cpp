//
// Created by Logan Morrison on 3/15/20.
//

#include "lanre/dm_models/higgs_portal.hpp"
#include <iostream>
#include <gtest/gtest.h>

using namespace lanre::dm_models;

TEST(HiggsPortal, BoltzmannTest) {
    HiggsPortal model{1e3, 1e3, 1.0, 1e-3};

    double rd = model.relic_density(1.0, 100.0, 1e-6, 1e-6);
    std::cout << "rd = " << rd << std::endl;
}