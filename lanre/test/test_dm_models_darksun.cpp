//
// Created by Logan Morrison on 3/19/20.
//
#include "lanre/dm_models/darksun.hpp"
#include <gtest/gtest.h>

using namespace lanre::dm_models;

class DarkSUNTest : public ::testing::Test {
protected:
    void SetUp() override {
    }

    DarkSUN model1{1e-3, 11, 1.0, 1.0, 1.0, 1.0, 1.0, 1e-2, false};
};

TEST_F(DarkSUNTest, BoltzmannTest) {

    //auto rds1 = model1.relic_densities(1e-2, 1e-4, "rodas");
    model1.solve_boltzmann(1e-5, 1e-5, "radau");
    //std::cout << "rodas: rds = " << rds1.first << ", " << rds1.second << std::endl;
    std::cout << "radau: rds = " << model1.get_rd_eta() << ", " << model1.get_rd_delta() << std::endl;
    std::cout << "xi_cmb = " << model1.get_xi_cmb() << std::endl;
    std::cout << "xi_bbn = " << model1.get_xi_bbn() << std::endl;
    std::cout << "xi_fo = " << model1.get_xi_fo() << std::endl;
}
