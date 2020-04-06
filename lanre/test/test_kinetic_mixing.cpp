//
// Created by Logan Morrison on 3/15/20.
//

#include "lanre/dm_models/kinetic_mixing.hpp"
#include <gtest/gtest.h>
#include <vector>

using namespace lanre::dm_models;

class KineticMixingTest : public ::testing::Test {
protected:
    void SetUp() override {
    }

    KineticMixing km1{1.0, 1e3, 1.0, 1e-3};
    KineticMixing km2{10.0, 1e3, 1.0, 1e-3};
    KineticMixing km3{100.0, 1e3, 1.0, 1e-3};
    KineticMixing km4{1e3, 1e3, 1.0, 1e-3};
    KineticMixing km5{1e4, 1e3, 1.0, 1e-3};
};

TEST_F(KineticMixingTest, Widths) {
    km1.vector_mediator_width("all");
}

TEST_F(KineticMixingTest, CrossSections) {

}

TEST_F(KineticMixingTest, RelicDensity) {
    std::vector<double> rds(4, 0.0);

    rds[0] = km1.relic_density(1.0, 100.0, 1e-5, 1e-5, "radau");
    rds[1] = km2.relic_density(1.0, 100.0, 1e-5, 1e-5, "radau");
    rds[2] = km3.relic_density(1.0, 100.0, 1e-5, 1e-5, "radau");
    rds[3] = km4.relic_density(1.0, 100.0, 1e-5, 1e-5, "radau");

    std::cout << "Radau" << std::endl;
    for (auto &rd : rds) {
        std::cout << "rd = " << rd << std::endl;
    }

    rds[0] = km1.relic_density(1.0, 100.0, 1e-5, 1e-5, "rodas");
    rds[1] = km2.relic_density(1.0, 100.0, 1e-5, 1e-5, "rodas");
    rds[2] = km3.relic_density(1.0, 100.0, 1e-5, 1e-5, "rodas");
    rds[3] = km4.relic_density(1.0, 100.0, 1e-5, 1e-5, "rodas");

    std::cout << "Rodas" << std::endl;
    for (auto &rd : rds) {
        std::cout << "rd = " << rd << std::endl;
    }

    std::cout << "Gondolo-Gelmini" << std::endl;
    rds[0] = km1.relic_density(1.0, 100.0, 1e-5, 1e-5, "gg");
    std::cout << "rd = " << rds[0] << std::endl;
    rds[1] = km2.relic_density(1.0, 100.0, 1e-5, 1e-5, "gg");
    std::cout << "rd = " << rds[1] << std::endl;
    rds[2] = km3.relic_density(1.0, 100.0, 1e-5, 1e-5, "gg");
    std::cout << "rd = " << rds[2] << std::endl;
    rds[3] = km4.relic_density(1.0, 100.0, 1e-5, 1e-5, "gg");
    std::cout << "rd = " << rds[3] << std::endl;

}

TEST_F(KineticMixingTest, ComputeXfMPU) {
    double xf = km1.compute_xf_mpu();
}

int main(int argc, char *argv[]) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}