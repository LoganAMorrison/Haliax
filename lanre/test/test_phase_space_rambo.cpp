//
// Created by Logan Morrison on 3/24/20.
//

#include "lanre/phase_space/rambo.hpp"
#include "lanre/constants.hpp"
#include <boost/math/special_functions/pow.hpp>
#include <iostream>
#include "gtest/gtest.h"

using namespace lanre;
using namespace lanre::phase_space;

double mat_squared_mu_to_e_nu_nu(const std::vector<FourVector> &momenta) {
    using namespace boost::math;
    double s = scalar_product(momenta[1] + momenta[2], momenta[1] + momenta[2]);
    double t = scalar_product(momenta[0] + momenta[2], momenta[0] + momenta[2]);
    return -16 * kG_FERMI * kG_FERMI * (s + t) * (s + t - kMUON_MASS * kMUON_MASS);
}

double mat_squared_ttbar_higgs_ttbar(const std::vector<FourVector> &momenta) {
    double Q = momenta[0].t + momenta[1].t;
    FourVector p1{Q / 2, 0.0, 0.0, sqrt(Q * Q / 4.0 - kTOP_QUARK_MASS * kTOP_QUARK_MASS)};

    double s = Q * Q;
    double t = scalar_product(p1 - momenta[0], p1 - momenta[0]);

    return (9 * pow(kALPHA_EM, 2) * pow(M_PI, 2) * pow(kTOP_QUARK_MASS, 4) * pow(s - 4 * pow(kTOP_QUARK_MASS, 2), 2)) /
            (pow(s - pow(kHIGGS_MASS, 2), 2) * pow(kW_BOSON_MASS, 4) * pow(kSIN_THETA_WEAK, 4)) +
            (9 * pow(kALPHA_EM, 2) * pow(M_PI, 2) * pow(kTOP_QUARK_MASS, 4) * pow(t - 4 * pow(kTOP_QUARK_MASS, 2), 2)) /
                    (pow(t - pow(kHIGGS_MASS, 2), 2) * pow(kW_BOSON_MASS, 4) * pow(kSIN_THETA_WEAK, 4)) -
            (3 * pow(kALPHA_EM, 2) * pow(M_PI, 2) * pow(kTOP_QUARK_MASS, 4) *
                    (s * t + 4 * (s + t) * pow(kTOP_QUARK_MASS, 2) - 16 * pow(kTOP_QUARK_MASS, 4))) /
                    ((s - pow(kHIGGS_MASS, 2)) * (-t + pow(kHIGGS_MASS, 2)) * pow(kW_BOSON_MASS, 4) *
                            pow(kSIN_THETA_WEAK, 4));
}

/**
 * Test that momentum is conserved and that the masses are correct.
 */
TEST(TestRambo, TestMomentumConservationAndMasses) {
    std::vector<double> isp_masses = {1.0, 2.0};
    std::vector<double> fsp_masses = {3.0, 4.0};
    double cme = 10.0;

    Rambo rambo{isp_masses, fsp_masses, cme};
    auto events = rambo.generate_events(100);

    for (auto &event:events) {
        FourVector sum = event.momenta[0] + event.momenta[1];
        ASSERT_NEAR(sum.t, cme, 1e-3);
        ASSERT_NEAR(sum.x, 0.0, 1e-3);
        ASSERT_NEAR(sum.y, 0.0, 1e-3);
        ASSERT_NEAR(sum.z, 0.0, 1e-3);

        ASSERT_NEAR(mass(event.momenta[0]), fsp_masses[0], 1e-3);
        ASSERT_NEAR(mass(event.momenta[1]), fsp_masses[1], 1e-3);
    }
}

/**
 * Check the 2->2 proccess with top quark annihilation into t tbar through the
 * SM Higgs boson.
 */
TEST(TestRambo, TestTopQuarkAnnihilationHiggs) {
    using namespace boost::math;

    double Q = 10.0 * kTOP_QUARK_MASS;
    std::vector<double> isp_masses = {kTOP_QUARK_MASS, kTOP_QUARK_MASS};
    std::vector<double> fsp_masses = {kTOP_QUARK_MASS, kTOP_QUARK_MASS};
    Rambo rambo{isp_masses, fsp_masses, Q, mat_squared_ttbar_higgs_ttbar};

    double s = Q * Q;
    double cs = (
            (3 * pow(kALPHA_EM, 2) * M_PI * pow(kTOP_QUARK_MASS, 4) *
                    ((log(pow(kHIGGS_MASS, 2) / (s + pow(kHIGGS_MASS, 2) - 4 * pow(kTOP_QUARK_MASS, 2))) *
                            (s - pow(kHIGGS_MASS, 2)) *
                            (-6 * pow(kHIGGS_MASS, 4) + 7 * pow(kHIGGS_MASS, 2) * (s + 4 * pow(kTOP_QUARK_MASS, 2)) -
                                    4 * (5 * s * pow(kTOP_QUARK_MASS, 2) + 4 * pow(kTOP_QUARK_MASS, 4)))) /
                            (s - 4 * pow(kTOP_QUARK_MASS, 2)) +
                            (6 * pow(kHIGGS_MASS, 8) + 48 * pow(s, 2) * pow(kTOP_QUARK_MASS, 4) -
                                    10 * pow(kHIGGS_MASS, 6) * (s + 4 * pow(kTOP_QUARK_MASS, 2)) + pow(kHIGGS_MASS, 4) *
                                    (3 * pow(s, 2) + 52 * s * pow(kTOP_QUARK_MASS, 2) + 112 * pow(kTOP_QUARK_MASS, 4)) +
                                    pow(kHIGGS_MASS, 2) * (7 * pow(s, 3) - 72 * pow(s, 2) * pow(kTOP_QUARK_MASS, 2) +
                                            32 * s * pow(kTOP_QUARK_MASS, 4) - 192 * pow(kTOP_QUARK_MASS, 6))) /
                                    (pow(kHIGGS_MASS, 2) * (s + pow(kHIGGS_MASS, 2) - 4 * pow(kTOP_QUARK_MASS, 2))))) /
                    (16. * s * pow(s - pow(kHIGGS_MASS, 2), 2) * pow(kW_BOSON_MASS, 4) * pow(kSIN_THETA_WEAK, 4))
    );

    auto res = rambo.compute_width_cross_section(1000000);
    std::cout << "avg = " << std::get<0>(res) << " +- " << std::get<1>(res) << std::endl;
    std::cout << "actual = " << cs << std::endl;

    ASSERT_NEAR(std::get<0>(res), cs, 3.0 * std::get<1>(res));
}

/**
 * Check 1->3 proccess with of muon decay into an electron and two neutrinos.
 */
TEST(TestRambo, TestMuonDecay) {
    using namespace boost::math;

    std::vector<double> isp_masses = {kMUON_MASS};
    std::vector<double> fsp_masses = {kELECTRON_MASS, 0.0, 0.0};
    Rambo rambo{isp_masses, fsp_masses, kMUON_MASS, mat_squared_mu_to_e_nu_nu};

    double width = (pow(kG_FERMI, 2) * pow(kMUON_MASS, 5)) / (192.0 * pow(M_PI, 3));

    auto res = rambo.compute_width_cross_section(100000);
    std::cout << "avg = " << std::get<0>(res) << " +- " << std::get<1>(res) << std::endl;
    std::cout << "actual = " << width << std::endl;

    ASSERT_NEAR(std::get<0>(res), width, 3.0 * std::get<1>(res));
}

