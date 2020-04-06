//
// Created by Logan Morrison on 3/20/20.
//

//
// Created by Logan Morrison on 2/4/20.
//

#include <lanre/dm_models/thermally_decoupled_model.hpp>
#include <lanre/cosmology/thermodynamic_particle.hpp>
#include <lanre/cosmology/standard_model.hpp>
#include "gtest/gtest.h"
#include <vector>
#include <cmath>

using namespace lanre;
using namespace lanre::cosmology;
using namespace lanre::dm_models;

class TestThermallyDecoupledModel : public ThermallyDecoupledModel {
public:
    ThermodynamicParticle eta;
    ThermodynamicParticle delta;
    ThermodynamicParticle dark_photon;

    double xi_inf{};
    double hd_inf{};
    double sum_g{};
    double gl{};
    double ml{};

    TestThermallyDecoupledModel(double t_lam, unsigned int t_N, double t_xi_inf, bool t_has_dp)
            : eta(ThermodynamicParticle{t_lam / sqrt(t_N), 1.0, 0}),
              delta(ThermodynamicParticle{t_lam * t_N, double(t_N + 1), t_N}),
              dark_photon(ThermodynamicParticle{0.0, (t_has_dp ? 2.0 : 0.0), 1}) {
        this->m_xi_inf = t_xi_inf;
        // quarks + gluons + dark photon
        m_hd_inf = (4.0 * t_N * 7.0 / 8.0) + (2.0 * (t_N * t_N - 1.0)) + dark_photon.get_g();
        m_sum_g = eta.get_g() + dark_photon.get_g() + delta.get_g() * (delta.get_spin2() % 2 == 0 ? 1.0 : 7.0 / 8.0);
        m_gl = (t_has_dp ? dark_photon.get_g() : eta.get_g());
        m_ml = (t_has_dp ? 0.0 : eta.get_mass());

        xi_inf = m_xi_inf;
        hd_inf = m_hd_inf;
        sum_g = m_sum_g;
        gl = m_gl;
        ml = m_ml;
    }

    double dark_heff(double Td) const override {
        return eta.heff(Td) + delta.heff(Td) + dark_photon.heff(Td);
    }

    double comp_xi_const_td(double Td) {
        return compute_xi_const_td(Td);
    }

    double comp_xi_const_tsm(double Td) {
        return compute_xi_const_tsm(Td);
    }

};


TEST(ComputeXiTest, TestConstantTDark) {

    constexpr double log_Td_min = -3.0;
    constexpr double log_Td_max = 3.0;
    constexpr int num_Tds = 1000;
    const double log_Td_step = (log_Td_max - log_Td_min) / double(num_Tds - 1);

    for (int i = 0; i < num_Tds; i++) {
        double Td = pow(10.0, log_Td_min + i * log_Td_step);
        auto model = TestThermallyDecoupledModel{1.0, 50, 1.0, true};

        double xi = model.comp_xi_const_td(Td);
        double Tsm = Td / xi;

        double hsm = sm_heff(Tsm);
        double hd = model.dark_heff(Td);

        double lhs = pow(xi, 3);
        double rhs = (hsm * model.hd_inf / hd / SM_HEFF_INF) * pow(model.xi_inf, 3);

        double residual = lhs - rhs;
        ASSERT_LE(fabs(residual), 1e-4);
    }

    for (int i = 0; i < num_Tds; i++) {
        double Td = pow(10.0, log_Td_min + i * log_Td_step);
        auto model = TestThermallyDecoupledModel{1e0, 50, 1e0, false};

        double xi = model.comp_xi_const_td(Td);
        double Tsm = Td / xi;

        double hsm = sm_heff(Tsm);
        double hd = model.dark_heff(Td);

        double rhs = (hsm * model.hd_inf / hd / SM_HEFF_INF) * pow(model.xi_inf, 3);
        double lhs = pow(xi, 3);

        double frac_diff = fabs(rhs - lhs) / rhs;

        ASSERT_LE(frac_diff, 1e-4);
    }
}

TEST(ComputeXiTest, TestConstantTDarkSmallXi) {

    constexpr double log_Td_min = -3.0;
    constexpr double log_Td_max = 3.0;
    constexpr int num_Tds = 1000;
    const double log_Td_step = (log_Td_max - log_Td_min) / double(num_Tds - 1);

    for (int i = 0; i < num_Tds; i++) {
        double Td = pow(10.0, log_Td_min + i * log_Td_step);
        auto model = TestThermallyDecoupledModel{1.0, 50, 1e-3, true};

        double xi = model.comp_xi_const_td(Td);
        double Tsm = Td / xi;

        double hsm = sm_heff(Tsm);
        double hd = model.dark_heff(Td);

        double residual = pow(xi, 3) - (hsm * model.hd_inf / hd / SM_HEFF_INF) * pow(model.xi_inf, 3);

        ASSERT_LE(fabs(residual), 1e-4);
    }

    for (int i = 0; i < num_Tds; i++) {
        double Td = pow(10.0, log_Td_min + i * log_Td_step);
        auto model = TestThermallyDecoupledModel{1e0, 50, 1e-3, false};

        double xi = model.comp_xi_const_td(Td);
        double Tsm = Td / xi;

        double hsm = sm_heff(Tsm);
        double hd = model.dark_heff(Td);

        double rhs = (hsm * model.hd_inf / hd / SM_HEFF_INF) * pow(model.xi_inf, 3);
        double lhs = pow(xi, 3);

        double frac_diff = fabs(rhs - lhs) / rhs;

        ASSERT_LE(frac_diff, 1e-4);
    }
}

TEST(ComputeXiTest, TestConstantTSM) {

    constexpr double log_Tsm_min = -3.0;
    constexpr double log_Tsm_max = 3.0;
    constexpr int num_Tsms = 1000;
    const double log_Tsm_step = (log_Tsm_max - log_Tsm_min) / double(num_Tsms - 1);

    for (int i = 0; i < num_Tsms; i++) {
        double Tsm = pow(10.0, log_Tsm_min + i * log_Tsm_step);
        auto model = TestThermallyDecoupledModel{1.0, 50, 1.0, true};

        double xi = model.comp_xi_const_tsm(Tsm);
        double Td = Tsm * xi;

        double hsm = sm_heff(Tsm);
        double hd = model.dark_heff(Td);

        double rhs = (hsm * model.hd_inf / hd / SM_HEFF_INF) * pow(model.xi_inf, 3);
        double lhs = pow(xi, 3);

        double frac_diff = fabs(rhs - lhs) / rhs;

        ASSERT_LE(frac_diff, 1e-4);
    }

    for (int i = 0; i < num_Tsms; i++) {
        double Tsm = pow(10.0, log_Tsm_min + i * log_Tsm_step);
        auto model = TestThermallyDecoupledModel{1e0, 50, 1e0, false};

        double xi = model.comp_xi_const_tsm(Tsm);
        double Td = Tsm * xi;

        double hsm = sm_heff(Tsm);
        double hd = model.dark_heff(Td);

        double rhs = (hsm * model.hd_inf / SM_HEFF_INF) * pow(model.xi_inf, 3);
        double lhs = hd * pow(xi, 3);

        double frac_diff = fabs(rhs - lhs) / rhs;

        ASSERT_LE(frac_diff, 1e-4);
    }
}

TEST(ComputeXiTest, TestConstantTSMSmallXi) {

    constexpr double log_Tsm_min = -3.0;
    constexpr double log_Tsm_max = 3.0;
    constexpr int num_Tsms = 1000;
    const double log_Tsm_step = (log_Tsm_max - log_Tsm_min) / double(num_Tsms - 1);

    for (int i = 0; i < num_Tsms; i++) {
        double Tsm = pow(10.0, log_Tsm_min + i * log_Tsm_step);
        auto model = TestThermallyDecoupledModel{1.0, 50, 1e-3, true};

        double xi = model.comp_xi_const_tsm(Tsm);
        double Td = Tsm * xi;

        double hsm = sm_heff(Tsm);
        double hd = model.dark_heff(Td);

        double rhs = (hsm * model.hd_inf / hd / SM_HEFF_INF) * pow(model.xi_inf, 3);
        double lhs = pow(xi, 3);

        double frac_diff = fabs(rhs - lhs) / rhs;

        ASSERT_LE(frac_diff, 1e-4);
    }

    for (int i = 0; i < num_Tsms; i++) {
        double Tsm = pow(10.0, log_Tsm_min + i * log_Tsm_step);
        auto model = TestThermallyDecoupledModel{1e0, 50, 1e-3, false};

        double xi = model.comp_xi_const_tsm(Tsm);
        double Td = Tsm * xi;

        double hsm = sm_heff(Tsm);
        double hd = model.dark_heff(Td);

        double rhs = (hsm * model.hd_inf / SM_HEFF_INF) * pow(model.xi_inf, 3);
        double lhs = hd * pow(xi, 3);

        double frac_diff = fabs(rhs - lhs) / rhs;

        ASSERT_LE(frac_diff, 1e-4);
    }
}


int main(int argc, char *argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
