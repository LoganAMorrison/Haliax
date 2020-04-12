//
// Created by Logan Morrison on 3/29/20.
//

#include "lanre/dm_models/dipole_dm.hpp"
#include <boost/math/special_functions/bessel.hpp>
#include "gtest/gtest.h"
#include <iostream>

using namespace lanre;
using namespace lanre::dm_models;


TEST(TestDipoleDM, TestTCS) {
    DipoleDM model{1e3, 1.5e3, 1e-3, 1e-3, 1e5};

    size_t num_xs = 100;
    double logxmin = -1;
    double logxmax = 3;
    double logxstep = (logxmax - logxmin) / double(num_xs - 1);
    for (int i = 0; i < num_xs; i++) {
        double x = pow(10.0, logxmin + i * logxstep);
        double tcs = model.thermal_cross_section(x);
        std::cout << "(x,tcs) = (" << x << "," << tcs << ")" << std::endl;
    }

}


TEST(TestDipoleDM, TestTCS2) {
    DipoleDM model{1e3, 1.5e3, 1e-3, 1e-3, 1e5};

    size_t num_xs = 100;
    double logxmin = -1;
    double logxmax = 3;
    double logxstep = (logxmax - logxmin) / double(num_xs - 1);
    for (int i = 0; i < num_xs; i++) {
        double x = pow(10.0, logxmin + i * logxstep);
        double tcs = model.thermal_cross_section2(x);
        std::cout << "(x,tcs) = (" << x << "," << tcs << ")" << std::endl;
    }

}


TEST(TestDipoleDM, TestRD) {
    double m1 = 1e3;
    double m2 = m1 + 1e-3;
    DipoleDM model{m1, m2, 0.0, 2e-3, 1e5};

    auto sol = model.solve_boltzmann(1.0, 5e3);
    if (sol.retcode == Success) {
        std::cout << "Success!" << std::endl;
    }

    for (size_t i = 0; i < sol.ts.size(); i++) {
        std::cout << sol.ts[i] << ", " << sol.us[i] << std::endl;
    }
    std::cout << "rd = " << model.relic_density(1.0, 1000.0) << std::endl;
}


TEST(TestDipoleDM, TestGamma) {
    double m1 = 100.0;
    double delta = 1.0 / 100.0;
    double ce = 1.0;
    double cm = 1.0;
    double lam = 1e4;
    double T = m1 / 1000.0;
    DipoleDM model{m1, m1 + delta, ce, cm, lam};

    std::vector<double> ws = {
            0.01, 0.1, 1.0, 10.0, 100.0
    };

    for (const auto &w: ws) {
        double gamma_integrand = model.gamma_integrand(w, T);
        std::cout << w << ", " << gamma_integrand << std::endl;
    }

    std::cout << model.gamma(T) << std::endl;
}

TEST(TestDipoleDM, TestTempSolve) {
    double m1 = 1e4;
    double m2 = m1 + 1e-3;
    DipoleDM model{m1, m2, 0.0, 2e-3, 1e6};
    double xstart = 20.0;
    double xend = m1 / kT_CMB;

    auto sol = model.solve_temperature(xstart, xend, 1e-3, 1e-9);
    if (sol.retcode == Success) {
        std::cout << "Success!" << std::endl;
    }

    std::cout << std::setprecision(15) << std::endl;
    for (size_t i = 0; i < sol.ts.size(); i++) {
        std::cout << sol.ts[i] << ", " << sol.us[i] << std::endl;
    }
}


int main(int argc, char *argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}