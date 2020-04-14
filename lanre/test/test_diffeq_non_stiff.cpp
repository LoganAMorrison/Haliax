//
// Created by Logan Morrison on 3/15/20.
//
#include "lanre/diffeq/dormand_prince.hpp"
#include "lanre/diffeq/problem.hpp"
#include "lanre/diffeq/integrator.hpp"
#include "gtest/gtest.h"

using namespace lanre;
using namespace lanre::diffeq;

class HarmonicOscillator : public ODEFunction {
private:
    double m_m;
    double m_k;
public:

    HarmonicOscillator(double t_m, double t_k) : ODEFunction(), m_m(t_m), m_k(t_k) {}

    ~HarmonicOscillator() = default;

    void dudt(Vector<double> &dx, const Vector<double> &x, const double /* t */) override {
        dx(0) = x(1);
        dx(1) = -m_k / m_m * x(0);
    }
};

TEST(DP5Test, HarmonicOscillatorTest) {
    double m = 1.0;
    double k = 1.0;

    HarmonicOscillator ho{m, k};
    double x0 = 1.0;
    double v0 = 0.0;
    double omega = sqrt(k / m);

    Vector<double> xinit(2);
    xinit(0) = x0;
    xinit(1) = v0;

    std::pair<double, double> tspan{0.0, 10.0};

    ODEProblem problem{std::make_shared<HarmonicOscillator>(ho), xinit, tspan};
    DormandPrince5 alg;
    ODEIntegratorOptions opts{};
    opts.reltol = 1e-3;
    opts.abstol = 1e-5;
    opts.dense = false;

    auto solution = solve(problem, alg, opts);

    auto analytic = [v0, x0, omega](double t_t) {
        return (v0 * sin(omega * t_t) + x0 * cos(omega * t_t)) / omega;
    };

    for (size_t i = 0; i < solution.ts.size(); i++) {
        double sol_analytic = analytic(solution.ts[i]);
        double sol_approx = solution.us[i](0);
        if (sol_analytic != 0.0 and sol_approx != 0) {
            double frac_diff = std::abs(sol_analytic - sol_approx) / sol_analytic;
            std::cout << "aprrox = " << sol_approx << ", analytic = " << sol_analytic << std::endl;
            ASSERT_LE(frac_diff, 2e-1);
        }
    }
}

TEST(DP8Test, HarmonicOscillatorTest) {
    double m = 1.0;
    double k = 1.0;

    HarmonicOscillator ho{m, k};
    double x0 = 1.0;
    double v0 = 0.0;
    double omega = sqrt(k / m);

    Vector<double> xinit(2);
    xinit(0) = x0;
    xinit(1) = v0;

    std::pair<double, double> tspan{0.0, 10.0};

    ODEProblem problem{std::make_shared<HarmonicOscillator>(ho), xinit, tspan};
    DormandPrince8 alg;
    ODEIntegratorOptions opts{};
    opts.reltol = 1e-5;
    opts.abstol = 1e-7;
    opts.dense = false;

    auto solution = solve(problem, alg, opts);

    auto analytic = [v0, x0, omega](double t_t) {
        return (v0 * sin(omega * t_t) + x0 * cos(omega * t_t)) / omega;
    };

    for (size_t i = 0; i < solution.ts.size(); i++) {
        double sol_analytic = analytic(solution.ts[i]);
        double sol_approx = solution.us[i](0);
        if (sol_analytic != 0.0 and sol_approx != 0) {
            double frac_diff = std::abs(sol_analytic - sol_approx) / sol_analytic;
            std::cout << "aprrox = " << sol_approx << ", analytic = " << sol_analytic << std::endl;
            ASSERT_LE(frac_diff, 2e-1);
        }
    }
}

int main(int argc, char *argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}