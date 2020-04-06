#include "lanre/diffeq/base.hpp"
#include "lanre/diffeq/decsol.hpp"
#include "lanre/diffeq/integrator.hpp"
#include "lanre/diffeq/problem.hpp"
#include "lanre/diffeq/algorithm.hpp"
#include "lanre/diffeq/function.hpp"
#include "lanre/diffeq/solution.hpp"
#include "lanre/diffeq/radau.hpp"
#include "lanre/diffeq/rodas.hpp"
#include "lanre/diffeq/seulex.hpp"
#include "lanre/diffeq/seulex2.hpp"
#include "lanre/cosmology/standard_model.hpp"
#include "lanre/cosmology/thermodynamic_particle.hpp"
#include "lanre/constants.hpp"
#include "gtest/gtest.h"

using namespace lanre;
using namespace lanre::diffeq;
using namespace lanre::cosmology;

static const std::vector<std::pair<double, double>> mass_sigmavs{
        std::make_pair(10.313897683787216, 1.966877938634266e-9),
        std::make_pair(104.74522360006331, 1.7597967261428258e-9),
        std::make_pair(1063.764854316313, 1.837766552668581e-9),
        std::make_pair(10000.0, 1.8795945459427076e-9),
};

struct DMModel : public ODEFunction {
    ThermodynamicParticle chi;
    double sigmav;

    DMModel(double t_mass, double t_sigmav) :
            ODEFunction(), chi(ThermodynamicParticle{t_mass, 2.0, 1}),
            sigmav(t_sigmav) {}

    ~DMModel() = default;

    void dudt(Vector<double> &dw, const Vector<double> &w, const double logx) override {
        double x = exp(logx);
        double T = chi.get_mass() / x;
        double s = sm_entropy_density(T);

        double weq = log(chi.neq(T) / s);
        double ww = w(0);

        double pf = -sqrt(M_PI / 45) * kM_PLANK * sm_sqrt_gstar(T) * T;

        // dW_e / dlogx
        dw(0) = pf * sigmav * (exp(ww) - exp(2.0 * weq - ww));
    }


    void dfdu(Matrix<double> &dfdw, const Vector<double> &w, const double logx) override {
        double x = exp(logx);
        double T = chi.get_mass() / x;
        double s = sm_entropy_density(T);

        double weq = log(chi.neq(T) / s);
        double ww = w(0);

        double pf = -sqrt(M_PI / 45) * kM_PLANK * sm_sqrt_gstar(T) * T;

        dfdw(0, 0) = pf * sigmav * (exp(ww) + exp(2.0 * weq - ww));
    }

};

struct VanDerPol : public ODEFunction {
    double mu;

    explicit VanDerPol(double t_mu) : ODEFunction(), mu(t_mu) {}

    void dudt(Vector<double> &du, const Vector<double> &u, const double) override {
        du(0) = u(1);
        du(1) = ((1 - u(0) * u(0)) * u(1) - u(0)) / mu;
    }

    void dfdu(Matrix<double> &df, const Vector<double> &u, const double) override {
        df(0, 0) = 0.0;
        df(0, 1) = 1.0;
        df(1, 0) = -(2 * u(0) * u(1) + 1.0) / mu;
        df(1, 1) = (1.0 - u(0) * u(0)) / mu;
    }

    void dfdt(Vector<double> &df, const Vector<double> &, const double) override {
        df(0) = 0.0;
        df(1) = 0.0;
    }
};

TEST(RadauTest, BoltzmannTest) {
    for (const auto &mass_sigmav: mass_sigmavs) {
        auto logx_span = std::make_pair(0.0, std::log(500.0));
        double abstol = 1e-6;
        double reltol = 1e-6;

        DMModel model{mass_sigmav.first, mass_sigmav.second};
        double Tinit = model.chi.get_mass() / exp(logx_span.first);

        Vector<double> winit{1};
        winit(0) = log(model.chi.neq(Tinit) / sm_entropy_density(Tinit));

        ODEProblem problem{std::make_shared<DMModel>(model), winit, logx_span};
        Radau5 alg{};
        ODEIntegratorOptions opts{};
        opts.abstol = 1e-6;
        opts.reltol = 1e-6;

        auto sol = solve(problem, alg);

        for (size_t i = 0; i < sol.ts.size(); i++) {
            std::cout << "t = " << sol.ts[i] << ", u = " << sol.us[i] << std::endl;
        }
        std::cout << "steps = " << sol.ts.size() << std::endl;

        double yinf = exp(sol.us.back()[0]);
        double rd = model.chi.get_mass() * yinf * kS_TODAY / kRHO_CRIT;
        double frac_diff = std::abs(rd - kOMEGA_H2_CDM) / kOMEGA_H2_CDM;

        std::cout << "rd = " << rd << std::endl;

        ASSERT_LE(frac_diff, 1e-1);
    }
}

TEST(RadauTest, VanDerPolTest) {
    auto tSpan = std::make_pair(0.0, 2.0);
    double reltol = 1e-7;
    double abstol = 1.0 * reltol;
    VanDerPol model{1e-6};

    Vector<double> uInit(2);
    uInit(0) = 2.0;
    uInit(1) = -0.66;

    ODEProblem problem{std::make_shared<VanDerPol>(model), uInit, tSpan};
    Radau5 alg;
    ODEIntegratorOptions opts{};
    opts.reltol = 1e-7;
    opts.abstol = opts.reltol;
    auto sol = solve(problem, alg, opts);

    for (size_t i = 0; i < sol.us.size(); i++) {
        std::cout << sol.ts[i] << ", " << sol.us[i][0] << ", " << sol.us[i][1] << std::endl;
    }

    double this_res1 = sol.us.back()[0];
    double this_res2 = sol.us.back()[1];

    double their_res1 = 0.1706157223e1;
    double their_res2 = -0.8928209476;

    double frac_err1 = std::abs((this_res1 - their_res1) / their_res1);
    double frac_err2 = std::abs((this_res2 - their_res2) / their_res2);


    ASSERT_LE(frac_err1, 1e-2);
    ASSERT_LE(frac_err2, 1e-2);

}

TEST(RodasTest, BoltzmannTest) {
    for (const auto &mass_sigmav: mass_sigmavs) {
        auto logx_span = std::make_pair(0.0, std::log(500.0));
        double abstol = 1e-6;
        double reltol = 1e-6;

        DMModel model{mass_sigmav.first, mass_sigmav.second};
        double Tinit = model.chi.get_mass() / exp(logx_span.first);

        Vector<double> winit{1};
        winit(0) = log(model.chi.neq(Tinit) / sm_entropy_density(Tinit));

        ODEProblem problem{std::make_shared<DMModel>(model), winit, logx_span};
        Rodas alg{};
        ODEIntegratorOptions opts{};
        opts.abstol = 1e-6;
        opts.reltol = 1e-6;

        auto sol = solve(problem, alg);

        for (size_t i = 0; i < sol.ts.size(); i++) {
            std::cout << "t = " << sol.ts[i] << ", u = " << sol.us[i] << std::endl;
        }
        std::cout << "steps = " << sol.ts.size() << std::endl;

        double yinf = exp(sol.us.back()[0]);
        double rd = model.chi.get_mass() * yinf * kS_TODAY / kRHO_CRIT;
        double frac_diff = std::abs(rd - kOMEGA_H2_CDM) / kOMEGA_H2_CDM;

        std::cout << "rd = " << rd << std::endl;

        ASSERT_LE(frac_diff, 1e-1);
    }
}

TEST(RodasTest, VanDerPolTest) {
    auto tSpan = std::make_pair(0.0, 2.0);
    double reltol = 1e-4;
    double abstol = 1e-6 * reltol;
    VanDerPol model{1e-6};

    Vector<double> uInit(2);
    uInit(0) = 2.0;
    uInit(1) = -0.66;

    ODEProblem problem{std::make_shared<VanDerPol>(model), uInit, tSpan};
    Rodas alg;
    ODEIntegratorOptions opts{};
    opts.reltol = 1e-4;
    opts.abstol = 1e-6 * opts.reltol;
    auto sol = solve(problem, alg, opts);

    for (size_t i = 0; i < sol.us.size(); i++) {
        std::cout << sol.ts[i] << ", " << sol.us[i][0] << ", " << sol.us[i][1] << std::endl;
    }

    double this_res1 = sol.us.back()[0];
    double this_res2 = sol.us.back()[1];

    double their_res1 = 0.1706157223e1;
    double their_res2 = -0.8928209476;

    double frac_err1 = std::abs((this_res1 - their_res1) / their_res1);
    double frac_err2 = std::abs((this_res2 - their_res2) / their_res2);


    ASSERT_LE(frac_err1, 1e-2);
    ASSERT_LE(frac_err2, 1e-2);

}

TEST(Seulex, BoltzmannTest) {
    for (const auto &mass_sigmav: mass_sigmavs) {
        auto logx_span = std::make_pair(0.0, std::log(500.0));
        double abstol = 1e-6;
        double reltol = 1e-6;

        DMModel model{mass_sigmav.first, mass_sigmav.second};
        double Tinit = model.chi.get_mass() / exp(logx_span.first);

        Vector<double> winit{1};
        winit(0) = log(model.chi.neq(Tinit) / sm_entropy_density(Tinit));

        ODEProblem problem{std::make_shared<DMModel>(model), winit, logx_span};
        Seulex2 alg{};
        ODEIntegratorOptions opts{};
        opts.abstol = 1e-6;
        opts.reltol = 1e-6;

        auto sol = solve(problem, alg);

        for (size_t i = 0; i < sol.ts.size(); i++) {
            std::cout << "t = " << sol.ts[i] << ", u = " << sol.us[i] << std::endl;
        }
        std::cout << "steps = " << sol.ts.size() << std::endl;

        double yinf = exp(sol.us.back()[0]);
        double rd = model.chi.get_mass() * yinf * kS_TODAY / kRHO_CRIT;
        double frac_diff = std::abs(rd - kOMEGA_H2_CDM) / kOMEGA_H2_CDM;

        std::cout << "rd = " << rd << std::endl;

        ASSERT_LE(frac_diff, 1e-1);
    }
}