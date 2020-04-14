//
// Created by Logan Morrison on 3/16/20.
//

#ifndef LANRE_DM_MODELS_CONSTANT_THERMAL_CROSS_SECTION_HPP
#define LANRE_DM_MODELS_CONSTANT_THERMAL_CROSS_SECTION_HPP

#include "lanre/diffeq/base.hpp"
#include "lanre/diffeq/radau.hpp"
#include "lanre/diffeq/function.hpp"
#include "lanre/cosmology/standard_model.hpp"
#include "lanre/cosmology/thermodynamic_particle.hpp"
#include "lanre/constants.hpp"
#include <memory>

namespace lanre {
namespace dm_models {

using namespace diffeq;
using namespace cosmology;

class ConstantThermalCrossSection {
private:
    double m_mx;
    double m_sigmav;
    int m_n;
public:
    ConstantThermalCrossSection(double mx, double sigmav, int n) : m_mx(mx), m_sigmav(sigmav), m_n(n) {}

    double get_mx() const {
        return m_mx;
    }

    double get_sigmav() const {
        return m_sigmav;
    }

    int get_n() const {
        return m_n;
    }

    void set_mx(double mx) {
        m_mx = mx;
    }

    void set_sigmav(double sigmav) {
        m_sigmav = sigmav;
    }

    void set_n(int n) {
        m_n = n;
    }

    double relic_density(double, double, double, double);
};

struct ConstantThermalCrossSectionBoltzmann : public ODEFunction {
    std::shared_ptr<ConstantThermalCrossSection> model;
    ThermodynamicParticle chi;

    ConstantThermalCrossSectionBoltzmann(
            std::shared_ptr<ConstantThermalCrossSection> t_model
    ) : model(t_model), chi{model->get_mx(), 2.0, 1} {}

    void dudt(Vector<double> &dw, const Vector<double> &w, double logx) override {
        double x = exp(logx);
        double T = chi.get_mass() / x;
        double s = cosmology::sm_entropy_density(T);

        double weq = log(chi.neq(T) / s);
        double ww = w(0);

        double pf = -sqrt(M_PI / 45) * kM_PLANK * cosmology::sm_sqrt_gstar(T) * T;
        double sigmav = model->get_sigmav();

        // dW_e / dlogx
        dw(0) = pf * sigmav * (exp(ww) - exp(2.0 * weq - ww));
    }

    void dfdu(Matrix<double> &J, const Vector<double> &w, double logx) override {
        double x = exp(logx);
        double T = chi.get_mass() / x;
        double s = cosmology::sm_entropy_density(T);

        double weq = log(chi.neq(T) / s);
        double ww = w(0);

        double pf = -sqrt(M_PI / 45) * kM_PLANK * cosmology::sm_sqrt_gstar(T) * T;
        double sigmav = model->get_sigmav();

        // dW_e / dlogx
        J(0, 0) = pf * sigmav * (exp(ww) + exp(2.0 * weq - ww));
    }
};

double ConstantThermalCrossSection::relic_density(
        const double xstart,
        const double xend,
        const double reltol,
        const double abstol
) {

    using namespace diffeq;
    using namespace cosmology;

    auto logx_span = std::make_pair(std::log(xstart), std::log(xend));

    ConstantThermalCrossSectionBoltzmann boltz{std::make_shared<ConstantThermalCrossSection>(*this)};
    double Tinit = boltz.chi.get_mass() / xstart;

    Vector<double> winit{1};
    winit(0) = log(boltz.chi.neq(Tinit) / sm_entropy_density(Tinit));

    ODEProblem problem{std::make_shared<ConstantThermalCrossSectionBoltzmann>(boltz), winit, logx_span};
    Radau5 alg{};
    ODEIntegratorOptions opts{};
    opts.abstol = abstol;
    opts.reltol = reltol;

    auto sol = solve(problem, alg);

    double yinf = exp(sol.us.back()[0]);
    return boltz.chi.get_mass() * yinf * kS_TODAY / kRHO_CRIT;
}

}
}

#endif //LANRE_DM_MODELS_CONSTANT_THERMAL_CROSS_SECTION_HPP
