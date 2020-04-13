//
// Created by Logan Morrison on 4/11/20.
//

#include "lanre/lanre.hpp"
#include "lanre/diffeq/base.hpp"
#include "lanre/diffeq/radau.hpp"
#include "lanre/diffeq/rodas.hpp"
#include "lanre/diffeq/function.hpp"
#include "lanre/cosmology/standard_model.hpp"
#include "lanre/cosmology/thermodynamic_particle.hpp"

#include "lanre/dm_models/kinetic_mixing/parameters.hpp"
#include "lanre/dm_models/kinetic_mixing/thermal_cross_section.hpp"

#ifndef LANRE_KINETIC_MIXING_BOLTZMANN_HPP
#define LANRE_KINETIC_MIXING_BOLTZMANN_HPP

namespace lanre {
namespace dm_models {
namespace kinetic_mixing {

struct KineticMixingBoltzman : public diffeq::ODEFunction {
    Parameters params;

    KineticMixingBoltzman(Parameters params) : params(params) {}

    void dudt(Vector<double> &dw, const Vector<double> &w, const double logx) override {
        double x = exp(logx);
        double T = params.mx / x;
        double s = cosmology::sm_entropy_density(T);

        double weq = log(cosmology::neq(T, params.mx, 2.0, 1) / s);
        double ww = w(0);

        double pf = -sqrt(M_PI / 45) * kM_PLANK * cosmology::sm_sqrt_gstar(T) * T;
        double sigmav = thermal_cross_section(params, x, "all", "all");

        // dW_e / dlogx
        dw(0) = pf * sigmav * (exp(ww) - exp(2.0 * weq - ww));
    }

    void dfdu(Matrix<double> &J, const Vector<double> &w, const double logx) override {
        double x = exp(logx);
        double T = params.mx / x;
        double s = cosmology::sm_entropy_density(T);

        double weq = log(cosmology::neq(T, params.mx, 2.0, 1) / s);
        double ww = w(0);

        double pf = -sqrt(M_PI / 45) * kM_PLANK * cosmology::sm_sqrt_gstar(T) * T;
        double sigmav = thermal_cross_section(params, x, "all", "all");

        // dW_e / dlogx
        J(0, 0) = pf * sigmav * (exp(ww) + exp(2.0 * weq - ww));
    }
};

/*
 * Compute the relic density of the dark matter.
 */
diffeq::ODESolution solve_boltzmann(
        const Parameters &params,
        double xstart = 1.0,
        double xend = 500.0,
        double reltol = 1e-3,
        double abstol = 1e-9,
        const std::string &t_alg = "radau"
) {
    using namespace diffeq;
    using namespace cosmology;

    auto logx_span = std::make_pair(std::log(xstart), std::log(xend));

    KineticMixingBoltzman boltz{params};
    double Tinit = params.mx / exp(logx_span.first);

    Vector<double> winit{1};
    winit(0) = log(neq(Tinit, params.mx, 2.0, 1) / sm_entropy_density(Tinit));

    ODEProblem problem{std::make_shared<KineticMixingBoltzman>(boltz), winit, logx_span};

    if (t_alg == "rodas") {
        Rodas alg{};
        ODEIntegratorOptions opts{};
        opts.abstol = abstol;
        opts.reltol = reltol;

        return solve(problem, alg, opts);
    } else {
        Radau5 alg{};
        ODEIntegratorOptions opts{};
        opts.abstol = abstol;
        opts.reltol = reltol;

        return solve(problem, alg, opts);
    }
}

double relic_density(
        const Parameters &params,
        double xstart = 1.0,
        double xend = 500.0,
        double reltol = 1e-3,
        double abstol = 1e-9,
        const std::string &t_alg = "radau"
) {
    auto sol = solve_boltzmann(params, xstart, xend, reltol, abstol, t_alg);
    double yinf = exp(sol.us.back()[0]);
    return params.mx * yinf * kS_TODAY / kRHO_CRIT;
}

}
}
}

#endif //LANRE_KINETIC_MIXING_BOLTZMANN_HPP
