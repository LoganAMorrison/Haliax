//
// Created by Logan Morrison on 4/12/20.
//

#ifndef LANRE_KINETIC_MIXING_RELIC_DENSITY_MPU_HPP
#define LANRE_KINETIC_MIXING_RELIC_DENSITY_MPU_HPP

#include "lanre/cosmology/standard_model.hpp"
#include "lanre/cosmology/thermodynamic_functions.hpp"
#include "lanre/integrate/qagi.hpp"

#include "lanre/dm_models/kinetic_mixing/parameters.hpp"
#include "lanre/dm_models/kinetic_mixing/thermal_cross_section.hpp"

#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/math/differentiation/finite_difference.hpp>

namespace lanre {
namespace dm_models {
namespace kinetic_mixing {

/**
 * Compute the residual equation for x_f using Morrison, Patel, Ulbricht method.
 * @param xf Freeze-out value for x = mass / T.
 * @return residual
 */
double xf_root_equation_mpu(const Parameters &params, double xf) {
    using cosmology::neq;
    using boost::math::differentiation::finite_difference_derivative;

    const double euler_gamma = 0.5772156649015328606065120;
    const double Tf = params.mx / xf;
    const double pf = sqrt(M_PI / 45.0) * params.mx * kM_PLANK;
    const double fx = (pf * cosmology::sm_sqrt_gstar(params.mx / xf) / (xf * xf) * thermal_cross_section(params, xf));
    const double yeq = neq(Tf, params.mx, 2.0, 1) / cosmology::sm_entropy_density(Tf);
    const double Qx = fx * yeq;

    // Compute f'(xf) / f(xf) = d/dx log(f(x))
    auto logf = [params](double x) {
        // Ignoring constant pf since it won't matter for derivative
        return log(cosmology::sm_sqrt_gstar(params.mx / x) / (x * x) * thermal_cross_section(params, x));
    };
    const double dlogfx = finite_difference_derivative(logf, xf);

    return Qx - (2.0 * exp(-(euler_gamma + dlogfx + 3.0 / (2.0 * xf))));
}

/**
 * Compute the freeze-out value of x = mass / T using the Morrison, Patel, Ulbricht
 * method.
 * @return xf
 */
double compute_xf_mpu(const Parameters &params) {
    using namespace boost::math::tools;
    auto f = [params](double xf) {
        return xf_root_equation_mpu(params, xf);
    };
    auto tol = [](double min, double max) { return abs(max - min) <= 1e-8; };
    std::pair<double, double> res = bisect(f, 1e-2, 100.0, tol);
    return (res.second + res.first) / 2.0;
}

/**
 * Compute the relic density using the Morrison, Patel, Ulbricht method
 * @return
 */
double relic_density_mpu(const Parameters &params) {
    using boost::math::quadrature::gauss_kronrod;

    const double xf = compute_xf_mpu(params);
    const double xinf = std::min(100.0 * xf, 350.0);
    const double pf = sqrt(M_PI / 45.0) * params.mx * kM_PLANK;

    auto f = [params, pf](double x) {
        return pf * cosmology::sm_sqrt_gstar(params.mx / x) / (x * x) * thermal_cross_section(params, x);
    };

    const double integal_f = gauss_kronrod<double, 15>::integrate(f, xf, std::numeric_limits<double>::infinity());
    const double Y0 = 1 / integal_f;

    return Y0 * params.mx * kS_TODAY / kRHO_CRIT;
}

}
}
}

#endif //LANRE_KINETIC_MIXING_RELIC_DENSITY_MPU_HPP
