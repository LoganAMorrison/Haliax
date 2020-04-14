//
// Created by Logan Morrison on 4/11/20.
//

#ifndef LANRE_KINETIC_MIXING_GONDOLO_GELMINI_HPP
#define LANRE_KINETIC_MIXING_GONDOLO_GELMINI_HPP

#include "lanre/autodiff/dual.hpp"
#include "lanre/cosmology/standard_model.hpp"
#include "lanre/cosmology/thermodynamic_functions.hpp"
#include "lanre/integrate/quad.hpp"

#include "lanre/dm_models/kinetic_mixing/parameters.hpp"
#include "lanre/dm_models/kinetic_mixing/thermal_cross_section.hpp"

#include <cmath>
#include <boost/math/tools/roots.hpp>

namespace lanre {
namespace dm_models {
namespace kinetic_mixing {

/**
 * Returns residual of root equation used to solve for x_star.
 *
 * See Eqn.(14) of arXiv:1204.3622v3 for similar expressions. Note that our
 * result is more exact since we do not assume `xstar` is large enough that
 * Yeq ~ x^{3/2} e^{-x} / h_sm. This may cause a bit of a slow down.
 *
 * @param params Model parameters
 * @param chi Dark matter particle
 * @param xstar Current value of x_star, the value of mass / temperature at
 *              which DM begins to freeze out.
 * @param delta Value of `delta` assumed for when DM begins to freeze out.
 *              Default value is the solution to
 *                  delta * (2 + delta) / (1 + delta) = 1,
 *              i.e., delta = (sqrt(5) - 1) / 2 = 0.618033988749895.
 *              See Eqn.(13) of arXiv:1204.3622v3 for details and other used
 *              values of delta. Value of xstar is logarithmically sensitive to
 *              this number.
 * @return res Residual of the root equation.
 */
double xstar_root_eqn_gondolo_gelmini(
        const Parameters &params,
        double xstar,
        double delta = 0.0
) {
    using cosmology::sm_entropy_density;
    using cosmology::sm_entropy_density_deriv;
    using cosmology::sm_sqrt_gstar;
    using cosmology::neq;
    using autodiff::Dual;

    double deltabar = delta <= 0.0 ? 1.0 : delta * (2.0 + delta) / (1.0 + delta);
    double T = params.mx / xstar;
    double lam = sqrt(M_PI / 45.0) * params.mx * kM_PLANK * sm_sqrt_gstar(T);
    double tcs = thermal_cross_section(params, xstar);

    double s = sm_entropy_density(T);
    double ds = sm_entropy_density_deriv(T);
    auto ndn = neq(Dual<double>{T, 1.0}, params.mx, 2.0, 1);
    double yeq = ndn.val / s;
    // This is dY/dT = -x^2/m dY/xd
    double dyeq = (s * ndn.eps - ndn.val * ds) / (s * s);
    dyeq *= -params.mx / (xstar * xstar);

    return (xstar * xstar * dyeq) + (lam * deltabar * tcs * yeq * yeq);
}

/**
 * Computes to value of `xstar`: the value of dm_mass / temperature such that
 * the DM begins to freeze out.
 *
 * @param params Model parameters
 * @param chi Dark matter particle
 * @param delta Value of `delta` assumed for when DM begins to freeze out.
 *              Default value is the solution to
 *                  delta * (2 + delta) / (1 + delta) = 1,
 *              i.e., delta = (sqrt(5) - 1) / 2 = 0.618033988749895.
 *              See Eqn.(13) of arXiv:1204.3622v3 for details and other used
 *              values of delta. Value of xstar is logarithmically sensitive to
 *              this number.
 * @return xstar Value of mass / temperature at which DM begins to freeze-out.
 */
double compute_xstar_gondolo_gelmini(
        const Parameters &params,
        double delta = 0.0
) {
    using namespace boost::math::tools;
    auto f = [params, delta](double xstar) {
        return xstar_root_eqn_gondolo_gelmini(params, xstar, delta);
    };
    auto tol = [](double min, double max) { return abs(max - min) <= 1e-8; };
    std::pair<double, double> res = bisect(f, 0.01, 100.0, tol);
    return (res.second + res.first) / 2.0;
}

/**
 * Computes the value of the integral of RHS of the Boltzmann equation with
 * Yeq set to zero from x_star to x_fo.
 *
 * @param xstar Value of mass / temperature at which DM begins to freeze-out.
 *              See `compute_xstar` for more details.
 * @return
 */
double compute_alpha_gondolo_gelmini(const Parameters &params, double xstar) {
    using integrate::Quad;
    using cosmology::sm_sqrt_gstar;

    auto integrand = [params](double x) {
        return sm_sqrt_gstar(params.mx / x) * thermal_cross_section(params, x) / (x * x);
    };

    double abstol = 1e-8;
    double reltol = 1e-5;
    double error;

    double integral = Quad<double>::integrate(
            integrand,
            xstar,
            std::numeric_limits<double>::infinity(),
            abstol,
            reltol,
            500,
            &error
    );
    double pf = sqrt(M_PI / 45.0) * params.mx * kM_PLANK;
    return pf * integral;
}

/**
 * Compute the relic density using method by Gondolo + Gelmini.
 * @return rd Relic density.
 */
double relic_density_gondolo_gelmini(const Parameters &params) {
    using cosmology::neq;
    using cosmology::sm_entropy_density;

    const double xstar = compute_xstar_gondolo_gelmini(params);
    const double alpha = compute_alpha_gondolo_gelmini(params, xstar);
    const double Tstar = params.mx / xstar;
    const double ystar = neq(Tstar, params.mx, 2.0, 1) / sm_entropy_density(Tstar);
    const double Y0 = ystar / (1.0 + ystar * alpha);

    return Y0 * params.mx * kS_TODAY / kRHO_CRIT;
}

}
}
}

#endif //LANRE_KINETIC_MIXING_GONDOLO_GELMINI_HPP
