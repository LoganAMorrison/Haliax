//
// Created by Logan Morrison on 4/11/20.
//

#ifndef LANRE_KINETIC_MIXING_THERMAL_CROSS_SECTION_HPP
#define LANRE_KINETIC_MIXING_THERMAL_CROSS_SECTION_HPP

#include "lanre/dm_models/kinetic_mixing/parameters.hpp"
#include "lanre/dm_models/kinetic_mixing/cross_sections.hpp"
#include "lanre/special_functions/besselk.hpp"
#include "lanre/integrate/qagi.hpp"
#include "lanre/integrate/qagp.hpp"

namespace lanre {
namespace dm_models {
namespace kinetic_mixing {

double thermal_cross_section(
        const Parameters &params,
        const double x,
        const std::string &state = "all",
        const std::string &channel = "all"
) {
    using special_functions::besselk1e;
    using special_functions::besselk2e;
    using integrate::qagi;
    using integrate::qagp;

    const double denom = 2.0 * besselk2e(x);
    const double pf = x / (denom * denom);
    auto integrand = [params, &x, &state, &channel](double z) {
        // z is equal to center of mass energy / dm mass
        const double z2 = z * z;
        const double sig = annihilation_cross_section(params, params.mx * z, state, channel);
        const double kernal = z2 * (z2 - 4.0) * besselk1e(x * z) * exp(-x * (z - 2.0));
        return sig * kernal;
    };
    double int_resonance = 0.0;
    double int_threshold = 0.0;
    double int_infinity;

    double epsabs = 1e-8;
    double epsrel = 1e-5;
    double abserr;
    int neval, ier;
    double bound;
    int inf = 1;

    double zmin = 2.0;
    double resonance_loc = params.mv / params.mx;
    double threshold_loc = 2.0 * params.mv / params.mx;

    // Resoance: sqrt(s) = mv = z * mx => z = mv / mx > 2 => mv > 2mx > mx
    // Resoance: sqrt(s) = 2mv = z * mx => z = 2mv / mx > 2 => mv > mx

    std::vector<double> sing_pts(2);
    if (params.mx < params.mv) { // Have resonance
        if (2.0 * params.mx < params.mv) { // Have threshold
            bound = threshold_loc;

            sing_pts[0] = zmin;
            sing_pts[1] = resonance_loc;
            int_resonance = qagp(
                    integrand,
                    sing_pts[0],
                    sing_pts[1],
                    sing_pts.size(),
                    sing_pts.data(),
                    epsabs,
                    epsrel,
                    &abserr,
                    &neval,
                    &ier
            );

            sing_pts[0] = resonance_loc;
            sing_pts[1] = threshold_loc;
            int_threshold = qagp(
                    integrand,
                    sing_pts[0],
                    sing_pts[1],
                    sing_pts.size(),
                    sing_pts.data(),
                    epsabs,
                    epsrel,
                    &abserr,
                    &neval,
                    &ier
            );
        } else { // No threshold
            bound = params.mv / params.mx;

            sing_pts[0] = zmin;
            sing_pts[1] = resonance_loc;

            int_resonance = qagp(
                    integrand,
                    sing_pts[0],
                    sing_pts[1],
                    sing_pts.size(),
                    sing_pts.data(),
                    epsabs,
                    epsrel,
                    &abserr,
                    &neval,
                    &ier
            );
        }
    } else {// No resonance or threshold
        bound = 2.0;
    }

    // Perform the integral from the bound to infinity
    // bound is either 2, resonance_loc or threshold_loc
    int_infinity = qagi(integrand, bound, inf, epsabs, epsrel, &abserr, &neval, &ier);

    return pf * (int_infinity + int_threshold + int_resonance);
}


}
}
}

#endif //LANRE_KINETIC_MIXING_THERMAL_CROSS_SECTION_HPP
