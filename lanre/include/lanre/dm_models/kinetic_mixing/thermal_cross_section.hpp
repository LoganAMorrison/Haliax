//
// Created by Logan Morrison on 4/11/20.
//

#ifndef LANRE_KINETIC_MIXING_THERMAL_CROSS_SECTION_HPP
#define LANRE_KINETIC_MIXING_THERMAL_CROSS_SECTION_HPP

#include "lanre/dm_models/kinetic_mixing/parameters.hpp"
#include "lanre/dm_models/kinetic_mixing/cross_sections.hpp"
#include "lanre/special_functions/besselk.hpp"
#include "lanre/integrate/quad.hpp"
#include <cmath>

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
    using integrate::Quad;

    const double denom = 2.0 * besselk2e(x);
    const double pf = x / (denom * denom);
    auto integrand = [params, &x, &state, &channel](double z) {
        // z is equal to center of mass energy / dm mass
        const double z2 = z * z;
        const double sig = annihilation_cross_section(params, params.mx * z, state, channel);
        const double kernal = z2 * (z2 - 4.0) * besselk1e(x * z) * exp(-x * (z - 2.0));
        if (isnan(sig) || isnan(kernal)) {
            // TODO: Fix this... Cross section returns NAN for large z.
            return 0.0;//sig * kernal;
        }
        return sig * kernal;
    };
    double int_resonance = 0.0;
    double int_threshold = 0.0;
    double int_infinity;

    double abstol = 1e-8;
    double reltol = 1e-3;
    double error;
    double lb;
    int ier;

    double zmin = 2.0;
    double resonance_loc = params.mv / params.mx;
    double threshold_loc = 2.0 * params.mv / params.mx;

    // Resonance: sqrt(s) = mv = z * mx => z = mv / mx > 2 => mv > 2mx > mx
    // Threshold: sqrt(s) = 2mv = z * mx => z = 2mv / mx > 2 => mv > mx

    if (params.mx < params.mv) { // Have resonance
        if (2.0 * params.mx < params.mv) { // Have threshold
            lb = threshold_loc;
            int_resonance = Quad<double>::integrate(integrand, zmin, resonance_loc, abstol, reltol, 500, &error);
            int_threshold = Quad<double>::integrate(integrand, resonance_loc, threshold_loc, abstol, reltol, 500,
                                                    &error);
        } else { // No threshold
            lb = resonance_loc;
            int_resonance = Quad<double>::integrate(integrand, zmin, resonance_loc, abstol, reltol, 500, &error);
        }
    } else {// No resonance or threshold
        lb = 2.0;
    }

    // Perform the integral from the bound to infinity
    // bound is either 2, resonance_loc or threshold_loc
    int_infinity = Quad<double>::integrate(
            integrand,
            lb,
            std::numeric_limits<double>::infinity(),
            abstol,
            reltol,
            500,
            &error,
            &ier
    );

    return pf * (int_infinity + int_threshold + int_resonance);
}


}
}
}

#endif //LANRE_KINETIC_MIXING_THERMAL_CROSS_SECTION_HPP
