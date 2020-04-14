//
// Created by Logan Morrison on 4/12/20.
//

#ifndef LANRE_KINETIC_MIXING_RELIC_DENSITY_BENDER_HPP
#define LANRE_KINETIC_MIXING_RELIC_DENSITY_BENDER_HPP

#include "lanre/dm_models/kinetic_mixing/parameters.hpp"
#include <cmath>

namespace lanre {
namespace dm_models {
namespace kinetic_mixing {

double compute_xf_bender(const Parameters &params) {
    double A = sqrt(M_PI / 2.0) * 45.0 / (4.0 * pow(M_PI, 4));
    double lam = sqrt(M_PI / 45.0) * params.mx * kM_PLANK;
}

double relic_density_bender(const Parameters &params) {
    return 0.0;
}

}
}
}

#endif //LANRE_KINETIC_MIXING_RELIC_DENSITY_BENDER_HPP
