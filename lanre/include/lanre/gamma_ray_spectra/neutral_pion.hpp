//
// Created by Logan Morrison on 4/5/20.
//

#ifndef LANRE_GAMMA_RAY_SPECTRA_NEUTRAL_PION_HPP
#define LANRE_GAMMA_RAY_SPECTRA_NEUTRAL_PION_HPP

#include "lanre/constants.hpp"
#include <vector>
#include <cmath>

namespace lanre {
namespace gamma_ray_spectra {

double decay_spectrum_neutral_pion(double egam, double epi) {
    if (epi < kNEUTRAL_PION_MASS) {
        return 0.0;
    }
    double gamma = epi / kNEUTRAL_PION_MASS;
    double beta = sqrt(1.0 - 1.0 / (gamma * gamma));
    double ret_val = 0.0;

    if (epi * (1 - beta) / 2.0 <= egam <= epi * (1 + beta) / 2.0) {
        ret_val = kBR_PI0_TO_GG * 2.0 / (epi * beta);
    }
    return ret_val;
}

std::vector<double> decay_spectrum_neutral_pion(std::vector<double> egams, double epi) {
    std::vector<double> spectra(egams.size(), 0.0);
    if (epi < kNEUTRAL_PION_MASS) {
        return egams;
    }
    for (size_t i = 0; i < egams.size(); i++) {
        spectra[i] = decay_spectrum_neutral_pion(egams[i], epi);
    }
    return spectra;
}

}
}

#endif //LANRE_GAMMA_RAY_SPECTRA_NEUTRAL_PION_HPP
