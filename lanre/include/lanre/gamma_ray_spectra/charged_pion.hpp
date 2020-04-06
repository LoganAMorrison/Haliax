//
// Created by Logan Morrison on 4/6/20.
//

#ifndef LANRE_GAMMA_RAY_SPECTRA_CHARGED_PION_HPP
#define LANRE_GAMMA_RAY_SPECTRA_CHARGED_PION_HPP

#include "lanre/constants.hpp"
#include "lanre/gamma_ray_spectra/muon.hpp"
#include "lanre/interpolate/univariate_spline.hpp"
#include <vector>

namespace lanre {
namespace gamma_ray_spectra {

/**
 * Returns dnde from pi-> l nu g.
 * @param egam
 * @param ml
 * @return
 */
static double dnde_pi_to_lnug(double egam, double ml) {
    static const double F_A_PI = 0.0119;
    static const double F_V_PI = 0.0254;
    static const double F_V_PI_SLOPE = 0.1;
    static const double fpi = 130.41 / sqrt(2.0);

    double mpi = kCHARGED_PION_MASS;
    double x = 2 * egam / mpi;
    double r = (ml / mpi) * (ml / mpi);

    if (0.0 <= x && x <= (1 - r)) {
        double F_V = F_V_PI * (1 + F_V_PI_SLOPE * (1 - x));
        double f = ((r + x - 1) * (
                mpi * mpi * x * x * x * x * (F_A_PI * F_A_PI + F_V * F_V) * (r * r - r * x + r - 2 * (x - 1) * (x - 1))
                        - 12 * sqrt(2) * fpi * mpi * r * (x - 1) * x * x * (F_A_PI * (r - 2 * x + 1) + F_V * x)
                        - 24 * fpi * fpi * r * (x - 1) * (4 * r * (x - 1) + (x - 2) * (x - 2))));
        double g = 12 * sqrt(2) * fpi * r * (x - 1) * (x - 1) * log(r / (1 - x)) * (
                mpi * x * x * (F_A_PI * (x - 2 * r) - F_V * x)
                        + sqrt(2) * fpi * (2 * r * r - 2 * r * x - x * x + 2 * x - 2));
        return kALPHA_EM * (f + g) / (24 * M_PI * mpi * fpi * fpi * (r - 1) * (r - 1)
                * (x - 1) * (x - 1) * r * x);
    } else {
        return 0.0;
    }
}

/**
 * Returns the maximum allowed gamma ray energy from a charged pion decay.
 *
 * @param epi
 * @return
 */
static double egam_max(double epi) {
    static const double mpi = kCHARGED_PION_MASS;
    static const double mmu = kMUON_MASS;
    static const double me = kELECTRON_MASS;
    static const double emu_pi_rf = (mpi * mpi + mmu * mmu) / (2.0 * mpi);
    static const double egam_max_mu_rf = (mmu * mmu - me * me) / (2.0 * mmu);

    double gamma_pi = epi / kCHARGED_PION_MASS;
    double beta_pi = sqrt(1.0 - 1.0 / (gamma_pi * gamma_pi));

    double gamma_mu = emu_pi_rf / mmu;
    double beta_mu = sqrt(1.0 - 1.0 / (gamma_mu * gamma_mu));

    return egam_max_mu_rf * gamma_pi * gamma_mu * (1.0 + beta_pi) * (1.0 + beta_mu);
}

}
}

#endif //LANRE_GAMMA_RAY_SPECTRA_CHARGED_PION_HPP
