//
// Created by Logan Morrison on 4/6/20.
//

#ifndef LANRE_GAMMA_RAY_SPECTRA_MUON_HPP
#define LANRE_GAMMA_RAY_SPECTRA_MUON_HPP

#include "lanre/constants.hpp"
#include "lanre/integrate/qagp.hpp"
#include <cmath>
#include <vector>

namespace lanre {
namespace gamma_ray_spectra {

/**
 * Form factor in differential branching fraction of radiative muon decay.
 * See p.18, eqn (54) of arXiv:hep-ph/9909265.
 * @param y 2 * (photon energy) / (muon mass)
 */
static double j_plus(double y) {
    double yconj = 1.0 - y;
    double r = (kELECTRON_MASS / kMUON_MASS);
    r *= r;
    double preFactor = kALPHA_EM * yconj / 6.0 / M_PI;
    double term1 = 3.0 * log(yconj / r) - (17.0 / 2.0);
    double term2 = -3.0 * log(yconj / r) + 7.0;
    double term3 = 2.0 * log(yconj / r) - (13.0 / 3.0);
    return preFactor * (term1 + term2 * yconj + term3 * yconj * yconj);
}

/**
 * Form factor in differential branching fraction of radiative muon decay.
 * See p.18, eqn (55) of arXiv:hep-ph/9909265.
 * @param y 2 * (photon energy) / (muon mass)
 */
static double j_minus(double y) {
    double yconj = 1.0 - y;
    double r = (kELECTRON_MASS / kMUON_MASS);
    r *= r;
    double preFactor = kALPHA_EM * yconj * yconj / 6.0 / M_PI;
    double term1 = 3.0 * log(yconj / r) - (93.0 / 12.0);
    double term2 = -4.0 * log(yconj / r) + (29.0 / 3.0);
    double term3 = 2.0 * log(yconj / r) - (55.0 / 12.0);
    return preFactor * (term1 + term2 * yconj + term3 * yconj * yconj);
}

/**
 * Differential branching fraction from: mu -> e nu nu gam.
 * See p.18, eqn (56) of arXiv:hep-ph/9909265.
 * @param y 2 * (photon energy) / (muon mass)
 */
static double dbdy(double y) {
    double result = 0.0;
    double r = (kELECTRON_MASS / kMUON_MASS);
    r *= r;
    if (0.0 <= y and y <= 1.0 - r) {
        result = (2.0 / y) * (j_plus(y) + j_minus(y));
    }
    return result;
}

/**
 * Compute integrand of dN_{\gamma}/dE_{\gamma} from mu -> e nu nu gamma
 * in laboratory frame. The integration variable is cl - the angle between
 * gamma ray and muon.
 *
 * @param cl Angle between gamma ray and muon in laboratory frame.
 * @param egam Gamma ray energy in laboratory frame.
 * @param emu Muon energy in laboratory frame.
 */
static double integrand(double cl, double egam, double emu) {
    double gamma = emu / kMUON_MASS;
    double beta = sqrt(1.0 - 1.0 / (gamma * gamma));
    double egam_mu_rf = gamma * egam * (1.0 - beta * cl);

    return (dbdy((2.0 / kMUON_MASS) * egam_mu_rf)
            / (emu * (1.0 - cl * beta)));
}

/**
 * Compute dN/dE from mu -> e nu nu gamma in the laborartory frame.
 *
 * @param egam Gamma ray energy in laboratory frame.
 * @param emu Muon energy in laboratory frame.
 * @return
 */
double decay_spectrum_muon(double egam, double emu) {
    using lanre::integrate::qagp;
    if (emu < kMUON_MASS) {
        return 0.0;
    } else if (emu == kMUON_MASS) {
        return 1.0;
    } else {
        double gamma = emu / kMUON_MASS;
        double beta = sqrt(1.0 - 1.0 / (gamma * gamma));
        double eng_gam_max = (0.5 * (kMUON_MASS - kELECTRON_MASS * kELECTRON_MASS / kMUON_MASS)
                * gamma * (1.0 + beta));
        if (0 <= egam and egam <= eng_gam_max) {
            auto f = [egam, emu](double cl) {
                return integrand(cl, egam, emu);
            };
            std::vector<double> pts = {-1.0, 1.0};
            double abserr;
            int neval;
            int ier;
            return qagp(f, -1.0, 1.0, pts.size(), pts.data(), 1e-10, 1e-4, &abserr, &neval, &ier);
        } else {
            return 0.0;
        }
    }
}

/**
 * Compute dN/dE from mu -> e nu nu gamma in the laborartory frame for a
 * vector of photon energies.
 *
 * @param egam Gamma ray energy in laboratory frame.
 * @param emu Muon energy in laboratory frame.
 * @return
 */
std::vector<double> decay_spectrum_muon(const std::vector<double> &egams, double emu) {
    std::vector<double> spectra(egams.size(), 0.0);
    if (emu >= kMUON_MASS) {
        for (size_t i = 0; i < egams.size(); i++) {
            spectra[i] = decay_spectrum_muon(egams[i], emu);
        }
    }
    return spectra;
}

}
}

#endif //LANRE_GAMMA_RAY_SPECTRA_MUON_HPP
