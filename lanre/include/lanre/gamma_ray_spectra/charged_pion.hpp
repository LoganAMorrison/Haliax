//
// Created by Logan Morrison on 4/6/20.
//

#ifndef LANRE_GAMMA_RAY_SPECTRA_CHARGED_PION_HPP
#define LANRE_GAMMA_RAY_SPECTRA_CHARGED_PION_HPP

#include "lanre/constants.hpp"
#include "lanre/gamma_ray_spectra/muon.hpp"
#include "lanre/integrate/qagp.hpp"
#include "lanre/interpolate/univariate_spline.hpp"
#include <vector>


namespace lanre {
namespace gamma_ray_spectra {

/*
class ChargedPionDecaySpectrum {
private:
    std::vector<double> m_interp_mu_egams;
    std::vector<double> m_interp_mu_dndes;
    interpolate::UnivariateSpline m_mu_spline;
public:
    ChargedPionDecaySpectrum() {
        auto log_egam_interval = std::make_pair(-5, 0.0);
        size_t ngams = 500;
        double step = (log_egam_interval.second - log_egam_interval.first) / double(ngams - 1);

        double logegam;
        m_interp_mu_egams.reserve(ngams);
        m_interp_mu_dndes.reserve(ngams);
        for (int i = 0; i < ngams; i++) {
            logegam = log_egam_interval.first + i * step;
            m_interp_mu_egams.push_back(
                    log_egam_interval.first + i * step
            );
            m_interp_mu_dndes.push_back(
                    decay_spectrum_muon(pow(10.0, logegam), 0.0)
            );

            std::vector<double> weights(ngams, 1.0);
            m_mu_spline = interpolate::UnivariateSpline(
                    m_interp_mu_egams, m_interp_mu_dndes, weights,
                    log_egam_interval, 1, 0.0, 0
            );
        }
    }
};
*/

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

/**
 * Returns the integrand of the differential radiative decay spectrum for
 * the charged pion.
 * @param cl Angle of photon w.r.t. charged pion in lab frame.
 * @param egam Energy of photon in laboratory frame.
 * @param epi Energy of pion in laboratory frame.
 */
double decay_spectrum_charged_pion_integrand(double cl, double egam, double epi, const std::string &mode) {
    // Make some aliases for readability
    const double &mpi = kCHARGED_PION_MASS;
    const double &mmu = kMUON_MASS;
    const double &me = kELECTRON_MASS;

    double gamma_pi = epi / kCHARGED_PION_MASS;
    double beta_pi = sqrt(1.0 - 1.0 / (gamma_pi * gamma_pi));

    double eng_mu_pi_rf = (mpi * mpi + mmu * mmu) / (2.0 * mpi);
    double gamma_mu = eng_mu_pi_rf / kMUON_MASS;
    double beta_mu = sqrt(1.0 - 1.0 / (gamma_mu * gamma_mu));

    double eng_gam_pi_rf = egam * gamma_pi * (1.0 - beta_pi * cl);

    double jac = 1. / (2.0 * gamma_pi * abs(1.0 - beta_pi * cl));

    double eng_gam_max_mu_rf = (mmu * mmu - me * me) / (2.0 * mmu);
    double eng_gam_max_pi_rf = eng_gam_max_mu_rf * gamma_mu * (1.0 + beta_mu);

    double dnde_munu = 0.0;
    if (0. < eng_gam_pi_rf and eng_gam_pi_rf < eng_gam_max_pi_rf) {
        dnde_munu = kBR_PI_TO_MUNU * jac * decay_spectrum_muon(eng_gam_pi_rf, eng_mu_pi_rf);
        //pow(10.0, loglog_mu_spline(log10(eng_gam_pi_rf)));
    }

    double dnde_munug = kBR_PI_TO_MUNU * jac * dnde_pi_to_lnug(eng_gam_pi_rf, mmu);
    double dnde_enug = kBR_PI_TO_ENU * jac * dnde_pi_to_lnug(eng_gam_pi_rf, me);

    if (mode == "total") {
        return dnde_munu + dnde_munug + dnde_enug;
    } else if (mode == "munu") {
        return dnde_munu;
    } else if (mode == "munug") {
        return dnde_munug;
    } else if (mode == "enug") {
        return dnde_enug;
    } else {
        return 0.0;
    }
}

double decay_spectrum_charged_pion(double egam, double epi, const std::string &mode) {
    using integrate::qagp;
    if (epi < kCHARGED_PION_MASS) {
        return 0.0;
    }

    auto f = [egam, epi, mode](double cl) {
        return decay_spectrum_charged_pion_integrand(cl, egam, epi, mode);
    };

    std::vector<double> pts = {-1.0, 1.0};
    double abserr;
    int neval;
    int ier;
    return qagp(f, -1.0, 1.0, pts.size(), pts.data(), 1e-10, 1e-4, &abserr, &neval, &ier);
}

std::vector<double> decay_spectrum_charged_pion(
        const std::vector<double> &egams, double epi, const std::string &mode) {
    std::vector<double> spectrum(egams.size(), 0.0);
    if (epi < kCHARGED_PION_MASS) {
        return spectrum;
    } else {
        for (size_t i = 0; i < egams.size(); i++) {
            spectrum[i] = decay_spectrum_charged_pion(egams[i], epi, mode);
        }
        return spectrum;
    }
}


}
}

#endif //LANRE_GAMMA_RAY_SPECTRA_CHARGED_PION_HPP
