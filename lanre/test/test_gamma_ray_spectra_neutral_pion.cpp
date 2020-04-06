//
// Created by Logan Morrison on 4/5/20.
//

#include "lanre/gamma_ray_spectra/neutral_pion.hpp"
#include "lanre/gamma_ray_spectra/muon.hpp"
#include "lanre/gamma_ray_spectra/charged_pion.hpp"
#include "lanre/constants.hpp"
#include <gtest/gtest.h>
#include <vector>
#include <iostream>

using namespace lanre;
using namespace lanre::gamma_ray_spectra;

TEST(TestDecaySpectra, TestMuonSpectrum) {
    double emu = 2.0 * kMUON_MASS;
    double logemu = log10(emu);
    size_t num_gams = 100;
    std::pair<double, double> loginterval{logemu - 3, logemu};
    double step = (loginterval.second - loginterval.first) / double(num_gams - 1);

    std::vector<double> egams(num_gams, 0.0);
    for (size_t i = 0; i < egams.size(); i++) {
        egams[i] = pow(10.0, i * step + loginterval.first);
    }

    auto spectrum = decay_spectrum_muon(egams, emu);

    for (size_t i = 0; i < egams.size(); i++) {
        double egam_mev = egams[i] * 1e3;
        double dnde_mev = spectrum[i] * 1e-3;

        std::cout << "(e, dNdE) = (" << egam_mev << ", " << dnde_mev << ")" << std::endl;
    }
}
