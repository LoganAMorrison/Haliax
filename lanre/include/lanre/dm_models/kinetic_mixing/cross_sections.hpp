//
// Created by Logan Morrison on 4/11/20.
//

#ifndef LANRE_KINETIC_MIXING_CROSS_SECTIONS_HPP
#define LANRE_KINETIC_MIXING_CROSS_SECTIONS_HPP

#include "lanre/dm_models/kinetic_mixing/parameters.hpp"
#include "lanre/constants.hpp"
#include <string>

namespace lanre {
namespace dm_models {
namespace kinetic_mixing {

/**
 * Compute the annihilation cross-section for dark matter to up-type quarks.
 *
 * @param params Model parameters.
 * @param Q Center of mass energy.
 * @param mf Mass of the up-type quark.
 * @param channel String specifying channel: "s" or "tu".
 * @return Cross section for x + xbar -> u + ubar
 */
double sigma_xx_to_ququ(
        const Parameters &params,
        double Q,
        double mf,
        const std::string &channel
) {
    if ((channel != "s" && channel != "all") || Q <= 2.0 * mf || Q <= 2.0 * params.mx) {
        return 0.0;
    }
    double temp1 = pow(params.mx, 2);
    double temp2 = pow(Q, 2);
    double temp3 = pow(mf, 2);
    double temp4 = pow(params.mv, 2);
    return ((kALPHA_EM * pow(params.eps, 2) * pow(params.gvxx, 2) * (2 * temp1 +
            temp2) * sqrt(temp2 - 4 * temp3) * (17 * temp2 +
            7 * temp3)) / (72. * pow(kCOS_THETA_WEAK, 2) *
            pow(Q, 2) * sqrt(-4 * temp1 + temp2) *
            (pow(-temp2 + temp4, 2) + temp4 * pow(params.widthv, 2))));
}

/**
 * Compute the annihilation cross-section for dark matter to down-type quarks.
 *
 * @param params Model parameters.
 * @param Q Center of mass energy.
 * @param mf Mass of the down-type quark.
 * @param channel String specifying channel: "s" or "tu".
 * @return Cross section for x + xbar -> d + dbar
 */
double sigma_xx_to_qdqd(
        const Parameters &params,
        double Q,
        double mf,
        const std::string &channel
) {
    if ((channel != "s" && channel != "all") || Q <= 2.0 * mf || Q <= 2.0 * params.mx) {
        return 0.0;
    }
    double temp1 = pow(Q, 2);
    double temp2 = pow(params.mx, 2);
    double temp3 = pow(mf, 2);
    double temp4 = pow(params.mv, 2);
    return ((kALPHA_EM * pow(params.eps, 2) * pow(params.gvxx, 2) * (temp1 +
            2 * temp2) * (5 * temp1 - 17 * temp3) * sqrt(temp1 - 4 * temp3)) /
            (72. * pow(kCOS_THETA_WEAK, 2) * pow(Q, 2) * sqrt(temp1 - 4 * temp2) *
                    (pow(-temp1 + temp4, 2) + temp4 * pow(params.widthv, 2))));
}

/**
 * Compute the annihilation cross-section for dark matter to leptons.
 *
 * @param params Model parameters.
 * @param Q Center of mass energy.
 * @param mf Mass of the lepton.
 * @param channel String specifying channel: "s" or "tu".
 * @return Cross section for x + xbar -> l + lbar
 */
double sigma_xx_to_ll(
        const Parameters &params,
        double Q,
        double mf,
        const std::string &channel
) {
    if ((channel != "s" && channel != "all") || Q <= 2.0 * mf || Q <= 2.0 * params.mx) {
        return 0.0;
    }
    double temp1 = pow(params.mx, 2);
    double temp2 = pow(Q, 2);
    double temp3 = pow(mf, 2);
    double temp4 = pow(params.mv, 2);
    return (kALPHA_EM * pow(params.eps, 2) * pow(params.gvxx, 2) * (2 * temp1 +
            temp2) * sqrt(temp2 - 4 * temp3) * (5 * temp2 + 7 * temp3)) /
            (24. * pow(kCOS_THETA_WEAK, 2) * pow(Q, 2) * sqrt(-4 * temp1 + temp2) *
                    (pow(-temp2 + temp4, 2) + temp4 * pow(params.widthv, 2)));
}

/**
 * Compute the annihilation cross-section for dark matter to neutrinos.
 *
 * @param params Model parameters.
 * @param Q Center of mass energy.
 * @param channel String specifying channel: "s" or "tu".
 * @return Cross section for x + xbar -> nu + nubar
 */
double sigma_xx_to_nunu(
        const Parameters &params,
        double Q,
        const std::string &channel
) {
    if (channel != "s" && channel != "all") {
        return 0.0;
    }
    double temp1 = pow(Q, 2);
    double temp2 = pow(params.mx, 2);
    double temp3 = pow(params.mv, 2);
    return (kALPHA_EM * pow(params.eps, 2) * pow(params.gvxx, 2) * sqrt(temp1) * (temp1 +
            2 * temp2)) / (24. * pow(kCOS_THETA_WEAK, 2) * sqrt(temp1 - 4 * temp2) *
            (pow(-temp1 + temp3, 2) + temp3 * pow(params.widthv, 2)));
}

/**
 * Compute the annihilation cross-section for dark matter to a higgs and z-boson.
 *
 * @param params Model parameters.
 * @param Q Center of mass energy.
 * @param channel String specifying channel: "s" or "tu".
 * @return Cross section for x + xbar -> H + Z.
 */
double sigma_xx_to_hz(
        const Parameters &params,
        double Q,
        const std::string &channel
) {
    if ((channel != "s" && channel != "all") || Q <= kHIGGS_MASS + kW_BOSON_MASS || Q <= 2.0 * params.mx) {
        return 0.0;
    }
    double temp1 = -Q;
    double temp2 = -kW_BOSON_MASS;
    double temp3 = pow(params.mx, 2);
    double temp4 = pow(Q, 2);
    double temp5 = pow(kW_BOSON_MASS, 2);
    double temp6 = pow(params.mv, 2);
    return (pow(kALPHA_EM, 2) * pow(params.eps, 2) * pow(params.gvxx, 2) * M_PI * sqrt(((
            kHIGGS_MASS + kW_BOSON_MASS + Q) * (kHIGGS_MASS + kW_BOSON_MASS +
            temp1) * (kHIGGS_MASS + Q + temp2) * (kHIGGS_MASS + temp1 +
            temp2)) / (-4 * temp3 + temp4)) * (2 * temp3 + temp4) * (pow(kHIGGS_MASS, 4) +
            pow(kW_BOSON_MASS, 4) + pow(Q, 4) + 10 * temp4 * temp5 -
            2 * pow(kHIGGS_MASS, 2) * (temp4 + temp5)) * pow(kHIGGS_VEV, 2)) /
            (48. * pow(kCOS_THETA_WEAK, 4) * pow(kW_BOSON_MASS, 2) * pow(Q, 5) * pow
                    (kSIN_THETA_WEAK, 2) * (pow(-temp4 + temp6, 2) +
                    temp6 * pow(params.widthv, 2)));
}

/**
 * Compute the annihilation cross-section for dark matter to a vector mediators.
 *
 * @param params Model parameters.
 * @param Q Center of mass energy.
 * @param channel String specifying channel: "s" or "tu".
 * @return Cross section for x + xbar -> v + v.
 */
double sigma_xx_to_vv(
        const Parameters &params,
        double Q,
        const std::string &channel
) {
    if (channel == "s" || Q <= 2.0 * params.mx || Q <= 2.0 * params.mv) {
        return 0.0;
    }
    double temp1 = pow(params.mx, 2);
    double temp2 = -4 * temp1;
    double temp3 = pow(Q, 2);
    double temp4 = temp2 + temp3;
    double temp5 = pow(params.mv, 4);
    double temp6 = pow(params.mv, 2);
    double temp7 = -4 * temp6;
    double temp8 = temp3 + temp7;
    double temp9 = 2 * temp6;
    double temp10 = -temp3;
    double temp11 = temp10 + temp9;
    double temp12 = -2 * temp6;
    double temp13 = sqrt(temp4);
    double temp14 = sqrt(temp8);
    double temp15 = temp13 * temp14;
    double temp16 = 2 * temp1;
    double temp17 = temp12 + temp15 + temp3;
    return (pow(params.gvxx, 4) * ((-48 * temp13 * temp14 * (4 * pow(params.mx, 4) +
            temp1 * temp3 + 2 * temp5)) / (temp5 + temp1 * temp8) + (48 * (temp1 * (temp12 +
            temp2 + temp3) * log(2) + temp3 * temp6 * log(4) + (2 * temp1 * (temp16 +
            temp6) - temp3 * (temp1 + temp9)) * log(-(pow(temp17, 2) / (-pow(Q, 4) +
            temp13 * temp14 * temp3 - 2 * temp5 + (-2 * temp13 * temp14 + 4 * temp3) * temp6 +
            2 * temp1 * temp8))) + temp11 * (temp12 + temp16 +
            temp3) * log(-(temp17 / (temp10 + temp15 +
            temp9))))) / temp11)) / (384. * M_PI * pow(Q, 2) * temp4);
}

/**
 * Compute the total or partial annihilation cross section of the dark matter.
 *
 * @param params Model parameters.
 * @param Q Center of mass energy.
 * @param state Final state to consider. Default is "all".
 * @param channel Chanel to consider, i.e. "s" or "tu".
 * @return Total or partial annihilation cross section.
 */
double annihilation_cross_section(
        const Parameters &params,
        double Q,
        const std::string &state,
        const std::string &channel
) {
    if (state == "all") {
        return (sigma_xx_to_ququ(params, Q, kUP_QUARK_MASS, channel) +
                sigma_xx_to_ququ(params, Q, kCHARM_QUARK_MASS, channel) +
                sigma_xx_to_ququ(params, Q, kTOP_QUARK_MASS, channel) +
                sigma_xx_to_qdqd(params, Q, kDOWN_QUARK_MASS, channel) +
                sigma_xx_to_qdqd(params, Q, kSTRANGE_QUARK_MASS, channel) +
                sigma_xx_to_qdqd(params, Q, kBOTTOM_QUARK_MASS, channel) +
                sigma_xx_to_ll(params, Q, kELECTRON_MASS, channel) +
                sigma_xx_to_ll(params, Q, kMUON_MASS, channel) +
                sigma_xx_to_ll(params, Q, kTAU_MASS, channel) +
                3 * sigma_xx_to_nunu(params, Q, channel) +
                sigma_xx_to_hz(params, Q, channel) +
                sigma_xx_to_vv(params, Q, channel));
    } else if (state == "u u" || state == "u") {
        return sigma_xx_to_ququ(params, Q, kUP_QUARK_MASS, channel);
    } else if (state == "c c" || state == "c") {
        return sigma_xx_to_ququ(params, Q, kCHARM_QUARK_MASS, channel);
    } else if (state == "t t" || state == "t") {
        return sigma_xx_to_ququ(params, Q, kTOP_QUARK_MASS, channel);
    } else if (state == "d d" || state == "d") {
        return sigma_xx_to_qdqd(params, Q, kDOWN_QUARK_MASS, channel);
    } else if (state == "s s" || state == "s") {
        return sigma_xx_to_qdqd(params, Q, kSTRANGE_QUARK_MASS, channel);
    } else if (state == "b b" || state == "b") {
        return sigma_xx_to_qdqd(params, Q, kBOTTOM_QUARK_MASS, channel);
    } else if (state == "e e" || state == "e") {
        return sigma_xx_to_ll(params, Q, kELECTRON_MASS, channel);
    } else if (state == "mu mu" || state == "mu") {
        return sigma_xx_to_ll(params, Q, kMUON_MASS, channel);
    } else if (state == "tau tau" || state == "tau") {
        return sigma_xx_to_ll(params, Q, kTAU_MASS, channel);
    } else if (state == "nue nue" || state == "nue" ||
            state == "num num" || state == "num" ||
            state == "nut nut" || state == "nut") {
        return sigma_xx_to_nunu(params, Q, channel);
    } else if (state == "h z" || state == "z h") {
        return sigma_xx_to_hz(params, Q, channel);
    } else if (state == "v v" || state == "v") {
        return sigma_xx_to_vv(params, Q, channel);
    }
    return 0.0;
}

}
}
}

#endif //LANRE_KINETIC_MIXING_CROSS_SECTIONS_HPP
