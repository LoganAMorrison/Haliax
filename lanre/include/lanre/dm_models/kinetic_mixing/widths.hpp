//
// Created by Logan Morrison on 4/11/20.
//

#ifndef LANRE_KINETIC_MIXING_WIDTHS_HPP
#define LANRE_KINETIC_MIXING_WIDTHS_HPP

#include "lanre/dm_models/kinetic_mixing/parameters.hpp"

namespace lanre {
namespace dm_models {
namespace kinetic_mixing {

/**
 * Compute the partial width for V -> up-type quarks.
 *
 * @param params Model parameters.
 * @param mf Mass of the up-type quark.
 * @return Partial width for V -> u + ubar
 */
double width_v_to_ququ(const Parameters &params, double mf) {
    if (params.mv <= 2.0 * mf) {
        return 0.0;
    }
    return (kALPHA_EM * pow(params.eps, 2) * sqrt(-4 * pow(mf, 2) + pow(params.mv, 2)) *
            (7 * pow(mf, 2) + 17 * pow(params.mv, 2))) /
            (72. * pow(kCOS_THETA_WEAK, 2) * pow(params.mv, 2));
}

/**
 * Compute the partial width for V -> down-type quarks.
 *
 * @param params Model parameters.
 * @param mf Mass of the down-type quark.
 * @return Partial width for V -> d + dbar
 */
double width_v_to_qdqd(const Parameters &params, double mf) {
    if (params.mv <= 2.0 * mf) {
        return 0.0;
    }
    return (kALPHA_EM * pow(params.eps, 2) * sqrt(-4 * pow(mf, 2) + pow(params.mv, 2)) *
            (-17 * pow(mf, 2) + 5 * pow(params.mv, 2))) /
            (72. * pow(kCOS_THETA_WEAK, 2) * pow(params.mv, 2));
}

/**
 * Compute the partial width for V -> leptons.
 *
 * @param params Model parameters.
 * @param mf Mass of the leptons.
 * @return Partial width for V -> l + lbar
 */
double width_v_to_ll(const Parameters &params, double mf) {
    if (params.mv <= 2.0 * mf) {
        return 0.0;
    }
    return (kALPHA_EM * pow(params.eps, 2) * sqrt(-4 * pow(mf, 2) + pow(params.mv, 2)) *
            (7 * pow(mf, 2) + 5 * pow(params.mv, 2))) /
            (24. * pow(kCOS_THETA_WEAK, 2) * pow(params.mv, 2));
}

/**
 * Compute the partial width for V -> neutrinos.
 *
 * @param params Model parameters.
 * @return Partial width for V -> nu + nubar
 */
double width_v_to_nunu(const Parameters &params) {
    return (kALPHA_EM * pow(params.eps, 2) * params.mv) / (24. * pow(kCOS_THETA_WEAK, 2));
}

/**
 * Compute the partial width for V -> higgs + z-boson.
 *
 * @param params Model parameters.
 * @return Partial width for V -> H + Z.
 */
double width_v_to_hz(const Parameters &params) {
    if (params.mv <= kHIGGS_MASS + kZ_BOSON_MASS) {
        return 0.0;
    }
    return (pow(kALPHA_EM, 2) * pow(params.eps, 2) *
            sqrt(-pow(kHIGGS_MASS, 2) + pow(pow(kHIGGS_MASS, 2) + pow(params.mv, 2) -
                                                    pow(kW_BOSON_MASS, 2), 2) /
                    (4. * pow(params.mv, 2))) * (2 * pow(params.mv, 2) * pow(kW_BOSON_MASS, 2) +
            pow(-pow(kHIGGS_MASS, 2) + pow(params.mv, 2) + pow(kW_BOSON_MASS, 2), 2) / 4.) * M_PI *
            pow(kHIGGS_VEV, 2)) /
            (6. * pow(kCOS_THETA_WEAK, 4) * pow(params.mv, 4) * pow(kW_BOSON_MASS, 2) * pow(kSIN_THETA_WEAK, 2));
}

/**
 * Compute the partial width for V -> dark matter.
 *
 * @param params Model parameters.
 * @return Partial width for V -> x + xbar.
 */
double width_v_to_xx(const Parameters &params) {
    if (params.mv <= 2.0 * params.mx) {
        return 0.0;
    }
    return (pow(params.gvxx, 2) * sqrt(-4 * pow(params.mx, 2) + pow(params.mv, 2)) *
            (2 * pow(params.mx, 2) + pow(params.mv, 2))) /
            (12. * pow(params.mv, 2) * M_PI);
}


/**
 * Compute the total width or parital of the vector meidator.
 *
 * @param params Model parameters
 * @param state Final state to consider. Default is "all".
 * @return Total or partial width of the vector mediator.
 */
double vector_mediator_width(const Parameters &params, const std::string &state = "all") {
    if (state == "all") {
        return (width_v_to_ququ(params, kUP_QUARK_MASS) +
                width_v_to_ququ(params, kCHARM_QUARK_MASS) +
                width_v_to_ququ(params, kTOP_QUARK_MASS) +
                width_v_to_qdqd(params, kDOWN_QUARK_MASS) +
                width_v_to_qdqd(params, kSTRANGE_QUARK_MASS) +
                width_v_to_qdqd(params, kBOTTOM_QUARK_MASS) +
                width_v_to_ll(params, kELECTRON_MASS) +
                width_v_to_ll(params, kMUON_MASS) +
                width_v_to_ll(params, kTAU_MASS) +
                3 * width_v_to_nunu(params) +
                width_v_to_hz(params) +
                width_v_to_xx(params));
    } else if (state == "u u" || state == "u") {
        return width_v_to_ququ(params, kUP_QUARK_MASS);
    } else if (state == "c c" || state == "c") {
        return width_v_to_ququ(params, kCHARM_QUARK_MASS);
    } else if (state == "t t" || state == "t") {
        return width_v_to_ququ(params, kTOP_QUARK_MASS);
    } else if (state == "d d" || state == "d") {
        return width_v_to_qdqd(params, kDOWN_QUARK_MASS);
    } else if (state == "s s" || state == "s") {
        return width_v_to_qdqd(params, kSTRANGE_QUARK_MASS);
    } else if (state == "b b" || state == "b") {
        return width_v_to_qdqd(params, kBOTTOM_QUARK_MASS);
    } else if (state == "e e" || state == "e") {
        return width_v_to_ll(params, kELECTRON_MASS);
    } else if (state == "mu mu" || state == "mu") {
        return width_v_to_ll(params, kMUON_MASS);
    } else if (state == "tau tau" || state == "tau") {
        return width_v_to_ll(params, kTAU_MASS);
    } else if (state == "nue nue" || state == "nue" ||
            state == "num num" || state == "num" ||
            state == "nut nut" || state == "nut") {
        return width_v_to_nunu(params);
    } else if (state == "x x" || state == "x") {
        return width_v_to_xx(params);
    } else if (state == "h z" || state == "z h") {
        return width_v_to_hz(params);
    }
}

}
}
}

#endif //LANRE_KINETIC_MIXING_WIDTHS_HPP
