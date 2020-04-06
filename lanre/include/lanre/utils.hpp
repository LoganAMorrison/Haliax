//
// Created by Logan Morrison on 3/15/20.
//

#ifndef LANRE_UTILS_HPP
#define LANRE_UTILS_HPP

#include "lanre/constants.hpp"
#include <string>

namespace lanre {
double string_to_fermion_mass(const std::string &f) {
    if (f == "e") {
        return kELECTRON_MASS;
    } else if (f == "mu") {
        return kMUON_MASS;
    } else if (f == "tau") {
        return kTAU_MASS;
    } else if (f == "u") {
        return kUP_QUARK_MASS;
    } else if (f == "c") {
        return kCHARM_QUARK_MASS;
    } else if (f == "t") {
        return kTOP_QUARK_MASS;
    } else if (f == "d") {
        return kDOWN_QUARK_MASS;
    } else if (f == "s") {
        return kSTRANGE_QUARK_MASS;
    } else if (f == "b") {
        return kBOTTOM_QUARK_MASS;
    } else {
        return 0.0;
    }
}

double string_to_num_colors(const std::string &f) {
    if (f == "u" || f == "c" || f == "t" ||
            f == "d" || f == "s" || f == "b") {
        return 3.0;
    } else {
        return 0.0;
    }
}
}

#endif //LANRE_UTILS_HPP
