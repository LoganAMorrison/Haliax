//
// Created by Logan Morrison on 3/20/20.
//

#include <lanre/constants.hpp>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace lanre;

PYBIND11_MODULE(constants, m) {
    // SM constants
    m.attr("g_fermi") = kG_FERMI;
    m.attr("higgs_vev") = kHIGGS_VEV;
    m.attr("alpha_em") = kALPHA_EM;
    m.attr("sin_theta_weak") = kSIN_THETA_WEAK;
    m.attr("cos_theta_weak") = kCOS_THETA_WEAK;
    // Cosmology constants
    m.attr("plank_mass") = kM_PLANK;
    m.attr("rho_crit") = kRHO_CRIT;
    m.attr("s_today") = kS_TODAY;
    m.attr("T_cmb") = kT_CMB;
    m.attr("T_bbn") = kT_BBN;
    m.attr("omega_h2_cdm") = kOMEGA_H2_CDM;
    // Masses
    m.attr("electron_mass") = kELECTRON_MASS;
    m.attr("muon_mass") = kMUON_MASS;
    m.attr("tau_mass") = kTAU_MASS;
    m.attr("up_quark_mass") = kUP_QUARK_MASS;
    m.attr("down_quark_mass") = kDOWN_QUARK_MASS;
    m.attr("strange_quark_mass") = kSTRANGE_QUARK_MASS;
    m.attr("charm_quark_mass") = kCHARM_QUARK_MASS;
    m.attr("bottom_quark_mass") = kBOTTOM_QUARK_MASS;
    m.attr("top_quark_mass") = kTOP_QUARK_MASS;
    m.attr("w_mass") = kW_BOSON_MASS;
    m.attr("z_mass") = kZ_BOSON_MASS;
    m.attr("higgs_mass") = kHIGGS_MASS;
    // Widths
    m.attr("width_w") = kW_BOSON_WIDTH;
    m.attr("width_z") = kZ_BOSON_WIDTH;
    m.attr("width_higgs") = kHIGGS_WIDTH;
}