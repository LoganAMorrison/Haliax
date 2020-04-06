//
// Created by Logan Morrison on 3/15/20.
//

#ifndef LANRE_CONSTANTS_HPP
#define LANRE_CONSTANTS_HPP

namespace lanre {
// Various physical constants
static constexpr double kG_FERMI = 1.1663787e-5;
static constexpr double kHIGGS_VEV = 246.21965;
static constexpr double kALPHA_EM = 1.0 / 137.0; //at p^2 = 0
static constexpr double kSIN_THETA_WEAK = 0.480853;
static constexpr double kSIN_THETA_WEAK_SQRD = 0.23122;
static constexpr double kCOS_THETA_WEAK = 0.876801;
static constexpr double kM_PLANK = 1.220910e19;
static constexpr double kRHO_CRIT = 1.05375e-5;
static constexpr double kS_TODAY = 2891.2;
static constexpr double kT_CMB = 2.56215e-10;
static constexpr double kT_BBN = 0.0001; // 0.1 MeV in GeV
static constexpr double kOMEGA_H2_CDM = 0.1198;

//Masses
static constexpr double kELECTRON_MASS = 0.5109989461e-3;
static constexpr double kMUON_MASS = 105.6583745e-3;
static constexpr double kTAU_MASS = 1776.86e-3;
static constexpr double kUP_QUARK_MASS = 2.16e-3;
static constexpr double kDOWN_QUARK_MASS = 4.67e-3;
static constexpr double kSTRANGE_QUARK_MASS = 93e-3;
static constexpr double kCHARM_QUARK_MASS = 1.27;
static constexpr double kBOTTOM_QUARK_MASS = 4.18;
static constexpr double kTOP_QUARK_MASS = 172.9;
static constexpr double kW_BOSON_MASS = 80.379;
static constexpr double kZ_BOSON_MASS = 91.1876;
static constexpr double kHIGGS_MASS = 125.10;
static constexpr double kNEUTRAL_PION_MASS = 134.9766e-3;
static constexpr double kCHARGED_PION_MASS = 139.57018e-3;
static constexpr double kNEUTRAL_KAON_MASS = 497.61e-3;
static constexpr double kLONG_KAON_MASS = 497.614e-3;
static constexpr double kCHARGED_KAON_MASS = 493.68e-3;

//Boson widths
static constexpr double kW_BOSON_WIDTH = 2.085;
static constexpr double kZ_BOSON_WIDTH = 2.4952;
static constexpr double kHIGGS_WIDTH = 4.07e-3;

// Branching ratios
static constexpr double kBR_PI0_TO_GG = 0.9882;         // Pi0   -> g   + g
static constexpr double kBR_PI_TO_MUNU = 0.9998;        // pi    -> mu  + nu
static constexpr double kBR_PI_TO_ENU = 0.000123;       // pi    -> e  + nu

static constexpr double kBR_KS_TO_PIPI = 0.6920;        // ks    -> pi  + pi
static constexpr double kBR_KS_TO_PI0PI0 = 0.3069;      // ks    -> pi0 + pi0

static constexpr double kBR_KL_TO_PIENU = 0.4055;       // kl    -> pi  + e   + nu
static constexpr double kBR_KL_TO_PIMUNU = 0.2704;      // kl    -> pi  + mu  + nu
static constexpr double kBR_KL_TO_3PI0 = 0.1952;        // kl    -> pi0 + pi0  + pi0
static constexpr double kBR_KL_TO_2PIPI0 = 0.1254;      // kl    -> pi  + pi  + pi0

static constexpr double kBR_K_TO_MUNU = 0.6356;         // k     -> mu  + nu
static constexpr double kBR_K_TO_PIPI0 = 0.2067;        // k     -> pi  + pi0
static constexpr double kBR_K_TO_3PI = 0.05583;         // k     -> pi  + pi  + pi
static constexpr double kBR_K_TO_PI0ENU = 0.0507;       // k     -> pi0 + e   + nu
static constexpr double kBR_K_TO_PI0MUNU = 0.03352;     // k     -> pi0 + mu  + nu
static constexpr double kBR_K_TO_PI2PI0 = 0.01760;      // k     -> pi  + pi0 + pi0

static constexpr double kBR_ETA_TO_GG = 0.3941;         // eta   -> g   + g
static constexpr double kBR_ETA_TO_3PI0 = 0.3268;       // eta   -> pi0 + pi0 + pi0
static constexpr double kBR_ETA_TO_2PIPI0 = 0.2292;     // eta   -> pi  + pi  + pi0
static constexpr double kBR_ETA_TO_2PIG = 0.0422;       // eta   -> pi  + pi  + g
static constexpr double kBR_ETAP_TO_2PIETA = 0.429;     // eta'  -> pi  + pi  + eta
static constexpr double kBR_ETAP_TO_RHOG = 0.291;       // eta'  -> rho + g
static constexpr double kBR_ETAP_2PI0ETA = 0.222;       // eta'  -> pi0 + pi0 + eta
static constexpr double kBR_ETAP_TO_OMEGAG = 0.0275;    // eta'  -> omega + g
static constexpr double kBR_ETAP_TO_GG = 0.0220;        // eta'  -> g   + g
static constexpr double kBR_ETAP_TO_3PI0 = 0.0214;      // eta'  -> pi0 + pi0 + pi-
static constexpr double kBR_ETAP_TO_MUMUG = 0.0108;     // eta'  -> mu  + mu  + g

static constexpr double kBR_OMEGA_TO_2PIPI0 = 0.892;    // omega -> pi + pi   + pi0
static constexpr double kBR_OMEGA_TO_PI0G = 0.0828;     // omega -> pi0 + g
static constexpr double kBR_OMEGA_TO_2PI = 0.0153;      // omega -> pi + pi
}

#endif //LANRE_CONSTANTS_HPP
