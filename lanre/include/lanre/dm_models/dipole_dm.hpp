//
// Created by Logan Morrison on 3/29/20.
//

#ifndef LANRE_DM_MODELS_KINETIC_RECOUPLING_HPP
#define LANRE_DM_MODELS_KINETIC_RECOUPLING_HPP

#include "lanre/constants.hpp"
#include "lanre/special_functions/besselk.hpp"
#include <string>
#include <cmath>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/special_functions/pow.hpp>

namespace lanre {
namespace dm_models {

class DipoleDM {
private:
    double m_m1; // Mass of chi1: the lightest DM fermion
    double m_m2; // Mass of chi2: the heavier DM fermion
    double m_ce; // Wilson coefficient for the electric dipole operator
    double m_cm; // Wilson coefficient for the magnetic dipole operator
    double m_lam; // Cut-off scale of the theory
    double m_width2; // Width of chi2

    double compute_width2() const;

    double sigma_x1_x1_to_g_g(double) const;

    double sigma_x1_x1_to_x2_x2(double) const;

    double sigma_x1_x2_to_w_w(double) const;

    double sigma_x1_x2_to_f_f(double, double) const;

    double weff(double) const;

    double thermal_cross_section_integrand(double, double) const;

public:
    DipoleDM(double m1, double m2, double ce, double cm, double lam)
            : m_m1(m1), m_m2(m2), m_ce(ce), m_cm(cm), m_lam(lam), m_width2(0.0) {
        m_width2 = compute_width2();
    }

    double get_m1() const { return m_m1; }

    double get_m2() const { return m_m2; }

    double get_ce() const { return m_ce; }

    double get_cm() const { return m_cm; }

    double get_lam() const { return m_lam; }

    double get_width2() const { return m_width2; }

    void set_m1(double m1) { m_m1 = m1; }

    void set_m2(double m2) { m_m2 = m2; }

    void set_ce(double ce) { m_ce = ce; }

    void set_cm(double cm) { m_cm = cm; }

    void set_lam(double lam) { m_lam = lam; }

    double gamma_integrand(double, double) const;

    double gamma(double) const;

    double annihilation_cross_section(double, int i, int j) const;

    double thermal_cross_section(double x) const;

};


/**
 * Compute the decay width of the heavier DM particle into the lighter and a
 * photon.
 * @return width
 */
double DipoleDM::compute_width2() const {
    return ((m_ce * m_ce + m_cm * m_cm) * pow(-m_m1 * m_m1 + m_m2 * m_m2, 3)) /
            (8. * m_lam * m_lam * m_m2 * m_m2 * m_m2 * M_PI);
}


/**
 * Compute the cross section for chi1 + chi1 -> photons.
 * @param Q Center of mass energy.
 * @return sigma Sigma(chi1 + chi1 -> photon + photon)
 */
double DipoleDM::sigma_x1_x1_to_g_g(double Q) const {
    if (Q <= 2 * m_m1) {
        return 0.0;
    }
    return (pow(pow(m_ce, 2) + pow(m_cm, 2), 2) * (-(((-4 * Q * sqrt(-4 * pow(m_m1, 2) + pow(Q, 2)) *
            (-24 * pow(pow(m_m1, 2) - pow(m_m2, 2), 4) +
                    2 * (pow(m_m1, 2) - 9 * pow(m_m2, 2)) * pow(pow(m_m1, 2) - pow(m_m2, 2), 2) * pow(Q, 2) +
                    (pow(m_m1, 4) + pow(m_m2, 4)) * pow(Q, 4) + pow(m_m2, 2) * pow(Q, 6))) /
            (pow(pow(m_m1, 2) - pow(m_m2, 2), 2) + pow(m_m2, 2) * pow(Q, 2)) +
            24 * pow(pow(m_m1, 2) - pow(m_m2, 2), 2) * (4 * pow(m_m1, 2) - 4 * pow(m_m2, 2) - pow(Q, 2)) *
                    log((2 * pow(m_m1, 2) - 2 * pow(m_m2, 2) - Q * (Q + sqrt(-4 * pow(m_m1, 2) + pow(Q, 2)))) /
                                (2 * pow(m_m1, 2) - 2 * pow(m_m2, 2) +
                                        Q * (-Q + sqrt(-4 * pow(m_m1, 2) + pow(Q, 2)))))) / pow(Q, 2)) +
            (24 * (pow(m_m1, 2) + pow(m_m2, 2)) *
                    (Q * sqrt(-4 * pow(m_m1, 2) + pow(Q, 2)) * (-2 * pow(m_m1, 2) + 2 * pow(m_m2, 2) + pow(Q, 2)) -
                            (pow(pow(m_m1, 2) - pow(m_m2, 2), 2) + pow(m_m2, 2) * pow(Q, 2)) *
                                    log(pow(-2 * pow(m_m1, 2) + 2 * pow(m_m2, 2) +
                                                    Q * (Q + sqrt(-4 * pow(m_m1, 2) + pow(Q, 2))), 2) /
                                                pow(2 * pow(m_m1, 2) - 2 * pow(m_m2, 2) +
                                                            Q * (-Q + sqrt(-4 * pow(m_m1, 2) + pow(Q, 2))), 2)))) /
                    (-2 * pow(m_m1, 2) + 2 * pow(m_m2, 2) + pow(Q, 2)))) /
            (192. * pow(m_lam, 4) * M_PI * (-4 * pow(m_m1, 2) + pow(Q, 2)));
}

/**
 * Compute the cross section for chi1 + chi1 -> chi2 + chi2.
 * @param Q Center of mass energy.
 * @return sigma Sigma(chi1 + chi1 -> chi2 + chi2)
 */
double DipoleDM::sigma_x1_x1_to_x2_x2(double Q) const {
    if (Q <= 2 * m_m1 || Q <= 2 * m_m2) {
        return 0.0;
    }
    return (16 * sqrt(-4 * pow(m_m1, 2) + pow(Q, 2)) * sqrt(-4 * pow(m_m2, 2) + pow(Q, 2)) * (4 *
            (3 * pow(pow(m_ce, 2) + pow(m_cm, 2), 2) * pow(m_m1, 4) +
                    6 * (-pow(m_ce, 4) + pow(m_cm, 4)) * pow(m_m1, 3) * m_m2 +
                    2 * (5 * pow(m_ce, 4) - 8 * pow(m_ce, 2) * pow(m_cm, 2) + 5 * pow(m_cm, 4)) * pow(m_m1, 2) *
                            pow(m_m2, 2) + 6 * (-pow(m_ce, 4) + pow(m_cm, 4)) * m_m1 * pow(m_m2, 3) +
                    3 * pow(pow(m_ce, 2) + pow(m_cm, 2), 2) * pow(m_m2, 4)) - (pow(m_ce, 2) + pow(m_cm, 2)) *
            (13 * (pow(m_ce, 2) + pow(m_cm, 2)) * pow(m_m1, 2) - 6 * (m_ce - m_cm) * (m_ce + m_cm) * m_m1 * m_m2 +
                    13 * (pow(m_ce, 2) + pow(m_cm, 2)) * pow(m_m2, 2)) * pow(Q, 2) +
            7 * pow(pow(m_ce, 2) + pow(m_cm, 2), 2) * pow(Q, 4)) -
            3 * pow(pow(m_ce, 2) + pow(m_cm, 2), 2) * (-2 * (pow(m_m1, 2) + pow(m_m2, 2)) + pow(Q, 2)) *
                    pow(2 * pow(m_m1, 2) + 2 * pow(m_m2, 2) - pow(Q, 2) +
                                sqrt(-4 * pow(m_m1, 2) + pow(Q, 2)) * sqrt(-4 * pow(m_m2, 2) + pow(Q, 2)), 2) -
            pow(pow(m_ce, 2) + pow(m_cm, 2), 2) * pow(2 * pow(m_m1, 2) + 2 * pow(m_m2, 2) - pow(Q, 2) +
                                                              sqrt(-4 * pow(m_m1, 2) + pow(Q, 2)) *
                                                                      sqrt(-4 * pow(m_m2, 2) + pow(Q, 2)), 3) +
            3 * pow(pow(m_ce, 2) + pow(m_cm, 2), 2) * (-2 * (pow(m_m1, 2) + pow(m_m2, 2)) + pow(Q, 2)) *
                    pow(-2 * pow(m_m1, 2) - 2 * pow(m_m2, 2) + pow(Q, 2) +
                                sqrt(-4 * pow(m_m1, 2) + pow(Q, 2)) * sqrt(-4 * pow(m_m2, 2) + pow(Q, 2)), 2) -
            pow(pow(m_ce, 2) + pow(m_cm, 2), 2) * pow(-2 * pow(m_m1, 2) - 2 * pow(m_m2, 2) + pow(Q, 2) +
                                                              sqrt(-4 * pow(m_m1, 2) + pow(Q, 2)) *
                                                                      sqrt(-4 * pow(m_m2, 2) + pow(Q, 2)), 3) + 12 *
            (-2 * pow(m_m1, 2) - 2 * pow(m_m2, 2) + pow(Q, 2) +
                    sqrt(-4 * pow(m_m1, 2) + pow(Q, 2)) * sqrt(-4 * pow(m_m2, 2) + pow(Q, 2))) * (pow(m_ce, 4) *
            (pow(m_m1 - m_m2, 2) * (3 * pow(m_m1, 2) - 2 * m_m1 * m_m2 + 3 * pow(m_m2, 2)) +
                    (-3 * pow(m_m1, 2) + 2 * m_m1 * m_m2 - 3 * pow(m_m2, 2)) * pow(Q, 2) + 2 * pow(Q, 4)) +
            2 * pow(m_ce, 2) * pow(m_cm, 2) * (3 * pow(m_m1, 4) - 14 * pow(m_m1, 2) * pow(m_m2, 2) + 3 * pow(m_m2, 4) -
                    3 * (pow(m_m1, 2) + pow(m_m2, 2)) * pow(Q, 2) + 2 * pow(Q, 4)) + pow(m_cm, 4) *
            (pow(m_m1 + m_m2, 2) * (3 * pow(m_m1, 2) + 2 * m_m1 * m_m2 + 3 * pow(m_m2, 2)) -
                    (3 * pow(m_m1, 2) + 2 * m_m1 * m_m2 + 3 * pow(m_m2, 2)) * pow(Q, 2) + 2 * pow(Q, 4))) + 24 *
            (pow(m_ce, 4) * (pow(m_m1 - m_m2, 2) * (3 * pow(m_m1, 2) - 2 * m_m1 * m_m2 + 3 * pow(m_m2, 2)) +
                    (-3 * pow(m_m1, 2) + 2 * m_m1 * m_m2 - 3 * pow(m_m2, 2)) * pow(Q, 2) + 2 * pow(Q, 4)) +
                    2 * pow(m_ce, 2) * pow(m_cm, 2) *
                            (3 * pow(m_m1, 4) - 14 * pow(m_m1, 2) * pow(m_m2, 2) + 3 * pow(m_m2, 4) -
                                    3 * (pow(m_m1, 2) + pow(m_m2, 2)) * pow(Q, 2) + 2 * pow(Q, 4)) + pow(m_cm, 4) *
                    (pow(m_m1 + m_m2, 2) * (3 * pow(m_m1, 2) + 2 * m_m1 * m_m2 + 3 * pow(m_m2, 2)) -
                            (3 * pow(m_m1, 2) + 2 * m_m1 * m_m2 + 3 * pow(m_m2, 2)) * pow(Q, 2) + 2 * pow(Q, 4))) *
            (pow(m_m1, 2) + pow(m_m2, 2) +
                    (-pow(Q, 2) + sqrt(-4 * pow(m_m1, 2) + pow(Q, 2)) * sqrt(-4 * pow(m_m2, 2) + pow(Q, 2))) / 2.) +
            (48 * (pow(m_ce, 2) + pow(m_cm, 2)) * pow(m_m1 - m_m2, 2) * pow(m_m1 + m_m2, 2) *
                    (2 * (pow(m_m1, 2) + pow(m_m2, 2)) *
                            (pow(m_ce, 2) * pow(m_m1 - m_m2, 2) + pow(m_cm, 2) * pow(m_m1 + m_m2, 2)) +
                            (pow(m_cm, 2) * pow(m_m1 - m_m2, 2) + pow(m_ce, 2) * pow(m_m1 + m_m2, 2)) * pow(Q, 2)) *
                    log(pow(m_m1, 2) + pow(m_m2, 2) - pow(Q, 2) / 2. -
                                (sqrt(-4 * pow(m_m1, 2) + pow(Q, 2)) * sqrt(-4 * pow(m_m2, 2) + pow(Q, 2))) / 2.)) /
                    (2 * (pow(m_m1, 2) + pow(m_m2, 2)) - pow(Q, 2)) -
            96 * (pow(m_ce, 2) + pow(m_cm, 2)) * pow(m_m1 - m_m2, 2) * pow(m_m1 + m_m2, 2) *
                    (pow(m_cm, 2) * (m_m1 + m_m2 - Q) * (m_m1 + m_m2 + Q) +
                            pow(m_ce, 2) * (pow(m_m1 - m_m2, 2) - pow(Q, 2))) *
                    log((-2 * pow(m_m1, 2) - 2 * pow(m_m2, 2) + pow(Q, 2) -
                            sqrt(-4 * pow(m_m1, 2) + pow(Q, 2)) * sqrt(-4 * pow(m_m2, 2) + pow(Q, 2))) /
                                (-2 * pow(m_m1, 2) - 2 * pow(m_m2, 2) + pow(Q, 2) +
                                        sqrt(-4 * pow(m_m1, 2) + pow(Q, 2)) * sqrt(-4 * pow(m_m2, 2) + pow(Q, 2)))) -
            (48 * (pow(m_ce, 2) + pow(m_cm, 2)) * pow(m_m1 - m_m2, 2) * pow(m_m1 + m_m2, 2) *
                    (2 * (pow(m_m1, 2) + pow(m_m2, 2)) *
                            (pow(m_ce, 2) * pow(m_m1 - m_m2, 2) + pow(m_cm, 2) * pow(m_m1 + m_m2, 2)) +
                            (pow(m_cm, 2) * pow(m_m1 - m_m2, 2) + pow(m_ce, 2) * pow(m_m1 + m_m2, 2)) * pow(Q, 2)) *
                    log(-pow(2 * pow(m_m1, 2) + 2 * pow(m_m2, 2) - pow(Q, 2) +
                                     sqrt(-4 * pow(m_m1, 2) + pow(Q, 2)) * sqrt(-4 * pow(m_m2, 2) + pow(Q, 2)), 2) /
                                (2. * (-2 * pow(m_m1, 2) - 2 * pow(m_m2, 2) + pow(Q, 2) +
                                        sqrt(-4 * pow(m_m1, 2) + pow(Q, 2)) * sqrt(-4 * pow(m_m2, 2) + pow(Q, 2)))))) /
                    (2 * (pow(m_m1, 2) + pow(m_m2, 2)) - pow(Q, 2)) +
            96 * (pow(m_ce, 2) + pow(m_cm, 2)) * pow(m_m1 - m_m2, 2) * pow(m_m1 + m_m2, 2) *
                    (pow(m_cm, 2) * (m_m1 + m_m2 - Q) * (m_m1 + m_m2 + Q) +
                            pow(m_ce, 2) * (pow(m_m1 - m_m2, 2) - pow(Q, 2))) *
                    log(-((-2 * pow(m_m1, 2) - 2 * pow(m_m2, 2) + pow(Q, 2) +
                            sqrt(-4 * pow(m_m1, 2) + pow(Q, 2)) * sqrt(-4 * pow(m_m2, 2) + pow(Q, 2))) /
                            (2 * pow(m_m1, 2) + 2 * pow(m_m2, 2) - pow(Q, 2) +
                                    sqrt(-4 * pow(m_m1, 2) + pow(Q, 2)) * sqrt(-4 * pow(m_m2, 2) + pow(Q, 2)))))) /
            (768. * pow(m_lam, 4) * M_PI * pow(Q, 2) * (-4 * pow(m_m1, 2) + pow(Q, 2)));
}

/**
 * Compute the cross section for chi1 + chi1 -> w^+ + w^-.
 * @param Q Center of mass energy.
 * @return sigma Sigma(chi1 + chi2 -> w^+ + w^-)
 */
double DipoleDM::sigma_x1_x2_to_w_w(double Q) const {
    if (Q <= m_m1 + m_m2 || Q <= 2 * kW_BOSON_MASS) {
        return 0.0;
    }
    return -(kALPHA_EM * sqrt(-4 * pow(kW_BOSON_MASS, 2) + pow(Q, 2)) *
            sqrt(-2 * (pow(m_m1, 2) + pow(m_m2, 2)) +
                         pow(pow(m_m1, 2) - pow(m_m2, 2), 2) / pow(Q, 2) + pow(Q, 2)) *
            (48 * pow(kW_BOSON_MASS, 6) + 68 * pow(kW_BOSON_MASS, 4) * pow(Q, 2) -
                    16 * pow(kW_BOSON_MASS, 2) * pow(Q, 4) - pow(Q, 6)) *
            (-(pow(m_ce, 2) * (m_m1 + m_m2 - Q) * (m_m1 + m_m2 + Q) * (2 * pow(m_m1 - m_m2, 2) + pow(Q, 2))) +
                    pow(m_cm, 2) * (-2 * pow(pow(m_m1, 2) - pow(m_m2, 2), 2) +
                            (pow(m_m1, 2) + 6 * m_m1 * m_m2 + pow(m_m2, 2)) * pow(Q, 2) + pow(Q, 4)))) /
            (96. * pow(m_lam, 2) * pow(kW_BOSON_MASS, 4) * (m_m1 - m_m2 - Q) * (m_m1 + m_m2 - Q) * pow(Q, 4) *
                    (m_m1 - m_m2 + Q) * (m_m1 + m_m2 + Q));
}

/**
 * Compute the cross section for chi1 + chi1 -> f + fbar.
 * @param Q Center of mass energy.
 * @param mf Mass of final state fermion.
 * @return sigma Sigma(chi1 + chi2 -> f + fbar)
 */
double DipoleDM::sigma_x1_x2_to_f_f(double Q, double mf) const {
    if (Q <= m_m1 + m_m2 || Q <= 2 * mf) {
        return 0.0;
    }
    return (kALPHA_EM * sqrt(-4 * pow(mf, 2) + pow(Q, 2)) * (2 * pow(mf, 2) + pow(Q, 2)) *
            sqrt(-2 * (pow(m_m1, 2) + pow(m_m2, 2)) + pow
                    (pow(m_m1, 2) - pow(m_m2, 2), 2) / pow(Q, 2) + pow(Q, 2)) *
            (-(pow(m_ce, 2) * (m_m1 + m_m2 - Q) * (m_m1 + m_m2 + Q) * (2 * pow(m_m1 - m_m2, 2) + pow(Q, 2))) +
                    pow(m_cm, 2) * (-2 * pow(pow(m_m1, 2) - pow(m_m2, 2), 2) +
                            (pow(m_m1, 2) + 6 * m_m1 * m_m2 + pow(m_m2, 2)) * pow(Q, 2) + pow(Q, 4)))) /
            (6. * pow(m_lam, 2) * (m_m1 - m_m2 - Q) * (m_m1 + m_m2 - Q) * pow(Q, 4) * (m_m1 - m_m2 + Q) *
                    (m_m1 + m_m2 + Q));
}

/**
 * Integrand of the transfer function
 * @param T Temperature
 * @param w energy of the SM photon
 * @return Transfer function integral
 */
double DipoleDM::gamma_integrand(double w, double T) const {
    const double gw = 1.0 / (exp(w / T) - 1.0);
    const double temp1 = pow(m_m2, 2);
    const double temp2 = pow(m_m1, 5);
    const double temp3 = pow(w, 2);
    const double temp4 = pow(m_m1, 4);
    const double temp5 = pow(m_width2, 2);
    const double temp6 = temp1 + temp5;
    const double temp7 = pow(m_m1, 2);
    const double temp8 = pow(m_m1, 3);
    const double temp9 = pow(m_m2, 4);
    const double temp10 = pow(m_ce, 2);
    const double temp11 = pow(m_cm, 2);
    const double temp12 = temp10 + temp11;
    const double temp13 = pow(temp12, 2);
    const double temp14 = pow(m_lam, -4);
    const double temp15 = 2 * w;
    const double temp16 = m_m1 + temp15;
    const double temp17 = 4 * temp8 * w;
    const double temp18 = -4 * m_m1 * temp1 * w;
    const double temp19 = -2 * temp3;
    const double temp20 = temp1 + temp19;
    const double temp21 = -2 * temp20 * temp7;
    const double temp22 = temp1 * temp6;
    const double temp23 = temp17 + temp18 + temp21 + temp22 + temp4;
    const double temp24 = pow(m_m1, 7);
    const double temp25 = pow(m_m2, 6);
    const double temp26 = pow(w, 3);
    const double temp27 = pow(m_m1, 6);
    const double temp28 = pow(w, 4);
    const double temp29 = pow(w, 5);
    const double temp30 = pow(w, 6);
    const double temp31 = pow(w, 7);
    const double temp32 = 2 * m_m1 * w;
    const double temp33 = pow(temp16, 3);
    const double temp34 = 6 * temp2 * w;
    const double temp35 = 3 * temp1 * w;
    const double temp36 = 4 * temp26;
    const double temp37 = temp35 + temp36;
    const double temp38 = -4 * temp37 * temp8;
    const double temp39 = -temp5;
    const double temp40 = temp1 + temp39;
    const double temp41 = 2 * temp40 * temp9;
    const double temp42 = -3 * temp1;
    const double temp43 = -12 * temp3;
    const double temp44 = temp42 + temp43 + temp5;
    const double temp45 = temp1 * temp44 * temp7;
    const double temp46 = 6 * temp9 * w;
    const double temp47 = -2 * temp1 * temp5 * w;
    const double temp48 = temp46 + temp47;
    const double temp49 = m_m1 * temp48;
    const double temp50 = temp27 + temp34 + temp38 + temp41 + temp45 +
            temp49;
    const double temp51 = -2 * temp1 * w;
    const double temp52 = temp1 * temp5;
    const double temp53 = pow(m_m1, 9);
    const double temp54 = -2 * temp1 * temp7;
    const double temp55 = -4 * temp8 * w;
    const double temp56 = 4 * m_m1 * temp1 * w;
    const double temp57 = 4 * temp3 * temp7;
    const double temp58 = temp4 + temp52 + temp54 + temp55 + temp56 +
            temp57 + temp9;
    const double temp59 = log(temp58);
    const double temp60 = pow(m_m1, 8);
    const double temp61 = pow(m_m2, 8);
    const double temp62 = pow(m_width2, 4);
    const double temp63 = pow(m_m1, 11);
    const double temp64 = pow(temp16, -2);
    const double temp65 = -2 * temp1 * temp4;
    const double temp66 = -4 * temp1 * temp8 * w;
    const double temp67 = temp1 * temp6 * temp7;
    const double temp68 = 4 * m_m1 * temp1 * temp6 * w;
    const double temp69 = 4 * temp1 * temp3 * temp6;
    const double temp70 = temp27 + temp65 + temp66 + temp67 + temp68 +
            temp69;
    const double temp71 = temp64 * temp70;
    const double temp72 = log(temp71);
    const double temp73 = pow(m_m1, 10);
    const double temp74 = -temp1;
    const double temp75 = temp32 + temp7 + temp74;
    const double temp76 = 1 / temp16;
    const double temp77 = -4 * temp3;
    const double temp78 = 1 / m_m2;
    const double temp79 = -temp7;
    const double temp80 = temp1 + temp32 + temp79;
    const double temp81 = 1 / m_width2;
    const double temp82 = temp78 * temp80 * temp81;
    const double temp83 = atan(temp82);
    const double temp84 = -temp24;
    const double temp85 = temp27 * w;
    const double temp86 = temp1 + temp3;
    const double temp87 = 3 * temp2 * temp86;
    const double temp88 = 2 * temp40 * temp9 * w;
    const double temp89 = temp42 + temp5 + temp77;
    const double temp90 = temp1 * temp7 * temp89 * w;
    const double temp91 = -3 * temp9;
    const double temp92 = -4 * temp28;
    const double temp93 = -6 * temp3;
    const double temp94 = temp5 + temp93;
    const double temp95 = temp1 * temp94;
    const double temp96 = temp91 + temp92 + temp95;
    const double temp97 = temp8 * temp96;
    const double temp98 = -(temp1 * temp3 * temp5);
    const double temp99 = -3 * temp3;
    const double temp100 = temp5 + temp99;
    const double temp101 = -(temp100 * temp9);
    const double temp102 = temp101 + temp25 + temp98;
    const double temp103 = m_m1 * temp102;
    const double temp104 = temp103 + temp84 + temp85 + temp87 + temp88 +
            temp90 + temp97;
    const double temp105 = -(m_m1 * temp1);
    const double temp106 = temp105 + temp51 + temp8;
    const double temp107 = m_m1 * m_m2 * m_width2;
    const double temp108 = 2 * m_m2 * w * m_width2;
    const double temp109 = temp107 + temp108;
    const double temp110 = 1 / temp109;
    const double temp111 = temp106 * temp110;
    const double temp112 = atan(temp111);
    const double temp113 = pow(temp16, 2);
    const double temp114 = pow(temp16, 4);
    const double temp115 = 5 * temp9;
    const double temp116 = -5 * temp1 * temp60;
    const double temp117 = -2 * temp53 * w;
    const double temp118 = 8 * temp1 * temp24 * w;
    const double temp119 = -12 * temp2 * temp9 * w;
    const double temp120 = 8 * temp25 * temp8 * w;
    const double temp121 = -3 * temp1 * temp5;
    const double temp122 = temp115 + temp121;
    const double temp123 = 2 * temp122 * temp27;
    const double temp124 = -10 * temp25;
    const double temp125 = 18 * temp5 * temp9;
    const double temp126 = 8 * temp1 * temp3 * temp5;
    const double temp127 = temp124 + temp125 + temp126;
    const double temp128 = temp127 * temp4;
    const double temp129 = -temp9;
    const double temp130 = temp129 + temp62;
    const double temp131 = 2 * m_m1 * temp130 * temp9 * w;
    const double temp132 = -6 * temp1 * temp5;
    const double temp133 = temp132 + temp62 + temp9;
    const double temp134 = -(temp133 * temp25);
    const double temp135 = -18 * temp1 * temp5;
    const double temp136 = -8 * temp3 * temp5;
    const double temp137 = temp115 + temp135 + temp136 + temp62;
    const double temp138 = temp137 * temp7 * temp9;
    const double temp139 = temp116 + temp117 + temp118 + temp119 +
            temp120 + temp123 + temp128 + temp131 + temp134 + temp138 + temp73;
    const double temp140 = -4 * temp1 * temp27;
    const double temp141 = 6 * temp4 * temp9;
    const double temp142 = -4 * temp25 * temp7;
    const double temp143 = -8 * temp1 * temp5 * temp8 * w;
    const double temp144 = 8 * m_m1 * temp5 * temp9 * w;
    const double temp145 = -(temp62 * temp9);
    const double temp146 = temp140 + temp141 + temp142 + temp143 +
            temp144 + temp145 + temp60 + temp61;
    const double temp147 = pow(temp16, -3);
    const double temp148 = temp5 + temp77;
    const double temp149 = pow(temp16, -4);
    const double temp150 = temp42 + temp5;
    const double temp151 = 57 * temp9;
    const double temp152 = -9 * temp5;
    const double temp153 = 9 * temp60;
    const double temp154 = 42 * temp24 * w;
    const double temp155 = -27 * temp1;
    const double temp156 = 48 * temp3;
    const double temp157 = temp155 + temp156;
    const double temp158 = temp157 * temp27;
    const double temp159 = 5 * temp1 * w;
    const double temp160 = temp159 + temp26;
    const double temp161 = -24 * temp160 * temp2;
    const double temp162 = 12 * temp150 * temp3 * temp9;
    const double temp163 = 3 * temp9;
    const double temp164 = -(temp1 * temp5);
    const double temp165 = 2 * temp3 * temp5;
    const double temp166 = temp163 + temp164 + temp165;
    const double temp167 = -12 * m_m1 * temp1 * temp166 * w;
    const double temp168 = -8 * temp28;
    const double temp169 = 12 * temp3;
    const double temp170 = temp152 + temp169;
    const double temp171 = temp1 * temp170;
    const double temp172 = temp151 + temp168 + temp171;
    const double temp173 = 2 * temp172 * temp8 * w;
    const double temp174 = -9 * temp9;
    const double temp175 = 32 * temp28;
    const double temp176 = -36 * temp3 * temp5;
    const double temp177 = 40 * temp3;
    const double temp178 = temp177 + temp5;
    const double temp179 = 3 * temp1 * temp178;
    const double temp180 = temp174 + temp175 + temp176 + temp179;
    const double temp181 = temp1 * temp180 * temp7;
    const double temp182 = 27 * temp9;
    const double temp183 = -56 * temp28;
    const double temp184 = 44 * temp3;
    const double temp185 = temp184 + temp5;
    const double temp186 = -3 * temp1 * temp185;
    const double temp187 = temp182 + temp183 + temp186;
    const double temp188 = temp187 * temp4;
    const double temp189 = temp153 + temp154 + temp158 + temp161 +
            temp162 + temp167 + temp173 + temp181 + temp188;
    const double temp190 = 2 * temp5;
    const double temp191 = 1 / temp70;
    const double temp192 = -(temp24 * w);
    const double temp193 = 3 * temp1 * temp2 * w;
    const double temp194 = 4 * temp1;
    const double temp195 = temp194 + temp3;
    const double temp196 = -(temp195 * temp27);
    const double temp197 = -(temp25 * temp5);
    const double temp198 = 6 * temp1;
    const double temp199 = 2 * temp3;
    const double temp200 = temp198 + temp199 + temp39;
    const double temp201 = temp1 * temp200 * temp4;
    const double temp202 = m_m1 * temp6 * temp9 * w;
    const double temp203 = 3 * temp1;
    const double temp204 = temp203 + temp5;
    const double temp205 = -(temp1 * temp204 * temp8 * w);
    const double temp206 = -4 * temp25;
    const double temp207 = temp1 * temp3 * temp5;
    const double temp208 = -temp3;
    const double temp209 = temp190 + temp208;
    const double temp210 = temp209 * temp9;
    const double temp211 = temp206 + temp207 + temp210;
    const double temp212 = temp211 * temp7;
    const double temp213 = temp192 + temp193 + temp196 + temp197 +
            temp201 + temp202 + temp205 + temp212 + temp60 + temp61;
    const double temp214 = -temp27;
    const double temp215 = 3 * temp1 * temp4;
    const double temp216 = -2 * temp2 * w;
    const double temp217 = 4 * temp1 * temp8 * w;
    const double temp218 = temp5 + temp74;
    const double temp219 = 2 * m_m1 * temp1 * temp218 * w;
    const double temp220 = temp6 * temp9;
    const double temp221 = -(temp1 * temp204 * temp7);
    const double temp222 = temp214 + temp215 + temp216 + temp217 +
            temp219 + temp220 + temp221;
    return (1024 * temp13 * temp14 * temp149 * temp2 * temp31 * (2 * temp24 +
            36 * temp2 * temp3 + temp4 * (40 * temp26 + temp51) + temp8 * (2 * temp1 * temp148
            + 16 * temp28 - 6 * temp9) + 2 * m_m1 * (2 * temp25 + 3 * temp1 * temp3 * temp5 +
            (temp190 - 5 * temp3) * temp9) + 14 * temp27 * w + temp1 * (-17 * temp1 - 8 * temp3
            + 7 * temp5) * temp7 * w + 5 * temp6 * temp9 * w)) / (3. * pow(temp23, 2)) -
            (4 * m_m1 * temp13 * temp14 * temp78 * temp81 * (-3 * temp112 * temp146 -
                    3 * temp146 * temp83 + (3 * temp139 * temp78 * temp81) / (1 +
                    pow(temp80, 2) / (pow(m_m2, 2) * pow(m_width2, 2))) -
                    3 * m_m2 * temp139 * temp191 * temp7 * m_width2 +
                    m_m2 * (-12 * temp149 * temp189 * temp3 + 3 * temp222 * temp59 - 3 * temp222 * temp72
                            + 12 * temp106 * temp191 * temp213 * temp7 * temp76 + 4 * temp147 * temp189 * w -
                            (12 * temp213 * (temp7 + temp74 - 2 * m_m1 * w)) / (temp21 + temp22 + temp4 +
                                    temp55 + temp56) + 4 * temp147 * temp3 * (21 * temp24 - 12 * temp2 * (5 * temp1 +
                            3 * temp3) + (-9 * temp1 * temp148 + temp151 - 40 * temp28) * temp8 -
                            6 * m_m1 * (3 * temp25 + 6 * temp1 * temp3 * temp5 - temp5 * temp9) + 48 * temp27 * w +
                            4 * temp1 * (30 * temp1 + temp152 + 16 * temp3) * temp7 * w + 12 * temp150 * temp9 * w
                            - 4 * temp4 * (28 * temp26 + 33 * temp1 * w))) * m_width2)) / 3. +
            (2 * m_m1 * temp13 * temp14 * (temp1 +
                    temp7) * (4 * m_m1 * temp75 * (-48 * temp24 * temp3 + 24 * m_m1 * temp25 * temp3 +
                    24 * temp16 * temp27 * temp3 + 48 * temp1 * temp16 * temp3 * temp4 +
                    9 * temp1 * temp113 * temp2 * temp59 + 12 * temp1 * temp24 * temp59 +
                    12 * temp113 * temp24 * temp59 - 3 * m_m1 * temp113 * temp25 * temp59 -
                    3 * temp27 * temp33 * temp59 - 12 * temp1 * temp33 * temp4 * temp59 +
                    12 * temp1 * temp2 * temp5 * temp59 + 18 * temp16 * temp25 * temp5 * temp59 +
                    9 * temp1 * temp16 * temp4 * temp5 * temp59 - 12 * temp16 * temp59 * temp60 -
                    3 * temp16 * temp59 * temp61 - 48 * temp1 * temp2 * temp28 * temp64 +
                    12 * temp16 * temp25 * temp59 * temp7 - 3 * temp1 * temp33 * temp5 * temp59 * temp7 -
                    9 * temp1 * temp113 * temp2 * temp72 - 12 * temp1 * temp24 * temp72 -
                    12 * temp113 * temp24 * temp72 + 3 * m_m1 * temp113 * temp25 * temp72 +
                    3 * temp27 * temp33 * temp72 + 12 * temp1 * temp33 * temp4 * temp72 -
                    12 * temp1 * temp2 * temp5 * temp72 - 18 * temp16 * temp25 * temp5 * temp72 -
                    9 * temp1 * temp16 * temp4 * temp5 * temp72 + 12 * temp16 * temp60 * temp72 +
                    3 * temp16 * temp61 * temp72 - 12 * temp16 * temp25 * temp7 * temp72 +
                    3 * temp1 * temp33 * temp5 * temp7 * temp72 + 48 * temp27 * temp28 * temp76 +
                    48 * temp1 * temp27 * temp3 * temp76 - 128 * temp30 * temp4 * temp76 +
                    24 * temp1 * temp3 * temp4 * temp5 * temp76 -
                    48 * temp1 * temp28 * temp5 * temp7 * temp76 - 48 * temp1 * temp28 * temp8 -
                    24 * temp1 * temp113 * temp3 * temp8 + 48 * temp1 * temp3 * temp5 * temp8 +
                    3 * temp1 * temp114 * temp59 * temp8 + 3 * temp25 * temp59 * temp8 +
                    128 * temp1 * temp30 * temp64 * temp8 - 3 * temp1 * temp114 * temp72 * temp8 -
                    3 * temp25 * temp72 * temp8 - 72 * m_m1 * temp3 * temp5 * temp9 -
                    12 * temp2 * temp59 * temp9 - 9 * temp16 * temp4 * temp59 * temp9 +
                    9 * m_m1 * temp113 * temp5 * temp59 * temp9 - 3 * temp16 * temp59 * temp62 * temp9 +
                    3 * temp33 * temp59 * temp7 * temp9 - 36 * temp16 * temp5 * temp59 * temp7 * temp9 +
                    12 * temp2 * temp72 * temp9 + 9 * temp16 * temp4 * temp72 * temp9 -
                    9 * m_m1 * temp113 * temp5 * temp72 * temp9 + 3 * temp16 * temp62 * temp72 * temp9 -
                    3 * temp33 * temp7 * temp72 * temp9 + 36 * temp16 * temp5 * temp7 * temp72 * temp9 -
                    24 * temp3 * temp4 * temp76 * temp9 + 48 * temp28 * temp7 * temp76 * temp9 -
                    48 * temp3 * temp8 * temp9 - 9 * temp5 * temp59 * temp8 * temp9 +
                    9 * temp5 * temp72 * temp8 * temp9 - 24 * m_m2 * temp104 * temp112 * m_width2 -
                    24 * m_m2 * temp104 * temp83 * m_width2) -
                    2 * temp147 * temp23 * (-288 * temp1 * temp26 * temp27 - 1248 * temp1 * temp2 * temp28
                            + 624 * temp24 * temp28 + 96 * m_m1 * temp25 * temp28 + 96 * temp27 * temp29 -
                            24 * temp1 * temp24 * temp3 - 1408 * temp2 * temp30 - 2304 * temp1 * temp29 * temp4 -
                            1280 * temp31 * temp4 + 24 * temp1 * temp2 * temp3 * temp5 +
                            48 * temp1 * temp26 * temp4 * temp5 + 24 * temp3 * temp53 + 6 * temp2 * temp25 * temp59
                            - 384 * temp1 * temp26 * temp27 * temp59 + 288 * temp1 * temp2 * temp28 * temp59 -
                            336 * temp24 * temp28 * temp59 - 96 * m_m1 * temp25 * temp28 * temp59 -
                            288 * temp27 * temp29 * temp59 - 288 * temp1 * temp24 * temp3 * temp59 +
                            1152 * temp1 * temp29 * temp4 * temp59 + 144 * temp25 * temp26 * temp5 * temp59 -
                            252 * temp1 * temp2 * temp3 * temp5 * temp59 +
                            216 * m_m1 * temp25 * temp3 * temp5 * temp59 -
                            648 * temp1 * temp26 * temp4 * temp5 * temp59 - 6 * temp1 * temp53 * temp59 +
                            72 * temp3 * temp53 * temp59 + 240 * temp26 * temp60 - 48 * temp26 * temp59 * temp60
                            - 24 * temp26 * temp59 * temp61 - 36 * m_m1 * temp3 * temp59 * temp61 +
                            3 * temp59 * temp63 + 96 * temp25 * temp26 * temp7 -
                            288 * temp1 * temp29 * temp5 * temp7 - 96 * temp25 * temp26 * temp59 * temp7 -
                            288 * temp1 * temp29 * temp5 * temp59 * temp7 - 6 * temp2 * temp25 * temp72 +
                            384 * temp1 * temp26 * temp27 * temp72 - 288 * temp1 * temp2 * temp28 * temp72 +
                            336 * temp24 * temp28 * temp72 + 96 * m_m1 * temp25 * temp28 * temp72 +
                            288 * temp27 * temp29 * temp72 + 288 * temp1 * temp24 * temp3 * temp72 -
                            1152 * temp1 * temp29 * temp4 * temp72 - 144 * temp25 * temp26 * temp5 * temp72 +
                            252 * temp1 * temp2 * temp3 * temp5 * temp72 -
                            216 * m_m1 * temp25 * temp3 * temp5 * temp72 +
                            648 * temp1 * temp26 * temp4 * temp5 * temp72 + 6 * temp1 * temp53 * temp72 -
                            72 * temp3 * temp53 * temp72 + 48 * temp26 * temp60 * temp72 +
                            24 * temp26 * temp61 * temp72 + 36 * m_m1 * temp3 * temp61 * temp72 -
                            3 * temp63 * temp72 + 96 * temp25 * temp26 * temp7 * temp72 +
                            288 * temp1 * temp29 * temp5 * temp7 * temp72 + 24 * temp25 * temp3 * temp8 -
                            1408 * temp1 * temp30 * temp8 - 144 * temp1 * temp28 * temp5 * temp8 +
                            768 * temp1 * temp30 * temp59 * temp8 + 18 * temp25 * temp5 * temp59 * temp8 -
                            720 * temp1 * temp28 * temp5 * temp59 * temp8 - 3 * temp59 * temp61 * temp8 -
                            768 * temp1 * temp30 * temp72 * temp8 - 18 * temp25 * temp5 * temp72 * temp8 +
                            720 * temp1 * temp28 * temp5 * temp72 * temp8 + 3 * temp61 * temp72 * temp8 -
                            24 * temp2 * temp3 * temp9 - 48 * temp26 * temp4 * temp9 -
                            288 * m_m1 * temp28 * temp5 * temp9 + 252 * temp2 * temp3 * temp59 * temp9 +
                            648 * temp26 * temp4 * temp59 * temp9 - 18 * temp2 * temp5 * temp59 * temp9 +
                            288 * m_m1 * temp28 * temp5 * temp59 * temp9 - 24 * temp26 * temp59 * temp62 * temp9 -
                            36 * m_m1 * temp3 * temp59 * temp62 * temp9 + 288 * temp29 * temp7 * temp9 -
                            288 * temp26 * temp5 * temp7 * temp9 + 288 * temp29 * temp59 * temp7 * temp9 +
                            288 * temp26 * temp5 * temp59 * temp7 * temp9 - 252 * temp2 * temp3 * temp72 * temp9 -
                            648 * temp26 * temp4 * temp72 * temp9 + 18 * temp2 * temp5 * temp72 * temp9 -
                            288 * m_m1 * temp28 * temp5 * temp72 * temp9 + 24 * temp26 * temp62 * temp72 * temp9 +
                            36 * m_m1 * temp3 * temp62 * temp72 * temp9 - 288 * temp29 * temp7 * temp72 * temp9 -
                            288 * temp26 * temp5 * temp7 * temp72 * temp9 + 144 * temp28 * temp8 * temp9 -
                            72 * temp3 * temp5 * temp8 * temp9 + 720 * temp28 * temp59 * temp8 * temp9 -
                            3 * temp59 * temp62 * temp8 * temp9 - 720 * temp28 * temp72 * temp8 * temp9 +
                            3 * temp62 * temp72 * temp8 * temp9 + 24 * temp25 * temp4 * temp59 * w -
                            36 * temp1 * temp27 * temp5 * temp59 * w - 72 * temp1 * temp59 * temp60 * w +
                            108 * temp25 * temp5 * temp59 * temp7 * w - 18 * temp59 * temp61 * temp7 * w -
                            24 * temp25 * temp4 * temp72 * w + 36 * temp1 * temp27 * temp5 * temp72 * w +
                            72 * temp1 * temp60 * temp72 * w - 108 * temp25 * temp5 * temp7 * temp72 * w +
                            18 * temp61 * temp7 * temp72 * w + 30 * temp59 * temp73 * w - 30 * temp72 * temp73 * w +
                            36 * temp27 * temp59 * temp9 * w - 72 * temp4 * temp5 * temp59 * temp9 * w -
                            18 * temp59 * temp62 * temp7 * temp9 * w - 36 * temp27 * temp72 * temp9 * w +
                            72 * temp4 * temp5 * temp72 * temp9 * w + 18 * temp62 * temp7 * temp72 * temp9 * w -
                            12 * m_m2 * temp112 * temp33 * temp50 * m_width2 -
                            12 * m_m2 * temp33 * temp50 * temp83 * m_width2))) / (3. * pow(temp52 +
                                                                                                   pow(temp75, 2), 2));
}

/**
 * Compute the thermal transfer function.
 * @param T
 * @return
 */
double DipoleDM::gamma(double T) const {
    using namespace boost::math::quadrature;
    using boost::math::pow;

    auto f = [T, this](double w) {
        return gamma_integrand(w, T);
    };
    double gam = gauss_kronrod<double, 15>::integrate(f, 0.0, std::numeric_limits<double>::infinity(), 5, 1e-9);

    return gam / (48.0 * pow<3>(M_PI * m_m1) * 2.0 * T);
}

/**
 * Returns the annihilations cross section for chi_i + chi_j -> X
 * @param Q Center of mass energy
 * @param i Index of chi_i. Equal to 1 or 2.
 * @param j Index of chi_j. Equal to 1 or 2.
 */
double DipoleDM::annihilation_cross_section(double Q, int i, int j) const {
    if (i == 1 && j == 1) {
        return (sigma_x1_x1_to_g_g(Q) +
                sigma_x1_x1_to_x2_x2(Q));
    } else if ((i == 1 && j == 2) || (i == 2 && j == 1)) {
        return (sigma_x1_x2_to_w_w(Q) +
                sigma_x1_x2_to_f_f(Q, kELECTRON_MASS) +
                sigma_x1_x2_to_f_f(Q, kMUON_MASS) +
                sigma_x1_x2_to_f_f(Q, kTAU_MASS) +
                sigma_x1_x2_to_f_f(Q, kUP_QUARK_MASS) +
                sigma_x1_x2_to_f_f(Q, kDOWN_QUARK_MASS) +
                sigma_x1_x2_to_f_f(Q, kSTRANGE_QUARK_MASS) +
                sigma_x1_x2_to_f_f(Q, kCHARM_QUARK_MASS) +
                sigma_x1_x2_to_f_f(Q, kBOTTOM_QUARK_MASS) +
                sigma_x1_x2_to_f_f(Q, kTOP_QUARK_MASS));
    } else {
        return 0.0;
    }

}

/**
 * Returns the effective W factor for computing the effective thermal cross
 * section.
 * @param Q Center of mass energy.
 */
double DipoleDM::weff(double Q) const {
    double sig11 = annihilation_cross_section(Q, 1, 1);
    double sig12 = annihilation_cross_section(Q, 1, 2);
    double sig22 = annihilation_cross_section(Q, 2, 2);
    double s = Q * Q;
    double p11 = sqrt(s - pow(m_m1 + m_m1, 2)) * sqrt(s - pow(m_m1 - m_m1, 2)) / (2 * Q);
    double p12 = sqrt(s - pow(m_m1 + m_m2, 2)) * sqrt(s - pow(m_m1 - m_m2, 2)) / (2 * Q);
    double p22 = sqrt(s - pow(m_m2 + m_m2, 2)) * sqrt(s - pow(m_m2 - m_m2, 2)) / (2 * Q);
    double peff = 0.5 * sqrt(s - 4 * m_m1 * m_m1);

    double w11 = 4 * p11 * Q * sig11;
    double w12 = 4 * p12 * Q * sig12;
    double w22 = 4 * p22 * Q * sig22;

    double wef = (p11 * w11 + 2 * p12 * w12 + p22 * w22) / peff;
}

/**
 * Returns the integrand of the effective thermally average cross section.
 * @param peff Effective momentum.
 * @param x Mass of DM over its temperature.
 */
double DipoleDM::thermal_cross_section_integrand(double peff, double x) const {
    using namespace special_functions;

    const double sqrts = 2.0 * sqrt(peff * peff + m_m1 * m_m1);
    const double T = m_m1 / x;
    return peff * peff * weff(sqrts) * besselk1e(sqrts / T) * exp(2.0 * x - sqrts / T);
}

/**
 * Compute the effective, thermally average annihilation cross section for
 * chi1 + chi2
 * @param x mass of chi1 / temperature.
 * @return <sigma_eff*v>
 */
double DipoleDM::thermal_cross_section(double x) const {
    using boost::math::quadrature::gauss_kronrod;
    using namespace special_functions;

    const double rat = m_m2 / m_m1;
    const double denom = m_m1 * m_m1 * (besselk2e(x) + rat * rat * besselk2e(rat * x) * exp(x * (1.0 - rat)));
    const double pf = 1.0 / (denom * denom * x / m_m1);
    auto integrand = [this, &x](double peff) {
        return thermal_cross_section_integrand(peff, x);
    };

    return pf * gauss_kronrod<double, 15>::integrate(integrand, 0.0, std::numeric_limits<double>::infinity(), 5, 1e-9);
}


}
}

#endif //LANRE_DM_MODELS_KINETIC_RECOUPLING_HPP
