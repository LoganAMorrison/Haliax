//
// Created by Logan Morrison on 3/15/20.
//

#ifndef LANRE_DM_MODELS_HIGGS_PORTAL_HPP
#define LANRE_DM_MODELS_HIGGS_PORTAL_HPP

#include "lanre/constants.hpp"
#include "lanre/utils.hpp"
#include "lanre/special_functions/besselk.hpp"
#include "lanre/diffeq/base.hpp"
#include "lanre/diffeq/radau.hpp"
#include "lanre/diffeq/rodas.hpp"
#include "lanre/diffeq/function.hpp"
#include "lanre/cosmology/standard_model.hpp"
#include "lanre/cosmology/thermodynamic_particle.hpp"
#include <string>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <cmath>
#include <Eigen/Dense>

namespace lanre {
namespace dm_models {

using namespace diffeq;
using namespace cosmology;

class HiggsPortal {
private:
    double m_mx;
    double m_ms;
    double m_gsxx;
    double m_smix;
    double m_width_s;
    ThermodynamicParticle m_chi;

    double width_s_to_ff(const std::string &) const;

    double width_s_to_ww() const;

    double width_s_to_zz() const;

    double width_s_to_hh() const;

    double width_s_to_xx() const;

    double sigma_xx_to_ff(double, const std::string &,
                          const std::string &) const;

    double sigma_xx_to_ww(double, const std::string &) const;

    double sigma_xx_to_zz(double, const std::string &) const;

    double sigma_xx_to_hh(double, const std::string &) const;

    double sigma_xx_to_ss(double, const std::string &) const;

    double thermal_cross_section_integrand(
            double, double, const std::string &,
            const std::string &) const;

public:
    HiggsPortal(double mx, double ms, double gsxx, double smix)
            : m_mx(mx), m_ms(ms), m_gsxx(gsxx), m_smix(smix),
              m_chi{mx, 2.0, 1} {
        m_width_s = scalar_partial_width("all");
    }

    double get_mx() const;

    double get_ms() const;

    double get_gsxx() const;

    double get_smix() const;

    double get_width_s() const;

    void set_mx(double);

    void set_ms(double);

    void set_gsxx(double);

    void set_smix(double);

    double scalar_partial_width(const std::string &) const;

    double annihilation_cross_section(double,
                                      const std::string &, const std::string &) const;

    double thermal_cross_section(double, const std::string &,
                                 const std::string &) const;

    double relic_density(double, double, double, double, const std::string &);
};

struct HiggsPortalBoltzmann : public ODEFunction {

    std::shared_ptr<HiggsPortal> model;
    cosmology::ThermodynamicParticle chi;

    HiggsPortalBoltzmann(std::shared_ptr<HiggsPortal> t_model)
            : model(t_model), chi{model->get_mx(), 2.0, 1} {}

    void dudt(Vector<double> &dw, const Vector<double> &w, const double logx) override {
        double x = exp(logx);
        double T = chi.get_mass() / x;
        double s = cosmology::sm_entropy_density(T);

        double weq = log(chi.neq(T) / s);
        double ww = w(0);

        double pf = -sqrt(M_PI / 45) * kM_PLANK * cosmology::sm_sqrt_gstar(T) * T;
        double sigmav = model->thermal_cross_section(x, "all", "all");

        // dW_e / dlogx
        dw(0) = pf * sigmav * (exp(ww) - exp(2.0 * weq - ww));
    }

    void dfdu(Matrix<double> &J, const Vector<double> &w, const double logx) override {
        double x = exp(logx);
        double T = chi.get_mass() / x;
        double s = cosmology::sm_entropy_density(T);

        double weq = log(chi.neq(T) / s);
        double ww = w(0);

        double pf = -sqrt(M_PI / 45) * kM_PLANK * cosmology::sm_sqrt_gstar(T) * T;
        double sigmav = model->thermal_cross_section(x, "all", "all");

        // dW_e / dlogx
        J(0, 0) = pf * sigmav * (exp(ww) + exp(2.0 * weq - ww));
    }
};

/*
 * Getters and setters
 */

/*
 * Get the value of dark matter mass
 */
double HiggsPortal::get_mx() const { return m_mx; }

/*
 * Get the value of the scalar mediator mass
 */
double HiggsPortal::get_ms() const { return m_ms; }

/*
 * Get the value of dark matter-scalar medaitor coupling
 */
double HiggsPortal::get_gsxx() const { return m_gsxx; }

/*
 * Get the value of the scalar-mediator-Higgs boson mixing angle
 */
double HiggsPortal::get_smix() const { return m_smix; }

/*
 * Get the scalar-mediator width
 */
double HiggsPortal::get_width_s() const { return m_width_s; }

/*
 * Set the value of the dark matter mass
 */
void HiggsPortal::set_mx(double mx) {
    m_mx = mx;
    m_width_s = scalar_partial_width("all");
}

/*
 * Set the value of the scalar-mediator mass
 */
void HiggsPortal::set_ms(double ms) {
    m_ms = ms;
    m_width_s = scalar_partial_width("all");
}

void HiggsPortal::set_gsxx(double gsxx) {
    m_gsxx = gsxx;
    m_width_s = scalar_partial_width("all");
}

void HiggsPortal::set_smix(double smix) {
    m_smix = smix;
    m_width_s = scalar_partial_width("all");
}


/*
 * Scalar-mediator partial decay widths
 */

/*
 * Partial width for scalar-mediator to dark matter
 */
double HiggsPortal::width_s_to_xx() const {
    if (2.0 * m_mx >= m_ms) {
        return 0.0;
    }
    return (-(pow(m_gsxx, 2) * pow(-4 * pow(m_mx, 2) + pow(m_ms, 2), 1.5) *
            (-1 + pow(m_smix, 2))) /
            (8. * pow(m_ms, 2) * M_PI));
}

/*
 * Partial width for scalar-mediator to SM fermions
 * @param f String representing the final state fermion
 */
double HiggsPortal::width_s_to_ff(const std::string &f) const {
    const double mf = string_to_fermion_mass(f);
    const double ncol = string_to_num_colors(f);

    if (2.0 * mf >= m_ms) {
        return 0.0;
    }
    double temp1 = mf * mf;
    return ((ncol * m_smix * m_smix *
            pow(m_ms * m_ms - 4 * temp1, 1.5) * temp1) /
            (8. * m_ms * m_ms * M_PI * kHIGGS_VEV * kHIGGS_VEV));
}

/*
 * Partial width for scalar-mediator to W bosons
 */
double HiggsPortal::width_s_to_ww() const {
    if (2.0 * kW_BOSON_MASS >= m_ms) {
        return 0.0;
    }
    double temp1 = m_ms * m_ms;
    double temp2 = kW_BOSON_MASS * kW_BOSON_MASS;
    return (kALPHA_EM * kALPHA_EM * M_PI * m_smix * m_smix * sqrt(temp1 - 4 * temp2) *
            (m_ms * m_ms * m_ms * m_ms + 12 * kW_BOSON_MASS * kW_BOSON_MASS * kW_BOSON_MASS * kW_BOSON_MASS -
                    4 * temp1 * temp2) * kHIGGS_VEV * kHIGGS_VEV) /
            (16. * m_ms * m_ms * kW_BOSON_MASS * kW_BOSON_MASS * kW_BOSON_MASS * kW_BOSON_MASS * kSIN_THETA_WEAK *
                    kSIN_THETA_WEAK * kSIN_THETA_WEAK * kSIN_THETA_WEAK);
}

/*
 * Partial width for scalar-mediator to Z bosons
 */
double HiggsPortal::width_s_to_zz() const {
    if (2.0 * kZ_BOSON_MASS >= m_ms) {
        return 0.0;
    }
    double temp1 = m_ms * m_ms;
    double temp2 = kZ_BOSON_MASS * kZ_BOSON_MASS;
    return (kALPHA_EM * kALPHA_EM * M_PI * m_smix * m_smix * sqrt(temp1 - 4 * temp2) *
            (m_ms * m_ms * m_ms * m_ms + 12 * pow(kZ_BOSON_MASS, 4) - 4 * temp1 * temp2) * kHIGGS_VEV * kHIGGS_VEV) /
            (32. * pow(kCOS_THETA_WEAK, 4) * m_ms * m_ms * pow(kZ_BOSON_MASS, 4) * pow(kSIN_THETA_WEAK, 4));
}

/*
 * Partial width for scalar-mediator to Higgs bosons
 */
double HiggsPortal::width_s_to_hh() const {
    if (2.0 * kHIGGS_MASS >= m_ms) {
        return 0.0;
    }

    const double temp1 = kHIGGS_MASS * kHIGGS_MASS;
    const double temp2 = m_ms * m_ms;
    const double temp3 = m_smix * m_smix;
    return (sqrt(-4 * temp1 + temp2) * pow(2 * temp1 + temp2, 2) * pow(-1 + temp3, 2) * temp3) /
            (32. * m_ms * m_ms * M_PI * kHIGGS_VEV * kHIGGS_VEV);
}

/*
 * Compute the total scalar-mediator width
 */
double HiggsPortal::scalar_partial_width(const std::string &state = "all") const {
    if (state == "all") {
        return (width_s_to_ff("e") +
                width_s_to_ff("mu") +
                width_s_to_ff("tau") +
                width_s_to_ff("u") +
                width_s_to_ff("c") +
                width_s_to_ff("t") +
                width_s_to_ff("d") +
                width_s_to_ff("s") +
                width_s_to_ff("b") +
                width_s_to_ww() +
                width_s_to_zz() +
                width_s_to_hh() +
                width_s_to_xx());
    } else if (state == "e" || state == "e e") {
        return width_s_to_ff("e");
    } else if (state == "mu" || state == "mu mu") {
        return width_s_to_ff("mu");
    } else if (state == "tau" || state == "tau tau") {
        return width_s_to_ff("tau");
    } else if (state == "u" || state == "u u") {
        return width_s_to_ff("u");
    } else if (state == "c" || state == "c c") {
        return width_s_to_ff("c");
    } else if (state == "t" || state == "t t") {
        return width_s_to_ff("t");
    } else if (state == "d" || state == "d d") {
        return width_s_to_ff("d");
    } else if (state == "s" || state == "s s") {
        return width_s_to_ff("s");
    } else if (state == "b" || state == "b b") {
        return width_s_to_ff("b");
    } else if (state == "w" || state == "w w") {
        return width_s_to_ww();
    } else if (state == "z" || state == "z z") {
        return width_s_to_zz();
    } else if (state == "h" || state == "h h") {
        return width_s_to_hh();
    } else {
        return 0.0;
    }
}


/*
 * Dark matter annihilation cross sections
 */

/*
 * Compute the dark matter annihilation cross section into SM fermions
 * @param Q Center of mass energy
 * @param channel String representing which channel to include, i.e. "s", "tu" or "all"
 * @param f String representing final state fermion
 */
double HiggsPortal::sigma_xx_to_ff(double Q, const std::string &channel, const std::string &f) const {
    if (channel != "s" && channel != "all") {
        return 0.0;
    }
    const double mf = string_to_fermion_mass(f);
    const double ncol = string_to_num_colors(f);

    if (Q < 2.0 * m_mx || Q < 2.0 * mf) {
        return 0.0;
    }
    const double temp1 = mf * mf;
    const double temp2 = Q * Q;
    const double temp3 = m_ms * m_ms;
    const double temp4 = m_width_s * m_width_s;
    const double temp5 = temp3 * temp4;
    return -(m_gsxx * m_gsxx * ncol * (-1 + m_smix) * m_smix * m_smix * (1 + m_smix) * temp1 *
            sqrt(-4 * m_mx * m_mx + temp2) * pow(-4 * temp1 + temp2, 1.5) *
            (pow(kHIGGS_MASS * kHIGGS_MASS - temp3, 2) + temp5)) /
            (16. * M_PI * pow(kHIGGS_MASS - Q, 2) * Q * Q * pow(kHIGGS_MASS + Q, 2) * (pow(-temp2 + temp3, 2) + temp5) *
                    kHIGGS_VEV * kHIGGS_VEV);
}

/*
 * Compute the dark matter annihilation cross section into W bosons
 * @param Q Center of mass energy
 * @param channel String representing which channel to include, i.e. "s", "tu" or "all"
 */
double HiggsPortal::sigma_xx_to_ww(double Q, const std::string &channel) const {
    if (channel != "s" && channel != "all") {
        return 0.0;
    }
    const double temp1 = Q * Q;
    const double temp2 = kW_BOSON_MASS * kW_BOSON_MASS;
    const double temp3 = kHIGGS_MASS * kHIGGS_MASS;
    const double temp4 = m_ms * m_ms;
    const double temp5 = m_ms * m_ms * m_ms * m_ms;
    const double temp6 = Q * Q * Q * Q;
    const double temp7 = m_width_s * m_width_s;
    return (-(kALPHA_EM * kALPHA_EM * m_gsxx * m_gsxx * M_PI * (-1 + m_smix) * m_smix * m_smix * (1 + m_smix) *
            sqrt(-4 * m_mx * m_mx + temp1) * sqrt(temp1 - 4 * temp2) *
            (12 * kW_BOSON_MASS * kW_BOSON_MASS * kW_BOSON_MASS * kW_BOSON_MASS - 4 * temp1 * temp2 + temp6) *
            (pow(kHIGGS_MASS, 4) - 2 * temp3 * temp4 + temp5 + temp4 * temp7) * kHIGGS_VEV * kHIGGS_VEV) /
            (32. * kW_BOSON_MASS * kW_BOSON_MASS * kW_BOSON_MASS * kW_BOSON_MASS * kSIN_THETA_WEAK * kSIN_THETA_WEAK *
                    kSIN_THETA_WEAK * kSIN_THETA_WEAK * pow(Q * Q * Q - Q * temp3, 2) *
                    (temp5 + temp6 + temp4 * (-2 * temp1 + temp7))));
}

/*
 * Compute the dark matter annihilation cross section into Z bosons
 * @param Q Center of mass energy
 * @param channel String representing which channel to include, i.e. "s", "tu" or "all"
 */
double HiggsPortal::sigma_xx_to_zz(double Q, const std::string &channel) const {
    if (channel != "s" && channel != "all") {
        return 0.0;
    }
    const double temp1 = Q * Q;
    const double temp2 = kZ_BOSON_MASS * kZ_BOSON_MASS;
    const double temp3 = kHIGGS_MASS * kHIGGS_MASS;
    const double temp4 = m_ms * m_ms;
    const double temp5 = m_ms * m_ms * m_ms * m_ms;
    const double temp6 = Q * Q * Q * Q;
    const double temp7 = m_width_s * m_width_s;
    return (-(kALPHA_EM * kALPHA_EM * m_gsxx * m_gsxx * M_PI * (-1 + m_smix) * m_smix * m_smix * (1 + m_smix) *
            sqrt(-4 * m_mx * m_mx + temp1) * sqrt(temp1 - 4 * temp2) *
            (12 * pow(kZ_BOSON_MASS, 4) - 4 * temp1 * temp2 + temp6) *
            (pow(kHIGGS_MASS, 4) - 2 * temp3 * temp4 + temp5 + temp4 * temp7) * kHIGGS_VEV * kHIGGS_VEV) /
            (64. * pow(kCOS_THETA_WEAK, 4) * pow(kZ_BOSON_MASS, 4) * kSIN_THETA_WEAK * kSIN_THETA_WEAK *
                    kSIN_THETA_WEAK * kSIN_THETA_WEAK * pow(Q * Q * Q - Q * temp3, 2) *
                    (temp5 + temp6 + temp4 * (-2 * temp1 + temp7))));
}

/*
 * Compute the dark matter annihilation cross section into Higgs bosons
 * @param Q Center of mass energy
 * @param channel String representing which channel to include, i.e. "s", "tu" or "all"
 */
double HiggsPortal::sigma_xx_to_hh(double Q, const std::string &channel) const {
    if (Q < 2.0 * m_mx && Q < 2.0 * kHIGGS_MASS) {
        return 0.0;
    }

    const double temp1 = pow(Q, 2);
    const double temp2 = pow(kHIGGS_MASS, 2);
    const double temp3 = pow(m_smix, 2);
    const double temp4 = pow(m_gsxx, 2);
    const double temp5 = 1 / M_PI;
    const double temp6 = pow(Q, -2);
    const double temp7 = pow(m_mx, 2);
    const double temp8 = -4 * temp7;
    const double temp9 = temp1 + temp8;
    const double temp10 = -4 * temp2;
    const double temp11 = temp1 + temp10;
    const double temp12 = temp11 * temp9;
    const double temp13 = sqrt(temp12);
    const double temp14 = -1 + temp3;
    const double temp15 = pow(temp14, 3);
    const double temp16 = pow(kHIGGS_VEV, -2);
    const double temp17 = pow(m_ms, 2);
    const double temp18 = -temp1;
    const double temp19 = pow(kHIGGS_MASS, 4);
    const double temp20 = -2 * temp2;
    const double temp21 = pow(m_gsxx, 4);
    const double temp22 = 1 / temp9;
    const double temp23 = pow(m_smix, 4);
    const double temp24 = pow(m_mx, 4);
    const double temp25 = 16 * temp24;
    const double temp26 = 2 * temp19;
    const double temp27 = -12 * temp2;
    const double temp28 = temp1 + temp27;
    const double temp29 = temp28 * temp7;
    const double temp30 = temp25 + temp26 + temp29;
    const double temp31 = temp11 * temp7;
    const double temp32 = temp19 + temp31;
    const double temp33 = 1 / temp32;
    const double temp34 = 2 * temp2;
    const double temp35 = temp1 + temp13 + temp20;
    const double temp36 = temp17 + temp34;
    const double temp37 = temp17 + temp18;
    const double temp38 = pow(temp37, 2);
    const double temp39 = pow(m_width_s, 2);
    const double temp40 = temp17 * temp39;
    const double temp41 = temp38 + temp40;
    const double temp42 = 1 / temp41;
    const double temp43 = -temp2;
    const double temp44 = temp1 + temp43;
    const double temp45 = 1 / temp44;
    const double temp46 = temp13 + temp18 + temp34;
    const double temp47 = 1 / temp46;
    const double temp48 = -(temp35 * temp47);
    const double temp49 = log(temp48);
    const double temp50 = pow(m_gsxx, 3);
    const double temp51 = temp18 + temp2;
    const double temp52 = pow(m_smix, 3);
    const double temp53 = -temp3;
    const double temp54 = 1 + temp53;
    const double temp55 = pow(temp54, 1.5);
    const double temp56 = 1 / kHIGGS_VEV;
    const double temp57 = 2 * temp13;
    const double temp58 = -temp13;
    const double temp59 = temp1 + temp20 + temp58;
    const double temp60 = 1 / temp35;
    const double temp61 = temp59 * temp60;
    const double temp62 = log(temp61);
    const double temp63 = -8 * temp7;
    const double temp64 = temp1 + temp34 + temp63;
    const double temp65 = temp49 * temp64;
    const double temp66 = temp57 + temp65;
    const double temp67 = -Q;
    const double temp68 = m_ms + temp67;
    const double temp69 = m_ms + Q;

    const double ampList[10] = {
            (-9 * temp13 * temp15 * temp16 * temp19 * temp3 * temp4 * temp5 * temp6) / (64. * pow(temp51, 2)),
            -(temp13 * temp15 * temp16 * temp3 * pow(temp36, 2) * temp4 * temp42 *
                    temp5 * temp6) /
                    64.,
            -(temp21 * temp22 * temp23 * temp5 * temp6 * (temp13 * temp30 * temp33 + temp62)) / 64.,
            (temp21 * temp22 * temp23 * (-(temp13 * temp30 * temp33) + temp49) * temp5 * temp6) / 64.,
            (3 * temp13 * temp15 * temp16 * (temp1 - temp17) * temp2 * temp3 * temp36 * temp4 * temp42 * temp45 *
                    temp5 * temp6) / 32.,
            (-3 * m_mx * temp2 * temp22 * temp45 * temp5 * temp50 * temp52 * temp55 * temp56 * temp6 * temp66) / 32.,
            (3 * m_mx * temp2 * temp22 * temp5 * temp50 * temp52 * temp55 * temp56 * temp6 *
                    (temp57 + temp62 * (temp18 + temp20 + 8 * temp7))) / (32. * temp51),
            -(
                    m_mx * temp22 * temp36 * temp42 * temp5 * temp50 * temp52 * temp55 * temp56 * temp6 * temp66 *
                            temp68 * temp69) /
                    32.,
            (m_mx * temp22 * temp36 * temp42 * temp5 * temp50 * temp52 * temp55 * temp56 * temp6 *
                    (-2 * temp13 + temp62 * temp64) * temp68 * temp69) / 32.,
            (temp21 * temp22 * temp23 * temp5 * temp6 * (temp58 +
                    ((temp19 - 16 * temp24 + 4 * temp1 * temp7) * log(16 * pow(temp32, 2))) / (temp1 + temp20))) / 32.};

    if (channel == "s") {
        return ampList[0] + ampList[1] + ampList[4];
    } else if (channel == "tu" || channel == "ut") {
        return ampList[2] + ampList[3] + ampList[9];
    } else {
        double sum = 0.0;
        for (int i = 0; i < 10; i++) {
            sum += ampList[i];
        }
        return sum;
    }
}

/*
 * Compute the dark matter annihilation cross section into scalar-mediators
 * @param Q Center of mass energy
 * @param channel String representing which channel to include, i.e. "s", "tu" or "all"
 */
double HiggsPortal::sigma_xx_to_ss(double Q, const std::string &channel) const {
    if (Q < 2.0 * m_mx || Q < 2.0 * m_ms) {
        return 0.0;
    }

    const double temp1 = pow(kHIGGS_MASS, 2);
    const double temp2 = pow(Q, 2);
    const double temp3 = pow(m_ms, 2);
    const double temp4 = pow(m_gsxx, 2);
    const double temp5 = 1 / M_PI;
    const double temp6 = pow(Q, -2);
    const double temp7 = pow(m_mx, 2);
    const double temp8 = -4 * temp7;
    const double temp9 = temp2 + temp8;
    const double temp10 = -4 * temp3;
    const double temp11 = temp10 + temp2;
    const double temp12 = temp11 * temp9;
    const double temp13 = sqrt(temp12);
    const double temp14 = -1 + m_smix;
    const double temp15 = pow(m_smix, 6);
    const double temp16 = 1 + m_smix;
    const double temp17 = pow(kHIGGS_VEV, -2);
    const double temp18 = -temp2;
    const double temp19 = pow(m_ms, 4);
    const double temp20 = 2 * temp3;
    const double temp21 = -2 * temp3;
    const double temp22 = pow(m_gsxx, 4);
    const double temp23 = 1 / temp9;
    const double temp24 = pow(m_smix, 2);
    const double temp25 = -1 + temp24;
    const double temp26 = pow(temp25, 2);
    const double temp27 = pow(m_mx, 4);
    const double temp28 = 16 * temp27;
    const double temp29 = 2 * temp19;
    const double temp30 = -12 * temp3;
    const double temp31 = temp2 + temp30;
    const double temp32 = temp31 * temp7;
    const double temp33 = temp28 + temp29 + temp32;
    const double temp34 = temp11 * temp7;
    const double temp35 = temp19 + temp34;
    const double temp36 = 1 / temp35;
    const double temp37 = temp13 * temp33 * temp36;
    const double temp38 = 8 * temp7;
    const double temp39 = temp2 + temp21 + temp38;
    const double temp40 = temp13 + temp18 + temp20;
    const double temp41 = 1 / temp40;
    const double temp42 = temp13 + temp2 + temp21;
    const double temp43 = -(temp41 * temp42);
    const double temp44 = log(temp43);
    const double temp45 = -(temp39 * temp44);
    const double temp46 = temp37 + temp45;
    const double temp47 = -(temp22 * temp23 * temp26 * temp46 * temp5 * temp6) / 64.;
    const double temp48 = temp1 + temp20;
    const double temp49 = temp18 + temp3;
    const double temp50 = pow(temp49, 2);
    const double temp51 = pow(m_width_s, 2);
    const double temp52 = temp3 * temp51;
    const double temp53 = temp50 + temp52;
    const double temp54 = 1 / temp53;
    const double temp55 = -temp1;
    const double temp56 = temp2 + temp55;
    const double temp57 = 1 / temp56;
    const double temp58 = pow(m_gsxx, 3);
    const double temp59 = temp1 + temp18;
    const double temp60 = pow(m_smix, 3);
    const double temp61 = -temp24;
    const double temp62 = 1 + temp61;
    const double temp63 = pow(temp62, 1.5);
    const double temp64 = 1 / kHIGGS_VEV;
    const double temp65 = temp18 + temp21 + temp38;
    const double temp66 = -Q;
    const double temp67 = m_ms + temp66;
    const double temp68 = m_ms + Q;
    const double temp69 = 4 * temp7;
    const double temp70 = temp18 + temp69;
    const double temp71 = 4 * temp3;
    const double temp72 = temp18 + temp71;
    const double temp73 = temp70 * temp72;
    const double temp74 = sqrt(temp73);
    const double temp75 = -8 * temp7;
    const double temp76 = temp2 + temp20 + temp75;
    const double temp77 = -temp13;
    const double temp78 = temp2 + temp21 + temp77;
    const double temp79 = 1 / temp42;
    const double temp80 = temp78 * temp79;
    const double temp81 = log(temp80);
    const double temp82 = temp2 + temp21;

    const double ampList[10] = {
            -(temp13 * temp14 * temp15 * temp16 * temp17 * temp4 * pow(temp48, 2) * temp5 * temp6) /
                    (64. * pow(temp59, 2)),
            (-9 * temp13 * temp14 * temp15 * temp16 * temp17 * temp19 * temp4 * temp5 * temp54 * temp6) / 64.,
            temp47,
            temp47,
            (-3 * temp15 * temp17 * temp25 * temp3 * temp4 * temp48 * temp5 * temp54 * temp57 * temp6 * temp67 *
                    temp68 * sqrt(-(temp11 * temp70))) / 32.,
            (m_mx * temp23 * temp48 * temp5 * temp57 * temp58 * temp6 * temp60 * temp63 * temp64 *
                    (-2 * temp13 + temp44 * temp65)) / 32.,
            (m_mx * temp23 * temp48 * temp5 * temp58 * temp6 * temp60 * temp63 * temp64 *
                    (2 * temp13 + temp65 * temp81)) / (32. * temp59),
            (-3 * m_mx * temp23 * temp3 * temp5 * temp54 * temp58 * temp6 * temp60 * temp63 * temp64 * temp67 * temp68 *
                    (2 * temp74 + temp44 * temp76)) / 32.,
            (3 * m_mx * temp23 * temp3 * temp5 * temp54 * temp58 * temp6 * temp60 * temp63 * temp64 * temp67 * temp68 *
                    (-2 * temp74 + temp76 * temp81)) / 32.,
            -(temp22 * temp23 * temp26 * temp5 * temp6 * (temp13 * temp82 +
                    (temp19 - 16 * temp27 + 4 * temp2 * temp7) * log(pow(temp40, 2) / pow(temp42, 2)))) /
                    (32. * temp82)};

    if (channel == "s") {
        return ampList[0] + ampList[1] + ampList[4];
    } else if (channel == "tu" || channel == "ut") {
        return ampList[2] + ampList[3] + ampList[9];
    } else {
        double sum = 0.0;
        for (int i = 0; i < 10; i++) {
            sum += ampList[i];
        }
        return sum;
    }
}

/*
 * Compute the total dark matter annihilation cross section
 * @param Q Center of mass energy
 * @param state String representing the final state to use.
 * @param channel String representing which channel to include, i.e. "s", "tu" or "all"
 */
double HiggsPortal::annihilation_cross_section(
        const double Q,
        const std::string &state = "all",
        const std::string &channel = "all"
) const {
    if (state == "all") {
        return (sigma_xx_to_ff(Q, channel, "e") +
                sigma_xx_to_ff(Q, channel, "mu") +
                sigma_xx_to_ff(Q, channel, "tau") +
                sigma_xx_to_ff(Q, channel, "u") +
                sigma_xx_to_ff(Q, channel, "c") +
                sigma_xx_to_ff(Q, channel, "t") +
                sigma_xx_to_ff(Q, channel, "d") +
                sigma_xx_to_ff(Q, channel, "s") +
                sigma_xx_to_ff(Q, channel, "b") +
                sigma_xx_to_ww(Q, channel) +
                sigma_xx_to_zz(Q, channel) +
                sigma_xx_to_hh(Q, channel) +
                sigma_xx_to_ss(Q, channel));
    } else if (state == "e" || state == "e e") {
        return sigma_xx_to_ff(Q, channel, "e");
    } else if (state == "mu" || state == "mu mu") {
        return sigma_xx_to_ff(Q, channel, "mu");
    } else if (state == "tau" || state == "tau tau") {
        return sigma_xx_to_ff(Q, channel, "tau");
    } else if (state == "u" || state == "u u") {
        return sigma_xx_to_ff(Q, channel, "u");
    } else if (state == "c" || state == "c c") {
        return sigma_xx_to_ff(Q, channel, "c");
    } else if (state == "t" || state == "t t") {
        return sigma_xx_to_ff(Q, channel, "t");
    } else if (state == "d" || state == "d d") {
        return sigma_xx_to_ff(Q, channel, "d");
    } else if (state == "s" || state == "s s") {
        return sigma_xx_to_ff(Q, channel, "s");
    } else if (state == "b" || state == "b b") {
        return sigma_xx_to_ff(Q, channel, "b");
    } else if (state == "w" || state == "w w") {
        return sigma_xx_to_ww(Q, channel);
    } else if (state == "z" || state == "z z") {
        return sigma_xx_to_zz(Q, channel);
    } else if (state == "h" || state == "h h") {
        return sigma_xx_to_hh(Q, channel);
    } else {
        return 0.0;
    }
}

/*
 * Compute the total dark matter thermal annihilation cross section integrand
 * @param z Center of mass energy divided by DM mass
 * @param x DM mass divided by the DM temperature
 * @param state String representing the final state to use.
 * @param channel String representing which channel to include, i.e. "s", "tu" or "all"
 */
double HiggsPortal::thermal_cross_section_integrand(
        const double z,
        const double x,
        const std::string &state,
        const std::string &channel
) const {
    using namespace special_functions;
    const double z2 = z * z;
    const double sig = annihilation_cross_section(m_mx * z, state, channel);
    const double kernal = z2 * (z2 - 4.0) * besselk1e(x * z) * exp(-x * (z - 2.0));
    return sig * kernal;
}

/*
 * Compute the total dark matter thermal annihilation cross section
 * @param x DM mass divided by the DM temperature
 * @param state String representing the final state to use.
 * @param channel String representing which channel to include, i.e. "s", "tu" or "all"
 */
double HiggsPortal::thermal_cross_section(
        const double x,
        const std::string &state = "all",
        const std::string &channel = "all"
) const {
    using boost::math::quadrature::gauss_kronrod;
    using namespace special_functions;

    const double denom = 2.0 * besselk2e(x);
    const double pf = x / (denom * denom);
    auto integrand = [this, &x, &state, &channel](double z) {
        return thermal_cross_section_integrand(z, x, state, channel);
    };

    return pf * gauss_kronrod<double, 15>::integrate(integrand, 2.0, std::numeric_limits<double>::infinity(), 5, 1e-9);
}

/*
 * Compute the relic density of the dark matter.
 */
double HiggsPortal::relic_density(
        double xstart = 1.0,
        double xend = 500.0,
        double reltol = 1e-3,
        double abstol = 1e-9,
        const std::string &t_alg = "radau"
) {
    using namespace diffeq;
    using namespace cosmology;

    auto logx_span = std::make_pair(std::log(xstart), std::log(xend));

    HiggsPortalBoltzmann boltz{std::make_shared<HiggsPortal>(*this)};
    double Tinit = boltz.chi.get_mass() / exp(logx_span.first);

    Vector<double> winit{1};
    winit(0) = log(boltz.chi.neq(Tinit) / sm_entropy_density(Tinit));

    ODEProblem problem{std::make_shared<HiggsPortalBoltzmann>(boltz), winit, logx_span};

    double yinf;
    if (t_alg == "rodas") {
        Rodas alg{};
        ODEIntegratorOptions opts{};
        opts.abstol = abstol;
        opts.reltol = reltol;

        auto sol = solve(problem, alg, opts);

        yinf = exp(sol.us.back()[0]);
    } else {
        Radau5 alg{};
        ODEIntegratorOptions opts{};
        opts.abstol = abstol;
        opts.reltol = reltol;

        auto sol = solve(problem, alg, opts);

        yinf = exp(sol.us.back()[0]);
    }
    return boltz.chi.get_mass() * yinf * kS_TODAY / kRHO_CRIT;


}

}
}

#endif //LANRE_DM_MODELS_HIGGS_PORTAL_HPP
