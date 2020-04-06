//
// Created by Logan Morrison on 3/15/20.
//

#ifndef LANRE_DM_MODELS_KINETIC_MIXING_HPP
#define LANRE_DM_MODELS_KINETIC_MIXING_HPP

#include "lanre/constants.hpp"
#include "lanre/special_functions/besselk.hpp"
#include "lanre/diffeq/base.hpp"
#include "lanre/diffeq/radau.hpp"
#include "lanre/diffeq/rodas.hpp"
#include "lanre/diffeq/function.hpp"
#include "lanre/cosmology/standard_model.hpp"
#include "lanre/cosmology/thermodynamic_particle.hpp"
#include <string>
#include <cmath>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/math/differentiation/finite_difference.hpp>


namespace lanre {
namespace dm_models {

using namespace diffeq;
using namespace cosmology;

class KineticMixing {
private:
    double m_mx;
    double m_mv;
    double m_gvxx;
    double m_eps;
    double m_width_v;

    ThermodynamicParticle m_chi;

    double width_v_to_ququ(double) const;

    double width_v_to_qdqd(double) const;

    double width_v_to_ll(double) const;

    double width_v_to_nunu() const;

    double width_v_to_hz() const;

    double width_v_to_xx() const;

    double sigma_xx_to_ququ(double, double, const std::string &) const;

    double sigma_xx_to_qdqd(double, double, const std::string &) const;

    double sigma_xx_to_ll(double, double, const std::string &) const;

    double sigma_xx_to_nunu(double, const std::string &) const;

    double sigma_xx_to_hz(double, const std::string &) const;

    double sigma_xx_to_vv(double, const std::string &) const;

    double thermal_cross_section_integrand(
            double, double, const std::string &,
            const std::string &) const;

    double xstar_root_eqn_gondolo_gelmini(double, double) const;

    double compute_xstar_gondolo_gelmini(double) const;

    double compute_alpha_gondolo_gelmini(double) const;

    double relic_density_gondolo_gelmini() const;

    double compute_xf_bender() const;

    double relic_density_bender() const;


public:
    KineticMixing(double mx, double mv, double gvxx, double eps)
            : m_mx(mx), m_mv(mv), m_gvxx(gvxx), m_eps(eps), m_chi{m_mx, 2.0, 1} {
        m_width_v = vector_mediator_width("all");
    }

    ~KineticMixing() = default;

    double get_mx() const;

    double get_mv() const;

    double get_gvxx() const;

    double get_eps() const;

    double get_width_v() const;

    void set_mx(double);

    void set_mv(double);

    void set_gvxx(double);

    void set_eps(double);

    double vector_mediator_width(const std::string &) const;

    double annihilation_cross_section(double,
                                      const std::string &, const std::string &) const;

    double thermal_cross_section(double, const std::string &,
                                 const std::string &) const;

    ODESolution solve_boltzmann(double, double, double, double, const std::string &);

    double relic_density(double, double, double, double, const std::string &);

    double xf_root_equation_mpu(double) const;

    double compute_xf_mpu() const;

    double relic_density_mpu() const;
};

struct KineticMixingBoltzman : public ODEFunction {

    std::shared_ptr<KineticMixing> model;
    cosmology::ThermodynamicParticle chi;

    KineticMixingBoltzman(std::shared_ptr<KineticMixing> t_model)
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

double KineticMixing::get_mx() const {
    return m_mx;
}

double KineticMixing::get_mv() const {
    return m_mv;
}

double KineticMixing::get_gvxx() const {
    return m_gvxx;
}

double KineticMixing::get_eps() const {
    return m_eps;
}

double KineticMixing::get_width_v() const {
    return m_width_v;
}

void KineticMixing::set_mx(double mx) {
    m_mx = mx;
    m_width_v = vector_mediator_width("all");
}

void KineticMixing::set_mv(double mv) {
    m_mv = mv;
    m_width_v = vector_mediator_width("all");
}

void KineticMixing::set_gvxx(double gvxx) {
    m_gvxx = gvxx;
    m_width_v = vector_mediator_width("all");
}

void KineticMixing::set_eps(double eps) {
    m_eps = eps;
    m_width_v = vector_mediator_width("all");
}


/*
 * Vector-mediator widths
 */

double KineticMixing::width_v_to_ququ(double mf) const {
    if (m_mv <= 2.0 * mf) {
        return 0.0;
    }
    return (kALPHA_EM * pow(m_eps, 2) * sqrt(-4 * pow(mf, 2) + pow(m_mv, 2)) * (7 * pow(mf, 2) + 17 * pow(m_mv, 2))) /
            (72. * pow(kCOS_THETA_WEAK, 2) * pow(m_mv, 2));
}

double KineticMixing::width_v_to_qdqd(double mf) const {
    if (m_mv <= 2.0 * mf) {
        return 0.0;
    }
    return (kALPHA_EM * pow(m_eps, 2) * sqrt(-4 * pow(mf, 2) + pow(m_mv, 2)) * (-17 * pow(mf, 2) + 5 * pow(m_mv, 2))) /
            (72. * pow(kCOS_THETA_WEAK, 2) * pow(m_mv, 2));
}

double KineticMixing::width_v_to_ll(double mf) const {
    if (m_mv <= 2.0 * mf) {
        return 0.0;
    }
    return (kALPHA_EM * pow(m_eps, 2) * sqrt(-4 * pow(mf, 2) + pow(m_mv, 2)) * (7 * pow(mf, 2) + 5 * pow(m_mv, 2))) /
            (24. * pow(kCOS_THETA_WEAK, 2) * pow(m_mv, 2));
}

double KineticMixing::width_v_to_nunu() const {
    return (kALPHA_EM * pow(m_eps, 2) * m_mv) / (24. * pow(kCOS_THETA_WEAK, 2));
}

double KineticMixing::width_v_to_hz() const {
    if (m_mv <= kHIGGS_MASS + kZ_BOSON_MASS) {
        return 0.0;
    }
    return (pow(kALPHA_EM, 2) * pow(m_eps, 2) * sqrt(-pow(kHIGGS_MASS, 2) + pow(pow(kHIGGS_MASS, 2) + pow(m_mv, 2) -
                                                                                        pow(kW_BOSON_MASS, 2), 2) /
            (4. * pow(m_mv, 2))) * (2 * pow(m_mv, 2) * pow(kW_BOSON_MASS, 2) +
            pow(-pow(kHIGGS_MASS, 2) + pow(m_mv, 2) + pow(kW_BOSON_MASS, 2), 2) / 4.) * M_PI * pow(kHIGGS_VEV, 2)) /
            (6. * pow(kCOS_THETA_WEAK, 4) * pow(m_mv, 4) * pow(kW_BOSON_MASS, 2) * pow(kSIN_THETA_WEAK, 2));
}

double KineticMixing::width_v_to_xx() const {
    if (m_mv <= 2.0 * m_mx) {
        return 0.0;
    }
    return (pow(m_gvxx, 2) * sqrt(-4 * pow(m_mx, 2) + pow(m_mv, 2)) * (2 * pow(m_mx, 2) + pow(m_mv, 2))) /
            (12. * pow(m_mv, 2) * M_PI);
}

double KineticMixing::vector_mediator_width(const std::string &state = "all") const {
    if (state == "all") {
        return (width_v_to_ququ(kUP_QUARK_MASS) +
                width_v_to_ququ(kCHARM_QUARK_MASS) +
                width_v_to_ququ(kTOP_QUARK_MASS) +
                width_v_to_qdqd(kDOWN_QUARK_MASS) +
                width_v_to_qdqd(kSTRANGE_QUARK_MASS) +
                width_v_to_qdqd(kBOTTOM_QUARK_MASS) +
                width_v_to_ll(kELECTRON_MASS) +
                width_v_to_ll(kMUON_MASS) +
                width_v_to_ll(kTAU_MASS) +
                3 * width_v_to_nunu() +
                width_v_to_hz() +
                width_v_to_xx());
    } else if (state == "u u" || state == "u") {
        return width_v_to_ququ(kUP_QUARK_MASS);
    } else if (state == "c c" || state == "c") {
        return width_v_to_ququ(kCHARM_QUARK_MASS);
    } else if (state == "t t" || state == "t") {
        return width_v_to_ququ(kTOP_QUARK_MASS);
    } else if (state == "d d" || state == "d") {
        return width_v_to_qdqd(kDOWN_QUARK_MASS);
    } else if (state == "s s" || state == "s") {
        return width_v_to_qdqd(kSTRANGE_QUARK_MASS);
    } else if (state == "b b" || state == "b") {
        return width_v_to_qdqd(kBOTTOM_QUARK_MASS);
    } else if (state == "e e" || state == "e") {
        return width_v_to_ll(kELECTRON_MASS);
    } else if (state == "mu mu" || state == "mu") {
        return width_v_to_ll(kMUON_MASS);
    } else if (state == "tau tau" || state == "tau") {
        return width_v_to_ll(kTAU_MASS);
    } else if (state == "nue nue" || state == "nue" ||
            state == "num num" || state == "num" ||
            state == "nut nut" || state == "nut") {
        return width_v_to_nunu();
    } else if (state == "x x" || state == "x") {
        return width_v_to_xx();
    } else if (state == "h z" || state == "z h") {
        return width_v_to_hz();
    }
}


/*
 * Dark matter annihilation cross sections
 */


double KineticMixing::sigma_xx_to_ququ(double Q, double mf, const std::string &channel) const {
    if ((channel != "s" && channel != "all") || Q <= 2.0 * mf || Q <= 2.0 * m_mx) {
        return 0.0;
    }
    double temp1 = pow(m_mx, 2);
    double temp2 = pow(Q, 2);
    double temp3 = pow(mf, 2);
    double temp4 = pow(m_mv, 2);
    return ((kALPHA_EM * pow(m_eps, 2) * pow(m_gvxx, 2) * (2 * temp1 +
            temp2) * sqrt(temp2 - 4 * temp3) * (17 * temp2 +
            7 * temp3)) / (72. * pow(kCOS_THETA_WEAK, 2) *
            pow(Q, 2) * sqrt(-4 * temp1 + temp2) *
            (pow(-temp2 + temp4, 2) + temp4 * pow(m_width_v, 2))));
}

double KineticMixing::sigma_xx_to_qdqd(double Q, double mf, const std::string &channel) const {
    if ((channel != "s" && channel != "all") || Q <= 2.0 * mf || Q <= 2.0 * m_mx) {
        return 0.0;
    }
    double temp1 = pow(Q, 2);
    double temp2 = pow(m_mx, 2);
    double temp3 = pow(mf, 2);
    double temp4 = pow(m_mv, 2);
    return ((kALPHA_EM * pow(m_eps, 2) * pow(m_gvxx, 2) * (temp1 +
            2 * temp2) * (5 * temp1 - 17 * temp3) * sqrt(temp1 - 4 * temp3)) /
            (72. * pow(kCOS_THETA_WEAK, 2) * pow(Q, 2) * sqrt(temp1 - 4 * temp2) *
                    (pow(-temp1 + temp4, 2) + temp4 * pow(m_width_v, 2))));
}

double KineticMixing::sigma_xx_to_ll(double Q, double mf, const std::string &channel) const {
    if ((channel != "s" && channel != "all") || Q <= 2.0 * mf || Q <= 2.0 * m_mx) {
        return 0.0;
    }
    double temp1 = pow(m_mx, 2);
    double temp2 = pow(Q, 2);
    double temp3 = pow(mf, 2);
    double temp4 = pow(m_mv, 2);
    return (kALPHA_EM * pow(m_eps, 2) * pow(m_gvxx, 2) * (2 * temp1 +
            temp2) * sqrt(temp2 - 4 * temp3) * (5 * temp2 + 7 * temp3)) /
            (24. * pow(kCOS_THETA_WEAK, 2) * pow(Q, 2) * sqrt(-4 * temp1 + temp2) *
                    (pow(-temp2 + temp4, 2) + temp4 * pow(m_width_v, 2)));
}

double KineticMixing::sigma_xx_to_nunu(double Q, const std::string &channel) const {
    if (channel != "s" && channel != "all") {
        return 0.0;
    }
    double temp1 = pow(Q, 2);
    double temp2 = pow(m_mx, 2);
    double temp3 = pow(m_mv, 2);
    return (kALPHA_EM * pow(m_eps, 2) * pow(m_gvxx, 2) * sqrt(temp1) * (temp1 +
            2 * temp2)) / (24. * pow(kCOS_THETA_WEAK, 2) * sqrt(temp1 - 4 * temp2) *
            (pow(-temp1 + temp3, 2) + temp3 * pow(m_width_v, 2)));
}

double KineticMixing::sigma_xx_to_hz(double Q, const std::string &channel) const {
    if ((channel != "s" && channel != "all") || Q <= kHIGGS_MASS + kW_BOSON_MASS || Q <= 2.0 * m_mx) {
        return 0.0;
    }
    double temp1 = -Q;
    double temp2 = -kW_BOSON_MASS;
    double temp3 = pow(m_mx, 2);
    double temp4 = pow(Q, 2);
    double temp5 = pow(kW_BOSON_MASS, 2);
    double temp6 = pow(m_mv, 2);
    return (pow(kALPHA_EM, 2) * pow(m_eps, 2) * pow(m_gvxx, 2) * M_PI * sqrt(((
            kHIGGS_MASS + kW_BOSON_MASS + Q) * (kHIGGS_MASS + kW_BOSON_MASS +
            temp1) * (kHIGGS_MASS + Q + temp2) * (kHIGGS_MASS + temp1 +
            temp2)) / (-4 * temp3 + temp4)) * (2 * temp3 + temp4) * (pow(kHIGGS_MASS, 4) +
            pow(kW_BOSON_MASS, 4) + pow(Q, 4) + 10 * temp4 * temp5 -
            2 * pow(kHIGGS_MASS, 2) * (temp4 + temp5)) * pow(kHIGGS_VEV, 2)) /
            (48. * pow(kCOS_THETA_WEAK, 4) * pow(kW_BOSON_MASS, 2) * pow(Q, 5) * pow
                    (kSIN_THETA_WEAK, 2) * (pow(-temp4 + temp6, 2) +
                    temp6 * pow(m_width_v, 2)));
}

double KineticMixing::sigma_xx_to_vv(double Q, const std::string &channel) const {
    if (channel == "s" || Q <= 2.0 * m_mx || Q <= 2.0 * m_mv) {
        return 0.0;
    }
    double temp1 = pow(m_mx, 2);
    double temp2 = -4 * temp1;
    double temp3 = pow(Q, 2);
    double temp4 = temp2 + temp3;
    double temp5 = pow(m_mv, 4);
    double temp6 = pow(m_mv, 2);
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
    return (pow(m_gvxx, 4) * ((-48 * temp13 * temp14 * (4 * pow(m_mx, 4) +
            temp1 * temp3 + 2 * temp5)) / (temp5 + temp1 * temp8) + (48 * (temp1 * (temp12 +
            temp2 + temp3) * log(2) + temp3 * temp6 * log(4) + (2 * temp1 * (temp16 +
            temp6) - temp3 * (temp1 + temp9)) * log(-(pow(temp17, 2) / (-pow(Q, 4) +
            temp13 * temp14 * temp3 - 2 * temp5 + (-2 * temp13 * temp14 + 4 * temp3) * temp6 +
            2 * temp1 * temp8))) + temp11 * (temp12 + temp16 +
            temp3) * log(-(temp17 / (temp10 + temp15 +
            temp9))))) / temp11)) / (384. * M_PI * pow(Q, 2) * temp4);
}

double KineticMixing::annihilation_cross_section(
        double Q,
        const std::string &state,
        const std::string &channel
) const {
    if (state == "all") {
        return (sigma_xx_to_ququ(Q, kUP_QUARK_MASS, channel) +
                sigma_xx_to_ququ(Q, kCHARM_QUARK_MASS, channel) +
                sigma_xx_to_ququ(Q, kTOP_QUARK_MASS, channel) +
                sigma_xx_to_qdqd(Q, kDOWN_QUARK_MASS, channel) +
                sigma_xx_to_qdqd(Q, kSTRANGE_QUARK_MASS, channel) +
                sigma_xx_to_qdqd(Q, kBOTTOM_QUARK_MASS, channel) +
                sigma_xx_to_ll(Q, kELECTRON_MASS, channel) +
                sigma_xx_to_ll(Q, kMUON_MASS, channel) +
                sigma_xx_to_ll(Q, kTAU_MASS, channel) +
                3 * sigma_xx_to_nunu(Q, channel) +
                sigma_xx_to_hz(Q, channel) +
                sigma_xx_to_vv(Q, channel));
    } else if (state == "u u" || state == "u") {
        return sigma_xx_to_ququ(Q, kUP_QUARK_MASS, channel);
    } else if (state == "c c" || state == "c") {
        return sigma_xx_to_ququ(Q, kCHARM_QUARK_MASS, channel);
    } else if (state == "t t" || state == "t") {
        return sigma_xx_to_ququ(Q, kTOP_QUARK_MASS, channel);
    } else if (state == "d d" || state == "d") {
        return sigma_xx_to_qdqd(Q, kDOWN_QUARK_MASS, channel);
    } else if (state == "s s" || state == "s") {
        return sigma_xx_to_qdqd(Q, kSTRANGE_QUARK_MASS, channel);
    } else if (state == "b b" || state == "b") {
        return sigma_xx_to_qdqd(Q, kBOTTOM_QUARK_MASS, channel);
    } else if (state == "e e" || state == "e") {
        return sigma_xx_to_ll(Q, kELECTRON_MASS, channel);
    } else if (state == "mu mu" || state == "mu") {
        return sigma_xx_to_ll(Q, kMUON_MASS, channel);
    } else if (state == "tau tau" || state == "tau") {
        return sigma_xx_to_ll(Q, kTAU_MASS, channel);
    } else if (state == "nue nue" || state == "nue" ||
            state == "num num" || state == "num" ||
            state == "nut nut" || state == "nut") {
        return sigma_xx_to_nunu(Q, channel);
    } else if (state == "h z" || state == "z h") {
        return sigma_xx_to_hz(Q, channel);
    } else if (state == "v v" || state == "v") {
        return sigma_xx_to_vv(Q, channel);
    }
    return 0.0;
}

double KineticMixing::thermal_cross_section_integrand(
        const double z,
        const double x,
        const std::string &state,
        const std::string &channel
) const {
    using special_functions::besselk1e;

    const double z2 = z * z;
    const double sig = annihilation_cross_section(m_mx * z, state, channel);
    const double kernal = z2 * (z2 - 4.0) * besselk1e(x * z) * exp(-x * (z - 2.0));
    return sig * kernal;
}


double KineticMixing::thermal_cross_section(
        const double x,
        const std::string &state = "all",
        const std::string &channel = "all"
) const {
    using boost::math::quadrature::gauss_kronrod;
    using special_functions::besselk2e;

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
ODESolution KineticMixing::solve_boltzmann(
        double xstart = 1.0,
        double xend = 500.0,
        double reltol = 1e-3,
        double abstol = 1e-9,
        const std::string &t_alg = "radau"
) {
    using namespace diffeq;
    using namespace cosmology;

    auto logx_span = std::make_pair(std::log(xstart), std::log(xend));

    KineticMixingBoltzman boltz{std::make_shared<KineticMixing>(*this)};
    double Tinit = boltz.chi.get_mass() / exp(logx_span.first);

    diffeq::Vector<double> winit{1};
    winit(0) = log(boltz.chi.neq(Tinit) / sm_entropy_density(Tinit));

    ODEProblem problem{std::make_shared<KineticMixingBoltzman>(boltz), winit, logx_span};

    if (t_alg == "rodas") {
        Rodas alg{};
        ODEIntegratorOptions opts{};
        opts.abstol = abstol;
        opts.reltol = reltol;

        return solve(problem, alg, opts);
    } else {
        Radau5 alg{};
        ODEIntegratorOptions opts{};
        opts.abstol = abstol;
        opts.reltol = reltol;

        return solve(problem, alg, opts);
    }
}

/*
 * Compute the relic density of the dark matter.
 */

/**
 * Compute the residual equation for x_f using Morrison, Patel, Ulbricht method.
 * @param xf Freeze-out value for x = mass / T.
 * @return residual
 */
double KineticMixing::xf_root_equation_mpu(double xf) const {
    using boost::math::differentiation::finite_difference_derivative;

    const double euler_gamma = 0.5772156649015328606065120;
    const double Tf = m_mx / xf;
    const double pf = sqrt(M_PI / 45.0) * m_mx * kM_PLANK;
    const double fx = (pf * sm_sqrt_gstar(m_mx / xf) / (xf * xf) * thermal_cross_section(xf));
    const double yeq = m_chi.neq(Tf) / sm_entropy_density(Tf);
    const double Qx = fx * yeq;

    // Compute f'(xf) / f(xf) = d/dx log(f(x))
    auto logf = [this](double x) {
        // Ignoring constant pf since it won't matter for derivative
        return log(sm_sqrt_gstar(m_mx / x) / (x * x) * thermal_cross_section(x));
    };
    const double dlogfx = finite_difference_derivative(logf, xf);

    return Qx - (2.0 * exp(-(euler_gamma + dlogfx + 3.0 / (2.0 * xf))));
}

/**
 * Compute the freeze-out value of x = mass / T using the Morrison, Patel, Ulbricht
 * method.
 * @return xf
 */
double KineticMixing::compute_xf_mpu() const {
    using namespace boost::math::tools;
    auto f = [this](double xf) {
        return xf_root_equation_mpu(xf);
    };
    auto tol = [](double min, double max) { return abs(max - min) <= 1e-8; };
    std::pair<double, double> res = bisect(f, 1e-2, 100.0, tol);
    return (res.second + res.first) / 2.0;
}

/**
 * Compute the relic density using the Morrison, Patel, Ulbricht method
 * @return
 */
double KineticMixing::relic_density_mpu() const {
    using boost::math::quadrature::gauss_kronrod;

    const double xf = compute_xf_mpu();
    const double xinf = std::min(100.0 * xf, 350.0);
    const double pf = sqrt(M_PI / 45.0) * m_mx * kM_PLANK;

    auto f = [this, pf](double x) {
        return pf * sm_sqrt_gstar(m_mx / x) / (x * x) * thermal_cross_section(x);
    };

    const double integal_f = gauss_kronrod<double, 15>::integrate(f, xf, std::numeric_limits<double>::infinity());
    const double Y0 = 1 / integal_f;

    return Y0 * m_mx * kS_TODAY / kRHO_CRIT;
}

/**
 * Returns residual of root equation used to solve for x_star.
 *
 * See Eqn.(14) of arXiv:1204.3622v3 for similar expressions. Note that our
 * result is more exact since we do not assume `xstar` is large enough that
 * Yeq ~ x^{3/2} e^{-x} / h_sm. This may cause a bit of a slow down.
 *
 * @param xstar Current value of x_star, the value of mass / temperature at
 *              which DM begins to freeze out.
 * @param delta Value of `delta` assumed for when DM begins to freeze out.
 *              Default value is the solution to
 *                  delta * (2 + delta) / (1 + delta) = 1,
 *              i.e., delta = (sqrt(5) - 1) / 2 = 0.618033988749895.
 *              See Eqn.(13) of arXiv:1204.3622v3 for details and other used
 *              values of delta. Value of xstar is logarithmically sensitive to
 *              this number.
 * @return res Residual of the root equation.
 */
double KineticMixing::xstar_root_eqn_gondolo_gelmini(double xstar, double delta = 0.0) const {
    double deltabar = delta <= 0.0 ? 1.0 : delta * (2.0 + delta) / (1.0 + delta);
    double T = m_mx / xstar;
    double lam = sqrt(M_PI / 45.0) * m_mx * kM_PLANK * sm_sqrt_gstar(T);
    double tcs = thermal_cross_section(xstar);

    double s = sm_entropy_density(T);
    double ds = sm_entropy_density_deriv(T);
    double neq = m_chi.neq(T);
    double dneq = m_chi.neq_deriv(T);
    double yeq = neq / s;
    // This is dY/dT = -x^2/m dY/xd
    double dyeq = (s * dneq - neq * ds) / (s * s);
    dyeq *= -m_mx / (xstar * xstar);

    return xstar * xstar * dyeq + lam * deltabar * tcs * yeq * yeq;
}

/**
 * Computes to value of `xstar`: the value of dm_mass / temperature such that
 * the DM begins to freeze out.
 * @param delta Value of `delta` assumed for when DM begins to freeze out.
 *              Default value is the solution to
 *                  delta * (2 + delta) / (1 + delta) = 1,
 *              i.e., delta = (sqrt(5) - 1) / 2 = 0.618033988749895.
 *              See Eqn.(13) of arXiv:1204.3622v3 for details and other used
 *              values of delta. Value of xstar is logarithmically sensitive to
 *              this number.
 * @return xstar Value of mass / temperature at which DM begins to freeze-out.
 */
double KineticMixing::compute_xstar_gondolo_gelmini(double delta = 0.0) const {
    using namespace boost::math::tools;
    auto f = [this, delta](double xstar) {
        return xstar_root_eqn_gondolo_gelmini(xstar, delta);
    };
    auto tol = [](double min, double max) { return abs(max - min) <= 1e-8; };
    std::pair<double, double> res = bisect(f, 0.01, 100.0, tol);
    return (res.second + res.first) / 2.0;
}

/**
 * Computes the value of the integral of RHS of the Boltzmann equation with
 * Yeq set to zero from x_star to x_fo.
 *
 * @param xstar Value of mass / temperature at which DM begins to freeze-out.
 *              See `compute_xstar` for more details.
 * @return
 */
double KineticMixing::compute_alpha_gondolo_gelmini(double xstar) const {
    using boost::math::quadrature::gauss_kronrod;
    double pf = sqrt(M_PI / 45.0) * m_mx * kM_PLANK;

    auto integrand = [this](double x) {
        return sm_sqrt_gstar(m_mx / x) * thermal_cross_section(x) / (x * x);
    };

    return pf * gauss_kronrod<double, 15>::integrate(integrand, xstar, 100 * xstar);
}

/**
 * Compute the relic density using method by Gondolo + Gelmini.
 * @return rd Relic density.
 */
double KineticMixing::relic_density_gondolo_gelmini() const {
    double xstar = compute_xstar_gondolo_gelmini();
    double alpha = compute_alpha_gondolo_gelmini(xstar);
    double Tstar = m_mx / xstar;
    double ystar = m_chi.neq(Tstar) / sm_entropy_density(Tstar);
    double Y0 = ystar / (1.0 + ystar * alpha);

    return Y0 * m_mx * kS_TODAY / kRHO_CRIT;
}

double KineticMixing::compute_xf_bender() const {
    double A = sqrt(M_PI / 2.0) * 45.0 / (4.0 * pow(M_PI, 4));
    double lam = sqrt(M_PI / 45.0) * m_mx * kM_PLANK;
}

double KineticMixing::relic_density_bender() const {

}


double KineticMixing::relic_density(
        double xstart = 1.0,
        double xend = 500.0,
        double reltol = 1e-3,
        double abstol = 1e-9,
        const std::string &t_alg = "radau"
) {
    using namespace diffeq;
    using namespace cosmology;

    auto logx_span = std::make_pair(std::log(xstart), std::log(xend));

    KineticMixingBoltzman boltz{std::make_shared<KineticMixing>(*this)};
    double Tinit = boltz.chi.get_mass() / exp(logx_span.first);

    diffeq::Vector<double> winit{1};
    winit(0) = log(boltz.chi.neq(Tinit) / sm_entropy_density(Tinit));

    ODEProblem problem{std::make_shared<KineticMixingBoltzman>(boltz), winit, logx_span};

    double yinf;
    if (t_alg == "rodas") {
        Rodas alg{};
        ODEIntegratorOptions opts{};
        opts.abstol = abstol;
        opts.reltol = reltol;

        auto sol = solve(problem, alg, opts);

        yinf = exp(sol.us.back()[0]);
    } else if (t_alg == "gg") {
        return relic_density_gondolo_gelmini();
    } else if (t_alg == "mpu") {
        return relic_density_mpu();
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

#endif //LANRE_DM_MODELS_KINETIC_MIXING_HPP
