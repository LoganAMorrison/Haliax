//
// Created by Logan Morrison on 3/29/20.
//

#ifndef LANRE_DM_MODELS_KINETIC_RECOUPLING_HPP
#define LANRE_DM_MODELS_KINETIC_RECOUPLING_HPP

#include "lanre/constants.hpp"
#include "lanre/cosmology/standard_model.hpp"
#include "lanre/cosmology/thermodynamic_particle.hpp"
#include "lanre/special_functions/besselk.hpp"
#include "lanre/diffeq/function.hpp"
#include "lanre/diffeq/problem.hpp"
#include "lanre/diffeq/radau.hpp"
#include "lanre/diffeq/rodas.hpp"
#include "lanre/integrate/qagi.hpp"
#include "lanre/integrate/qagp.hpp"
#include "lanre/interpolate/univariate_spline.hpp"
#include <string>
#include <cmath>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/special_functions/pow.hpp>
#include <utility>

namespace lanre {
namespace dm_models {

using namespace diffeq;
using namespace cosmology;

const std::vector<double> gausslaguerre_nodes = {
        0.01438614699541903, 0.0758036120233568, 0.18631410205718493, 0.34596918099142754, 0.5548109375809135,
        0.8128912841156709, 1.1202738350075414, 1.4770343299238304, 1.8832608263423967, 2.339053849646035,
        2.84452654275536, 3.3998048274457138, 4.005027581758653, 4.660346835568911, 5.365927985585115,
        6.121950030804021, 6.9286058293761705, 7.786102377862518, 8.694661113922168, 9.65451824355508,
        10.665925094121674, 11.729148494472227, 12.844471183641032, 14.012192249694277, 15.2326276004667,
        16.50611046808199, 17.832991949326388, 19.213641584136067, 20.64844797466835, 22.1378194476567,
        23.68218476300238, 25.281993871834043, 26.937718727574264, 28.64985415389129, 30.418918773790942,
        32.24545600452067, 34.130035123421514, 36.07325241037997, 38.075732373107094, 40.138129062115546,
        42.2611274829848, 44.4454451143118, 46.69183354065154, 49.001080210772436, 51.374010332703996,
        53.81148891835566, 56.314422991961706, 58.88376397828209, 61.520510288396146, 64.22571012310156,
        67.00046451641931, 69.84593064455838, 72.76332542897458, 75.75392946593993, 78.81909131941147,
        81.96023221906013, 85.1788512112189, 88.47653081739462, 91.85494326304931, 95.31585734883173,
        98.86114604761353,
        102.49279492391656, 106.21291148804681, 110.02373561603092, 113.92765118897162, 117.92719913257322,
        122.02509207044162, 126.2242308447504, 130.5277232067994, 134.9389050402274, 139.46136455424016,
        144.09896997721273, 148.85590139775826, 153.73668754797302, 158.7462485117131, 163.88994558258725,
        169.17363981000304, 174.60376118237662, 180.1873909402457, 185.93236023966696, 191.84736937224832,
        197.94213310214326, 204.2275595670305, 210.71597286157694, 217.4213932720015, 224.3598947888746,
        231.5500680251725, 239.01362975131494, 246.7762409672485, 254.86862925704742, 263.3281684691579,
        272.20117002409256, 281.54632828389737, 291.4401336163771, 301.9858552516392, 313.3295340040755,
        325.6912634370265, 339.4351019234496, 355.2613118885341, 374.9841128343427
};

const std::vector<double> gausslaguerre_weights = {
        0.03639260588336492, 0.07967674621266874, 0.11211510334236749, 0.1303566129751494, 0.13404333972844779,
        0.1254070907806824, 0.10831411209725411, 0.08709663846998496, 0.0655510093123086, 0.04634013358264724,
        0.03084630862767996, 0.019367828113979857, 0.011485442360179551, 0.006438951001610651, 0.003414979989692477,
        0.0017143197401822045, 0.0008148715915878162, 0.00036685483659949086, 0.00015645207417810467,
        6.321087052888555e-5, 2.41957522651884e-5, 8.774309763755412e-6, 3.014267486000932e-6, 9.808335899345412e-7,
        3.0226387435322415e-7, 8.820058395295924e-8, 2.4364258562006706e-8, 6.369711373901789e-9, 1.5756003204596974e-9,
        3.686329201346082e-10, 8.154798924224659e-11, 1.705062556826075e-11, 3.368214170866727e-12,
        6.283524955366568e-13, 1.1064980159833095e-13, 1.8383501754549047e-14, 2.8801150572305977e-15,
        4.2526128973372586e-16, 5.914432467252218e-17, 7.74308910250013e-18, 9.53624944907476e-19,
        1.1040945169869212e-19, 1.2008469704304912e-20, 1.226007007490482e-21, 1.1740217146467282e-22,
        1.0535940646087198e-23, 8.853231768474294e-25, 6.95919877548344e-26, 5.112369815465843e-27,
        3.506272148817131e-28, 2.2426462726939212e-29, 1.3362094234457466e-30, 7.407428957575754e-32,
        3.815860137136348e-33, 1.8242001976828574e-34, 8.081632336064125e-36, 3.3130638018623833e-37,
        1.254832997687672e-38, 4.3837909801617215e-40, 1.410139229217339e-41, 4.168863999668573e-43,
        1.13048190447714e-44, 2.8060374967870913e-46, 6.361270570982087e-48, 1.3139811959781107e-49,
        2.4668106921703373e-51, 4.197747022871623e-53, 6.456269104911767e-55, 8.947337248621289e-57,
        1.113566974699702e-58, 1.2402350422551838e-60, 1.2313765726081835e-62, 1.0853674489049016e-64,
        8.454956490952747e-67, 5.792630756606999e-69, 3.4718547763480512e-71, 1.8098595335601085e-73,
        8.153752608410954e-76, 3.1524725487932437e-78, 1.0379005899152711e-80, 2.884875523307375e-83,
        6.704801225968573e-86, 1.2889637464360883e-88, 2.0248483451529067e-91, 2.5633922199660193e-94,
        2.573961507515936e-97, 2.0126547874800456e-100, 1.1994844193823018e-103, 5.3121327315397424e-107,
        1.695969259446738e-110, 3.76209729863441e-114, 5.539641754449383e-118, 5.110640477091653e-122,
        2.739965469400114e-126, 7.713611492631752e-131, 9.882494600958801e-136, 4.646863007294177e-141,
        5.62603729501982e-147, 8.90503140588921e-154, 3.2465651634358123e-162
};

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

    double sigma_x2_x2_to_g_g(double) const;

    double sigma_x1_x2_to_w_w(double) const;

    double sigma_x1_x2_to_f_f(double, double) const;

    double weff(double) const;

    double thermal_cross_section_integrand(double, double) const;


public:
    ThermodynamicParticle x1, x2;

    DipoleDM(double m1, double m2, double ce, double cm, double lam)
            : m_m1(m1), m_m2(m2), m_ce(ce), m_cm(cm), m_lam(lam), m_width2(0.0),
              x1{m1, 2.0, 1}, x2{m2, 2.0, 1} {
        m_width2 = compute_width2();
    }

    double get_m1() const { return m_m1; }

    double get_m2() const { return m_m2; }

    double get_ce() const { return m_ce; }

    double get_cm() const { return m_cm; }

    double get_lam() const { return m_lam; }

    double get_width2() const { return m_width2; }

    void set_m1(double m1) {
        m_m1 = m1;
        x1.set_mass(m1);
        m_width2 = compute_width2();
    }

    void set_m2(double m2) {
        m_m2 = m2;
        x2.set_mass(m2);
        m_width2 = compute_width2();
    }

    void set_ce(double ce) {
        m_ce = ce;
        m_width2 = compute_width2();
    }

    void set_cm(double cm) {
        m_cm = cm;
        m_width2 = compute_width2();
    }

    void set_lam(double lam) {
        m_lam = lam;
        m_width2 = compute_width2();
    }

    double gamma_integrand(double, double) const;

    double gamma(double) const;

    double average_p6_E3(double T, double y) const;

    double annihilation_cross_section(double, int i, int j) const;

    double thermal_cross_section(double x) const;

    double thermal_cross_section2(double x) const;

    ODESolution solve_boltzmann(double, double, double, double, const std::string &) const;

    double relic_density(double, double, double, double, const std::string &) const;

    ODESolution solve_temperature(double, double, double, double, const std::string &);

};

/**
 * Struct for solving the Boltzmann equation for the DipoleDM model.
 */
struct DipoleDMBoltzmann : public diffeq::ODEFunction {
    const std::shared_ptr<DipoleDM> model;
    const ThermodynamicParticle x1;
    const ThermodynamicParticle x2;

    explicit DipoleDMBoltzmann(std::shared_ptr<DipoleDM> mod)
            : model(std::move(mod)),
              x1{model->get_m1(), 2.0, 1},
              x2{model->get_m2(), 2.0, 1} {
    }

    void dudt(Vector<double> &dw, const Vector<double> &w, double logx) override {
        double x = exp(logx);
        double T = model->get_m1() / x;

        double sigeff = model->thermal_cross_section(x);
        double pf = -sqrt(M_PI / 45.0) * T * kM_PLANK * sm_sqrt_gstar(T) * sigeff;
        double ww = w(0);
        double weq = log((x1.neq(T) + x2.neq(T)) / sm_entropy_density(T));
        double rhs = pf * (exp(ww) - exp(2 * weq - ww));
        //std::cout << "x, T, rhs = " << x << ", " << T << ", " << rhs << std::endl;
        dw(0) = pf * (exp(ww) - exp(2 * weq - ww));
    }

    void dfdu(Matrix<double> &df, const Vector<double> &w, double logx) override {
        double x = exp(logx);
        double T = model->get_m1() / x;

        double sigeff = model->thermal_cross_section(x);
        double pf = -sqrt(M_PI / 45.0) * T * kM_PLANK * sm_sqrt_gstar(T) * sigeff;

        double ww = w(0);
        double weq = log((x1.neq(T) + x2.neq(T)) / sm_entropy_density(T));

        df(0, 0) = pf * (exp(ww) + exp(2 * weq - ww));
    }
};

struct DipoleDMBoltzmannTemp : public diffeq::ODEFunction {
    const std::shared_ptr<DipoleDM> model;

    explicit DipoleDMBoltzmannTemp(std::shared_ptr<DipoleDM> mod) : model(std::move(mod)) {}

    void dudt(Vector<double> &dy, const Vector<double> &y, double logx) override {
        double x = exp(logx);
        double T = model->get_m1() / x;
        double s = sm_entropy_density(T);
        double yy = y(0);

        double Tx = yy * pow(s, 2.0 / 3.0) / model->get_m1();
        double gam = model->gamma(Tx);
        double H = sqrt(4.0 * pow(M_PI, 3) * sm_geff(T) / 45.0) * T * T / kM_PLANK;
        double gt = sm_sqrt_gstar(T) * sqrt(sm_geff(T)) / sm_heff(T) - 1.0;
        double Ht = H / (1.0 + gt);

        double ddy = yy * (gam / (x * Ht) * (T / Tx - 1.0) +
                H / (3.0 * Tx * x * Ht) * model->average_p6_E3(T, yy)
        );

        dy(0) = ddy;
    }
};


/**
 * Compute the decay width of the heavier DM particle into the lighter and a
 * photon.
 * @return width
 */
double DipoleDM::compute_width2() const {
    double ce2 = m_ce * m_ce;
    double cm2 = m_cm * m_cm;
    double m22 = m_m2 * m_m2;
    double m12 = m_m1 * m_m1;
    double m23 = m22 * m_m2;
    double lam2 = m_lam * m_lam;

    return (ce2 + cm2) * pow(m22 - m12, 3) / (8.0 * M_PI * lam2 * m23);
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
    double temp1 = pow(m_m1, 2);
    double temp2 = pow(Q, 2);
    double temp3 = -2 * temp1;
    double temp4 = pow(m_m2, 2);
    double temp5 = 2 * temp4;
    double temp6 = temp2 + temp3 + temp5;
    double temp7 = -4 * temp1;
    double temp8 = temp2 + temp7;
    double temp9 = -Q;
    double temp10 = sqrt(temp8);
    double temp11 = temp10 + temp9;
    double temp12 = pow(m_m1, 4);
    double temp13 = 36 * temp12;
    double temp14 = pow(m_m2, 4);
    double temp15 = 36 * temp14;
    double temp16 = 9 * temp4;
    double temp17 = 2 * temp2;
    double temp18 = temp16 + temp17;
    double temp19 = -8 * temp1 * temp18;
    double temp20 = pow(Q, 3);
    double temp21 = Q + temp10;
    double temp22 = pow(m_width2, 2);
    double temp23 = temp2 + temp22;
    double temp24 = -12 * temp23 * temp4;
    double temp25 = pow(m_m1, 8);
    double temp26 = pow(m_m1, 6);
    double temp27 = -4 * temp26 * temp4;
    double temp28 = pow(m_m2, 6);
    double temp29 = -6 * temp22 * temp4;
    double temp30 = pow(m_width2, 4);
    double temp31 = temp14 + temp29 + temp30;
    double temp32 = temp14 * temp31;
    double temp33 = -temp10;
    double temp34 = Q + temp33;
    double temp35 = Q * temp34;
    double temp36 = temp3 + temp35 + temp5;
    double temp37 = 1 / temp36;
    double temp38 = 2 * m_m2 * temp37 * m_width2;
    double temp39 = atan(temp38);
    double temp40 = pow(Q, 4);
    double temp41 = 3 * temp14;
    double temp42 = -7 * temp22 * temp4;
    double temp43 = temp41 + temp42;
    double temp44 = 2 * temp12 * temp43;
    double temp45 = -5 * temp14 * temp22;
    double temp46 = temp28 + temp45;
    double temp47 = -4 * temp1 * temp46;
    double temp48 = temp25 + temp27 + temp32 + temp44 + temp47;
    double temp49 = Q * temp21;
    double temp50 = temp3 + temp49 + temp5;
    double temp51 = 1 / temp50;
    double temp52 = 2 * m_m2 * temp51 * m_width2;
    double temp53 = atan(temp52);
    double temp54 = -temp4;
    double temp55 = temp1 + temp54;
    double temp56 = -(temp22 * temp4);
    double temp57 = temp14 + temp56;
    double temp58 = 6 * temp12 * temp57;
    double temp59 = -3 * temp14 * temp22;
    double temp60 = temp28 + temp59;
    double temp61 = -4 * temp1 * temp60;
    double temp62 = temp25 + temp27 + temp32 + temp58 + temp61;
    double temp63 = 2 * temp1;
    double temp64 = -2 * temp1 * temp4;
    double temp65 = temp12 + temp14 + temp56 + temp64;
    double temp66 = -2 * temp4;
    double temp67 = Q * temp10;
    return (pow(pow(m_ce, 2) + pow(m_cm, 2), 2) * (-12 * temp2 * temp39 * temp48 +
            48 * temp1 * temp22 * temp4 * temp40 * temp53 + 12 * temp2 * temp48 * temp53 +
            24 * temp39 * temp55 * temp62 - 24 * temp53 * temp55 * temp62 -
            (m_m2 * Q * temp11 * (temp13 + temp15 + temp19 + temp11 * temp20 +
                    temp24) * temp6 * m_width2) / 2. - (m_m2 * Q * temp21 * (temp13 + temp15 + temp19
            - temp20 * temp21 + temp24) * temp6 * m_width2) / 2. +
            48 * temp1 * temp22 * temp4 * temp40 * atan((2 * m_m2 * m_width2) / (-temp2 +
                    temp63 + temp66 + temp67)) + 6 * m_m2 * ((temp12 - 4 * temp1 * temp4 -
            temp4 * (temp22 + temp4)) * temp40 + 8 * pow(temp55, 2) * temp65 -
            4 * temp2 * (temp54 + temp63) * temp65) * m_width2 * log((temp12 + temp14 +
            (temp20 * temp21) / 2. - temp1 * (Q * (2 * Q + temp10) + temp5) + temp4 * (temp2
            + temp22 + temp67)) / (temp12 + temp14 + (temp20 * temp34) / 2. +
            (-(Q * temp10) + temp2 + temp22) * temp4 + temp1 * (Q * (-2 * Q + temp10) +
            temp66))))) / (96. * pow(m_lam, 4) * m_m2 * M_PI * pow(Q, 2) * temp6 * temp8 * m_width2);
}

/**
 * Compute the cross section for chi2 + chi2 -> photons.
 * @param Q Center of mass energy.
 * @return sigma Sigma(chi2 + chi2 -> photon + photon)
 */
double DipoleDM::sigma_x2_x2_to_g_g(double Q) const {
    if (Q <= 2 * m_m1) {
        return 0.0;
    }
    double temp1 = pow(m_m1, 2);
    double temp2 = pow(Q, 2);
    double temp3 = pow(m_m2, 2);
    double temp4 = -temp3;
    double temp5 = -4 * temp1;
    double temp6 = temp2 + temp5;
    double temp7 = temp1 + temp4;
    double temp8 = -Q;
    double temp9 = sqrt(temp6);
    double temp10 = temp8 + temp9;
    double temp11 = 2 * temp1;
    double temp12 = -2 * temp3;
    double temp13 = Q * temp10;
    double temp14 = temp11 + temp12 + temp13;
    double temp15 = pow(temp7, 4);
    double temp16 = pow(temp7, 2);
    double temp17 = -2 * temp16;
    double temp18 = temp1 * temp2;
    double temp19 = temp17 + temp18;
    double temp20 = Q + temp9;
    double temp21 = -temp2;
    double temp22 = temp11 + temp12 + temp21;
    double temp23 = -2 * temp1;
    double temp24 = 2 * temp3;
    double temp25 = Q * temp20;
    double temp26 = temp23 + temp24 + temp25;
    double temp27 = pow(m_m2, 4);
    double temp28 = pow(Q, 4);
    double temp29 = pow(temp26, 2);
    return (pow(pow(m_ce, 2) + pow(m_cm, 2), 2) * (-pow(temp14, 3) / 12. +
            (4 * temp15) / temp14 - (4 * temp15) / (temp11 + temp12 - Q * temp20) +
            (pow(temp14, 2) * temp22) / 4. + 2 * temp19 * temp26 - pow(temp26, 3) / 12. -
            (temp22 * temp29) / 4. + 4 * temp19 * (temp1 + (Q * temp10) / 2. + temp4) -
            ((8 * pow(m_m1, 8) + 8 * pow(m_m2, 8) + 4 * pow(m_m2, 6) * temp2 -
                    temp27 * temp28 - 8 * pow(m_m1, 6) * (temp2 + 4 * temp3) -
                    4 * temp1 * temp3 * (8 * temp27 + temp28 + 4 * temp2 * temp3) +
                    pow(m_m1, 4) * (48 * temp27 + temp28 +
                            20 * temp2 * temp3)) * log(temp29 / pow(temp14, 2))) / temp22)) / (16. * pow(m_lam, 4)
            * M_PI * pow(Q, 2) * temp6);
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
    double temp1 = -Q;
    double temp2 = -m_m2;
    double temp3 = pow(m_m1, 2);
    double temp4 = pow(m_m2, 2);
    double temp5 = pow(Q, 2);
    double temp6 = pow(kW_BOSON_MASS, 2);
    double temp7 = m_m1 + m_m2 + temp1;
    double temp8 = m_m1 + m_m2 + Q;
    double temp9 = -temp4;
    double temp10 = temp3 + temp9;
    double temp11 = pow(temp10, 2);
    double temp12 = pow(Q, 4);
    return -(kALPHA_EM * sqrt(temp11 / pow(Q, 2) - 2 * (temp3 + temp4) +
                                      temp5) * sqrt(temp5 - 4 * temp6) * (48 * temp6 * temp6 * temp6 - pow(Q, 6) +
            68 * temp6 * temp6 * temp5 - 16 * temp12 * temp6) * (pow(m_cm, 2) * (-2 * temp11 +
            temp12 + (6 * m_m1 * m_m2 + temp3 + temp4) * temp5) -
            pow(m_ce, 2) * (2 * pow(m_m1 + temp2, 2) +
                    temp5) * temp7 * temp8)) / (96. * pow(m_lam, 2) * temp6 * temp6 * pow(Q, 4) * (m_m1 + Q +
            temp2) * (m_m1 + temp1 + temp2) * temp7 * temp8);
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
    double temp1 = -Q;
    double temp2 = -m_m2;
    double temp3 = pow(mf, 2);
    double temp4 = pow(Q, 2);
    double temp5 = pow(m_m1, 2);
    double temp6 = pow(m_m2, 2);
    double temp7 = m_m1 + m_m2 + temp1;
    double temp8 = m_m1 + m_m2 + Q;
    double temp9 = -temp6;
    double temp10 = temp5 + temp9;
    double temp11 = pow(temp10, 2);
    return (kALPHA_EM * sqrt(-4 * temp3 + temp4) * (2 * temp3 +
            temp4) * sqrt(temp11 / pow(Q, 2) + temp4 - 2 * (temp5 +
            temp6)) * (pow(m_cm, 2) * (pow(Q, 4) - 2 * temp11 + temp4 * (6 * m_m1 * m_m2 +
            temp5 + temp6)) - pow(m_ce, 2) * (2 * pow(m_m1 + temp2, 2) +
            temp4) * temp7 * temp8)) / (6. * pow(m_lam, 2) * pow(Q, 4) * (m_m1 + Q +
            temp2) * (m_m1 + temp1 + temp2) * temp7 * temp8);
}

/**
 * Integrand of the transfer function
 * @param T Temperature
 * @param w energy of the SM photon
 * @return Transfer function integral
 */
double DipoleDM::gamma_integrand(double w, double T) const {
    double temp1 = pow(m_m2, 2);
    double temp2 = pow(m_m1, 4);
    double temp3 = pow(m_width2, 2);
    double temp4 = pow(m_m1, 3);
    double temp5 = pow(m_m1, 2);
    double temp6 = pow(m_m2, 4);
    double temp7 = temp1 + temp3;
    double temp8 = pow(m_width2, 4);
    double temp9 = pow(m_m2, 6);
    double temp10 = pow(m_width2, 6);
    double temp11 = pow(m_m2, 8);
    double temp12 = pow(m_width2, 8);
    double temp13 = pow(m_m2, 10);
    double temp14 = pow(m_m2, 18);
    double temp15 = pow(temp7, 5);
    double temp16 = pow(w, 2);
    double temp17 = pow(m_m2, 16);
    double temp18 = pow(temp7, 4);
    double temp19 = pow(m_m2, 14);
    double temp20 = pow(temp7, 3);
    double temp21 = pow(m_m2, 12);
    double temp22 = pow(temp7, 2);
    double temp23 = pow(temp7, 6);
    double temp24 = pow(m_m1, 29);
    double temp25 = pow(m_m1, 28);
    double temp26 = pow(m_m1, 30);
    double temp27 = pow(w, 3);
    double temp28 = pow(m_m1, 26);
    double temp29 = pow(m_m1, 27);
    double temp30 = pow(m_m1, 24);
    double temp31 = pow(m_m1, 25);
    double temp32 = pow(m_m1, 5);
    double temp33 = pow(m_m1, 21);
    double temp34 = pow(m_m1, 22);
    double temp35 = pow(m_m1, 23);
    double temp36 = pow(m_m1, 20);
    double temp37 = pow(m_m1, 17);
    double temp38 = pow(m_m1, 9);
    double temp39 = pow(m_m1, 18);
    double temp40 = pow(m_m1, 8);
    double temp41 = pow(m_m1, 6);
    double temp42 = pow(m_m1, 16);
    double temp43 = pow(m_m1, 19);
    double temp44 = pow(m_m1, 7);
    double temp45 = pow(m_m1, 14);
    double temp46 = pow(m_m1, 10);
    double temp47 = pow(m_m1, 15);
    double temp48 = pow(m_m1, 12);
    double temp49 = pow(m_m1, 11);
    double temp50 = pow(m_m1, 13);
    return (8 * pow(pow(m_ce, 2) + pow(m_cm, 2), 2) * temp4 * (temp2 +
            temp1 * temp3 - 2 * temp1 * temp5 +
            temp6) * pow(w, 8) * (365472 * temp14 * temp15 * temp16 + (740745 * temp1 +
            413042 * temp16) * temp25 + 63595 * temp26 + temp1 * temp28 * (436745 * temp1 +
            5992788 * temp16 + 91185 * temp3) +
            1728 * temp17 * temp18 * temp5 * (569 * temp16 * temp3 + 6 * temp1 * (71 * temp16 +
                    315 * temp3) + 1890 * temp6) + 2 * temp34 * temp6 * (4014688 * temp16 * temp3 +
            temp1 * (-38345168 * temp16 + 2655410 * temp3) + 8860695 * temp6 +
            15275 * temp8) + temp30 * (245938 * temp1 * temp16 * temp3 + (13450346 * temp16 +
            1706015 * temp3) * temp6 - 10777345 * temp9) + 2 * temp39 * (1205 * temp10 +
            5 * (4946156 * temp16 + 154811 * temp3) * temp6 + 2106324 * temp16 * temp8 +
            temp1 * (-69172784 * temp16 * temp3 + 4515595 * temp8) -
            75607855 * temp9) * temp9 - 36 * temp19 * temp2 * temp20 * (2 * (73658 * temp16 +
            292485 * temp3) * temp6 + temp1 * (50304 * temp16 * temp3 - 7905 * temp8) -
            4692 * temp16 * temp8 + 592875 * temp9) + 2 * temp36 * temp6 * ((66429754 * temp16
            - 13453320 * temp3) * temp6 + 34262 * temp16 * temp8 +
            temp1 * (22383072 * temp16 * temp3 + 613105 * temp8) + 13609415 * temp9) +
            temp13 * temp46 * temp7 * (50973925 * temp11 - 9 * (35 * temp12 +
                    16092 * temp10 * temp16) + temp6 * (658540268 * temp16 * temp3 -
                    35793090 * temp8) + 20 * temp1 * (47789 * temp10 + 4473339 * temp16 * temp8) +
                    (724794260 * temp16 - 27247260 * temp3) * temp9) +
            48 * temp21 * temp22 * temp41 * (1311465 * temp11 + 189 * temp10 * temp16 +
                    temp6 * (1704619 * temp16 * temp3 + 151365 * temp8) + temp1 * (-6795 * temp10 +
                    535251 * temp16 * temp8) + (1474357 * temp16 + 1469625 * temp3) * temp9) +
            temp11 * temp45 * (-247599705 * temp11 - 865 * temp12 + 280464 * temp10 * temp16
                    + temp6 * (782252656 * temp16 * temp3 - 53538830 * temp8) +
                    20 * temp1 * (298201 * temp10 + 562252 * temp16 * temp8) - 20 * (-54220084 * temp16
                    + 15398371 * temp3) * temp9) + 2 * temp42 * temp9 * (134103655 * temp11 -
            2634 * temp10 * temp16 - temp6 * (40111446 * temp16 * temp3 + 8049355 * temp8) +
            temp1 * (175475 * temp10 + 25973742 * temp16 * temp8) + (-290753742 * temp16 +
            74643625 * temp3) * temp9) - temp13 * temp40 * temp7 * (96619485 * temp13 -
            162 * temp12 * temp16 + 15 * temp1 * (1355 * temp12 + 61176 * temp10 * temp16) +
            2 * temp11 * (144511991 * temp16 + 59013150 * temp3) + temp6 * (-7996020 * temp10
            + 162459748 * temp16 * temp8) + 2 * (218370964 * temp16 * temp3 +
            6695235 * temp8) * temp9) + temp11 * temp48 * (95976965 * temp13 -
            3294 * temp12 * temp16 + temp1 * (-8715 * temp12 + 16443408 * temp10 * temp16) +
            2 * temp11 * (-565359371 * temp16 + 120083490 * temp3) +
            4 * temp6 * (2022905 * temp10 - 94198711 * temp16 * temp8) +
            2 * (-732604144 * temp16 * temp3 + 69665175 * temp8) * temp9) -
            272160 * m_m1 * temp14 * temp15 * w + 90740 * temp24 * w + 765880 * temp1 * temp29 * w
            - 20 * temp1 * (140423 * temp1 - 7337 * temp3) * temp31 * w +
            6480 * temp17 * temp18 * (1697 * temp1 + 897 * temp3) * temp4 * w -
            40 * (265467 * temp1 - 59581 * temp3) * temp35 * temp6 * w -
            1080 * temp19 * temp20 * temp32 * (33130 * temp1 * temp3 + 46547 * temp6 +
                    455 * temp8) * w + 40 * temp33 * temp6 * (-23725 * temp1 * temp3 + 2631056 * temp6 +
            2007 * temp8) * w - 160 * temp43 * (425630 * temp1 * temp3 + 2302858 * temp6 -
            16691 * temp8) * temp9 * w + 120 * temp21 * temp22 * temp44 * (-467 * temp10 +
            854067 * temp3 * temp6 + 129333 * temp1 * temp8 + 794107 * temp9) * w -
            80 * temp11 * temp47 * (-15682 * temp10 + 9043644 * temp3 * temp6 +
                    1216791 * temp1 * temp8 + 11445371 * temp9) * w + 40 * temp37 * temp9 * (587 * temp10
            + 8141273 * temp3 * temp6 + 272263 * temp1 * temp8 + 18423717 * temp9) * w -
            40 * temp13 * temp49 * (673190 * temp1 * temp10 + 6529157 * temp11 - 705 * temp12 +
                    5677332 * temp6 * temp8 + 11532594 * temp3 * temp9) * w +
            20 * temp11 * temp50 * (679300 * temp1 * temp10 + 34668601 * temp11 - 103 * temp12
                    + 11939142 * temp6 * temp8 + 43033940 * temp3 * temp9) * w -
            20 * temp13 * temp38 * (-162809 * temp1 * temp12 + 1674283 * temp13 +
                    4496683 * temp11 * temp3 + 519910 * temp10 * temp6 + 3505182 * temp8 * temp9 +
                    63 * pow(m_width2, 10)) * w) * pow(sinh(w / (2.0 * T)), -2)) / (15. * pow(m_lam, 4) * (
            temp2 - 2 * (temp1 - 2 * temp16) * temp5 + temp1 * temp7 - 4 * m_m1 * temp1 * w +
                    4 * temp4 * w) * (12719 * pow(m_m1, 33) + pow(m_m1, 31) * (84554 * temp1 -
            64932 * temp16) + 4756320 * m_m1 * temp14 * temp16 * temp23 +
            4665600 * temp14 * temp23 * temp27 - 4 * temp1 * temp24 * (112473 * temp1 +
            135906 * temp16 - 7739 * temp3) +
            288 * temp15 * temp17 * temp4 * (-1925 * temp16 * temp3 + 9 * temp1 * (-13565 * temp16
                    + 84 * temp3) + 756 * temp6) + temp31 * temp6 * (temp1 * (4578448 * temp16 -
            594686 * temp3) - 4304 * temp16 * temp3 + 10439591 * temp6 + 24347 * temp8) +
            temp29 * (28452 * temp1 * temp16 * temp3 + 2 * (829418 * temp16 +
                    180005 * temp3) * temp6 - 832342 * temp9) -
            12 * temp18 * temp19 * temp32 * (2 * (-5198203 * temp16 + 82689 * temp3) * temp6 +
                    4554 * temp16 * temp8 + temp1 * (-1686796 * temp16 * temp3 + 4467 * temp8) +
                    160911 * temp9) - 8 * temp35 * temp6 * ((3554699 * temp16 + 941983 * temp3) * temp6
            - 20059 * temp16 * temp8 - 8 * temp1 * (-12658 * temp16 * temp3 + 7831 * temp8) +
            4389879 * temp9) + 4 * temp33 * temp9 * (1648 * temp10 + 2 * (6973099 * temp16 +
            4785789 * temp3) * temp6 + 634750 * temp16 * temp8 +
            3 * temp1 * (-95044 * temp16 * temp3 + 73373 * temp8) + 16942675 * temp9) -
            temp21 * temp22 * temp38 * (18892195 * temp11 + 9 * temp10 * (-3984 * temp16 +
                    7 * temp3) + 2 * temp6 * (-98132888 * temp16 * temp3 + 3388995 * temp8) -
                    16 * temp1 * (9311 * temp10 - 717953 * temp16 * temp8) + 24 * (-15452110 * temp16 +
                    1075801 * temp3) * temp9) + 4 * temp20 * temp21 * temp44 * (1975965 * temp11 -
            288 * temp10 * temp16 + 3 * temp6 * (-7849976 * temp16 * temp3 + 111961 * temp8) +
            temp1 * (-969 * temp10 + 914452 * temp16 * temp8) + (-67088284 * temp16 +
            2312817 * temp3) * temp9) + 3 * temp11 * temp37 * (21391919 * temp11 + 103 * temp12
            + 977728 * temp10 * temp16 + 2 * temp6 * (-24927072 * temp16 * temp3 +
            6946129 * temp8) + 12 * temp1 * (54863 * temp10 - 1624368 * temp16 * temp8) +
            4 * (-126800 * temp16 + 8786081 * temp3) * temp9) -
            4 * temp43 * temp9 * (20935147 * temp11 - 22098 * temp10 * temp16 +
                    temp6 * (-9588566 * temp16 * temp3 + 3652433 * temp8) - temp1 * (73673 * temp10 +
                    644470 * temp16 * temp8) + (11186894 * temp16 + 21389573 * temp3) * temp9) -
            4 * temp13 * temp50 * (3941085 * temp13 + temp1 * (-295903 * temp12 +
                    20119728 * temp10 * temp16) + temp12 * (-312130 * temp16 + 59 * temp3) +
                    temp11 * (-37962698 * temp16 + 4423991 * temp3) + temp6 * (-2451762 * temp10 +
                    56285788 * temp16 * temp8) + (3075232 * temp16 * temp3 -
                    1672894 * temp8) * temp9) + 2 * temp13 * temp49 * temp7 * (13387473 * temp13 -
            1014 * temp12 * temp16 + temp1 * (-2487 * temp12 + 7561184 * temp10 * temp16) +
            2 * temp11 * (-80749775 * temp16 + 10280478 * temp3) - 4 * temp6 * (89541 * temp10
            - 7683553 * temp16 * temp8) + (-93975536 * temp16 * temp3 +
            6817806 * temp8) * temp9) - 2 * temp11 * temp47 * (10325267 * temp13 -
            10214 * temp12 * temp16 - 3 * temp1 * (11311 * temp12 + 2435800 * temp10 * temp16)
            + 10 * temp11 * (620433 * temp16 + 3078982 * temp3) + 36 * temp6 * (136695 * temp10
            - 2653093 * temp16 * temp8) + 2 * (-56556916 * temp16 * temp3 +
            12277753 * temp8) * temp9) + 18148 * pow(m_m1, 32) * w -
            4 * temp1 * temp25 * (217383 * temp1 + 209512 * temp16 - 24593 * temp3) * w +
            12960 * temp15 * temp17 * temp5 * (-8 * temp16 * temp3 + temp1 * (-2332 * temp16 +
                    133 * temp3) + 133 * temp6) * w - 4 * temp1 * temp28 * (9960 * temp16 * temp3 -
            temp1 * (288992 * temp16 + 291783 * temp3) + 848521 * temp6) * w +
            16 * temp30 * temp6 * (temp1 * (1650800 * temp16 - 112449 * temp3) -
                    80172 * temp16 * temp3 + 2118728 * temp6 + 7397 * temp8) * w -
            48 * temp18 * temp19 * temp2 * (6 * (-341545 * temp16 + 50059 * temp3) * temp6 +
                    490 * temp16 * temp8 + temp1 * (-326460 * temp16 * temp3 + 6783 * temp8) +
                    293571 * temp9) * w - 16 * temp34 * temp6 * ((8974454 * temp16 +
            1758253 * temp3) * temp6 - 2344 * temp16 * temp8 -
            2 * temp1 * (-50531 * temp16 * temp3 + 70504 * temp8) + 7320933 * temp9) * w +
            8 * temp36 * temp9 * (5649 * temp10 + (47900620 * temp16 +
                    18393251 * temp3) * temp6 - 4036 * temp16 * temp8 +
                    temp1 * (8539848 * temp16 * temp3 + 434261 * temp8) + 29774943 * temp9) * w -
            8 * temp21 * temp22 * temp40 * (15289633 * temp11 + 63 * temp12 -
                    9720 * temp10 * temp16 + 8 * temp6 * (-2742283 * temp16 * temp3 + 661902 * temp8) +
                    temp1 * (-166058 * temp10 + 45296 * temp16 * temp8) + 10 * (-3715936 * temp16 +
                    2075097 * temp3) * temp9) * w + 8 * temp20 * temp21 * temp41 * (6744399 * temp11 -
            72 * temp10 * temp16 + temp6 * (-8945600 * temp16 * temp3 + 1117329 * temp8) +
            temp1 * (-3339 * temp10 + 383548 * temp16 * temp8) + 3 * (-8441708 * temp16 +
            2621689 * temp3) * temp9) * w + 4 * temp11 * temp42 * (60622233 * temp11 +
            1553 * temp12 + 299584 * temp10 * temp16 + 2 * temp6 * (71798656 * temp16 * temp3 +
            24204759 * temp8) + 20 * temp1 * (122325 * temp10 + 571376 * temp16 * temp8) +
            28 * (6851552 * temp16 + 3964043 * temp3) * temp9) * w -
            8 * temp39 * temp9 * (38774933 * temp11 - 4680 * temp10 * temp16 +
                    5 * temp6 * (7015880 * temp16 * temp3 + 1613419 * temp8) -
                    7 * temp1 * (30325 * temp10 - 45328 * temp16 * temp8) + (81526704 * temp16 +
                    43157663 * temp3) * temp9) * w + 4 * temp13 * temp46 * temp7 * (42235227 * temp13 -
            40 * temp12 * temp16 + temp1 * (-2637 * temp12 + 2420200 * temp10 * temp16) +
            temp11 * (-93276704 * temp16 + 64235556 * temp3) - 4 * temp6 * (363483 * temp10 +
            2178230 * temp16 * temp8) + 18 * (-4878948 * temp16 * temp3 +
            1141613 * temp8) * temp9) * w - 4 * temp13 * temp48 * (28490915 * temp13 +
            temp1 * (-1844233 * temp12 + 2864320 * temp10 * temp16) +
            temp12 * (-199528 * temp16 + 339 * temp3) + temp11 * (-123107656 * temp16 +
            36517223 * temp3) - 2 * temp6 * (6965045 * temp10 + 29484616 * temp16 * temp8) -
            2 * (90842368 * temp16 * temp3 + 2029605 * temp8) * temp9) * w -
            4 * temp11 * temp45 * (14355025 * temp13 - 3560 * temp12 * temp16 -
                    temp1 * (130975 * temp12 + 1415248 * temp10 * temp16) +
                    132 * temp11 * (1265610 * temp16 + 474829 * temp3) + 4 * temp6 * (3316517 * temp10
                    + 8677264 * temp16 * temp8) + 2 * (92259432 * temp16 * temp3 +
                    29563723 * temp8) * temp9) * w + 4 * temp26 * (-17016 * temp27 + 41047 * temp1 * w)));
}

/**
 * Compute the thermal transfer function.
 * @param T
 * @return
 */
double DipoleDM::gamma(double T) const {
    using integrate::qagp;
    using integrate::qagi;
    // Integration parameters
    double epsabs = 1e-10;
    double epsrel = 1e-5;
    double abserr;
    int neval, ier;
    int inf = 1;

    auto f = [T, this](double w) {
        return gamma_integrand(w, T);
    };

    // Set break pt at m2 - m1
    double brkpt = m_m2 - m_m1;
    std::vector<double> pts = {0.0, brkpt};
    // Integrate from brkpt -> infinity
    double gam = qagi(f, brkpt, inf, epsabs, epsrel, &abserr, &neval, &ier);
    // Integrate from 0 -> brkpt
    gam += qagp(f, 0.0, brkpt, pts.size(), pts.data(), epsabs, epsrel, &abserr, &neval, &ier);

    return gam / (48.0 * pow(M_PI * m_m1, 3) * 2.0 * T);
}


double DipoleDM::average_p6_E3(double T, double y) const {
    using integrate::qagp;
    using integrate::qagi;
    // Integration parameters
    double epsabs = 1e-10;
    double epsrel = 1e-5;
    double abserr;
    int neval, ier;
    int inf = 1;

    auto f = [this, T, y](double p) {
        double E = sqrt(p * p + m_m1 * m_m1);
        double s = sm_entropy_density(T);
        double Tx = y * pow(s, 2.0 / 3.0) / m_m1;
        double pf = x1.get_g() / (2.0 * M_PI * M_PI * x1.neq(Tx));
        return pow(p * p / E, 3) * exp(-E / Tx);
    };

    return qagi(f, 0.0, inf, epsabs, epsrel, &abserr, &neval, &ier);
}

/**
 * Returns the annihilations cross section for chi_i + chi_j -> X
 * @param Q Center of mass energy
 * @param i Index of chi_i. Equal to 1 or 2.
 * @param j Index of chi_j. Equal to 1 or 2.
 */
double DipoleDM::annihilation_cross_section(double Q, int i, int j) const {
    if (i == 1 && j == 1) {
        return sigma_x1_x1_to_g_g(Q);
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
        return sigma_x2_x2_to_g_g(Q);
    }

}

/**
 * Returns the effective W factor for computing the effective thermal cross
 * section.
 * @param Q Center of mass energy.
 */
double DipoleDM::weff(double Q) const {
    double sig11 = 0, sig12 = 0, sig22 = 0;
    double p11 = 0, p12 = 0, p22 = 0;
    double w11 = 0, w12 = 0, w22 = 0;

    double s = Q * Q;
    if (Q > 2 * m_m1) {
        sig11 = annihilation_cross_section(Q, 1, 1);
        p11 = sqrt(s - pow(m_m1 + m_m1, 2)) * sqrt(s - pow(m_m1 - m_m1, 2)) / (2 * Q);
        w11 = 4 * p11 * Q * sig11;
    } else {
        return 0.0;
    }
    if (Q > m_m1 + m_m2) {
        sig12 = annihilation_cross_section(Q, 1, 2);
        p12 = sqrt(s - pow(m_m1 + m_m2, 2)) * sqrt(s - pow(m_m1 - m_m2, 2)) / (2 * Q);
        w12 = 4 * p12 * Q * sig12;
    }
    if (Q > 2 * m_m2) {
        sig22 = annihilation_cross_section(Q, 2, 2);
        p22 = sqrt(s - pow(m_m2 + m_m2, 2)) * sqrt(s - pow(m_m2 - m_m2, 2)) / (2 * Q);
        w22 = 4 * p22 * Q * sig22;
    }

    return (p11 * w11 + 2 * p12 * w12 + p22 * w22) / p11;
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
    const double wef = weff(sqrts);
    return peff * peff * wef * besselk1e(sqrts / T) * exp(2.0 * x - sqrts / T);
}

/**
 * Compute the effective, thermally average annihilation cross section for
 * chi1 + chi2
 * @param x mass of chi1 / temperature.
 * @return <sigma_eff*v>
 */
double DipoleDM::thermal_cross_section(double x) const {
    //using boost::math::quadrature::gauss_kronrod;
    using namespace integrate;
    using namespace special_functions;

    const double rat = m_m2 / m_m1;
    const double denom = m_m1 * m_m1 * (
            besselk2e(x) + rat * rat * besselk2e(rat * x) * exp(x * (1.0 - rat))
    );
    const double pf = 1.0 / (denom * denom * x / m_m1);
    auto integrand = [this, &x](double peff) {
        const double res = thermal_cross_section_integrand(peff, x);
        return res;
    };

    double epsabs = 1e-10;
    double epsrel = 1e-5;
    double abserr;
    int neval, ier;

    return pf * qagi(integrand, 0.0, 1, epsabs, epsrel, &abserr, &neval, &ier);
}

double DipoleDM::thermal_cross_section2(double x) const {
    using special_functions::besselk2e;
    using namespace integrate;

    auto f = [this, x](double st) -> double {
        // Perform inner integral using gauss-laguerre quadrature
        double inner = 0.0;
        for (size_t i = 0; i < gausslaguerre_nodes.size(); i++) {
            double w = gausslaguerre_nodes[i];
            double wt = gausslaguerre_weights[i];
            double sqrt_fac = sqrt(((st - 1.0) * ((1.0 + w) * (1.0 + w) - 4.0 * st * x * x)) / st);
            inner += wt * (
                    (1.0 + w) * sqrt_fac +
                            2.0 * x * x * log((1 + w - sqrt_fac) / (1 + w + sqrt_fac))
            );
        }
        inner *= exp(2.0 * x - 1) / (8.0 * st * x * x * x);
        double s = st * (4.0 * m_m1 * m_m1);
        // TODO: Only using x1 + x1 -> ... Should x1 + x2 -> ... and x2 + x2 -> .. be included?
        double sig = annihilation_cross_section(sqrt(s), 1, 1);
        double res = inner * sig * 4.0 * st * (2.0 * st - 1.0) * x * x * x / (3.0 * pow(besselk2e(x), 2));
        return res;
    };

    double epsabs = 1e-10;
    double epsrel = 1e-5;
    double abserr;
    int neval, ier;
    // Perform outter integral
    return qagi(f, 1.0, 1, epsabs, epsrel, &abserr, &neval, &ier);
}

/*
 * Compute the relic density of the dark matter.
 */
ODESolution DipoleDM::solve_boltzmann(
        double xstart = 1.0,
        double xend = 500.0,
        double reltol = 1e-3,
        double abstol = 1e-9,
        const std::string &t_alg = "radau"
) const {
    using namespace diffeq;
    using namespace cosmology;

    auto logx_span = std::make_pair(std::log(xstart), std::log(xend));

    DipoleDMBoltzmann boltz{std::make_shared<DipoleDM>(*this)};
    double Tinit = m_m1 / exp(logx_span.first);

    diffeq::Vector<double> winit{1};
    winit(0) = log(boltz.x1.neq(Tinit) / sm_entropy_density(Tinit));

    ODEProblem problem{std::make_shared<DipoleDMBoltzmann>(boltz), winit, logx_span};

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

double DipoleDM::relic_density(
        double xstart = 1.0,
        double xend = 500.0,
        double reltol = 1e-3,
        double abstol = 1e-9,
        const std::string &t_alg = "radau"
) const {
    auto sol = solve_boltzmann(xstart, xend, reltol, abstol, t_alg);
    double yinf = exp(sol.us.back()[0]);
    return m_m1 * yinf * kS_TODAY / kRHO_CRIT;
}


ODESolution DipoleDM::solve_temperature(
        double xstart = 50.0,
        double xend = 1e4,
        double reltol = 1e-5,
        double abstol = 1e-9,
        const std::string &t_alg = "radau"
) {
    using namespace diffeq;
    using namespace cosmology;

    auto logx_span = std::make_pair(std::log(xstart), std::log(xend));

    DipoleDMBoltzmannTemp boltz{std::make_shared<DipoleDM>(*this)};
    double Tinit = m_m1 / exp(logx_span.first);

    diffeq::Vector<double> yinit{1};
    // Assume that DM is in kinetic equillibrium to start
    yinit(0) = 1.73216203306263 * m_m1 * pow(sm_heff(Tinit), -2.0 / 3.0) / Tinit;

    ODEProblem problem{std::make_shared<DipoleDMBoltzmannTemp>(boltz), yinit, logx_span};

    if (t_alg == "rodas") {
        Rodas alg{};
        ODEIntegratorOptions opts{};
        opts.abstol = abstol;
        opts.reltol = reltol;
        opts.dtmax = (logx_span.second - logx_span.first) / 1000.0;

        return solve(problem, alg, opts);
    } else {
        Radau5 alg{};
        ODEIntegratorOptions opts{};
        opts.abstol = abstol;
        opts.reltol = reltol;
        opts.dtmax = (logx_span.second - logx_span.first) / 1000.0;

        return solve(problem, alg, opts);
    }
}


}
}

#endif //LANRE_DM_MODELS_KINETIC_RECOUPLING_HPP
