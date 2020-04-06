/* Created by Logan Morrison on 1/27/20.
 *
 * The "c" abstract base class defines functions to compute the
 * temperature ratio between the standard model bath, and a secluded (decoupled)
 * sector. The temperature ratio, xi, is defined as the ratio of decoupled to
 * standard model temperature: xi = T_d / T_sm.
 *
 * Classes deriving from ThermallyDecoupledModel must define the following:
 *      double m_xi_inf: ratio of dark to SM temperature at above EW scale
 *      double m_hd_inf: value of entropic dof in dark sector above EW scale
 *      double m_sum_g: sum of internal dof in the dark sector
 *      double m_gl: internal dof of the lightest particle in dark sector
 *      double m_ml: mass of the lightest dof in dark sector
 *
 *      double dark_heff(double Td): function to compute the entropic
 *                                          dof in the dark sector
 *
 * The temperature ratio, xi = Td / Tsm, can be computed in two ways. The
 * first, if T_d is known, then we can determine xi using:
 *      compute_xi_const_td(Td)
 * second, if T_sm is known, then:
 *      compute_xi_const_tsm(Tsm)
 */

#ifndef LANRE_DM_MODELS_THERMALLY_DECOUPLED_MODEL_HPP
#define LANRE_DM_MODELS_THERMALLY_DECOUPLED_MODEL_HPP

#include "lanre/cosmology/standard_model.hpp"
#include <boost/math/tools/roots.hpp>
#include <boost/math/special_functions/lambert_w.hpp>

using namespace lanre::cosmology;

namespace lanre {
namespace dm_models {

/**
 * Abstract base class for computing the ratio of dark to SM temperature for a
 * model in which the dark sector is thermally decoupled from the SM.
 */
class ThermallyDecoupledModel {
protected:

    /**
     * Compute the upper bound on xi assuming Td is known.
     * @param Td Temperature of the dark bath
     * @return Upper bound on xi
     */
    double xi_upper_bound_const_td(double Td) const;

    /**
     * Compute the lower bound on xi assuming Td is known.
     * @param Td Temperature of the dark bath
     * @return Lower bound on xi
     */
    double xi_lower_bound_const_td(double Td) const;

    /**
     * Compute the upper bound on xi assuming Tsm is known.
     * @param Tsm Temperature of the standard model bath
     * @return Upper bound on xi
     */
    double xi_upper_bound_const_tsm(double Tsm) const;

    /**
     * Compute the lower bound on xi assuming Tsm is known.
     * @param Tsm Temperature of the standard model bath
     * @return Lower bound on xi
     */
    double xi_lower_bound_const_tsm(double Tsm) const;

    double m_xi_inf = 1.0; // ratio of dark to SM temperature at above EW scale
    double m_hd_inf = 1.0; // value of entropic dof in dark sector above EW scale
    double m_sum_g = 1.0; // sum of internal dof in the dark sector
    double m_gl = 1.0; // internal dof of the lightest particle in dark sector
    double m_ml = 0.0; // mass of the lightest dof in dark sector

    ThermallyDecoupledModel() = default;

    ~ThermallyDecoupledModel() = default;

    /**
    * Compute the entropic dof of the dark sector
    * @param Td Temperature of the dark bath
    * @return h_eff entropic dof of the dark sector
    */
    virtual double dark_heff(double Td) const = 0;

    /**
     * Compute xi = Td/Tsm assuming Td is known
     * @param Td Temperature of the dark bath
     * @return Value of xi = Td / Tsm
     */
    double compute_xi_const_td(double Td) const;

    /**
     * Compute xi assuming Tsm is known
     * @param Tsm Temperature of the standard model bath
     * @return Value of xi = Td / Tsm
     */
    double compute_xi_const_tsm(double Tsm) const;
};

/**
 * Lower bound on xi = Tsm / Td assuming Td is a constant.
 * @param Td Dark sector temperature.
 * @return Lower bound on xi.
 */
double ThermallyDecoupledModel::xi_lower_bound_const_td(double Td) const {
    double hd = dark_heff(Td);
    return cbrt(m_hd_inf / hd * SM_HEFF_0 / SM_HEFF_INF) * m_xi_inf;
}

/**
 * Upper bound on xi = Tsm / Td assuming Td is a constant.
 * @param Td Dark sector temperature.
 * @return Upper bound on xi.
 */
double ThermallyDecoupledModel::xi_upper_bound_const_td(double Td) const {
    double hd = dark_heff(Td);
    return cbrt(m_hd_inf / hd) * m_xi_inf;
}

/**
 * Lower bound on xi = Tsm / Td assuming Tsm is a constant.
 * @param Tsm Temperature of the SM bath.
 * @return Lower bound.
 */
double ThermallyDecoupledModel::xi_lower_bound_const_tsm(double Tsm) const {
    return cbrt(sm_heff(Tsm) * m_hd_inf / m_sum_g / SM_HEFF_INF) * m_xi_inf;
}

/**
 * Upper bound on xi = Tsm / Td assuming Tsm is a constant
 * @param Tsm Temperature of the SM bath.
 * @return Lower bound.
 */
double ThermallyDecoupledModel::xi_upper_bound_const_tsm(double Tsm) const {
    if (m_ml <= 0.0)
        return cbrt(sm_heff(Tsm) * m_hd_inf / m_gl / SM_HEFF_INF) * m_xi_inf;
    else {
        using namespace boost::math;

        double xl = m_ml / Tsm;
        double hsm = sm_heff(Tsm);

        double lw_arg_num = pow(45.0 * m_gl * SM_HEFF_INF * xl * xl * xl, 2);
        double lw_arg_den = pow(4.0 * m_hd_inf * hsm * m_xi_inf, 2) * pow(M_PI, 7);
        return 2.0 * xl / lambert_w0(lw_arg_num / lw_arg_den);
    }
}

/**
 * Compute the value of xi = Td / Tsm assuming Td is a constant.
 * @param Td Temperature of the dark bath.
 * @return Ratio of dark to SM temperatures: xi = Td / Tsm.
 */
double ThermallyDecoupledModel::compute_xi_const_td(double Td) const {
    using namespace boost::math::tools;
    double hd = dark_heff(Td);

    double lb = xi_lower_bound_const_td(Td);
    double ub = xi_upper_bound_const_td(Td);

    double c1 = m_hd_inf * pow(m_xi_inf, 3) / SM_HEFF_INF;
    auto f = [hd, Td, c1](double xi) {
        return hd * pow(xi, 3) - sm_heff(Td / xi) * c1;
    };
    auto tol = [](double min, double max) { return abs(max - min) <= 1e-8; };

    std::pair<double, double> res = bisect(f, 0.8 * lb, 1.2 * ub, tol);
    return (res.second + res.first) / 2.0;

    //return brent(f, 0.8 * lb, 1.2 * ub, 1e-4);
}

/**
 * Compute the value of xi = Td / Tsm assuming Tsm is a constant.
 * @param Tsm Temperature of the SM bath.
 * @return Ratio of dark to SM temperatures: xi = Td / Tsm.
 */
double ThermallyDecoupledModel::compute_xi_const_tsm(double Tsm) const {
    using namespace boost::math::tools;
    double ub = xi_upper_bound_const_tsm(Tsm);
    double lb = xi_lower_bound_const_tsm(Tsm);
    double hsm = sm_heff(Tsm);

    double c1 = hsm * m_hd_inf * pow(m_xi_inf, 3) / SM_HEFF_INF;
    auto f = [this, Tsm, c1](double xi) {
        return dark_heff(xi * Tsm) * pow(xi, 3) - c1;
    };
    auto tol = [](double min, double max) { return abs(max - min) <= 1e-8; };
    std::pair<double, double> res = bisect(f, 0.8 * lb, 1.2 * ub, tol);
    return (res.second + res.first) / 2.0;
    //return brent(f, 0.8 * lb, 1.2 * ub, 1e-4);
}
}
}

#endif //LANRE_DM_MODELS_THERMALLY_DECOUPLED_MODEL_HPP
