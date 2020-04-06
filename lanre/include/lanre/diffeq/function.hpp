//
// Created by Logan Morrison on 3/15/20.
//

#ifndef LANRE_DIFFEQ_FUNCTION_HPP
#define LANRE_DIFFEQ_FUNCTION_HPP

#include "lanre/diffeq/base.hpp"

namespace lanre {
namespace diffeq {

struct ODEFunction {

    ODEFunction() = default;

    ~ODEFunction() = default;

    /**
     * Computes the RHS of ODE: partial u(i) / partial t
     * @param df vector to store time derivatives of u: du_i / dt
     * @param u Current value of solution vector
     * @param t Current time value
     */
    virtual void dudt(Vector<double> &du, const Vector<double> &u, double t) = 0;

    /**
     * Computes the gradient of RHS of ODE wrt time: df_i / dt. By default,
     * gradient is computed by finite differences.
     *
     * @param df vector to store gradient
     * @param u Current value of solution vector
     * @param t Current time value
     */
    virtual void dfdt(Vector<double> &df, const Vector<double> &u, double t) {
        double delt = std::sqrt(std::numeric_limits<double>::epsilon() * std::max(1.e-5, std::abs(t)));
        Vector<double> du(u.size());
        dudt(du, u, t);
        dudt(df, u, t + delt);
        for (int j = 0; j < u.size(); j++) {
            df(j) = (df(j) - du(j)) / delt;
        }
    }

    /**
     * Computes the jacobian matrix: df_i / du_j. By default, jacobian is
     * computed by finite differences.
     *
     * @param df matrix to store jacobian
     * @param u Current value of solution vector
     * @param t Current time value
     */
    virtual void dfdu(Matrix<double> &df, const Vector<double> &u, double t) {
        static double eps = std::numeric_limits<double>::epsilon();

        int n = u.size();
        Vector<double> utemp(u);
        Vector<double> du0(n);
        Vector<double> du1(n);
        Vector<double> du2(n);
        Vector<double> du3(n);
        dudt(du0, u, t);

        double delt, ysafe;
        for (int i = 0; i < n; i++) {
            ysafe = utemp(i);
            delt = sqrt(eps * std::max(1.0e-5, fabs(ysafe)));
            utemp(i) = ysafe + delt;
            dudt(du1, utemp, t);
            for (int j = 0; j < n; j++)
                df(j, i) = (du1(j) - du0(j)) / delt;
            utemp(i) = ysafe;
        }
    }

    Matrix<double> mass_matrix{};
};

}
}

#endif //LANRE_DIFFEQ_FUNCTION_HPP
