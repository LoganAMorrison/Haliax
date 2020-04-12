//
// Created by Logan Morrison on 4/11/20.
//

#ifndef LANRE_OPTIMIZE_CHKDER_HPP
#define LANRE_OPTIMIZE_CHKDER_HPP

#include "lanre/lanre.hpp"

namespace lanre {
namespace optimize {

/**
 * Checks the gradients of M functions in N variables.
 *
 * This function checks the gradients of M nonlinear functions
 * in N variables, evaluated at a point x, for consistency with
 * the functions themselves.  The user must call chkder twice,
 * first with mode = 1 and then with mode = 2.
 *
 * mode = 1: on input, x must contain the point of evaluation.
 *           on output, xp is set to a neighboring point.
 *
 * mode = 2: on input, fvec must contain the functions and the
 *           rows of fjac must contain the gradients
 *           of the respective functions each evaluated
 *           at x, and fvecp must contain the functions
 *           evaluated at xp.
 *           on output, err contains measures of correctness of
 *           the respective gradients.
 *
 * The function does not perform reliably if cancellation or
 * rounding errors cause a severe loss of significance in the
 * evaluation of a function. therefore, none of the components
 * of x should be unusually small (in particular, zero) or any
 * other value which may cause loss of significance.
 *
 * @param x the point at which the jacobian is to be checked.
 * @param fvec is an array of length m. on input when mode = 2,
 *        fvec must contain the functions evaluated at x.
 * @param fjac is an m by n array. on input when mode = 2,
 *        the rows of fjac must contain the gradients of
 *        the respective functions evaluated at x.
 * @param xp is an array of length n. on output when mode = 1,
 *        xp is set to a neighboring point of x.
 * @param fvecp is an array of length m. on input when mode = 2,
 *        fvecp must contain the functions evaluated at xp.
 * @param mode is an integer input variable set to 1 on the first call
 *        and 2 on the second. other values of mode are equivalent
 *        to mode = 1.
 * @param err  is an array of length m. on output when mode = 2,
 *        err contains measures of correctness of the respective
 *        gradients. if there is no severe loss of significance,
 *        then if err(i) is 1.0 the i-th gradient is correct,
 *        while if err(i) is 0.0 the i-th gradient is incorrect.
 *        for values of err between 0.0 and 1.0, the categorization
 *        is less certain. in general, a value of err(i) greater
 *        than 0.5 indicates that the i-th gradient is probably
 *        correct, while a value of err(i) less than 0.5 indicates
 *        that the i-th gradient is probably incorrect.
 */
void chkder(
        Vector<double> &x,
        Vector<double> &fvec,
        Matrix<double> &fjac,
        Vector<double> &xp,
        Vector<double> &fvecp,
        int mode,
        Vector<double> &err
) {
    const double factor = 100.0;
    const double epsmch = std::numeric_limits<double>::epsilon();
    const double eps = sqrt(epsmch);

    const int n = x.size();
    const int m = fjac.rows();

    if (mode == 1) {
        for (int j = 0; j < n; j++) {
            double temp = eps * std::abs(x(j));
            if (temp == 0.0) temp = eps;
            xp(j) = x(j) + temp;
        }
    } else {
        double epsf = factor * epsmch;
        double epslog = log10(eps);
        for (int i = 0; i < m; i++) {
            err(i) = 0.0;
        }
        for (int j = 0; j < n; j++) {
            double temp = std::abs(x(j));
            if (temp == 0.0) temp = 1.0;
            for (int i = 0; i < m; i++) {
                err(i) = err(i) + temp * fjac(i, j);
            }
        }
        for (int i = 0; i < m; i++) {
            double temp = 1.0;
            if (fvec(i) != 0.0 && fvecp(i) != 0.0 &&
                    std::abs(fvecp(i) - fvec(i)) >= epsf * std::abs(fvec(i))) {
                temp = eps * std::abs((fvecp(i) - fvec(i)) / eps - err(i)) / (std::abs(fvec(i)) + std::abs(fvecp(i)));
            }
            err(i) = 1.0;
            if (temp > epsmch && temp < eps) {
                err(i) = (log10(temp) - epslog) / epslog;
            }
            if (temp >= eps) {
                err(i) = 0.0;
            }
        }
    }
}

}
}
#endif //LANRE_OPTIMIZE_CHKDER_HPP
