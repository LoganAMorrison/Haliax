//
// Created by Logan Morrison on 4/10/20.
//

#ifndef LANRE_OPTIMIZE_FDJAC1_HPP
#define LANRE_OPTIMIZE_FDJAC1_HPP

#include "lanre/lanre.hpp"

namespace lanre {
namespace optimize {

/**
 * Computes a forward-difference approximation
 * to the n by n jacobian matrix associated with a specified
 * problem of n functions in n variables.
 *
 * If the jacobian has a banded form, then function evaluations are saved
 * by only approximating the nonzero terms.
 *
 * @tparam Function Function type.
 * @param fcn Function specifying the root equations (inplace).
 * @param n The number of functions and variables.
 * @param x The evaluation point.
 * @param fvec The functions evaluated at x.
 * @param fjac The approximate jacobian matrix at x.
 * @param ldfjac is a positive integer input variable not less than n
 *               which specifies the leading dimension of the array fjac
 * @param iflag is an integer variable which can be used to terminate
 *              the execution of fdjac1. see description of fcn.
 * @param ml is a nonnegative integer input variable which specifies
 *           the number of subdiagonals within the band of the
 *           jacobian matrix. if the jacobian is not banded, set
 *           ml to at least n - 1.
 * @param mu is a nonnegative integer input variable which specifies
 *           the number of superdiagonals within the band of the
 *           jacobian matrix. if the jacobian is not banded, set
 *           mu to at least n - 1.
 * @param epsfcn is an input variable used in determining a suitable
 *               step length for the forward-difference approximation. this
 *               approximation assumes that the relative errors in the
 *               functions are of the order of epsfcn. if epsfcn is less
 *               than the machine precision, it is assumed that the relative
 *               errors in the functions are of the order of the machine
 *               precision.
 * @param wa1 Work array of length n.
 * @param wa2 Work array of length n. if ml + mu + 1 is at
 *            least n, then the jacobian is considered dense, and wa2 is
 *            not referenced.
 */
template<class Function>
void fdjac1(Function fcn,
            Vector<double> &x,
            Vector<double> &fvec,
            Matrix<double> &fjac,
            int &iflag,
            int ml,
            int mu,
            double epsfcn,
            Vector<double> &wa1,
            Vector<double> &wa2
) {
    const double epsmch = std::numeric_limits<double>::epsilon();
    const double eps = sqrt(std::max(epsfcn, epsmch));
    const int msum = ml + mu + 1;
    const int n = x.size();

    //
    //  Computation of dense approximate jacobian.
    //
    if (n <= msum) {
        for (int j = 0; j < n; j++) {
            double temp = x[j];
            double h = eps * fabs(temp);
            if (h == 0.0) {
                h = eps;
            }
            x[j] = temp + h;
            fcn(x, wa1, iflag);
            if (iflag < 0) {
                break;
            }
            x[j] = temp;
            for (int i = 0; i < n; i++) {
                fjac(i, j) = (wa1[i] - fvec[i]) / h;
            }
        }
    }
//
//  Computation of a banded approximate jacobian.
//
    else {
        for (int k = 0; k < msum; k++) {
            for (int j = k; j < n; j = j + msum) {
                wa2[j] = x[j];
                double h = eps * fabs(wa2[j]);
                if (h == 0.0) {
                    h = eps;
                }
                x[j] = wa2[j] + h;
            }
            fcn(x, wa1, iflag);
            if (iflag < 0) {
                break;
            }
            for (int j = k; j < n; j = j + msum) {
                x[j] = wa2[j];
                double h = eps * fabs(wa2[j]);
                if (h == 0.0) {
                    h = eps;
                }
                for (int i = 0; i < n; i++) {
                    if (j - mu <= i && i <= j + ml) {
                        fjac(i, j) = (wa1[i] - fvec[i]) / h;
                    } else {
                        fjac(i, j) = 0.0;
                    }
                }
            }
        }
    }
}
}
}

#endif //LANRE_OPTIMIZE_FDJAC1_HPP
