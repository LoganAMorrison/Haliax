//
// Created by Logan Morrison on 4/10/20.
//

#ifndef LANRE_OPTIMIZE_R1MPYQ_HPP
#define LANRE_OPTIMIZE_R1MPYQ_HPP

#include "lanre/lanre.hpp"

namespace lanre {
namespace optimize {

/**
 *  Given an m by n matrix a, this subroutine computes a*q where
 *  q is the product of 2*(n - 1) transformations
 *      gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1)
 *  and gv(i), gw(i) are givens rotations in the (i,n) plane which
 *  eliminate elements in the i-th and n-th planes, respectively.
 *  q itself is not given, rather the information to recover the
 *  gv, gw rotations is supplied.
 *
 * @param a is an m by n array.
 *          on input a must contain the matrix
 *          to be postmultiplied by the orthogonal matrix q
 *          described above.
 *          on output a*q has replaced a.
 * @param v is an input array of length n. v(i) must contain the
 *          information necessary to recover the givens rotation gv(i)
 *          described above.
 * @param w is an input array of length n. w(i) must contain the
 *          information necessary to recover the givens rotation gw(i)
 *          described above.
 */
template<int n, int m>
void r1mpyq(
        Eigen::Matrix<double, n, m> &a,
        Vector<double> &v,
        Vector<double> &w
) {
    for (int j = n - 2; 0 <= j; j--) {
        double c, s;
        if (1.0 < fabs(v(j))) {
            c = 1.0 / v[j];
            s = sqrt(1.0 - c * c);
        } else {
            s = v(j);
            c = sqrt(1.0 - s * s);
        }
        for (int i = 0; i < m; i++) {
            double temp = c * a(i, j) - s * a(i, n - 1);
            a(i, n - 1) = s * a(i, j) + c * a(i, n - 1);
            a(i, j) = temp;
        }
    }
    //
    //  Apply the second set of Givens rotations to A.
    //
    for (int j = 0; j < n - 1; j++) {
        double s, c;
        if (1.0 < fabs(w[j])) {
            c = 1.0 / w[j];
            s = sqrt(1.0 - c * c);
        } else {
            s = w[j];
            c = sqrt(1.0 - s * s);
        }
        for (int i = 0; i < m; i++) {
            double temp = c * a(i, j) + s * a(i, n - 1);
            a(i, n - 1) = -s * a(i, j) + c * a(i, n - 1);
            a(i, j) = temp;
        }
    }

}

}
}

#endif //LANRE_OPTIMIZE_R1MPYQ_HPP
