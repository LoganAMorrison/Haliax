/**
 * Created by Logan Morrison on 4/10/20.
 *
 * Adapted from:
 * Original FORTRAN77 version by Jorge More, Burt Garbow, Ken Hillstrom.
 * C++ version by John Burkardt.
 *
 * Licensing:
 *   This code may freely be copied, modified, and used for any purpose.
 *
 * Reference:
 *  Jorge More, Burton Garbow, Kenneth Hillstrom,
 *  User Guide for MINPACK-1,
 *  Technical Report ANL-80-74,
 *  Argonne National Laboratory, 1980.
 */

#ifndef LANRE_OPTIMIZE_R1UPDT_HPP
#define LANRE_OPTIMIZE_R1UPDT_HPP

#include "lanre/lanre.hpp"

namespace lanre {
namespace optimize {

/**
 * Updates the Q factor after a rank one update of the matrix.
 *
 * Given an M by N lower trapezoidal matrix S, an M-vector U,
 * and an N-vector V, the problem is to determine an
 * orthogonal matrix Q such that:
 *      (S + U*V') * Q
 * is again lower trapezoidal.
 *
 * This function determines q as the product of 2*(n - 1) transformations
 *       gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1)
 * where gv(i), gw(i) are givens rotations in the (i,n) plane
 * which eliminate elements in the i-th and n-th planes,
 * respectively.
 *
 * Q itself is not accumulated, rather the
 * information to recover the gv, gw rotations is returned.
 *
 *
 * @param m Number of rows of s.
 * @param n Number of columns of s. n must not exceed m.
 * @param s On input, the lower trapezoidal matrix s stored by columns.
 *          On output, the lower trapezoidal matrix produced as described above.
 * @param ls A positive value not less than (n*(2*m-n+1))/2
 * @param u Contains the vector U
 * @param v On input v must contain the vector v.
 *          On output v(i) contains the information necessary to
 *          recover the givens rotation gv(i) described above
 * @param w An output array of length m. w(i) contains information
 *          necessary to recover the givens rotation gw(i) described
 *          above
 * @return sing is set true if any of the diagonal elements of the output s are
 *         zero. Otherwise sing is set false.
 */
bool r1updt(
        Vector<double> &s,
        Vector<double> &u,
        Vector<double> &v,
        Vector<double> &w
) {
    double cotan;
    double cs;
    double giant = std::numeric_limits<double>::max();
    int i;
    int j;
    int jj;
    int l;
    int nm1;
    const double p25 = 0.25;
    const double p5 = 0.5;
    double sn;
    bool sing;
    double tan;
    double tau;
    double temp;

    const int m = u.size();
    const int n = v.size();

    jj = (n * (2 * m - n + 1)) / 2 - (m - n);
    //
    //  Move the nontrivial part of the last column of S into W.
    //
    l = jj;
    for (i = n; i <= m; i++) {
        w[i - 1] = s[l - 1];
        l = l + 1;
    }
    //
    //  Rotate the vector V into a multiple of the N-th unit vector
    //  in such a way that a spike is introduced into W.
    //
    nm1 = n - 1;

    for (j = n - 1; 1 <= j; j--) {
        jj = jj - (m - j + 1);
        w[j - 1] = 0.0;

        if (v[j - 1] != 0.0) {
            //
            //  Determine a Givens rotation which eliminates the J-th element of V.
            //
            if (fabs(v[n - 1]) < fabs(v[j - 1])) {
                cotan = v[n - 1] / v[j - 1];
                sn = p5 / sqrt(p25 + p25 * cotan * cotan);
                cs = sn * cotan;
                tau = 1.0;
                if (1.0 < fabs(cs) * giant) {
                    tau = 1.0 / cs;
                }
            } else {
                tan = v[j - 1] / v[n - 1];
                cs = p5 / sqrt(p25 + p25 * tan * tan);
                sn = cs * tan;
                tau = sn;
            }
            //
            //  Apply the transformation to V and store the information
            //  necessary to recover the Givens rotation.
            //
            v[n - 1] = sn * v[j - 1] + cs * v[n - 1];
            v[j - 1] = tau;
            //
            //  Apply the transformation to S and extend the spike in W.
            //
            l = jj;
            for (i = j; i <= m; i++) {
                temp = cs * s[l - 1] - sn * w[i - 1];
                w[i - 1] = sn * s[l - 1] + cs * w[i - 1];
                s[l - 1] = temp;
                l = l + 1;
            }
        }
    }
    //
    //  Add the spike from the rank 1 update to W.
    //
    for (i = 1; i <= m; i++) {
        w[i - 1] = w[i - 1] + v[n - 1] * u[i - 1];
    }
    //
    //  Eliminate the spike.
    //
    sing = false;

    for (j = 1; j <= nm1; j++) {
        //
        //  Determine a Givens rotation which eliminates the
        //  J-th element of the spike.
        //
        if (w[j - 1] != 0.0) {

            if (fabs(s[jj - 1]) < fabs(w[j - 1])) {
                cotan = s[jj - 1] / w[j - 1];
                sn = p5 / sqrt(p25 + p25 * cotan * cotan);
                cs = sn * cotan;
                tau = 1.0;
                if (1.0 < fabs(cs) * giant) {
                    tau = 1.0 / cs;
                }
            } else {
                tan = w[j - 1] / s[jj - 1];
                cs = p5 / sqrt(p25 + p25 * tan * tan);
                sn = cs * tan;
                tau = sn;
            }
            //
            //  Apply the transformation to s and reduce the spike in w.
            //
            l = jj;

            for (i = j; i <= m; i++) {
                temp = cs * s[l - 1] + sn * w[i - 1];
                w[i - 1] = -sn * s[l - 1] + cs * w[i - 1];
                s[l - 1] = temp;
                l = l + 1;
            }
            //
            //  Store the information necessary to recover the givens rotation.
            //
            w[j - 1] = tau;
        }
        //
        //  Test for zero diagonal elements in the output s.
        //
        if (s[jj - 1] == 0.0) {
            sing = true;
        }
        jj = jj + (m - j + 1);
    }
    //
    //  Move W back into the last column of the output S.
    //
    l = jj;
    for (i = n; i <= m; i++) {
        s[l - 1] = w[i - 1];
        l = l + 1;
    }
    if (s[jj - 1] == 0.0) {
        sing = true;
    }
    return sing;
}

}
}

#endif //LANRE_OPTIMIZE_R1UPDT_HPP
