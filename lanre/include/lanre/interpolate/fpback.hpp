// Created by Logan Morrison on 4/4/20.
//
// Adapted from the FORTRAN pacakge Dierckx written by:
//      p.dierckx
//          dept. computer science, k.u. leuven
//          celestijnenlaan 200a, b-3001 heverlee, belgium.
//          e-mail : Paul.Dierckx@cs.kuleuven.ac.be
//

#ifndef LANRE_INTERPOLATE_FPBACK_HPP
#define LANRE_INTERPOLATE_FPBACK_HPP

namespace lanre {
namespace interpolate {
namespace dierckx {

/**
 * calculates the solution of the system of
 * equations a*c = z with a a n x n upper triangular matrix
 * of bandwidth k.
 */
void fpback(double **a, const double *z, int n, int k, double *c, int nest) {
    double store;
    int i, i1, j, k1, l, m;

    k1 = k - 1;
    c[n - 1] = z[n - 1] / a[n - 1][0];
    i = n - 1;
    if (i == 0) return;

    for (j = 2; j <= n; j++) {
        store = z[i - 1];
        i1 = k1;
        if (j <= k1) i1 = j - 1;
        m = i;
        for (l = 1; l <= i1; l++) {
            m = m + 1;
            store = store - c[m - 1] * a[i - 1][l + 1 - 1];
        }
        c[i - 1] = store / a[i - 1][0];
        i = i - 1;
    }
}


} // namespace dierckx
} // namespace interpolate
} // namespace lanre

#endif //LANRE_INTERPOLATE_FPBACK_HPP