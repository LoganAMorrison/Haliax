// Created by Logan Morrison on 4/4/20.
//
// Adapted from the FORTRAN pacakge Dierckx written by:
//      p.dierckx
//          dept. computer science, k.u. leuven
//          celestijnenlaan 200a, b-3001 heverlee, belgium.
//          e-mail : Paul.Dierckx@cs.kuleuven.ac.be
//

#ifndef LANRE_INTERPOLATE_FPCHEC_HPP
#define LANRE_INTERPOLATE_FPCHEC_HPP

namespace lanre {
namespace interpolate {
namespace dierckx {

/**
 * fpchec verifies the number and the position of the knots
 * t(j),j=1,2,...,n of a spline of degree k, in relation to the number
 * and the position of the data points x(i),i=1,2,...,m. if all of the
 * following conditions are fulfilled, the error parameter ier is set
 * to zero. if one of the conditions is violated ier is set to ten.
 *
 * 1) k+1 <= n-k-1 <= m
 * 2) t(1) <= t(2) <= ... <= t(k+1)
 *    t(n-k) <= t(n-k+1) <= ... <= t(n)
 * 3) t(k+1) < t(k+2) < ... < t(n-k)
 * 4) t(k+1) <= x(i) <= t(n-k)
 * 5) the conditions specified by schoenberg and whitney must hold
 *    for at least one subset of data points, i.e. there must be a
 *    subset of data points y(j) such that
 *    t(j) < y(j) < t(j+k+1), j=1,2,...,n-k-1
 */
void fpchec(const double *x, int m, double *t, int n, int k, int &ier) {
    int i, j, k1, k2, l, nk1, nk2, nk3;
    double tj, tl;

    k1 = k + 1;
    k2 = k1 + 1;
    nk1 = n - k1;
    nk2 = nk1 + 1;
    ier = 10;

    //  check condition no 1:
    // k+1 <= n-k-1 <= m
    if (nk1 < k1 || nk1 > m) return;
    // check condition no 2:
    // t(1) <= t(2) <= ... <= t(k+1)
    // t(n-k) <= t(n-k+1) <= ... <= t(n)
    j = n;
    for (i = 1; i <= k; i++) {
        if (t[i - 1] > t[i]) return;
        if (t[j - 1] < t[j - 2]) return;
        j = j - 1;
    }
    // check condition no 3:
    // t(k+1) < t(k+2) < ... < t(n-k)
    for (i = k2; i <= nk2; i++) {
        if (t[i - 1] <= t[i - 2]) return;
    }
    // check condition no 4:
    // t(k+1) <= x(i) <= t(n-k)
    if (x[0] < t[k1 - 1] || x[m - 1] > t[nk2 - 1]) return;
    // check condition no 5:
    // the conditions specified by schoenberg and whitney must hold
    // for at least one subset of data points, i.e. there must be a
    // subset of data points y(j) such that
    // t(j) < y(j) < t(j+k+1), j=1,2,...,n-k-1
    if (x[0] >= t[k2 - 1] || x[m - 1] <= t[nk1 - 1]) return;
    i = 1;
    l = k2;
    nk3 = nk1 - 1;
    if (nk3 < 2) goto _70;
    for (j = 2; j <= nk3; j++) {
        tj = t[j - 1];
        l = l + 1;
        tl = t[l - 1];
        _40:
        i = i + 1;
        if (i >= m) return;
        if (x[i - 1] <= tj) goto _40;
        if (x[i - 1] >= tl) return;
    }

    _70:
    ier = 0;
}


} // namespace dierckx
} // namespace interpolate
} // namespace lanre

#endif //LANRE_INTERPOLATE_FPCHEC_HPP