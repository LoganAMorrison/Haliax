// Created by Logan Morrison on 4/4/20.
//
// Adapted from the FORTRAN pacakge Dierckx written by:
//      p.dierckx
//          dept. computer science, k.u. leuven
//          celestijnenlaan 200a, b-3001 heverlee, belgium.
//          e-mail : Paul.Dierckx@cs.kuleuven.ac.be
//

#ifndef LANRE_INTERPOLATE_FPADER_HPP
#define LANRE_INTERPOLATE_FPADER_HPP

namespace lanre {
namespace interpolate {
namespace dierckx {

/**
 * fpader calculates the derivatives
 * d(j) = s(x) , j=1,2,...,k1
 * of a spline of order k1 at the point t(l)<=x<t(l+1), using the
 * stable recurrence scheme of de boor
 */
void fpader(
        const double *t,
        int n,
        const double *c,
        int k1,
        double x,
        int l,
        double *d
) {
    int i, ik, j, jj, j1, j2, ki, kj, li, lj, lk;
    double ak, fac, one;
    double h[20];

    one = 0.1e+01;
    lk = l - k1;
    for (i = 1; i <= k1; i++) {
        ik = i + lk;
        h[i - 1] = c[ik - 1];
    }
    kj = k1;
    fac = one;
    for (j = 1; j <= k1; j++) {
        ki = kj;
        j1 = j + 1;
        if (j == 1) goto _300;
        i = k1;
        for (jj = j; jj <= k1; jj++) {
            li = i + lk;
            lj = li + kj;
            h[i - 1] = (h[i - 1] - h[i - 1 - 1]) / (t[lj - 1] - t[li - 1]);
            i = i - 1;
        }
        _300:
        for (i = j; i <= k1; i++) {
            d[i - 1] = h[i - 1];
        }
        if (j == k1) goto _600;
        for (jj = j1; jj <= k1; jj++) {
            ki = ki - 1;
            i = k1;
            for (j2 = jj; j2 <= k1; j2++) {
                li = i + lk;
                lj = li + ki;
                d[i - 1] = ((x - t[li - 1]) * d[i - 1] + (t[lj - 1] - x) * d[i - 1 - 1]) / (t[lj - 1] - t[li - 1]);
                i = i - 1;
            }
        }
        _600:
        d[j - 1] = d[k1 - 1] * fac;
        ak = k1 - j;
        fac = fac * ak;
        kj = kj - 1;
    }
}

} // namespace dierckx
} // namespace interpolate
} // namespace lanre

#endif //LANRE_INTERPOLATE_FPADER_HPP