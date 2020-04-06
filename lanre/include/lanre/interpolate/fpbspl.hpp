// Created by Logan Morrison on 4/4/20.
//
// Adapted from the FORTRAN pacakge Dierckx written by:
//      p.dierckx
//          dept. computer science, k.u. leuven
//          celestijnenlaan 200a, b-3001 heverlee, belgium.
//          e-mail : Paul.Dierckx@cs.kuleuven.ac.be
//

#ifndef LANRE_INTERPOLATE_FPBSPL_HPP
#define LANRE_INTERPOLATE_FPBSPL_HPP

namespace lanre {
namespace interpolate {
namespace dierckx {

/**
 * evaluates the (k+1) non-zero b-splines of
 * degree k at t(l) <= x < t(l+1) using the stable recurrence
 * relation of de boor and cox.
 */
void fpbspl(const double *t, int &n, int &k, double &x, int &l, double h[20]) {
    double f, one;
    int i, j, li, lj;
    double hh[19];

    one = 1.0;
    h[0] = one;
    for (j = 1; j <= k; j++) {
        for (i = 1; i <= j; i++) {
            hh[i - 1] = h[i - 1];
        }
        h[0] = 0.0;
        for (i = 1; i <= j; i++) {
            li = l + i;
            lj = li - j;
            if (t[li - 1] != t[lj - 1]) goto _15;
            h[i + 1 - 1] = 0.0;
            continue;
            _15:
            f = hh[i - 1] / (t[li - 1] - t[lj - 1]);
            h[i - 1] = h[i - 1] + f * (t[li - 1] - x);
            h[i + 1 - 1] = f * (x - t[lj - 1]);
        }
    }
}


} // namespace dierckx
} // namespace interpolate
} // namespace lanre

#endif //LANRE_INTERPOLATE_FPBSPL_HPP
