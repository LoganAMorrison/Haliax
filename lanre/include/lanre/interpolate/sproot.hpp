//
// Created by Logan Morrison on 4/4/20.
//

#ifndef LANRE_INTERPOLATE_SPROOT_HPP
#define LANRE_INTERPOLATE_SPROOT_HPP

#include "lanre/interpolate/fpcuro.hpp"

namespace lanre {
namespace interpolate {
namespace dierckx {

void sproot(
        const double *t,
        int n,
        const double *c,
        double *zero,
        const int &mest,
        int &m,
        int &ier
) {
    int i, j, j1, l, n4;
    double ah, a0 = 0, a1 = 0, a2 = 0, a3 = 0, bh, b0 = 0, b1 = 0, c1 = 0, c2 = 0, c3, c4, c5, d4, d5, h1, h2,
            three, two, t1, t2, t3, t4, t5, zz;
    bool z0, z1, z2, z3, z4, nz0, nz1, nz2, nz3, nz4;
    double y[3];

    // set some constants
    two = 2.0;
    three = 3.0;
    // before starting computations a data check is made. if the input data
    // are invalid, control is immediately repassed to the calling program.
    n4 = n - 4;
    ier = 10;
    if (n < 8) goto _800;
    j = n;
    for (i = 1; i <= 3; i++) {
        if (t[i - 1] > t[i + 1 - 1]) goto _800;
        if (t[j - 1] < t[j - 1 - 1]) goto _800;
        j = j - 1;
    }
    for (i = 4; i <= n4; i++) {
        if (t[i - 1] >= t[i + 1 - 1]) goto _800;
    }
    // the problem considered reduces to finding the zeros of the cubic
    // polynomials pl(x) which define the cubic spline in each knot
    // interval t(l)<=x<=t(l+1). a zero of pl(x) is also a zero of s(x) on
    // the condition that it belongs to the knot interval.
    // the cubic polynomial pl(x) is determined by computing s(t(l)),
    // s'(t(l)),s(t(l+1)) and s'(t(l+1)). in fact we only have to compute
    // s(t(l+1)) and s'(t(l+1)); because of the continuity conditions of
    // splines and their derivatives, the value of s(t(l)) and s'(t(l))
    // is already known from the foregoing knot interval.
    ier = 0;
    // evaluate some constants for the first knot interval
    h1 = t[4 - 1] - t[3 - 1];
    h2 = t[5 - 1] - t[4 - 1];
    t1 = t[4 - 1] - t[2 - 1];
    t2 = t[5 - 1] - t[3 - 1];
    t3 = t[6 - 1] - t[4 - 1];
    t4 = t[5 - 1] - t[2 - 1];
    t5 = t[6 - 1] - t[3 - 1];
    // calculate a0 = s(t(4)) and ah = s'(t(4)).
    c1 = c[1 - 1];
    c2 = c[2 - 1];
    c3 = c[3 - 1];
    c4 = (c2 - c1) / t4;
    c5 = (c3 - c2) / t5;
    d4 = (h2 * c1 + t1 * c2) / t4;
    d5 = (t3 * c2 + h1 * c3) / t5;
    a0 = (h2 * d4 + h1 * d5) / t2;
    ah = three * (h2 * c4 + h1 * c5) / t2;
    z1 = true;
    if (ah < 0.0) z1 = false;
    nz1 = !z1;
    m = 0;
    // main loop for the different knot intervals.
    for (l = 4; l <= n4; l++) {
        // evaluate some constants for the knot interval t(l) <= x <= t(l+1).
        h1 = h2;
        h2 = t[l + 2 - 1] - t[l + 1 - 1];
        t1 = t2;
        t2 = t3;
        t3 = t[l + 3 - 1] - t[l + 1 - 1];
        t4 = t5;
        t5 = t[l + 3 - 1] - t[l - 1];
        // find a0 = s(t(l)), ah = s'(t(l)), b0 = s(t(l+1)) and bh = s'(t(l+1)).
        c1 = c2;
        c2 = c3;
        c3 = c[l - 1];
        c4 = c5;
        c5 = (c3 - c2) / t5;
        d4 = (h2 * c1 + t1 * c2) / t4;
        d5 = (h1 * c3 + t3 * c2) / t5;
        b0 = (h2 * d4 + h1 * d5) / t2;
        bh = three * (h2 * c4 + h1 * c5) / t2;
        // test whether pl(x) could have a zero in the range
        // t(l) <= x <= t(l+1).
        z3 = true;
        if (b1 < 0.00) z3 = false;
        nz3 = !z3;
        if (a0 * b0 <= 0.0) goto _100;
        z0 = true;
        if (a0 < 0.0) z0 = false;
        nz0 = !z0;
        z2 = true;
        if (a2 < 0.) z2 = false;
        nz2 = !z2;
        z4 = true;
        if ((3.0 * a3 + a2) < 0.0) z4 = false;
        nz4 = !z4;
        if (!(((z0 && ((nz1 && (z3 || (z2 && nz4))) || (nz2 && z3 && z4))) ||
                (nz0 && ((z1 && (nz3 || (nz2 && z4))) || (z2 && nz3 && nz4)))))) {
            goto _200;
        }
        // find the zeros of ql(y).
        _100:
        fpcuro(a3, a2, a1, a0, y, j);
        if (j == 0) goto _200;
        // find which zeros of pl(x) are zeros of s(x).
        for (i = 1; i <= j; i++) {
            if (y[i - 1] < 0.0 || y[i - 1] > 1.0) continue;
            //  test whether the number of zeros of s(x) exceeds mest.
            if (m >= mest) goto _700;
            m = m + 1;
            zero[m - 1] = t[l - 1] + h1 * y[i - 1];
        }
        _200:
        a0 = b0;
        ah = bh;
        z1 = z3;
        nz1 = nz3;
    }
    // the zeros of s(x) are arranged in increasing order.
    if (m < 2) goto _800;
    for (i = 2; i <= m; i++) {
        j = i;
        _350:
        j1 = j - 1;
        if (j1 == 0) continue;
        if (zero[j - 1] >= zero[j1 - 1]) continue;
        zz = zero[j - 1];
        zero[j - 1] = zero[j1 - 1];
        zero[j1 - 1] = zz;
        j = j1;
        goto _350;
    }
    j = m;
    m = 1;
    for (i = 2; i <= j; i++) {
        if (zero[i - 1] == zero[m - 1]) continue;
        m = m + 1;
        zero[m - 1] = zero[i - 1];
    }
    goto _800;

    _700:
    ier = 1;

    _800:;
}

} // namespace dierckx
} // namespace interpolate
} // namespace lanre

#endif //LANRE_INTERPOLATE_SPROOT_HPP
