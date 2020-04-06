//
// Created by Logan Morrison on 4/4/20.
//

#ifndef LANRE_INTERPOLATE_SPLDER_HPP
#define LANRE_INTERPOLATE_SPLDER_HPP

#include "lanre/interpolate/fpbspl.hpp"
#include <cmath>

namespace lanre {
namespace interpolate {
namespace dierckx {

void splder(
        const double *t,
        int n,
        const double *c,
        int k,
        int nu,
        const double *x,
        double *y,
        int m,
        int e,
        double *wrk,
        int &ier
) {
    int i, j, kk, k1, k2, l, ll, l1, l2, nk1, nk2, nn;
    double ak, arg, fac, sp, tb, te;
    int k3;
    double h[6];

    // before starting computations a data check is made. if the input data
    // are invalid control is immediately repassed to the calling program.
    ier = 10;
    if (nu < 0 || nu > k) goto _200;
    if (m < 1) goto _200;
    ier = 0;
    // fetch tb and te, the boundaries of the approximation interval.
    k1 = k + 1;
    k3 = k1 + 1;
    nk1 = n - k1;
    tb = t[k1 - 1];
    te = t[nk1];
    // the derivative of order nu of a spline of degree k is a spline of
    // degree k-nu,the b-spline coefficients wrk(i) of which can be found
    // using the recurrence scheme of de boor.
    l = 1;
    kk = k;
    nn = n;
    for (i = 1; i <= nk1; i++) {
        wrk[i - 1] = c[i - 1];
    }
    if (nu == 0) goto _100;
    nk2 = nk1;
    for (j = 1; j <= nu; j++) {
        ak = kk;
        nk2 = nk2 - 1;
        l1 = l;
        for (i = 1; i <= nk2; i++) {
            l1 = l1 + 1;
            l2 = l1 + kk;
            fac = t[l2 - 1] - t[l1 - 1];
            if (fac <= 0.0) continue;
            wrk[i - 1] = ak * (wrk[i] - wrk[i - 1]) / fac;
        }
        l = l + 1;
        kk = kk - 1;
    }
    if (kk != 0) goto _100;
    // if nu=k the derivative is a piecewise constant function
    j = 1;
    for (i = 1; i <= m; i++) {
        arg = x[i - 1];
        // check if arg is in the support
        if (arg < tb || arg > te) {
            if (e == 0) {
                goto _65;
            } else if (e == 1) {
                y[i - 1] = 0.0;
                goto _90;
            } else if (e == 2) {
                ier = 1;
                goto _200;
            }
        }
        //  search for knot interval t(l) <= arg < t(l+1)
        _65:
        if (arg >= t[l - 1] || l + 1 == k3) goto _70;
        l1 = l;
        l = l - 1;
        j = j - 1;
        goto _65;

        _70:
        if (arg < t[l] || l == nk1) goto _80;
        l = l + 1;
        j = j + 1;
        goto _70;

        _80:
        y[i - 1] = wrk[j - 1];
        _90:;
    }
    goto _200;

    _100:
    l = k1;
    l1 = l + 1;
    k2 = k1 - nu;
    //  main loop for the different points.
    for (i = 1; i <= m; i++) {
        // fetch a new x-value arg.
        arg = x[i - 1];
        // check if arg is in the support
        if (arg < tb || arg > te) {
            if (e == 0) {
                goto _135;
            } else if (e == 1) {
                y[i - 1] = 0.0;
                goto _180;
            } else if (e == 2) {
                ier = 1;
                goto _200;
            }
        }
        // search for knot interval t(l) <= arg < t(l+1)
        _135:
        if (arg >= t[l - 1] || l1 == k3) goto _140;
        l1 = l;
        l = l - 1;
        goto _135;

        _140:
        if (arg < t[l1 - 1] || l == nk1) goto _150;
        l = l1;
        l1 = l + 1;
        goto _140;

        // evaluate the non-zero b-splines of degree k-nu at arg.
        _150:
        fpbspl(t, n, kk, arg, l, h);
        // find the value of the derivative at x=arg.
        sp = 0.0;
        ll = l - k1;
        for (j = 1; j <= k2; j++) {
            ll = ll + 1;
            sp = sp + wrk[ll - 1] * h[j - 1];
        }
        y[i - 1] = sp;
        _180:;
    }
    _200:;
}

} // namespace dierckx
} // namespace interpolate
} // namespace lanre

#endif //LANRE_INTERPOLATE_SPLDER_HPP
