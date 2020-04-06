//
// Created by Logan Morrison on 4/4/20.
//


#ifndef LANRE_INTERPOLATE_SPLEV_HPP
#define LANRE_INTERPOLATE_SPLEV_HPP

#include "lanre/interpolate/fpbspl.hpp"
#include <cmath>

namespace lanre {
namespace interpolate {
namespace dierckx {

/**
 * Evaluates in a number of points x(i),i=1,2,...,m
 * a spline s(x) of degree k, given in its b-spline representation.
 *
 * @param t Array,length n, which contains the position of the knots.
 * @param n Integer, giving the total number of knots of s(x).
 * @param c Array,length n, which contains the b-spline coefficients.
 * @param k Integer, giving the degree of s(x).
 * @param x Array,length m, which contains the points where s(x) must
 *          be evaluated.
 * @param y Array,length m, giving the value of s(x) at the different
 *          points.
 * @param m Integer, giving the number of points where s(x) must be
 *          evaluated.
 * @param e Integer, if 0 the spline is extrapolated from the end
 *          spans for points not in the support, if 1 the spline
 *          evaluates to zero for those points, if 2 ier is set to
 *          1 and the subroutine returns, and if 3 the spline evaluates
 *          to the value of the nearest boundary point.
 * @param ier Error flag
 *               ier = 0 : normal return
 *               ier = 1 : argument out of bounds and e == 2
 *               ier =10 : invalid input data (see restrictions)
 */
void splev(const double *t, int &n, const double *c, int &k, const double *x,
           double *y, int m, int e, int &ier) {
    int i, j, k1, l, ll, l1, nk1;
    int k2;
    double arg, sp, tb, te;
    double h[20];

    // before starting computations a data check is made. if the input data
    // are invalid control is immediately repassed to the calling program.
    ier = 10;
    if (m < 1) return;
    ier = 0;
    // fetch tb and te, the boundaries of the approximation interval.
    k1 = k + 1;
    k2 = k1 + 1;

    nk1 = n - k1;
    tb = t[k1 - 1];
    te = t[nk1];
    l = k1;
    l1 = l + 1;

    // main loop for the different point
    for (i = 1; i <= m; i++) {
        // fetch a new x-value arg.
        arg = x[i - 1];
        if (arg < tb || arg > te) {
            if (e == 0) {
                goto _35;
            } else if (e == 1) {
                y[i - 1] = 0;
                goto _80;
            } else if (e == 2) {
                ier = 1;
                return;
            } else if (e == 3) {
                if (arg < tb) {
                    arg = tb;
                } else {
                    arg = te;
                }
            }
        }
        // search for knot interval t(l) <= arg < t(l+1)
        _35:
        if (arg >= t[l - 1] || l1 == k2) goto _40;
        l1 = l;
        l = l - 1;
        goto _35;

        _40:
        if (arg < t[l1 - 1] || l == nk1) goto _50;
        l = l1;
        l1 = l + 1;
        goto _40;

        //evaluate the non-zero b-splines at arg.
        _50:
        fpbspl(t, n, k, arg, l, h);
        //  find the value of s(x) at x=arg.
        sp = 0.0;
        ll = l - k1;
        for (j = 1; j <= k1; j++) {
            ll = ll + 1;
            sp = sp + c[ll - 1] * h[j - 1];
        }
        y[i - 1] = sp;
        _80:;
    }
}

double splev(const double *t, int &n, const double *c, int &k, double x, int e, int &ier) {
    int i, j, k1, l, ll, l1, nk1;
    int k2;
    double arg, sp, tb, te, y = 0;
    double h[20];

    // before starting computations a data check is made. if the input data
    // are invalid control is immediately repassed to the calling program.
    ier = 0;
    // fetch tb and te, the boundaries of the approximation interval.
    k1 = k + 1;
    k2 = k1 + 1;

    nk1 = n - k1;
    tb = t[k1 - 1];
    te = t[nk1 + 1 - 1];
    l = k1;
    l1 = l + 1;

    // main loop for the different point
    // fetch a new x-value arg.
    arg = x;
    if (arg < tb || arg > te) {
        if (e == 0) {
            goto _35;
        } else if (e == 1) {
            y = 0;
            goto _80;
        } else if (e == 2) {
            ier = 1;
            return y;
        } else if (e == 3) {
            if (arg < tb) {
                arg = tb;
            } else {
                arg = te;
            }
        }
    }
    // search for knot interval t(l) <= arg < t(l+1)
    _35:
    if (arg >= t[l - 1] || l1 == k2) goto _40;
    l1 = l;
    l = l - 1;
    goto _35;

    _40:
    if (arg < t[l1 - 1] || l == nk1) goto _50;
    l = l1;
    l1 = l + 1;
    goto _40;

    //evaluate the non-zero b-splines at arg.
    _50:
    fpbspl(t, n, k, arg, l, h);
    //  find the value of s(x) at x=arg.
    sp = 0.0;
    ll = l - k1;
    for (j = 1; j <= k1; j++) {
        ll = ll + 1;
        sp = sp + c[ll - 1] * h[j - 1];
    }
    y = sp;
    _80:
    return y;
}

} // namespace dierckx
} // namespace interpolate
} // namespace lanre

#endif //LANRE_INTERPOLATE_SPLEV_HPP
