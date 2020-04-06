//
// Created by Logan Morrison on 4/4/20.
//

#ifndef LANRE_INTERPOLATE_FPINTB_HPP
#define LANRE_INTERPOLATE_FPINTB_HPP

namespace lanre {
namespace interpolate {
namespace dierckx {

/**
 * calculates integrals of the normalized b-splines
 * nj,k+1(x) of degree k, defined on the set of knots t(j),j=1,2,...n.
 * it makes use of the formulae of gaffney for the calculation of
 * indefinite integrals of b-splines.
 *
 * @param t Length n, containing the position of the knots.
 * @param n integer value, giving the number of knots.
 * @param bint (out) array,length nk1, containing the integrals of the b-splines.
 * @param nk1 integer value, giving the number of b-splines of degree k,
 *            defined on the set of knots ,i.e. nk1 = n-k-1.
 * @param x lower bound
 * @param y upper bound
 */
void fpintb(
        const double *t,
        int n,
        double *bint,
        int nk1,
        double x,
        double y
) {
    int i, ia, ib, it, j, j1, k, k1, l, li, lj, lk, l0, min;
    double a, ak, arg, b, f, one;
    double aint[6], h[6], h1[6];

    // initialization.
    one = 1.0;
    k1 = n - nk1;
    ak = k1;
    k = k1 - 1;
    for (i = 1; i <= nk1; i++) {
        bint[i - 1] = 0.0;
    }
    // the integration limits are arranged in increasing order.
    a = x;
    b = y;
    min = 0;
    if (a < b) goto _30;
    if (a == b) return;
    goto _20;

    _20:
    a = y;
    b = x;
    min = 1;

    _30:
    if (a < t[k1 - 1]) a = t[k1 - 1];
    if (b > t[nk1]) b = t[nk1];
    if (a > b) return;

    // using the expression of gaffney for the indefinite integral of a
    // b-spline we find that
    // bint(j) = (t(j+k+1)-t(j))*(res(j,b)-res(j,a))/(k+1)
    //     where for t(l) <= x < t(l+1)
    //     res(j,x) = 0, j=1,2,...,l-k-1
    //              = 1, j=l+1,l+2,...,nk1
    //              = aint(j+k-l+1), j=l-k,l-k+1,...,l
    //                = sumi((x-t(j+i))*nj+i,k+1-i(x)/(t(j+k+1)-t(j+i)))
    //                  i=0,1,...,k
    l = k1;
    l0 = l + 1;
    arg = a;
    for (it = 1; it <= 2; it++) {
        //  search for the knot interval t(l) <= arg < t(l+1).
        _40:
        if (arg < t[l0 - 1] || l == nk1) goto _50;
        l = l0;
        l0 = l + 1;
        goto _40;
        // calculation of aint(j), j=1,2,...,k+1.
        // initialization.
        _50:
        for (j = 1; j <= k1; j++) {
            aint[j - 1] = 0.0;
        }
        aint[0] = (arg - t[l - 1]) / (t[l] - t[l - 1]);
        h1[0] = one;
        for (j = 1; j <= k; j++) {
            // evaluation of the non-zero b-splines of degree j at arg,i.e.
            // h(i+1) = nl-j+i,j(arg), i=0,1,...,j.
            h[0] = 0.0;
            for (i = 1; i <= j; i++) {
                li = l + i;
                lj = li - j;
                f = h1[i - 1] / (t[li - 1] - t[lj - 1]);
                h[i - 1] = h[i - 1] + f * (t[li - 1] - arg);
                h[i] = f * (arg - t[lj - 1]);
            }
            // updating of the integrals aint.
            j1 = j + 1;
            for (i = 1; i <= j1; i++) {
                li = l + i;
                lj = li - j1;
                aint[i - 1] = aint[i - 1] + h[i - 1] * (arg - t[lj - 1]) / (t[li - 1] - t[lj - 1]);
                h1[i - 1] = h[i - 1];
            }
        }
        if (it == 2) goto _100;
        //  updating of the integrals bint
        lk = l - k;
        ia = lk;
        for (i = 1; i <= k1; i++) {
            bint[lk - 1] = -aint[i - 1];
            lk = lk + 1;
        }
        arg = b;
    }
    //  updating of the integrals bint.
    _100:
    lk = l - k;
    ib = lk - 1;
    for (i = 1; i <= k1; i++) {
        bint[lk - 1] = bint[lk - 1] + aint[i - 1];
        lk = lk + 1;
    }
    if (ib < ia) goto _130;
    for (i = ia; i <= ib; i++) {
        bint[i - 1] = bint[i - 1] + one;
    }
    // the scaling factors are taken into account.
    _130:
    f = one / ak;
    for (i = 1; i <= nk1; i++) {
        j = i + k1;
        bint[i - 1] = bint[i - 1] * (t[j - 1] - t[i - 1]) * f;
    }
    // the order of the integration limits is taken into account.
    if (min == 0) return;
    for (i = 1; i <= nk1; i++) {
        bint[i - 1] = -bint[i - 1];
    }
}

} // namespace dierckx
} // namespace interpolate
} // namespace lanre

#endif //LANRE_INTERPOLATE_FPINTB_HPP
