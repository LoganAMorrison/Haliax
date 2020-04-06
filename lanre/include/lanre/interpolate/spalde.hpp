//
// Created by Logan Morrison on 4/4/20.
//

#ifndef LANRE_INTERPOLATE_SPALDE_HPP
#define LANRE_INTERPOLATE_SPALDE_HPP

#include "lanre/interpolate/fpader.hpp"

namespace lanre {
namespace interpolate {
namespace dierckx {

/**
 * evaluates at a point x all the derivatives
 * (j-1)
 * d(j) = s     (x) , j=1,2,...,k1
 * of a spline s(x) of order k1 (degree k=k1-1), given in its b-spline
 * representation.
 * @param t array,length n, which contains the position of the knots.
 * @param n integer, giving the total number of knots of s(x).
 * @param c array,length n, which contains the b-spline coefficients.
 * @param k1 integer, giving the order of s(x) (order=degree+1)
 * @param x real, which contains the point where the derivatives must
 *          be evaluated.
 * @param d (out) array,length k1, containing the derivative values of s(x).
 * @param ier  error flag
 */
void spalde(
        const double *t,
        int n,
        const double *c,
        int k1,
        double x,
        double *d,
        int &ier
) {
    int l, nk1;

    // before starting computations a data check is made. if the input data
    // are invalid control is immediately repassed to the calling program.
    ier = 10;
    nk1 = n - k1;
    if (x < t[k1 - 1] || x > t[nk1 + 1 - 1]) return;
    // search for knot interval t(l) <= x < t(l+1)
    l = k1;
    _100:
    if (x < t[l + 1 - 1] || l == nk1) goto _200;
    l = l + 1;
    goto _100;
    _200:
    if (t[l - 1] >= t[l + 1 - 1]) return;
    ier = 0;
    //  calculate the derivatives.
    fpader(t, n, c, k1, x, l, d);
}

} // namespace dierckx
} // namespace interpolate
} // namespace lanre

#endif //LANRE_INTERPOLATE_SPALDE_HPP
