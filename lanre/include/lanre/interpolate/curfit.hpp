// Created by Logan Morrison on 4/4/20.
//
// Adapted from the FORTRAN pacakge Dierckx written by:
//      p.dierckx
//          dept. computer science, k.u. leuven
//          celestijnenlaan 200a, b-3001 heverlee, belgium.
//          e-mail : Paul.Dierckx@cs.kuleuven.ac.be
//

#ifndef LANRE_INTERPOLATE_CURFIT_HPP
#define LANRE_INTERPOLATE_CURFIT_HPP

#include "lanre/interpolate/fpchec.hpp"
#include "lanre/interpolate/fpcurf.hpp"

namespace lanre {
namespace interpolate {
namespace dierckx {

/**
 * given the set of data points (x(i),y(i)) and the set of positive
 * numbers w(i),i=1,2,...,m,subroutine curfit determines a smooth spline
 * approximation of degree k on the interval xb <= x <= xe.
 * if iopt=-1 curfit calculates the weighted least-squares spline
 * according to a given set of knots.
 * if iopt>=0 the number of knots of the spline s(x) and the position
 * t(j),j=1,2,...,n is chosen automatically by the routine. the smooth-
 * ness of s(x) is then achieved by minimalizing the discontinuity
 * jumps of the k-th derivative of s(x) at the knots t(j),j=k+2,k+3,...,
 * n-k-1. the amount of smoothness is determined by the condition that
 * f(p)=sum((w(i)*(y(i)-s(x(i))))**2) be <= s, with s a given non-
 * negative constant, called the smoothing factor.
 * the fit s(x) is given in the b-spline representation (b-spline coef-
 * ficients c(j),j=1,2,...,n-k-1) and can be evaluated by means of
 * subroutine splev.
 *
 *
 */
void curfit(
        int &iopt,
        int &m,
        double *x,
        double *y,
        double *w,
        double &xb,
        double &xe,
        int &k,
        double &s,
        int &nest,
        int &n,
        double *t,
        double *c,
        double &fp,
        double *fpint,
        double *z,
        double **a,
        double **b,
        double **g,
        double **q,
        int *iwrk,
        int &ier
) {
    double tol;
    int i, j, k1, k2, maxit, nmin;

    // we set up the parameters tol and maxit
    maxit = 20;
    tol = 0.1e-02;
    // before starting computations a data check is made. if the input data
    // are invalid, control is immediately repassed to the calling program.
    ier = 10;
    if (k <= 0 || k > 5) return;
    k1 = k + 1;
    k2 = k1 + 1;
    if ((iopt < (-1)) || iopt > 1) return;
    nmin = 2 * k1;
    if (m < k1 || nest < nmin) return;
    //lwest = m * k1 + nest * (7 + 3 * k);
    //if (lwrk < lwest) return;
    if (xb > x[0] || xe < x[m - 1]) return;
    for (i = 2; i <= m; i++) {
        if (x[i - 2] > x[i - 1]) return;
    }
    if (iopt >= 0) goto _30;
    if (n < nmin || n > nest) return;
    j = n;
    for (i = 1; i <= k1; i++) {
        t[i - 1] = xb;
        t[j - 1] = xe;
        j = j - 1;
    }
    fpchec(x, m, t, n, k, ier);
    if (ier == 0) goto _40;
    return;
    _30:
    if (s < 0.0) return;
    if (s == 0.0 && nest < (m + k1)) return;
    // we partition the working space and determine the spline approximation.
    _40:

    fpcurf(iopt, x, y, w, m, xb, xe, k, s, nest, tol, maxit, k1, k2, n, t, c, fp,
           fpint, z, a, b, g, q, iwrk, ier);
}

} // namespace dierckx
} // namespace interpolate
} // namespace lanre

#endif //LANRE_INTERPOLATE_CURFIT_HPP