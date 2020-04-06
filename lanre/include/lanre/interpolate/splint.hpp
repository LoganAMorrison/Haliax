//
// Created by Logan Morrison on 4/4/20.
//

#ifndef LANRE_INTERPOLATE_SPLINT_HPP
#define LANRE_INTERPOLATE_SPLINT_HPP

#include "lanre/interpolate/fpintb.hpp"

namespace lanre {
namespace interpolate {
namespace dierckx {

double splint(
        const double *t,
        int n,
        const double *c,
        int k, double a,
        double b,
        double *wrk
) {
    int i, nk1;
    double result;

    nk1 = n - k - 1;
    // calculate the integrals wrk(i) of the normalized b-splines
    // ni,k+1(x), i=1,2,...nk1.
    fpintb(t, n, wrk, nk1, a, b);
    // calculate the integral of s(x).
    result = 0.0;
    for (i = 1; i <= nk1; i++) {
        result = result + c[i - 1] * wrk[i - 1];
    }
    return result;
}

} // namespace dierckx
} // namespace interpolate
} // namespace lanre

#endif //LANRE_INTERPOLATE_SPLINT_HPP