// Created by Logan Morrison on 4/4/20.
//
// Adapted from the FORTRAN pacakge Dierckx written by:
//      p.dierckx
//          dept. computer science, k.u. leuven
//          celestijnenlaan 200a, b-3001 heverlee, belgium.
//          e-mail : Paul.Dierckx@cs.kuleuven.ac.be
//

#ifndef LANRE_INTERPOLATE_FPCURO_HPP
#define LANRE_INTERPOLATE_FPCURO_HPP

namespace lanre {
namespace interpolate {
namespace dierckx {

/**
 * finds the real zeros of a cubic polynomial
 * p(x) = a*x**3+b*x**2+c*x+d.
 */
void fpcuro(const double &a, const double &b, const double &c, const double &d, double *x, int &n) {
    int i;
    double a1, b1, c1, df, disc, d1, e3, f, four, half, ovfl, pi3, p3, q, r,
            step, tent, three, two, u, u1, u2, y;

    // set constants
    two = 0.2e+01;
    three = 0.3e+01;
    four = 0.4e+01;
    ovfl = 0.1e+05;
    half = 0.5e+0;
    tent = 0.1e+0;
    e3 = tent / 0.3;
    pi3 = atan(0.1e+01) / 0.75;
    a1 = fabs(a);
    b1 = fabs(b);
    c1 = fabs(c);
    d1 = fabs(d);
    // test whether p(x) is a third degree polynomial.
    if (std::max({b1, c1, d1}) < a1 * ovfl) goto _300;
    // test whether p(x) is a second degree polynomial.
    if (fmax(c1, d1) < b1 * ovfl) goto _200;
    // test whether p(x) is a first degree polynomial.
    if (d1 < c1 * ovfl) goto _100;
    // p(x) is a constant function.
    n = 0;
    return;
    //  p(x) is a first degree polynomial.
    _100:
    n = 1;
    x[0] = -d / c;
    goto _500;

    //  p(x) is a second degree polynomial.
    _200:
    disc = c * c - four * b * d;
    n = 0;
    if (disc < 0.0) return;
    n = 2;
    u = sqrt(disc);
    b1 = b + b;
    x[0] = (-c + u) / b1;
    x[1] = (-c - u) / b1;
    goto _500;
    // p(x) is a third degree polynomial.

    _300:
    b1 = b / a * e3;
    c1 = c / a;
    d1 = d / a;
    q = c1 * e3 - b1 * b1;
    r = b1 * b1 * b1 + (d1 - b1 * c1) * half;
    disc = q * q * q + r * r;
    if (disc > 0.0) goto _400;
    u = sqrt(abs(q));
    if (r < 0.0) u = -u;
    p3 = atan2(sqrt(-disc), abs(r)) * e3;
    u2 = u + u;
    n = 3;
    x[0] = -u2 * cos(p3) - b1;
    x[1] = u2 * cos(pi3 - p3) - b1;
    x[2] = u2 * cos(pi3 + p3) - b1;
    goto _500;

    _400:
    u = sqrt(disc);
    u1 = -r + u;
    u2 = -r - u;
    n = 1;
    x[0] = std::copysign(pow(abs(u1), e3), u1) + std::copysign(pow(abs(u2), e3), u2) - b1;

    // apply a newton iteration to improve the accuracy of the roots.
    _500:
    for (i = 1; i <= n; i++) {
        y = x[i - 1];
        f = ((a * y + b) * y + c) * y + d;
        df = (three * a * y + two * b) * y + c;
        step = 0.0;
        if (fabs(f) < fabs(df) * tent) step = f / df;
        x[i - 1] = y - step;
    }
}


} // namespace dierckx
} // namespace interpolate
} // namespace lanre

#endif //LANRE_INTERPOLATE_FPCURO_HPP
