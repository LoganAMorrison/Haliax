//
// Created by Logan Morrison on 3/20/20.
//

#ifndef LANRE_ROOTS_HPP
#define LANRE_ROOTS_HPP

#include <cmath>
#include <stdexcept>

namespace lanre {

/**
 * Find the root of a function given a bracketing interval using
 * Van Wijngaarden-Dekker-Brent method.
 * @tparam T Function type
 * @param func Function to find root of
 * @param x1 Lower bracket
 * @param x2 Upper bracket
 * @param tol Tolerance
 * @param max_iters Maximum number of iterations
 * @return Value of the root
 */
template<class T>
double brent(T &func, const double x1, const double x2, const double tol, const int max_iters = 100) {
    const double EPS = std::numeric_limits<double>::epsilon();
    double a = x1;
    double b = x2;
    double c = x2;
    double d{};
    double e{};
    double fa = func(a);
    double fb = func(b);
    double fc{};
    double p{};
    double q{};
    double r{};
    double s{};
    double tol1{};
    double xm{};


    if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
        throw std::runtime_error("Function values at endpoints don't have oppositive signs.");
    }
    fc = fb;
    for (int iter = 0; iter < max_iters; iter++) {
        // Rename a, b, c and adjust the bounding interval d.
        if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
            c = a;
            fc = fa;
            e = d = b - a;
        }
        if (abs(fc) < abs(fb)) {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }
        // Convergence check
        tol1 = 2.0 * EPS * abs(b) + 0.5 * tol;
        xm = 0.5 * (c - b);
        if (abs(xm) <= tol1 || fb == 0.0) {
            return b;
        }
        if (abs(e) >= tol1 && abs(fa) > abs(fb)) {
            //Attempt inverse quadratic interpolation.
            s = fb / fa;
            if (a == c) {
                p = 2.0 * xm * s;
                q = 1.0 - s;
            } else {
                q = fa / fc;
                r = fb / fc;
                p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
                q = (q - 1.0) * (r - 1.0) * (s - 1.0);
            }
            // Check whether in bounds.
            if (p > 0.0) {
                q = -q;
            }
            p = abs(p);
            double min1 = 3.0 * xm * q - abs(tol1 * q);
            double min2 = abs(e * q);
            if (2.0 * p < (min1 < min2 ? min1 : min2)) {
                e = d;
                d = p / q;
            } else {
                d = xm;
                e = d;
            }
        } else {
            d = xm;
            e = d;
        }
        a = b;
        fa = fb;
        if (abs(d) > tol1)
            b += d;
        else
            b += std::copysign(tol1, xm);
        fb = func(b);
    }
    throw std::runtime_error("Brent: Maximum number of iterations exceeded.");
}
}

#endif //LANRE_ROOTS_HPP
