//
// Created by Logan Morrison on 4/4/20.
//

#ifndef LANRE_INTERPOLATE_FPRATI_HPP
#define LANRE_INTERPOLATE_FPRATI_HPP

namespace lanre {
namespace interpolate {
namespace dierckx {

static double fprati(
        double &p1,
        double &f1,
        double &p2,
        double &f2,
        double &p3,
        double &f3
) {
    double h1, h2, h3, p;

    if (p3 > 0.) goto _10;
    //  value of p in case p3 = infinity.
    p = (p1 * (f1 - f3) * f2 - p2 * (f2 - f3) * f1) / ((f1 - f2) * f3);
    goto _20;
    //  value of p in case p3 ^= infinity.
    _10:
    h1 = f1 * (f2 - f3);
    h2 = f2 * (f3 - f1);
    h3 = f3 * (f1 - f2);
    p = -(p1 * p2 * h3 + p2 * p3 * h1 + p3 * p1 * h2) / (p1 * h1 + p2 * h2 + p3 * h3);
    // adjust the value of p1,f1,p3 and f3 such that f1 > 0 and f3 < 0.
    _20:
    if (f2 < 0.) goto _30;
    p1 = p2;
    f1 = f2;
    goto _40;
    _30:
    p3 = p2;
    f3 = f2;
    _40:
    return p;
}

} // namespace dierckx
} // namespace interpolate
} // namespace lanre

#endif //LANRE_INTERPOLATE_FPRATI_HPP
