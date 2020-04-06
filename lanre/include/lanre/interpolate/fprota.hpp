//
// Created by Logan Morrison on 4/4/20.
//

#ifndef LANRE_INTERPOLATE_FPROTA_HPP
#define LANRE_INTERPOLATE_FPROTA_HPP

namespace lanre {
namespace interpolate {
namespace dierckx {

/**
 * applies a givens rotation to a and b
 */
static void fprota(const double cos, const double sin, double &a, double &b) {
    double stor1, stor2;
    stor1 = a;
    stor2 = b;
    b = cos * stor2 + sin * stor1;
    a = cos * stor1 - sin * stor2;
}

} // namespace dierckx
} // namespace interpolate
} // namespace lanre

#endif //LANRE_INTERPOLATE_FPROTA_HPP
