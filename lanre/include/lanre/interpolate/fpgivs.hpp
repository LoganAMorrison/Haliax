// Created by Logan Morrison on 4/4/20.
//
// Adapted from the FORTRAN pacakge Dierckx written by:
//      p.dierckx
//          dept. computer science, k.u. leuven
//          celestijnenlaan 200a, b-3001 heverlee, belgium.
//          e-mail : Paul.Dierckx@cs.kuleuven.ac.be
//

#ifndef LANRE_INTERPOLATE_FPGIVS_HPP
#define LANRE_INTERPOLATE_FPGIVS_HPP

#include <cmath>

namespace lanre {
namespace interpolate {
namespace dierckx {

/**
 * Compute a givens rotation angles given piv and ww.
 */
void fpgivs(const double piv, double &ww, double &cos, double &sin) {
    double dd, one, store;
    one = 1.0;
    store = abs(piv);
    if (store >= ww) {
        dd = store * sqrt(one + pow(ww / piv, 2));
    } else {
        dd = ww * sqrt(one + pow(piv / ww, 2));
    }
    cos = ww / dd;
    sin = piv / dd;
    ww = dd;
}

} // namespace dierckx
} // namespace interpolate
} // namespace lanre

#endif //LANRE_INTERPOLATE_FPGIVS_HPP

