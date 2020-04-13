//
// Created by Logan Morrison on 4/3/20.
//

#ifndef LANRE_INTEGRATE_QEXT_HPP
#define LANRE_INTEGRATE_QEXT_HPP

#include "lanre/integrate/base.hpp"
#include <cmath>
#include <array>

namespace lanre {
namespace integrate {

double qext(
        int *n,
        std::array<double, 52> &epstab,
        double *abserr,
        std::array<double, 3> &res3la,
        int *nres
) {
    static const double epmach = std::numeric_limits<double>::epsilon();
    static const double oflow = std::numeric_limits<double>::max();

    double delta1, delta2, delta3,
            epsinf, error, err1, err2, err3, e0, e1, e1abs, e2, e3,
            res, result, ss, tol1, tol2, tol3;
    int i, ib, ib2, ie, indx, k1, k2, k3, limexp, newelm, num;

    (*nres)++;
    *abserr = oflow;
    result = epstab[*n - 1];
    if (*n < 3) goto _100;
    limexp = 50;
    epstab[*n + 1] = epstab[*n - 1];
    newelm = (*n - 1) / 2;
    epstab[*n - 1] = oflow;
    num = *n;
    k1 = *n;

    for (i = 1; i <= newelm; i++) {
        k2 = k1 - 1;
        k3 = k1 - 2;
        res = epstab[k1 + 1];
        e0 = epstab[k3 - 1];
        e1 = epstab[k2 - 1];
        e2 = res;
        e1abs = std::abs(e1);
        delta2 = e2 - e1;
        err2 = std::abs(delta2);
        tol2 = std::max(std::abs(e2), e1abs) * epmach;
        delta3 = e1 - e0;
        err3 = std::abs(delta3);
        tol3 = std::max(e1abs, std::abs(e0)) * epmach;

        if (err2 > tol2 || err3 > tol3) goto _10;
        // if e0, e1 and e2 are equal to within machine
        // accuracy, convergence is assumed.
        // result = e2
        // abserr = abs(e1-e0)+abs(e2-e1)
        result = res;
        *abserr = err2 + err3;
        goto _100;

        _10:
        e3 = epstab[k1 - 1];
        epstab[k1 - 1] = e1;
        delta1 = e1 - e3;
        err1 = std::abs(delta1);
        tol1 = std::max(e1abs, std::abs(e3)) * epmach;
        // if two elements are very close to each other, omit
        // a part of the table by adjusting the value of n
        if (err1 <= tol1 || err2 <= tol2 || err3 <= tol3) goto _20;
        ss = 1.0 / delta1 + 1.0 / delta2 - 1.0 / delta3;
        epsinf = std::abs(ss * e1);

        // test to detect irregular behaviour in the table, and
        // eventually omit a part of the table adjusting the value
        // of n.
        if (epsinf > 0.1e-3) goto _30;

        _20:
        *n = i + i - 1;
        goto _50;

        // compute a new element and eventually adjust
        // the value of result.
        _30:
        res = e1 + 0.1e1 / ss;
        epstab[k1 - 1] = res;
        k1 = k1 - 2;
        error = err2 + std::abs(res - e2) + err3;
        if (error > *abserr) continue;
        *abserr = error;
        result = res;
    }

    //shift the table
    _50:
    if (*n == limexp) *n = 2 * (limexp / 2) - 1;
    ib = 1;
    if ((num / 2) * 2 == num) ib = 2;
    ie = newelm + 1;
    for (i = 1; i <= ie; i++) {
        ib2 = ib + 2;
        epstab[ib - 1] = epstab[ib2 - 1];
        ib = ib2;
    }
    if (num == *n) goto _80;
    indx = num - *n + 1;
    for (i = 1; i <= *n; i++) {
        epstab[i - 1] = epstab[indx - 1];
        indx = indx + 1;
    }

    _80:
    if (*nres < 4) goto _90;
    res3la[*nres - 1] = result;
    *abserr = oflow;
    goto _100;

    // compute error estimate
    _90:
    *abserr = std::abs(result - res3la[2]) + std::abs(result - res3la[1])
            + std::abs(result - res3la[0]);
    res3la[0] = res3la[1];
    res3la[1] = res3la[2];
    res3la[2] = result;

    _100:
    *abserr = std::max(*abserr, 5.0 * epmach * std::abs(result));
    return result;
}

}
}

#endif //LANRE_INTEGRATE_QEXT_HPP
