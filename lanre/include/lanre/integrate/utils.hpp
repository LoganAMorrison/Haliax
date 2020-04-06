//
// Created by Logan Morrison on 4/3/20.
//

#ifndef LANRE_INTEGRATE_UTILS_HPP
#define LANRE_INTEGRATE_UTILS_HPP

#include <cmath>
#include <iostream>

namespace lanre {
namespace integrate {
static void xerror(const char *mess, int nmess, int nerr, int level) {
    if (1 <= level) {
        std::cout << mess << std::endl;
        std::cout << "Error number: " << nerr << ", message lvl = " << level << std::endl;
    }
}

/**
 * Solve a general tridiagonal linear system.
 * @param n Order of the tridiagonal matrix
 * @param c Contains the subdiagonal of the tridiagonal matrix in
 *          entries C(2:N).  On output, C is destroyed
 * @param d On input, the diagonal of the matrix. On output, D is destroyed.
 * @param e Contains the superdiagonal of the tridiagonal matrix
 *          in entries E(1:N-1).  On output E is destroyed.
 * @param b On input, the right hand side. On output, the solution.
 * @param info Error flag
 *               0, normal value.
 *               k, the k-th element of the diagonal becomes exactly zero.
 *                  The routine returns if this error condition is detected.
 */
static void gtsl(int n, double *c, double *d, double *e, double *b, int &info) {
    int k;
    double t;

    info = 0;
    c[0] = d[0];

    if (2 <= n) {
        d[0] = e[0];
        e[0] = 0.0e+00;
        e[n - 1] = 0.0e+00;

        for (k = 1; k <= n - 1; k++) {
            // Find the larger of the two rows.
            if (abs(c[k - 1]) <= abs(c[k + 1 - 1])) {
                // Interchange rows.
                t = c[k + 1 - 1];
                c[k + 1 - 1] = c[k - 1];
                c[k - 1] = t;

                t = d[k + 1 - 1];
                d[k + 1 - 1] = d[k - 1];
                d[k - 1] = t;

                t = e[k + 1 - 1];
                e[k + 1 - 1] = e[k - 1];
                e[k - 1] = t;

                t = b[k + 1 - 1];
                b[k + 1 - 1] = b[k - 1];
                b[k - 1] = t;
            }

            // Zero elements.
            if (c[k - 1] == 0.0) {
                info = k;
                return;
            }

            t = -c[k + 1 - 1] / c[k - 1];
            c[k + 1 - 1] = d[k + 1 - 1] + t * d[k - 1];
            d[k + 1 - 1] = e[k + 1 - 1] + t * e[k - 1];
            e[k + 1 - 1] = 0.0e+00;
            b[k + 1 - 1] = b[k + 1 - 1] + t * b[k - 1];

        }
    }

    if (c[n - 1] == 0.0) {
        info = n;
        return;
    }

    // Back solve.
    b[n - 1] = b[n - 1] / c[n - 1];

    if (1 < n) {
        b[n - 1 - 1] = (b[n - 1 - 1] - d[n - 1 - 1] * b[n - 1]) / c[n - 1 - 1];

        for (k = n - 2; n >= 1; k--) {
            b[k - 1] = (b[k - 1] - d[k - 1] * b[k + 1 - 1] - e[k - 1] * b[k + 2 - 1]) / c[k - 1];
        }
    }
}

static void qpsrt(int limit, int last, int maxerr, double ermax, double *elist, int *iord, int nrmax) {
    double errmax, errmin;
    int i, ibeg, ido, isucc, j, jbnd, jupbn, k, lim;

    // check whether the list contains more than
    // two error estimates.
    if (last > 2) {
        goto _10;
    }
    iord[0] = 1;
    iord[1] = 2;
    goto _90;

    // this part of the routine is only executed if, due to a
    // difficult integrand, subdivision increased the error
    // estimate. in the normal case the insert procedure should
    // start after the nrmax-th largest error estimate.
    _10:
    errmax = elist[maxerr - 1];
    ido = nrmax - 1;
    for (i = 1; i <= ido; i++) {
        isucc = iord[nrmax - 1 - 1];
        if (errmax <= elist[isucc - 1]) goto _30;
        iord[nrmax - 1] = isucc;
        nrmax = nrmax - 1;
    }

    // compute the number of elements in the list to be maintained
    // in descending order. this number depends on the number of
    // subdivisions still allowed.
    _30:
    jupbn = last;
    if (last > (limit / 2 + 2)) {
        jupbn = limit + 3 - last;
    }
    errmin = elist[last - 1];

    // insert errmax by traversing the list top-down,
    // starting comparison from the element elist(iord(nrmax+1)).
    jbnd = jupbn - 1;
    ibeg = nrmax + 1;
    for (i = ibeg; i <= jbnd; i++) {
        isucc = iord[i - 1];
        if (errmax >= elist[isucc - 1]) goto _60;
        iord[i - 1 - 1] = isucc;
    }
    iord[jbnd - 1] = maxerr;
    iord[jupbn - 1] = last;
    goto _90;

    // insert errmin by traversing the list bottom-up.
    _60:
    iord[i - 1 - 1] = maxerr;
    k = jbnd;
    for (j = i; j <= jbnd; j++) {
        isucc = iord[k - 1];
        if (errmin < elist[isucc - 1]) goto _80;
        iord[k + 1 - 1] = isucc;
        k = k - 1;
    }
    iord[i - 1] = last;
    goto _90;

    _80:
    iord[k + 1 - 1] = last;

    // set maxerr and ermax.
    _90:
    maxerr = iord[nrmax - 1];
    ermax = elist[maxerr - 1];
}

static void qelg(int n, double *epstab, double &result, double &abserr, double *res3la, int &nres) {
    double delta1, delta2, delta3,
            epmach, epsinf, error, err1, err2, err3, e0, e1, e1abs, e2, e3,
            oflow, res, ss, tol1, tol2, tol3;
    int i, ib, ib2, ie, indx, k1, k2, k3, limexp, newelm, num;

    epmach = std::numeric_limits<double>::epsilon();
    oflow = std::numeric_limits<double>::max();
    nres = nres + 1;
    abserr = oflow;
    result = epstab[n - 1];
    if (n < 3) {
        goto _100;
    }
    limexp = 50;
    epstab[n + 2 - 1] = epstab[n - 1];
    newelm = (n - 1) / 2;
    epstab[n - 1] = oflow;
    num = n;
    k1 = n;

    for (i = 1; i <= newelm; i++) {
        k2 = k1 - 1;
        k3 = k1 - 2;
        res = epstab[k1 + 2 - 1];
        e0 = epstab[k3 - 1];
        e1 = epstab[k2 - 1];
        e2 = res;
        e1abs = fabs(e1);
        delta2 = e2 - e1;
        err2 = fabs(delta2);
        tol2 = fmax(abs(e2), e1abs) * epmach;
        delta3 = e1 - e0;
        err3 = fabs(delta3);
        tol3 = fmax(e1abs, fabs(e0)) * epmach;
        if (err2 > tol2 || err3 > tol3) {
            goto _10;
        }

        // if e0, e1 and e2 are equal to machine accuracy, convergence is assumed.
        result = res;
        abserr = err2 + err3;
        goto _100;

        _10:
        e3 = epstab[k1 - 1];
        epstab[k1 - 1] = e1;
        delta1 = e1 - e3;
        err1 = fabs(delta1);
        tol1 = fmax(e1abs, fabs(e3)) * epmach;

        // if two elements are very close to each other, omit
        // a part of the table by adjusting the value of n
        if (err1 <= tol1 || err2 <= tol2 || err3 <= tol3) goto _20;
        ss = 0.1e+01 / delta1 + 0.1e+01 / delta2 - 0.1e+01 / delta3;
        epsinf = abs(ss * e1);

        // test to detect irregular behaviour in the table, and
        // eventually omit a part of the table adjusting the value
        // of n.
        if (epsinf > 0.1e-03) goto _30;

        _20:
        n = i + i - 1;
        goto _50;

        // compute a new element and eventually adjust
        // the value of result.
        _30:
        res = e1 + 0.1e+01 / ss;
        epstab[k1 - 1] = res;
        k1 = k1 - 2;
        error = err2 + fabs(res - e2) + err3;
        if (error <= abserr) {
            abserr = error;
            result = res;
        }
    }
    _50:
    if (n == limexp) {
        n = 2 * (limexp / 2) - 1;
    }
    ib = 1;
    if ((num / 2) * 2 == num) {
        ib = 2;
    }
    ie = newelm + 1;
    for (i = 1; i <= ie; i++) {
        ib2 = ib + 2;
        epstab[ib - 1] = epstab[ib2 - 1];
        ib = ib2;
    }
    if (num == n) {
        goto _80;
    }
    indx = num - n + 1;
    for (i = 1; i <= n; i++) {
        epstab[i - 1] = epstab[indx - 1];
        indx = indx + 1;
    }

    _80:
    if (nres >= 4) {
        goto _90;
    }
    res3la[nres - 1] = result;
    abserr = oflow;
    goto _100;

    // compute error estimate
    _90:
    abserr = fabs(result - res3la[2]) + abs(result - res3la[1])
            + fabs(result - res3la[0]);
    res3la[0] = res3la[1];
    res3la[1] = res3la[2];
    res3la[2] = result;

    _100:
    abserr = fmax(abserr, 0.5e+01 * epmach * abs(result));
}

}
}

#endif //LANRE_INTEGRATE_UTILS_HPP
