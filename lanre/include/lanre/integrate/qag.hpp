//
// Created by Logan Morrison on 3/30/20.
//

#ifndef LANRE_INTEGRATE_QAG_HPP
#define LANRE_INTEGRATE_QAG_HPP

#include "lanre/integrate/qsrt.hpp"
#include "lanre/integrate/qext.hpp"
#include "lanre/integrate/qk.hpp"
#include <cmath>
#include <vector>
#include <array>
#include <numeric>

namespace lanre {
namespace integrate {

template<class Integrand, size_t Limit>
double qage(Integrand f, double a, double b, double epsabs, double epsrel,
            int irule, double *abserr, int *neval, int *ier, int *last) {
    double area, area1, area2, area12, a1, a2, b1, b2, c, defabs;
    double defab1, defab2, errbnd, errmax, error1, error2;
    double erro12, errsum, resabs, result;
    double alist[Limit], blist[Limit], rlist[Limit], elist[Limit];
    int iroff1, iroff2, k, keyf, maxerr, nrmax, iord[Limit], limit;

    const double epmach = std::numeric_limits<double>::epsilon();
    const double uflow = std::numeric_limits<double>::min();

    limit = Limit - 1;
    *ier = 0;
    *neval = 0;
    *last = 0;
    result = 0.0;
    *abserr = 0.0;
    alist[0] = a;
    blist[0] = b;
    rlist[0] = 0.0;
    elist[0] = 0.0;
    iord[0] = 0;
    defabs = 0.0;
    resabs = 0.0;
    if ((epsabs < 0.0) && (epsrel < 0.0))
        *ier = 6;
    if (*ier == 6) return result;

    /* First approximation to the integral. */
    keyf = irule;
    if (irule <= 0) keyf = 1;
    if (irule >= 7) keyf = 6;
    c = keyf;
    *neval = 0;
    switch (keyf) {
        case 1:
            result = qk15(f, a, b, abserr, &defabs, &resabs);
            break;
        case 2:
            result = qk21(f, a, b, abserr, &defabs, &resabs);
            break;
        case 3:
            result = qk31(f, a, b, abserr, &defabs, &resabs);
            break;
        case 4:
            result = qk41(f, a, b, abserr, &defabs, &resabs);
            break;
        case 5:
            result = qk51(f, a, b, abserr, &defabs, &resabs);
            break;
        case 6:
            result = qk61(f, a, b, abserr, &defabs, &resabs);
            break;
        default:
            result = 0.0;
            break;
    }
    *last = 0;
    rlist[0] = result;
    elist[0] = *abserr;
    iord[0] = 0;

    /* Test on accuracy. */
    errbnd = std::max(epsabs, epsrel * fabs(result));
    if ((*abserr <= 50.0 * epmach * defabs) && (*abserr > errbnd))
        *ier = 2;
    if (limit == 0) *ier = 1;
    if ((*ier != 0) || ((*abserr <= errbnd) && (*abserr != resabs)) ||
            (*abserr == 0.0))
        goto _60;

    /* Initialization. */
    errmax = *abserr;
    maxerr = 0;
    area = result;
    errsum = *abserr;
    nrmax = 0;
    iroff1 = 0;
    iroff2 = 0;

    /* Main Loop. */
    for (*last = 1; *last <= limit; (*last)++) {
        /* Bisect the subinterval with the largest error estimate. */
        a1 = alist[maxerr];
        b1 = 0.5 * (alist[maxerr] + blist[maxerr]);
        a2 = b1;
        b2 = blist[maxerr];
        switch (keyf) {
            case 1:
                area1 = qk15(f, a1, b1, &error1, &resabs, &defab1);
                area2 = qk15(f, a2, b2, &error2, &resabs, &defab2);
                break;
            case 2:
                area1 = qk21(f, a1, b1, &error1, &resabs, &defab1);
                area2 = qk21(f, a2, b2, &error2, &resabs, &defab2);
                break;
            case 3:
                area1 = qk31(f, a1, b1, &error1, &resabs, &defab1);
                area2 = qk31(f, a2, b2, &error2, &resabs, &defab2);
                break;
            case 4:
                area1 = qk41(f, a1, b1, &error1, &resabs, &defab1);
                area2 = qk41(f, a2, b2, &error2, &resabs, &defab2);
                break;
            case 5:
                area1 = qk51(f, a1, b1, &error1, &resabs, &defab1);
                area2 = qk51(f, a2, b2, &error2, &resabs, &defab2);
                break;
            case 6:
                area1 = qk61(f, a1, b1, &error1, &resabs, &defab1);
                area2 = qk61(f, a2, b2, &error2, &resabs, &defab2);
                break;
        }

        /* Improve previous approximations to integral and error,
        and test for accuracy. */
        (*neval) += 1;
        area12 = area1 + area2;
        erro12 = error1 + error2;
        errsum = errsum + erro12 - errmax;
        area = area + area12 - rlist[maxerr];
        if ((defab1 != error1) && (defab2 != error2)) {
            if ((fabs(rlist[maxerr] - area12) <= 1.0e-5 * fabs(area12)) &&
                    (erro12 >= .99 * errmax))
                iroff1++;
            if ((*last > 9) && (erro12 > errmax))
                iroff2++;
        }
        rlist[maxerr] = area1;
        rlist[*last] = area2;
        errbnd = std::max(epsabs, epsrel * fabs(area));
        if (errsum > errbnd) {

            /* Test for roundoff error and eventually set error flag. */
            if ((iroff1 > 6) || (iroff2 > 20))
                *ier = 2;

            /* Set error flag in the case that the number of subintervals
                equals the limit. */
            if (*last == limit)
                *ier = 1;

            /* Set error flag in the case of bad integrand behavior at a
                point of the integration range. */
            if (std::max(fabs(a1), fabs(b2)) <= (1.0 + c * 1000.0 * epmach) *
                    (fabs(a2) + 1.0e4 * uflow))
                *ier = 3;
        }
        /* Append the newly-created intervals to the list. */

        if (error2 <= error1) {
            alist[*last] = a2;
            blist[maxerr] = b1;
            blist[*last] = b2;
            elist[maxerr] = error1;
            elist[*last] = error2;
        } else {
            alist[maxerr] = a2;
            alist[*last] = a1;
            blist[*last] = b1;
            rlist[maxerr] = area2;
            rlist[*last] = area1;
            elist[maxerr] = error2;
            elist[*last] = error1;
        }

        /* Call DQSORT to maintain the descending ordering in the list of
            error estimates and select the subinterval with the
            largest error estimate (to be bisected next). */

        qsrt(limit, *last, &maxerr, &errmax, elist, iord, &nrmax);
        if ((*ier != 0) || (errsum <= errbnd)) break;
    }

    /* Compute final result. */

    result = 0.0;
    for (k = 0; k <= *last; k++) {
        result += rlist[k];
    }
    *abserr = errsum;
    _60:
    if (keyf != 1)
        *neval = (10 * keyf + 1) * (2 * (*neval) + 1);
    else
        *neval = 30 * (*neval) + 15;

    return result;
}


template<class Integrand, size_t Limit = 500>
double qag(
        Integrand f,
        double a,
        double b,
        double epsabs,
        double epsrel,
        int key,
        double *abserr,
        int *neval,
        int *ier,
        int *last,
        std::array<double, Limit> &alist,
        std::array<double, Limit> &blist,
        std::array<double, Limit> &rlist,
        std::array<double, Limit> &elist,
        std::array<int, Limit> &iord
) {
    // check validity of lenw.
    *ier = 6;
    *neval = 0;
    *last = 0;
    double result = 0.0;
    *abserr = 0.0;
    if (Limit >= 1) {
        // prepare call for dqage.
        result = qage(f, a, b, epsabs, epsrel, key, abserr, neval,
                      ier, alist, blist, rlist, elist, iord, last);
    }
    return result;
}


template<class Integrand, size_t Limit = 500>
double qag(
        Integrand f,
        double a,
        double b,
        double epsabs,
        double epsrel,
        int key,
        double *abserr
) {
    // check validity of lenw.
    int ier = 6;
    int neval = 0;
    int last = 0;
    double result = 0.0;
    *abserr = 0.0;
    if (Limit >= 1) {
        // prepare call for dqage.
        std::array<double, Limit> alist = {0};
        std::array<double, Limit> blist = {0};
        std::array<double, Limit> rlist = {0};
        std::array<double, Limit> elist = {0};
        std::array<int, Limit> iord = {0};

        result = qage(f, a, b, epsabs, epsrel, key, abserr, &neval,
                      &ier, alist, blist, rlist, elist, iord, &last);
    }
    return result;
}

}

}

#endif //LANRE_INTEGRATE_QAG_HPP
